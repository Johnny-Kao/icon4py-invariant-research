"""IC-01 Research: generate and execute the actual grid_init FFI on CPU.
Uses current ICON4Py Research source; horizontal/MPI setup is stubbed.
Fortran caller is a small driver, not a complete ICON simulation.
"""
import argparse
import os
import pathlib
import subprocess
import sys
import tempfile

from icon4py.bindings import grid_wrapper
from icon4py.tools import py2fgen
from icon4py.tools.py2fgen import _codegen, _generator
import ic01_bridge_fixture as fixture

parser = argparse.ArgumentParser()
parser.add_argument("mode", choices=("old", "fixed"))
mode = parser.parse_args().mode
fixed = mode == "fixed"
name = "ic01_bindings"
args = grid_wrapper.grid_init.param_descriptors
assert ("htop_moist_proc" in args) == fixed
plugin = py2fgen.get_cffi_description(
    [grid_wrapper.grid_init, fixture.prepare_fixture, fixture.read_result], name
)
sources = py2fgen.render(plugin)
types = {
    py2fgen.INT32: "integer(c_int)",
    py2fgen.INT64: "integer(c_long)",
    py2fgen.FLOAT64: "real(c_double)",
    py2fgen.FLOAT32: "real(c_float)",
    py2fgen.BOOL: "logical(c_bool)",
}
decls, assignments, keywords = [], [], []
scalars = {
    "lowest_layer_thickness": "20.0_c_double",
    "model_top_height": "60000.0_c_double",
    "stretch_factor": "0.6_c_double",
    "flat_height": "20000.0_c_double",
    "rayleigh_damping_height": "35000.0_c_double",
    "mean_cell_area": "1.0_c_double",
    "vertical_size": "90_c_int",
    "limited_area": ".false._c_bool",
    "backend": "0_c_int",
    "comm_id": "0_c_int",
    "htop_moist_proc": "htop_moist_proc",
}
for key, arg in args.items():
    if isinstance(arg, py2fgen.ArrayParamDescriptor):
        shape = "91" if key == "vct_a" else ",".join(["2"] * arg.rank)
        decls.append(f"{types[arg.dtype]} :: {key}({shape})")
        if key != "vct_a":
            value = ".true._c_bool" if arg.dtype == py2fgen.BOOL else "1"
            assignments.append(f"{key} = {value}")
        keywords.append(f"{key}={key}")
    else:
        value = scalars.get(key, "2_c_int")
        keywords.append(f"{key}={value}")
call = "call grid_init( &\n" + ", &\n".join(keywords + ["rc=rc"]) + ")"
driver = "\n".join([
    "program research", "use iso_c_binding", f"use {name}",
    "implicit none", "integer(c_int) :: rc, i, k, expected, observed(1)",
    "real(c_double) :: htop_moist_proc, requested(3), out(90), reference(90)",
    *decls, *assignments,
    "call prepare_fixture(vct_a, rc)", "if (rc /= 1) stop 1",
    "requested = [22500.0_c_double, 30000.0_c_double, 15000.0_c_double]",
    "do i=1,3", "htop_moist_proc=requested(i)",
    call, "if (rc /= 1) stop 2",
    "call read_result(observed, out, rc)", "if (rc /= 1) stop 3",
    "expected=0",
    "do k=1,90",
    "  if (0.5_c_double*(vct_a(k)+vct_a(k+1)) < htop_moist_proc) then",
    "    expected=k-1", "    exit", "  endif", "enddo",
    "reference=0.0_c_double",
    "reference(expected+1:90)=0.00017_c_double",
    "print *, 'requested=', htop_moist_proc, 'expected=', expected, 'observed=', observed(1)",
    "print *, 'different increments=', count(abs(out-reference)>1.0e-10_c_double)",
    "if (i == 1 .or. " + (".true." if fixed else ".false.") + ") then",
    "  if (observed(1) /= expected) stop 4",
    "  if (any(abs(out-reference)>1.0e-10_c_double)) stop 5",
    "else",
    "  if (observed(1) == expected) stop 6",
    "  if (count(abs(out-reference)>1.0e-10_c_double) /= abs(observed(1)-expected)) stop 7",
    "endif", "enddo", "end program research", "",
])
with tempfile.TemporaryDirectory() as d:
    p = pathlib.Path(d)
    (p / "bindings.f90").write_text(sources.f90)
    (p / "driver.f90").write_text(driver)
    _generator.generate_and_compile_cffi_plugin(
        name, _codegen.generate_c_header(plugin), sources.py, p
    )
    subprocess.run(["gfortran", "-cpp", "-ffree-line-length-none", "bindings.f90",
                    "driver.f90", "-L.", f"-l{name}", "-Wl,-rpath="+str(p),
                    "-o", "research"], cwd=p, check=True)
    env = os.environ.copy()
    env["PYTHONPATH"] = os.pathsep.join(filter(None, [*sys.path, env.get("PYTHONPATH", "")]))
    subprocess.run([str(p / "research")], cwd=p, env=env, check=True)
print(f"PASS {mode}: generated Fortran/C/Python bridge and native Exner arithmetic")
