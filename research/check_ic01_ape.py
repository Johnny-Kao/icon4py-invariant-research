"""IC-01: APE-configured grid, reference level selection and native Exner stencil.
Not a coupled ICON/ICON4Py forecast. The Fortran reference is the published
level-selection loop, not a build of ICON itself.
"""
import dataclasses
import pathlib
import subprocess
import tempfile

import gt4py.next as gtx
import numpy as np
from icon4py.model.common import dimension as dims
from icon4py.model.common.grid import vertical as v
from icon4py.model.common.type_alias import vpfloat, wpfloat
from icon4py.model.atmosphere.dycore.stencils.update_dynamical_exner_time_increment import (
    update_dynamical_exner_time_increment,
)

# Public ICON run/run_ICONAM_R2B4_APE settings; other grid values are
# ICON4Py defaults, not an archived Fortran savepoint.
cfg = v.VerticalGridConfig(
    num_levels=90, model_top_height=60000., lowest_layer_thickness=20.,
    flat_height=20000., stretch_factor=0.6, rayleigh_damping_height=35000.,
    htop_moist_proc=30000.,
)
a,b = v.get_vct_a_and_vct_b(cfg, None)
z = a.asnumpy()
mid = (z[:-1] + z[1:]) * .5
binding_cfg = dataclasses.replace(cfg, htop_moist_proc=22500.)
expected = int(v.VerticalGrid(config=cfg, vct_a=a, vct_b=b).kstart_moist)
actual = int(v.VerticalGrid(config=binding_cfg, vct_a=a, vct_b=b).kstart_moist)

# Independently compiled equivalent of the Fortran selection loop, no nesting.
src = """program reference
use iso_fortran_env, only: real64
implicit none
real(real64) :: z(91), h
integer :: jk, kstart, i
read(*,*) z
do i=1,2
  read(*,*) h
  kstart=1
  do jk=1,90
    if (0.5_real64*(z(jk)+z(jk+1)) < h) then
      kstart=jk
      exit
    endif
  enddo
  print *, kstart
enddo
end program
"""
with tempfile.TemporaryDirectory() as directory:
    p = pathlib.Path(directory)
    (p/"reference.f90").write_text(src)
    subprocess.run(["gfortran", str(p/"reference.f90"), "-o", str(p/"reference")], check=True)
    result = subprocess.run([str(p/"reference")], input=" ".join(map(str,z))+"\n30000\n22500\n", text=True, capture_output=True, check=True)
ref30,ref225 = map(int,result.stdout.split())
print(f"Fortran reference (1-based) 30km={ref30}, 22.5km={ref225}; Python (0-based) correct={expected}, binding-default={actual}")
assert (expected,actual)==(ref30-1,ref225-1)
assert expected < actual
print(f"Different full levels={actual-expected}; midpoint range={mid[actual-1]:.6f}..{mid[expected]:.6f} m")
# Identical synthetic state in both paths. Predictor initialization mimics
# its published vertical mask; actual corrector arithmetic uses native stencil.
shape=(2,90)
initial=np.broadcast_to(np.linspace(.25,.99,90),shape).astype(wpfloat).copy()
next_exner=initial+wpfloat(.0002)
physics=np.full(shape,1e-7,dtype=vpfloat)
dt=wpfloat(60.); nsub=wpfloat(5.)
def execute(k):
    accum=np.zeros(shape,dtype=vpfloat)
    accum[:,k:]=initial[:,k:].astype(vpfloat)
    make=lambda x:gtx.as_field((dims.CellDim,dims.KDim),x)
    out=make(accum.copy())
    update_dynamical_exner_time_increment.with_backend(None)(
        exner=make(next_exner.copy()),ddt_exner_phy=make(physics.copy()),
        exner_dyn_incr=out,ndyn_substeps_var=nsub,dtime=dt,
        horizontal_start=gtx.int32(0),horizontal_end=gtx.int32(2),
        vertical_start=gtx.int32(k),vertical_end=gtx.int32(90),
        offset_provider={},
    )
    oracle=accum.copy()
    oracle[:,k:]=(next_exner[:,k:]-(accum[:,k:].astype(wpfloat)+nsub*dt*physics[:,k:].astype(wpfloat))).astype(vpfloat)
    np.testing.assert_allclose(out.asnumpy(),oracle,rtol=1e-6,atol=1e-8)
    return out.asnumpy()
good=execute(expected)
old=execute(actual)
delta=np.abs(good-old)
print(f"Native stencil differing cells={np.count_nonzero(delta)}; max abs Exner increment difference={delta.max():.12g}")
assert np.count_nonzero(delta)==2*(actual-expected)
np.testing.assert_array_equal(good[:,actual:],old[:,actual:])
np.testing.assert_array_equal(execute(actual),old)
print("PASS: default control; nondefault selection mismatch; native stencil differential; common layers identical")
