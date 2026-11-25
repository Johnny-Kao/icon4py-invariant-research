"""Research-only check of the pinned ICON caller and official style checker.

The compiled bridge differential lives on the Research main branch.
This does not execute the full ICON VERIFY mode.
"""
import pathlib
import re
import subprocess
import sys

root = pathlib.Path(__file__).resolve().parents[2]
caller = root / "src/atm_dyn_iconam/mo_icon4py_interfaces.f90"
source = caller.read_text()
imports = re.findall(r"(?im)^\s*USE mo_nonhydrostatic_config,\s*ONLY:[^\n]*(?:\n[ \t]*&[^\n]*)*", source)
assert len(imports) == 1
assert re.search(r"\bhtop_moist_proc\b", imports[0])
call = re.search(r"(?ims)^\s*CALL grid_init\s*\(.*?\brc\s*=\s*rc\s*\)", source)
assert call is not None
assert len(re.findall(r"\bhtop_moist_proc\s*=\s*htop_moist_proc\b", call.group(), re.I)) == 1
subprocess.run(
    [sys.executable, str(root / "scripts/ci_scripts/check_fortran.py"), str(caller)],
    check=True,
)
print("PASS: ICON caller forwards configured height; official Fortran checker passes")
