# IC-01: coupled ICON validation

This branch is a minimal source snapshot of the official public ICON commit
`7de853356c97db40d9f57105c58e72c3bf822cb0`. It is not a full ICON build.
The separate `ic01/icon4py` branch provides the companion binding change.
The generated ABI must be rebuilt with both changes present.

## Research checks

Run `python3 test/ic01/check_caller.py` against the caller and upstream Fortran
checker. The compiled Fortran/C/Python differential is recorded on Research
main (Actions run 45589262645).

## Native ICON regression (not run)

Use an existing supported single-domain, single-block ICON4Py experiment on
the project's CSCS Santis runner (or equivalent configured ICON environment).
Build ICON with `--enable-icon4py` using the companion Python bindings, regenerate
the interface, and source the generated `externals/icon4py/run.env`.
Use the documented VERIFY mode in the experiment's existing &run_nml.
Run the same experiment for 22500, 30000 and 15000 m in
`&nonhydrostatic_nml` (maintaining other configuration constraints, including
the water-vapour substepping height). Compare the existing VERIFY diagnostic
`exner_dyn_incr` for both cores, including default-control results.
Do not change the existing operational experiment defaults.

The full ICON build, official experiment execution, and tolerance comparison
are not yet verified. This directory is Research evidence, not a proposal to
add a new upstream ICON test framework.
