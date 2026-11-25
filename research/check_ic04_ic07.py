"""Focused IC-04/07 controls, not a complete coupled simulation."""
from types import SimpleNamespace
import numpy as np
import gt4py.next as gtx
from icon4py.model.atmosphere.dycore import solve_nonhydro, dycore_utils
from icon4py.model.common import dimension as dims, model_options

for v,r in [(0.25,-0.1),(0.125,0.2),(0.0,-0.25)]:
    config=solve_nonhydro.NonHydrostaticConfig(veladv_offctr=v,rhotheta_offctr=r)
    params=solve_nonhydro.NonHydrostaticParams(config)
    np.testing.assert_allclose(params.advection_explicit_weight_parameter,0.5-v)
    np.testing.assert_allclose(params.advection_implicit_weight_parameter,0.5+v)
    np.testing.assert_allclose(params.rhotheta_implicit_weight_parameter,0.5+r)
try:
    solve_nonhydro.NonHydrostaticConfig(extra_diffu=False)
except NotImplementedError as exc:
    assert "extra_diffu" in str(exc)
else:
    raise AssertionError("unsupported extra_diffu silently accepted")
print("IC-04 PASS: 3 nondefault coefficient configurations and explicit unsupported flag")

w=np.array([0.,.05,.25,1.],dtype=np.float64)
field=gtx.as_field((dims.KHalfDim,),w)
target=gtx.as_field((dims.KHalfDim,),np.zeros_like(w))
native=model_options.setup_program(program=dycore_utils.compute_rayleigh_damping_factor,
                                  backend=None,constant_args={"rayleigh_w":field})
calls=[]
def compute(**kwargs):
    calls.append(float(kwargs["dtime"]))
    native(**kwargs)
state=SimpleNamespace(_dtime_previous_substep=0.,rayleigh_damping_factor=target,
                      _compute_rayleigh_damping_factor=compute)
for dt in [1.,2.,2.,.5,1.]:
    out=solve_nonhydro.SolveNonhydro._get_rayleigh_damping_factor(state,dt)
    np.testing.assert_allclose(out.asnumpy(),1/(1+dt*w))
assert calls==[1.,2.,.5,1.],calls
print("IC-07 PASS: five alternating substeps; four recalculations, one correct cache hit")
