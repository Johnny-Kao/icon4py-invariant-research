"""Compare native transport density stencils with ICON's density floor.

Reference: ICON 7de853356c97db40d9f57105c58e72c3bf822cb0,
src/advection/mo_advection_stepping.f90, lines 379-384.
Synthetic input; atmospheric reachability is not established.
"""
import gt4py.next as gtx
import numpy as np

from icon4py.model.atmosphere.tracer_advection.stencils.apply_density_increment import (
    apply_density_increment,
)
from icon4py.model.atmosphere.tracer_advection.stencils.apply_horizontal_density_increment import (
    apply_horizontal_density_increment,
)
from icon4py.model.common import dimension as dims

mass = np.ones((6, 2), dtype=np.float64)
transfer = np.array([0., .25, .9, .95, 1., 1.1])
flux = np.zeros((6, 3), dtype=np.float64)
flux[:, 1] = transfer
increment = np.diff(flux, axis=1)
def cell(x):
    return gtx.as_field((dims.CellDim, dims.KDim), x.copy())
def half(x):
    return gtx.as_field((dims.CellDim, dims.KHalfDim), x.copy())
metric = gtx.as_field((dims.KDim,), np.ones(2, dtype=np.float64))
kwargs = dict(
    p_mflx_contra_v=half(flux),
    deepatmo_divzl=metric,
    deepatmo_divzu=metric,
    p_dtime=np.float64(1.),
    horizontal_start=gtx.int32(0),
    horizontal_end=gtx.int32(6),
    vertical_start=gtx.int32(0),
    vertical_end=gtx.int32(2),
    offset_provider={},
)
odd = cell(np.zeros_like(mass))
even = cell(np.zeros_like(mass))
horizontal = cell(np.zeros_like(mass))
apply_density_increment.with_backend(None)(rhodz_in=cell(mass), rhodz_out=odd, even_timestep=False, **kwargs)
apply_density_increment.with_backend(None)(rhodz_in=cell(mass), rhodz_out=even, even_timestep=True, **kwargs)
apply_horizontal_density_increment.with_backend(None)(p_rhodz_new=cell(mass), rhodz_ast2=horizontal, **kwargs)
ref = np.maximum(.1 * mass, mass - increment)
old = odd.asnumpy()
np.testing.assert_allclose(even.asnumpy(), mass + increment)
np.testing.assert_allclose(horizontal.asnumpy(), old)
np.testing.assert_allclose(old, mass - increment)
np.testing.assert_allclose(old[:3], ref[:3], atol=1e-15)
assert np.count_nonzero(old < .1 - 1e-15) == 3
assert np.count_nonzero(old <= 0) == 2
print("transfer_fraction, native_odd, ICON_reference")
for i,x in enumerate(transfer):
    print(f"{x:.2f}, {old[i,0]:.6f}, {ref[i,0]:.6f}")
print("PASS: even branch and normal flux controls; odd branch violates density floor in 3 cases")
print("Native horizontal density stencil reproduces the same discrepancy.")
print("No full transport run; frequency and reachability with real model inputs remain unverified.")
