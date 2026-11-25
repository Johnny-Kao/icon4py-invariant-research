"""IC-06/08: native transport arithmetic and bound-program isolation.

Synthetic fields only; not a complete atmospheric trajectory.
"""
import argparse
import numpy as np
import gt4py.next as gtx

from icon4py.model.atmosphere.tracer_advection.stencils.apply_density_increment import apply_density_increment
from icon4py.model.atmosphere.tracer_advection.stencils.integrate_tracer_horizontally import integrate_tracer_horizontally
from icon4py.model.common import dimension as dims, model_backends, model_options
from icon4py.model.common.grid.simple import simple_grid

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--backend", choices=("embedded", "gtfn_cpu"), default="embedded")
    args = parser.parse_args()
    backend = model_backends.BACKENDS[args.backend]
    grid = simple_grid(num_levels=2)
    n = grid.num_cells
    mass = np.ones((n, 2), dtype=np.float64)
    transfers = np.resize(np.array([0., .25, .9, .95, 1., 1.1]), n)
    flux = np.zeros((n, 3))
    flux[:, 1] = transfers
    inc = np.diff(flux, axis=1)
    def cell(x):
        return gtx.as_field((dims.CellDim, dims.KDim), x.copy())
    def half(x):
        return gtx.as_field((dims.CellDim, dims.KHalfDim), x.copy())
    def level(x):
        return gtx.as_field((dims.KDim,), np.asarray(x, dtype=np.float64))
    domain = dict(horizontal_start=gtx.int32(0), horizontal_end=gtx.int32(n),
                  vertical_start=gtx.int32(0), vertical_end=gtx.int32(2))
    metric = level([1., 1.])
    def bound(dt, even):
        return model_options.setup_program(
            program=apply_density_increment, backend=backend,
            constant_args=dict(p_dtime=np.float64(dt), even_timestep=even,
                               deepatmo_divzl=metric, deepatmo_divzu=metric),
            horizontal_sizes={k:v for k,v in domain.items() if k.startswith("horizontal")},
            vertical_sizes={k:v for k,v in domain.items() if k.startswith("vertical")},
        )
    odd = bound(1., False)
    even = bound(.5, True)
    other = bound(.25, False)
    def density(program):
        out = cell(np.zeros_like(mass))
        program(rhodz_in=cell(mass), p_mflx_contra_v=half(flux), rhodz_out=out)
        return out
    # A/B/A/C/A checks independent specializations on the same native program.
    for prog, expected in [(odd,mass-inc),(even,mass+.5*inc),(odd,mass-inc),
                           (other,mass-.25*inc),(odd,mass-inc)]:
        np.testing.assert_allclose(density(prog).asnumpy(), expected, atol=1e-14)
    print("IC-08 PASS: 5 alternating bound-program calls, different dt and parity")
    old = density(odd)
    ref = cell(np.maximum(.1*mass,mass-inc))
    # Zero horizontal tracer flux isolates the downstream denominator effect.
    # Independent synthetic fields do not claim a dynamically consistent air-mass trajectory.
    horizontal = model_options.setup_program(
        program=integrate_tracer_horizontally, backend=backend,
        constant_args=dict(deepatmo_divh=metric, p_dtime=np.float64(1.),
                           geofac_div=gtx.as_field((dims.CellDim,dims.C2EDim),np.ones((n,3)))),
        horizontal_sizes={k:v for k,v in domain.items() if k.startswith("horizontal")},
        vertical_sizes={k:v for k,v in domain.items() if k.startswith("vertical")},
        offset_provider=grid.connectivities,
    )
    def tracer(denom):
        out = cell(np.zeros_like(mass))
        with np.errstate(divide="ignore", invalid="ignore"):
            horizontal(p_mflx_tracer_h=gtx.as_field((dims.EdgeDim,dims.KDim),np.zeros((grid.num_edges,2))),
                       tracer_now=cell(.01*mass), rhodz_now=cell(mass),
                       rhodz_new=denom, tracer_new_hor=out)
        return out.asnumpy()
    a,b=tracer(old),tracer(ref)
    np.testing.assert_allclose(a[transfers<.9],b[transfers<.9])
    np.testing.assert_allclose(a[transfers==.95,0],.2)
    np.testing.assert_allclose(b[transfers==.95,0],.1)
    assert np.all(~np.isfinite(a[transfers==1.,0]))
    assert np.all(a[transfers>1.,0]<0)
    assert np.all(np.isfinite(b)) and np.all(b>0)
    for i in range(6):
        print(f"IC-06 {transfers[i]:.2f} density={old.asnumpy()[i,0]:.6g} tracer={a[i,0]:.6g} floored_ref={b[i,0]:.6g}")
    print("IC-06 PASS: reproduced downstream numerical discrepancy; full-model reachability UNKNOWN")

if __name__ == "__main__":
    main()
