"""Test-only CPU fixture for the genuine grid_init binding."""
from types import SimpleNamespace
import gt4py.next as gtx
import numpy as np
from icon4py.bindings import grid_wrapper, common, debug_utils
from icon4py.model.common import dimension as dims, model_backends
from icon4py.model.common.grid import vertical as v
from icon4py.model.common.decomposition import definitions
from icon4py.model.common.type_alias import vpfloat, wpfloat
from icon4py.model.atmosphere.dycore.stencils.update_dynamical_exner_time_increment import update_dynamical_exner_time_increment
from icon4py.tools import py2fgen

def array(dtype):
    return py2fgen.ArrayParamDescriptor(rank=1, dtype=dtype, memory_space=py2fgen.MemorySpace.HOST, is_optional=False)

@py2fgen.export(param_descriptors={"z": array(py2fgen.FLOAT64)})
def prepare_fixture(z):
    # Only the unrelated horizontal geometry and MPI machinery is replaced.
    # grid_init itself, vertical config, grid and stencil are unmodified.
    common.construct_icon_grid = lambda **kwargs: SimpleNamespace()
    common.construct_decomposition = lambda **kwargs: (
        definitions.SingleNodeProcessProperties(), None, definitions.SingleNodeExchange()
    )
    debug_utils.print_grid_decomp_info = lambda **kwargs: None
    common.select_backend = lambda *args, **kwargs: None
    model_backends.get_allocator = lambda *args, **kwargs: None
    cfg = v.VerticalGridConfig(num_levels=90, model_top_height=60000., lowest_layer_thickness=20., flat_height=20000., stretch_factor=.6, rayleigh_damping_height=35000.)
    a, _ = v.get_vct_a_and_vct_b(cfg, None)
    z[:] = a.asnumpy()

@py2fgen.export(param_descriptors={"index": array(py2fgen.INT32), "output": array(py2fgen.FLOAT64)})
def read_result(index, output):
    k = int(grid_wrapper.grid_state.vertical_grid.kstart_moist)
    index[0] = k
    initial = np.linspace(.25, .99, 90, dtype=wpfloat)[None,:]
    new = initial + wpfloat(.0002)
    phys = np.full((1,90), 1.e-7, dtype=vpfloat)
    accum = np.zeros((1,90), dtype=vpfloat)
    accum[:,k:] = initial[:,k:].astype(vpfloat)
    make = lambda a: gtx.as_field((dims.CellDim,dims.KDim), a)
    out = make(accum)
    update_dynamical_exner_time_increment.with_backend(None)(
        exner=make(new), ddt_exner_phy=make(phys), exner_dyn_incr=out,
        ndyn_substeps_var=wpfloat(5.), dtime=wpfloat(60.),
        horizontal_start=gtx.int32(0), horizontal_end=gtx.int32(1),
        vertical_start=gtx.int32(k), vertical_end=gtx.int32(90), offset_provider={}
    )
    output[:] = out.asnumpy()[0]
