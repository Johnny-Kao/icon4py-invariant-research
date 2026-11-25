"""Two-rank synthetic exchange through the actual MPI/GHEX implementation."""
import numpy as np
import gt4py.next as gtx
from mpi4py import MPI
from icon4py.model.common import dimension as dims
from icon4py.model.common.decomposition import definitions as d, mpi_decomposition as m

comm=MPI.COMM_WORLD
assert comm.size==2
rank=comm.rank
info=d.DecompositionInfo()
for i,dim in enumerate(dims.horizontal_dims(),start=1):
    owned=np.arange(rank*i,(rank+1)*i,dtype=np.int32)
    halo=np.arange((1-rank)*i,(2-rank)*i,dtype=np.int32)
    gids=np.concatenate([owned,halo])
    mask=np.arange(2*i)<i
    levels=np.where(mask,d.DecompositionFlag.OWNED.value,d.DecompositionFlag.FIRST_HALO_LEVEL.value).astype(np.int32)
    info.set_dimension(dim,gids,mask,levels)
props=m.MPICommProcessProperties(comm)
exchange=d.create_exchange(props,info)
for dim in dims.horizontal_dims():
    gids=info.global_index(dim)
    mask=info.owner_mask(dim)
    n=len(gids)
    # Include padding like Fortran nproma; only actual decomposition rows may change.
    array=np.full((n+2,3),-999.,dtype=np.float64)
    field=gtx.as_field((dim,dims.KDim),array)
    for step in range(3):
        expected=gids[:,None]*100+step*10+np.arange(3)
        field.ndarray[:n][mask]=expected[mask]
        field.ndarray[:n][~mask]=-123.
        exchange.exchange(dim,field,stream=d.BLOCK)
        np.testing.assert_array_equal(field.asnumpy()[:n],expected)
        np.testing.assert_array_equal(field.asnumpy()[n:],-999.)
    print(f"IC-05 PASS rank={rank} dim={dim.value} three exchanges including cached field and padding",flush=True)
comm.Barrier()
