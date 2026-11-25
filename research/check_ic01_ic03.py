"""Research-only checks; synthetic inputs, not operational forecast validation."""
import ast
from pathlib import Path
import numpy as np
import gt4py.next as gtx
from icon4py.model.common import dimension as dims
from icon4py.model.common.decomposition.definitions import DecompositionInfo, DecompositionFlag
from icon4py.model.common.grid import vertical as v, horizontal as h, grid_refinement as r

def ic01():
    path=Path("bindings/src/icon4py/bindings/grid_wrapper.py")
    tree=ast.parse(path.read_text())
    grid_init=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=="grid_init")
    calls=[n for n in ast.walk(grid_init) if isinstance(n,ast.Call) and isinstance(n.func,ast.Attribute) and n.func.attr=="VerticalGridConfig"]
    keys={k.arg for k in calls[0].keywords}
    assert "htop_moist_proc" not in keys
    z=np.linspace(30000.,0.,61)
    default=v.VerticalGridConfig(num_levels=60).htop_moist_proc
    a=int(v.VerticalGrid._determine_start_level_of_moist_physics(z,default))
    b=int(v.VerticalGrid._determine_start_level_of_moist_physics(z,15000.))
    print(f"IC01 default={default} configured=15000 indices={a},{b}; differing layers={abs(a-b)}; hypothetical configuration, not Fortran integration")
    assert a!=b

def ic02():
    dim=dims.CellDim
    def bounds(values, levels):
        values=np.asarray(values,dtype=np.int32)
        levels=np.asarray(levels,dtype=np.int32)
        d=DecompositionInfo()
        d.set_dimension(dim,np.arange(len(values),dtype=np.int32),levels==DecompositionFlag.OWNED.value,levels)
        s,e=r.compute_domain_bounds(dim,{dim:gtx.as_field((dim,),values)},d)
        k=h.domain(dim)(h.Zone.INTERIOR)
        return int(s[k]),int(e[k])
    O=DecompositionFlag.OWNED.value
    H=DecompositionFlag.FIRST_HALO_LEVEL.value
    control=bounds([1,2,3,4,5,0,0,0],[O]*7+[H])
    empty=bounds([1,2,3,4,0,0,0],[O]*6+[H])
    print(f"IC02 control={control} expected=(5,7); absent nudging synthetic rank={empty} expected by zone semantics=(4,6)")
    assert control==(5,7)
    print("IC02", "counterexample" if empty!=(4,6) else "no discrepancy", "(rank reachability unverified)")

def ic03():
    for levels in (60,90):
        cfg=v.VerticalGridConfig(num_levels=levels)
        a,b=v.get_vct_a_and_vct_b(cfg,None)
        z=a.asnumpy()
        grid=v.VerticalGrid(config=cfg,vct_a=a,vct_b=b)
        nflat=int(grid.nflatlev)
        # Analytic SLEVE columns, using a range of full and smoothed topography heights.
        topo=np.arange(0.,5501.,25.)
        full=np.repeat(topo, len(topo))
        smooth=np.tile(topo,len(topo))
        keep=np.abs(full-smooth)<=3000
        full,smooth=full[keep],smooth[keep]
        zz=np.broadcast_to(z,(len(full),len(z))).copy()
        def decay(scale):
            exp=cfg.SLEVE_decay_exponent
            top=cfg.model_top_height
            return np.sinh((top/scale)**exp-(z[nflat+1:-1]/scale)**exp)/np.sinh((top/scale)**exp)
        zz[:,nflat+1:-1]+=smooth[:,None]*decay(cfg.SLEVE_decay_scale_1)+(full-smooth)[:,None]*decay(cfg.SLEVE_decay_scale_2)
        zz[:,-1]=full
        dz=z[:-1]-z[1:]
        ratio=((500-dz)/400)**2
        floor=np.maximum(np.where(dz<100,dz/3,np.where(dz<500,(ratio/3+(1-ratio)/2)*dz,250*(dz/500)**(1/3))),min(50,cfg.lowest_layer_thickness))
        try:
            result=v._check_and_correct_layer_thickness(vertical_coordinate=zz,vct_a=z,SLEVE_minimum_layer_thickness_1=100.,SLEVE_minimum_relative_layer_thickness_1=1/3,SLEVE_minimum_layer_thickness_2=500.,SLEVE_minimum_relative_layer_thickness_2=0.5,lowest_layer_thickness=cfg.lowest_layer_thickness)
        except Exception as exc:
            print(f"IC03 levels={levels}: blocked by {type(exc).__name__}: {exc}")
            continue
        thickness=result[:,:-1]-result[:,1:]
        fail=thickness<floor-1e-8
        rows=np.flatnonzero(fail.any(axis=1))
        print(f"IC03 levels={levels} columns={len(full)} violating_columns={len(rows)} min_thickness={thickness.min():.9g} nonpositive={np.count_nonzero(thickness<=0)}")
        if len(rows):
            i=rows[0]; k=int(np.flatnonzero(fail[i])[0])
            print(f"IC03 example terrain={full[i]} smoothed={smooth[i]} layer={k} thickness={thickness[i,k]} floor={floor[k]}")

if __name__=="__main__":
    ic01()
    ic02()
    ic03()
