! Copyright (C) 2004- DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
! SPDX-License-Identifier: BSD-3-Clause
! Project URL: https://icon-model.org

MODULE mo_icon4py_interfaces

#ifdef __USE_ICON4PY
  USE, INTRINSIC :: ISO_C_BINDING
  USE mo_kind,                 ONLY: wp, vp
  USE mo_nonhydrostatic_config,ONLY: itime_scheme, divdamp_order, ndyn_substeps, &
       &                             ndyn_substeps_var, divdamp_fac_o2, l_zdiffu_t, damp_height, &
       &                             veladv_offctr, rhotheta_offctr, &
       &                             divdamp_fac, divdamp_fac2, divdamp_fac3, divdamp_fac4, &
       &                             divdamp_z, divdamp_z2, divdamp_z3, divdamp_z4, &
       &                             divdamp_type, damp_height, rayleigh_type, &
       &                             igradp_method, iadv_rhotheta, lextra_diffu
  USE mo_parallel_config,      ONLY: nproma
  USE mo_initicon_config,      ONLY: init_mode, iau_wgt_dyn, is_iau_active
  USE mo_diffusion_config,     ONLY: diffusion_config
  USE mo_dynamics_config,      ONLY: ldeepatmo
  USE mo_grid_config,          ONLY: n_dom, l_limited_area
  USE mo_interpol_config,      ONLY: nudge_max_coeff
  USE mo_gridref_config,       ONLY: denom_diffu_v
  USE mo_turbdiff_config,      ONLY: turbdiff_config
  USE mo_run_config,           ONLY: icon4py_backend, lvert_nest, iforcing, &
       &                             icon4py_mode, ICON4PY_MODE_OFF,        &
       &                             ICON4PY_MODE_SUBSTITUTE, ICON4PY_MODE_VERIFY
  USE mo_sleve_config,         ONLY: stretch_fac, min_lay_thckn, top_height, flat_height
  USE mo_vertical_grid,        ONLY: nflat_gradp
  USE mo_vertical_coord_table, ONLY: vct_a, vct_b
  USE mo_model_domain,         ONLY: p_patch, t_patch
  USE mo_nonhydro_types,       ONLY: t_nh_state, t_nh_prog, t_nh_diag, t_nh_metrics
  USE mo_nonhydro_state,       ONLY: p_nh_state, p_nh_state_lists
  !$ser verbatim USE mo_nonhydro_state, ONLY: jstep_ptr, nstep_ptr
  USE mo_intp_data_strc,       ONLY: p_int_state, t_int_state
  USE mo_loopindices,          ONLY: get_indices_c, get_indices_e
  USE mo_impl_constants,       ONLY: min_rlcell, min_rledge, min_rlcell_int, MODE_IAU
  USE mo_impl_constants_grf,   ONLY: grf_bdywidth_c
  USE mo_fortran_tools,        ONLY: assert_acc_device_only
  USE mo_exception,            ONLY: message, message_text, finish
  USE mo_prepadv_types,        ONLY: t_prepare_adv
  USE mo_mpi,                  ONLY: get_my_mpi_work_id, get_my_mpi_work_communicator
  USE icon4py_bindings,        ONLY: solve_nh_run, diffusion_run, diffusion_init, grid_init, solve_nh_init
  USE mo_icon4py_helpers,      ONLY: reorder_geographical_coordinates, reorder_tangent_vectors, &
       &                             reorder_tangent_vectors_sparse

  !! ICON4Py verification
  !$ser verbatim USE mo_ser_all, ONLY: serialize_all
  !$ser verbatim USE mo_icon4py_verification, ONLY: serialize_diagnostics, serialize_prognostics, serialize_metrics, serialize_patch, serialize_interpolation_state

  USE mo_solve_nonhydro,       ONLY: solve_nh
  USE mo_nh_diffusion,         ONLY: diffusion
  USE mo_fortran_tools,        ONLY: copy
  USE mo_nh_verify_utils,      ONLY: verify_field_wp, print_verify_header, print_verify_footer

  !$ser verbatim USE mo_nonhydro_state,       ONLY: jstep_ptr, nstep_ptr, mo_solve_nonhydro_ctr, diffctr

  !! ICON4Py verification
  !$ser verbatim USE mo_icon4py_verification, ONLY: serialize_solve_nonhydro_init, serialize_solve_nonhydro_exit, serialize_solve_nonhydro_final, serialize_solve_nonhydro_30to38_exit, serialize_solve_nonhydro_30to38_init, serialize_solve_nonhydro_14to28_init, serialize_solve_nonhydro_14to28_exit, serialize_solve_nonhydro_41to60_init, serialize_solve_nonhydro_41to60_exit
  !$ser verbatim USE mo_icon4py_verification, ONLY:serialize_diffusion_init, serialize_diffusion_exit
  !$ser verbatim USE mtime,                   ONLY: datetime
  !$ser verbatim USE mo_time_config,          ONLY: time_config


  IMPLICIT NONE

  PRIVATE

  !> module name string
  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_icon4py_interfaces'

  ! Geographical coordinates and tangent vectors for Python granule processing
  REAL(wp), ALLOCATABLE :: primal_normal_x(:,:)
  REAL(wp), ALLOCATABLE :: primal_normal_y(:,:)
  REAL(wp), ALLOCATABLE :: primal_normal_vert_x(:,:,:)
  REAL(wp), ALLOCATABLE :: primal_normal_vert_y(:,:,:)
  REAL(wp), ALLOCATABLE :: dual_normal_vert_x(:,:,:)
  REAL(wp), ALLOCATABLE :: dual_normal_vert_y(:,:,:)
  REAL(wp), ALLOCATABLE :: primal_normal_cell_x(:,:,:)
  REAL(wp), ALLOCATABLE :: primal_normal_cell_y(:,:,:)
  REAL(wp), ALLOCATABLE :: dual_normal_cell_x(:,:,:)
  REAL(wp), ALLOCATABLE :: dual_normal_cell_y(:,:,:)
  REAL(wp), ALLOCATABLE :: cell_center_lat(:,:)
  REAL(wp), ALLOCATABLE :: cell_center_lon(:,:)
  REAL(wp), ALLOCATABLE :: edges_center_lat(:,:)
  REAL(wp), ALLOCATABLE :: edges_center_lon(:,:)

  ! Buffers passed to the granule via host_data use_device; module scope to
  ! survive past *_init, contiguous to bypass the binding's contiguous-dummy
  ! host-temp copy, logical(c_bool) for the py2fgen interface.
  LOGICAL(c_bool), ALLOCATABLE :: c_owner_mask_cb(:), e_owner_mask_cb(:), v_owner_mask_cb(:)
  LOGICAL(c_bool), ALLOCATABLE :: mask_prog_halo_c_cb(:)
  INTEGER,         ALLOCATABLE :: c2e_cc(:,:), e2c_cc(:,:), c2e2c_cc(:,:), e2c2e_cc(:,:), &
                                  e2v_cc(:,:), v2e_cc(:,:), v2c_cc(:,:), e2c2v_cc(:,:), c2v_cc(:,:)

  !! ICON4Py verification
  !$ser verbatim TYPE(datetime)::date

  PUBLIC :: icon4py_special_fields_create, icon4py_special_fields_delete, icon4py_init
  PUBLIC :: icon4py_solve_nh, icon4py_diffusion

  CONTAINS

  !-------------------------------------------------------------------------
  !>
  !! Organizes nonhydrostatic time stepping
  !! Currently we assume to have only one grid level.
  !!

  SUBROUTINE icon4py_special_fields_create

    ! Parallel processing variables
    INTEGER :: jg

    IF (icon4py_mode == ICON4PY_MODE_OFF) RETURN

    ! Reorder coordinates and tangent vectors
    jg = n_dom
    IF( n_dom /= 1 ) THEN
      CALL finish('perform_nh_stepping','Diffusion Python granule does not support nesting!')
    ENDIF
    CALL reorder_geographical_coordinates(p_patch(jg)%cells%center, cell_center_lat, cell_center_lon)
    CALL reorder_geographical_coordinates(p_patch(jg)%edges%center, edges_center_lat, edges_center_lon)
    CALL reorder_tangent_vectors(p_patch(jg)%edges%primal_normal, primal_normal_x, primal_normal_y)
    CALL reorder_tangent_vectors_sparse(p_patch(jg)%edges%primal_normal_vert, primal_normal_vert_x, primal_normal_vert_y)
    CALL reorder_tangent_vectors_sparse(p_patch(jg)%edges%dual_normal_vert, dual_normal_vert_x, dual_normal_vert_y)
    CALL reorder_tangent_vectors_sparse(p_patch(jg)%edges%primal_normal_cell, primal_normal_cell_x, primal_normal_cell_y)
    CALL reorder_tangent_vectors_sparse(p_patch(jg)%edges%dual_normal_cell, dual_normal_cell_x, dual_normal_cell_y)

  END SUBROUTINE icon4py_special_fields_create

  SUBROUTINE icon4py_init

  INTEGER :: jg, this_pe, p_comm
  ! Return code for Python functions called from ICON
  INTEGER(c_int) :: rc

  CHARACTER(len=*), PARAMETER :: routine = modname//':icon4py_init'

  IF (icon4py_mode == ICON4PY_MODE_OFF) RETURN

#ifdef _OPENACC
  ! GPU-enabled, but possible with python dycore granule
  WRITE(message_text,'(a)') 'This executable was built for:  BUILD_GPU2PY'
  CALL message(routine, message_text)
#else
  ! CPU-only, but possible with python dycore granule
  WRITE(message_text,'(a)') 'This executable was built for:  BUILD_CPU2PY'
  CALL message(routine, message_text)
#endif

  ! Copying data from host to device
  !$ACC ENTER DATA COPYIN(cell_center_lat, cell_center_lon, edges_center_lat, edges_center_lon, primal_normal_x, primal_normal_y) &
  !$ACC   COPYIN(primal_normal_vert_x, primal_normal_vert_y, dual_normal_vert_x, dual_normal_vert_y) &
  !$ACC   COPYIN(primal_normal_cell_x, primal_normal_cell_y, dual_normal_cell_x, dual_normal_cell_y, vct_b)

  ! get mpi communicator for granule
  this_pe = get_my_mpi_work_id()
  p_comm  = get_my_mpi_work_communicator()

  jg = n_dom
  IF( n_dom /= 1 ) THEN
    CALL finish(routine,'Diffusion Python granule does not support nesting!')
  ENDIF

  IF (p_patch(jg)%nblks_e > 1) THEN
    WRITE(message_text,'(a,i0,a,i0,a,i0,a)') &
      'ICON4Py requires a single block per patch: got nblks_e=', &
      p_patch(jg)%nblks_e, ' with nproma=', nproma, &
      ' and n_patch_edges=', p_patch(jg)%n_patch_edges, &
      '. Use nblocks_e=1 or increase nproma (>= n_patch_edges) in &parallel_nml.'
    CALL finish(routine, message_text)
  ENDIF

  WRITE(message_text,'(a,i6,a,i6)') 'Diffusion Granule MPI Communicator:', this_pe,' Communicator: ', p_comm
  CALL message(routine, message_text)

  WRITE(message_text,'(a,i6,a,i6)') 'diffusion initialized for ', n_dom,' domains, zd_diffcoef size ', SIZE(p_nh_state(1)%metrics%zd_diffcoef,1)
  CALL message(routine, message_text)

  CALL message('Icon4Py', "Grid initialization")

  c_owner_mask_cb = LOGICAL(p_patch(jg)%cells%decomp_info%owner_mask(:, 1), c_bool)
  e_owner_mask_cb = LOGICAL(p_patch(jg)%edges%decomp_info%owner_mask(:, 1), c_bool)
  v_owner_mask_cb = LOGICAL(p_patch(jg)%verts%decomp_info%owner_mask(:, 1), c_bool)

  c2e_cc   = p_patch(jg)%cells%edge_idx(:,1,:)
  e2c_cc   = p_patch(jg)%edges%cell_idx(:,1,:)
  c2e2c_cc = p_patch(jg)%cells%neighbor_idx(:,1,:)
  e2c2e_cc = p_patch(jg)%edges%quad_idx(:,1,:)
  e2v_cc   = p_patch(jg)%edges%vertex_idx(:,1,:)
  v2e_cc   = p_patch(jg)%verts%edge_idx(:,1,:)
  v2c_cc   = p_patch(jg)%verts%cell_idx(:,1,:)
  e2c2v_cc = p_patch(jg)%edges%vertex_idx(:,1,:)
  c2v_cc   = p_patch(jg)%cells%vertex_idx(:,1,:)

  !$ACC ENTER DATA COPYIN(c2e_cc, e2c_cc, c2e2c_cc, e2c2e_cc, e2v_cc, v2e_cc, v2c_cc, e2c2v_cc, c2v_cc)

  CALL grid_init( &
    cell_starts                 = p_patch(jg)%cells%start_index, &
    cell_ends                   = p_patch(jg)%cells%end_index, &
    vertex_starts               = p_patch(jg)%verts%start_index, &
    vertex_ends                 = p_patch(jg)%verts%end_index, &
    edge_starts                 = p_patch(jg)%edges%start_index, &
    edge_ends                   = p_patch(jg)%edges%end_index, &
    c2e                         = c2e_cc, &
    e2c                         = e2c_cc, &
    c2e2c                       = c2e2c_cc, &
    e2c2e                       = e2c2e_cc, &
    e2v                         = e2v_cc, &
    v2e                         = v2e_cc, &
    v2c                         = v2c_cc, &
    e2c2v                       = e2c2v_cc, &
    c2v                         = c2v_cc, &
    c_owner_mask                = c_owner_mask_cb, &
    e_owner_mask                = e_owner_mask_cb, &
    v_owner_mask                = v_owner_mask_cb, &
    c_glb_index                 = p_patch(jg)%cells%decomp_info%glb_index(:), &
    e_glb_index                 = p_patch(jg)%edges%decomp_info%glb_index(:), &
    v_glb_index                 = p_patch(jg)%verts%decomp_info%glb_index(:), &
    tangent_orientation         = p_patch(jg)%edges%tangent_orientation(:,1), &
    inverse_primal_edge_lengths = p_patch(jg)%edges%inv_primal_edge_length(:,1), &
    inv_dual_edge_length        = p_patch(jg)%edges%inv_dual_edge_length(:,1), &
    inv_vert_vert_length        = p_patch(jg)%edges%inv_vert_vert_length(:,1), &
    edge_areas                  = p_patch(jg)%edges%area_edge(:,1), &
    f_e                         = p_patch(jg)%edges%f_e(:,1), &
    cell_center_lat             = cell_center_lat(:, 1), &
    cell_center_lon             = cell_center_lon(:, 1), &
    cell_areas                  = p_patch(jg)%cells%area(:,1), &
    primal_normal_vert_x        = primal_normal_vert_x(:,:,1), &
    primal_normal_vert_y        = primal_normal_vert_y(:,:,1), &
    dual_normal_vert_x          = dual_normal_vert_x(:,:,1), &
    dual_normal_vert_y          = dual_normal_vert_y(:,:,1), &
    primal_normal_cell_x        = primal_normal_cell_x(:,:,1), &
    primal_normal_cell_y        = primal_normal_cell_y(:,:,1), &
    dual_normal_cell_x          = dual_normal_cell_x(:,:,1), &
    dual_normal_cell_y          = dual_normal_cell_y(:,:,1), &
    edge_center_lat             = edges_center_lat(:, 1), &
    edge_center_lon             = edges_center_lon(:, 1), &
    primal_normal_x             = primal_normal_x(:, 1), &
    primal_normal_y             = primal_normal_y(:, 1), &
    vct_a                       = vct_a, &
    lowest_layer_thickness      = min_lay_thckn, &
    model_top_height            = top_height, &
    stretch_factor              = stretch_fac, &
    flat_height                 = flat_height, &
    rayleigh_damping_height     = damp_height(jg), &
    mean_cell_area              = p_patch(jg)%geometry_info%mean_cell_area, &
    comm_id                     = p_comm, &
    num_vertices                = p_patch(jg)%n_patch_verts, &
    num_cells                   = p_patch(jg)%n_patch_cells, &
    num_edges                   = p_patch(jg)%n_patch_edges, &
    vertical_size               = p_patch(jg)%nlev, &
    limited_area                = LOGICAL(l_limited_area, c_bool), &
    backend                     = icon4py_backend, &
    rc                          = rc )

  CALL icon4py_check_rc(rc, routine, 'grid_init')

  CALL message('Icon4py', "Diffusion granule initialization")

  CALL diffusion_init( &
    theta_ref_mc                = p_nh_state(jg)%metrics%theta_ref_mc(:,:,1), &
    wgtfac_c                    = p_nh_state(jg)%metrics%wgtfac_c(:,:,1), &
    e_bln_c_s                   = p_int_state(jg)%e_bln_c_s(:,:,1), &
    geofac_div                  = p_int_state(jg)%geofac_div(:,:,1), &
    geofac_grg_x                = p_int_state(jg)%geofac_grg(:,:,1,1), &
    geofac_grg_y                = p_int_state(jg)%geofac_grg(:,:,1,2), &
    geofac_n2s                  = p_int_state(jg)%geofac_n2s(:,:,1), &
    nudgecoeff_e                = p_int_state(jg)%nudgecoeff_e(:,1), &
    rbf_vec_coeff_v             = p_int_state(jg)%rbf_vec_coeff_v(:,:,:,1), &
    zd_cellidx                  = p_nh_state(jg)%metrics%zd_indlist, &
    zd_vertidx                  = p_nh_state(jg)%metrics%zd_vertidx, &
    zd_intcoef                  = p_nh_state(jg)%metrics%zd_intcoef, &
    zd_diffcoef                 = p_nh_state(jg)%metrics%zd_diffcoef, &
    ndyn_substeps               = ndyn_substeps, &
    diffusion_type              = diffusion_config(jg)%hdiff_order, &
    hdiff_w                     = LOGICAL(diffusion_config(jg)%lhdiff_w, c_bool), &
    hdiff_vn                    = LOGICAL(diffusion_config(jg)%lhdiff_vn, c_bool), &
    hdiff_smag_w                = LOGICAL(diffusion_config(jg)%lhdiff_smag_w, c_bool), &
    zdiffu_t                    = LOGICAL(l_zdiffu_t, c_bool), &
    type_t_diffu                = diffusion_config(jg)%itype_t_diffu, &
    type_vn_diffu               = diffusion_config(jg)%itype_vn_diffu, &
    hdiff_efdt_ratio            = diffusion_config(jg)%hdiff_efdt_ratio, &
    hdiff_w_efdt_ratio          = diffusion_config(jg)%hdiff_w_efdt_ratio, &
    smagorinski_scaling_factor  = diffusion_config(jg)%hdiff_smag_fac, &
    smagorinski_scaling_factor2 = diffusion_config(jg)%hdiff_smag_fac2, &
    smagorinski_scaling_factor3 = diffusion_config(jg)%hdiff_smag_fac3, &
    smagorinski_scaling_factor4 = diffusion_config(jg)%hdiff_smag_fac4, &
    smagorinski_scaling_height  = diffusion_config(jg)%hdiff_smag_z, &
    smagorinski_scaling_height2 = diffusion_config(jg)%hdiff_smag_z2, &
    smagorinski_scaling_height3 = diffusion_config(jg)%hdiff_smag_z3, &
    smagorinski_scaling_height4 = diffusion_config(jg)%hdiff_smag_z4, &
    hdiff_temp                  = LOGICAL(diffusion_config(jg)%lhdiff_temp, c_bool), &
    denom_diffu_v               = denom_diffu_v, &
    nudge_max_coeff             = nudge_max_coeff, &
    itype_sher                  = turbdiff_config(jg)%itype_sher, &
    iforcing                    = iforcing,&
    a_hshr                      = turbdiff_config(jg)%a_hshr,&
    loutshs                     = LOGICAL(turbdiff_config(jg)%loutshs, c_bool), &
    backend                     = icon4py_backend, &
    rc                          = rc)

  CALL icon4py_check_rc(rc, routine, 'diffusion_init')

  WRITE(message_text,'(a,i6,a,i6)') 'Dycore Granule MPI Communicator:', this_pe,' Communicator: ', p_comm
  CALL message(routine, message_text)

  WRITE(message_text,'(a,i6,a,i6)') 'dycore initialized for ', n_dom,' domains, zd_diffcoef size ', SIZE(p_nh_state(1)%metrics%zd_diffcoef,1)
  CALL message(routine, message_text)

  ! Array arguments
  WRITE(message_text,'(a,i6,a,f10.3,a,f10.3)') 'vct_a shape: ', SIZE(vct_a, 1), ' min: ', MINVAL(vct_a), ' max: ', MAXVAL(vct_a)
  CALL message(routine, message_text)

  WRITE(message_text,'(a,i6,a,f10.3,a,f10.3)') 'vct_b shape: ', SIZE(vct_b, 1), ' min: ', MINVAL(vct_b), ' max: ', MAXVAL(vct_b)
  CALL message(routine, message_text)

  WRITE(message_text,'(a,i6,a,f10.3,a,f10.3)') 'cell_areas shape: ', SIZE(p_patch(jg)%cells%area, 1), ' min: ', MINVAL(p_patch(jg)%cells%area), ' max: ', MAXVAL(p_patch(jg)%cells%area)
  CALL message(routine, message_text)

  WRITE(message_text,'(a,i6,a,i6,a,f10.3,a,f10.3)') 'primal_normal_cell_x shape: ', SIZE(primal_normal_cell_x, 1), ' x ', SIZE(primal_normal_cell_x, 2), ' min: ', MINVAL(primal_normal_cell_x), ' max: ', MAXVAL(primal_normal_cell_x)
  CALL message(routine, message_text)

  WRITE(message_text,'(a,i6,a,i6,a,f10.3,a,f10.3)') 'primal_normal_cell_y shape: ', SIZE(primal_normal_cell_y, 1), ' x ', SIZE(primal_normal_cell_y, 2), ' min: ', MINVAL(primal_normal_cell_y), ' max: ', MAXVAL(primal_normal_cell_y)
  CALL message(routine, message_text)

  WRITE(message_text,'(a,i6,a,i6,a,f10.3,a,f10.3)') 'dual_normal_cell_x shape: ', SIZE(dual_normal_cell_x, 1), ' x ', SIZE(dual_normal_cell_x, 2), ' min: ', MINVAL(dual_normal_cell_x), ' max: ', MAXVAL(dual_normal_cell_x)
  CALL message(routine, message_text)

  WRITE(message_text,'(a,i6,a,i6,a,f10.3,a,f10.3)') 'dual_normal_cell_y shape: ', SIZE(dual_normal_cell_y, 1), ' x ', SIZE(dual_normal_cell_y, 2), ' min: ', MINVAL(dual_normal_cell_y), ' max: ', MAXVAL(dual_normal_cell_y)
  CALL message(routine, message_text)

  WRITE(message_text,'(a,i6,a,f10.3,a,f10.3)') 'edge_areas shape: ', SIZE(p_patch(jg)%edges%area_edge, 1), ' min: ', MINVAL(p_patch(jg)%edges%area_edge), ' max: ', MAXVAL(p_patch(jg)%edges%area_edge)
  CALL message(routine, message_text)

  WRITE(message_text,'(a,i6,a,f10.3,a,f10.3)') 'tangent_orientation shape: ', SIZE(p_patch(jg)%edges%tangent_orientation, 1), ' min: ', MINVAL(p_patch(jg)%edges%tangent_orientation), ' max: ', MAXVAL(p_patch(jg)%edges%tangent_orientation)
  CALL message(routine, message_text)

  WRITE(message_text,'(a,i6,a,f10.3,a,f10.3)') 'inverse_primal_edge_lengths shape: ', SIZE(p_patch(jg)%edges%inv_primal_edge_length, 1), ' min: ', MINVAL(p_patch(jg)%edges%inv_primal_edge_length), ' max: ', MAXVAL(p_patch(jg)%edges%inv_primal_edge_length)
  CALL message(routine, message_text)

  WRITE(message_text,'(a,i6,a,f10.3,a,f10.3)') 'inverse_dual_edge_lengths shape: ', SIZE(p_patch(jg)%edges%inv_dual_edge_length, 1), ' min: ', MINVAL(p_patch(jg)%edges%inv_dual_edge_length), ' max: ', MAXVAL(p_patch(jg)%edges%inv_dual_edge_length)
  CALL message(routine, message_text)

  WRITE(message_text,'(a,i6,a,f10.3,a,f10.3)') 'inverse_vertex_vertex_lengths shape: ', SIZE(p_patch(jg)%edges%inv_vert_vert_length, 1), ' min: ', MINVAL(p_patch(jg)%edges%inv_vert_vert_length), ' max: ', MAXVAL(p_patch(jg)%edges%inv_vert_vert_length)
  CALL message(routine, message_text)

  WRITE(message_text,'(a,i6,a,i6,a,f10.3,a,f10.3)') 'primal_normal_vert_x shape: ', SIZE(primal_normal_vert_x, 1), ' x ', SIZE(primal_normal_vert_x, 2), ' min: ', MINVAL(primal_normal_vert_x), ' max: ', MAXVAL(primal_normal_vert_x)
  CALL message(routine, message_text)

  WRITE(message_text,'(a,i6,a,i6,a,f10.3,a,f10.3)') 'primal_normal_vert_y shape: ', SIZE(primal_normal_vert_y, 1), ' x ', SIZE(primal_normal_vert_y, 2), ' min: ', MINVAL(primal_normal_vert_y), ' max: ', MAXVAL(primal_normal_vert_y)
  CALL message(routine, message_text)

  WRITE(message_text,'(a,i6,a,i6,a,f10.3,a,f10.3)') 'dual_normal_vert_x shape: ', SIZE(dual_normal_vert_x, 1), ' x ', SIZE(dual_normal_vert_x, 2), ' min: ', MINVAL(dual_normal_vert_x), ' max: ', MAXVAL(dual_normal_vert_x)
  CALL message(routine, message_text)

  WRITE(message_text,'(a,i6,a,i6,a,f10.3,a,f10.3)') 'dual_normal_vert_y shape: ', SIZE(dual_normal_vert_y, 1), ' x ', SIZE(dual_normal_vert_y, 2), ' min: ', MINVAL(dual_normal_vert_y), ' max: ', MAXVAL(dual_normal_vert_y)
  CALL message(routine, message_text)

  WRITE(message_text,'(a,i6,a,f10.3,a,f10.3)') 'f_e shape: ', SIZE(p_patch(jg)%edges%f_e, 1), ' min: ', MINVAL(p_patch(jg)%edges%f_e), ' max: ', MAXVAL(p_patch(jg)%edges%f_e)
  CALL message(routine, message_text)

  WRITE(message_text,'(a,i6,a,i6,a,f10.3,a,f10.3)') 'c_lin_e shape: ', SIZE(p_int_state(jg)%c_lin_e, 1), ' x ', SIZE(p_int_state(jg)%c_lin_e, 2), ' min: ', MINVAL(p_int_state(jg)%c_lin_e), ' max: ', MAXVAL(p_int_state(jg)%c_lin_e)
  CALL message(routine, message_text)

  WRITE(message_text,'(a,i6,a,i6,a,f10.3,a,f10.3)') 'c_intp shape: ', SIZE(p_int_state(jg)%cells_aw_verts, 1), ' x ', SIZE(p_int_state(jg)%cells_aw_verts, 2), ' min: ', MINVAL(p_int_state(jg)%cells_aw_verts), ' max: ', MAXVAL(p_int_state(jg)%cells_aw_verts)
  CALL message(routine, message_text)

  WRITE(message_text,'(a,i6,a,i6,a,f10.3,a,f10.3)') 'e_flx_avg shape: ', SIZE(p_int_state(jg)%e_flx_avg, 1), ' x ', SIZE(p_int_state(jg)%e_flx_avg, 2), ' min: ', MINVAL(p_int_state(jg)%e_flx_avg), ' max: ', MAXVAL(p_int_state(jg)%e_flx_avg)
  CALL message(routine, message_text)

  WRITE(message_text,'(a,i6,a,i6,a,f10.3,a,f10.3)') 'geofac_grdiv shape: ', SIZE(p_int_state(jg)%geofac_grdiv, 1), ' x ', SIZE(p_int_state(jg)%geofac_grdiv, 2), ' min: ', MINVAL(p_int_state(jg)%geofac_grdiv), ' max: ', MAXVAL(p_int_state(jg)%geofac_grdiv)
  CALL message(routine, message_text)

  WRITE(message_text,'(a,i6,a,i6,a,f10.3,a,f10.3)') 'geofac_rot shape: ', SIZE(p_int_state(jg)%geofac_rot, 1), ' x ', SIZE(p_int_state(jg)%geofac_rot, 2), ' min: ', MINVAL(p_int_state(jg)%geofac_rot), ' max: ', MAXVAL(p_int_state(jg)%geofac_rot)
  CALL message(routine, message_text)

  WRITE(message_text,'(a,i6,a,i6,a,i6,a,f10.3,a,f10.3)') 'pos_on_tplane_e_1 shape: ', SIZE(p_int_state(jg)%pos_on_tplane_e, 1), ' x ', SIZE(p_int_state(jg)%pos_on_tplane_e, 2), ' x ', SIZE(p_int_state(jg)%pos_on_tplane_e, 3), ' min: ', MINVAL(p_int_state(jg)%pos_on_tplane_e(:,:,1,1)), ' max: ', MAXVAL(p_int_state(jg)%pos_on_tplane_e(:,:,1,1))
  CALL message(routine, message_text)

  WRITE(message_text,'(a,i6,a,i6,a,i6,a,f10.3,a,f10.3)') 'pos_on_tplane_e_2 shape: ', SIZE(p_int_state(jg)%pos_on_tplane_e, 1), ' x ', SIZE(p_int_state(jg)%pos_on_tplane_e, 2), ' x ', SIZE(p_int_state(jg)%pos_on_tplane_e, 3), ' min: ', MINVAL(p_int_state(jg)%pos_on_tplane_e(:,:,2,1)), ' max: ', MAXVAL(p_int_state(jg)%pos_on_tplane_e(:,:,2,1))
  CALL message(routine, message_text)

  WRITE(message_text,'(a,i6,a,f10.3,a,f10.3)') 'cell_center_lat shape: ', SIZE(cell_center_lat, 1), ' min: ', MINVAL(cell_center_lat), ' max: ', MAXVAL(cell_center_lat)
  CALL message(routine, message_text)

  WRITE(message_text,'(a,i6,a,f10.3,a,f10.3)') 'cell_center_lon shape: ', SIZE(cell_center_lon, 1), ' min: ', MINVAL(cell_center_lon), ' max: ', MAXVAL(cell_center_lon)
  CALL message(routine, message_text)

  WRITE(message_text,'(a,i6,a,f10.3,a,f10.3)') 'edge_center_lat shape: ', SIZE(edges_center_lat, 1), ' min: ', MINVAL(edges_center_lat), ' max: ', MAXVAL(edges_center_lat)
  CALL message(routine, message_text)

  WRITE(message_text,'(a,i6,a,f10.3,a,f10.3)') 'edge_center_lon shape: ', SIZE(edges_center_lon, 1), ' min: ', MINVAL(edges_center_lon), ' max: ', MAXVAL(edges_center_lon)
  CALL message(routine, message_text)

  WRITE(message_text,'(a,i6,a,f10.3,a,f10.3)') 'primal_normal_x shape: ', SIZE(primal_normal_x, 1), ' min: ', MINVAL(primal_normal_x), ' max: ', MAXVAL(primal_normal_x)
  CALL message(routine, message_text)

  WRITE(message_text,'(a,i6,a,f10.3,a,f10.3)') 'primal_normal_y shape: ', SIZE(primal_normal_y, 1), ' min: ', MINVAL(primal_normal_y), ' max: ', MAXVAL(primal_normal_y)
  CALL message(routine, message_text)

  ! Scalar arguments
  WRITE(message_text,'(a,f10.3)') 'rayleigh_damping_height: ', damp_height(jg)
  CALL message(routine, message_text)

  WRITE(message_text,'(a,i6)') 'itime_scheme: ', itime_scheme
  CALL message(routine, message_text)

  WRITE(message_text,'(a,i6)') 'iadv_rhotheta: ', iadv_rhotheta
  CALL message(routine, message_text)

  WRITE(message_text,'(a,i6)') 'igradp_method: ', igradp_method
  CALL message(routine, message_text)

  WRITE(message_text,'(a,f10.3)') 'ndyn_substeps: ', REAL(ndyn_substeps, kind=wp)
  CALL message(routine, message_text)

  WRITE(message_text,'(a,i6)') 'rayleigh_type: ', rayleigh_type
  CALL message(routine, message_text)

  WRITE(message_text,'(a,i6)') 'divdamp_order: ', divdamp_order
  CALL message(routine, message_text)

  WRITE(message_text,'(a,l1)') 'is_iau_active: ', is_iau_active
  CALL message(routine, message_text)

  WRITE(message_text,'(a,f10.3)') 'iau_wgt_dyn: ', iau_wgt_dyn
  CALL message(routine, message_text)

  WRITE(message_text,'(a,i6)') 'divdamp_type: ', divdamp_type
  CALL message(routine, message_text)

  WRITE(message_text,'(a,l1)') 'l_vert_nested: ', lvert_nest
  CALL message(routine, message_text)

  WRITE(message_text,'(a,f10.3)') 'rhotheta_offctr: ', rhotheta_offctr
  CALL message(routine, message_text)

  WRITE(message_text,'(a,f10.3)') 'veladv_offctr: ', veladv_offctr
  CALL message(routine, message_text)

  WRITE(message_text,'(a,f10.3)') 'nudge_max_coeff: ', nudge_max_coeff
  CALL message(routine, message_text)

  WRITE(message_text,'(a,f10.3)') 'divdamp_fac: ', divdamp_fac
  CALL message(routine, message_text)

  WRITE(message_text,'(a,f10.3)') 'divdamp_fac2: ', divdamp_fac2
  CALL message(routine, message_text)

  WRITE(message_text,'(a,f10.3)') 'divdamp_fac3: ', divdamp_fac3
  CALL message(routine, message_text)

  WRITE(message_text,'(a,f10.3)') 'divdamp_fac4: ', divdamp_fac4
  CALL message(routine, message_text)

  WRITE(message_text,'(a,f10.3)') 'divdamp_z: ', divdamp_z
  CALL message(routine, message_text)

  WRITE(message_text,'(a,f10.3)') 'divdamp_z2: ', divdamp_z2
  CALL message(routine, message_text)

  WRITE(message_text,'(a,f10.3)') 'divdamp_z3: ', divdamp_z3
  CALL message(routine, message_text)

  WRITE(message_text,'(a,f10.3)') 'divdamp_z4: ', divdamp_z4
  CALL message(routine, message_text)

  WRITE(message_text,'(a,f10.3)') 'lowest_layer_thickness: ', min_lay_thckn
  CALL message(routine, message_text)

  WRITE(message_text,'(a,f10.3)') 'model_top_height: ', top_height
  CALL message(routine, message_text)

  WRITE(message_text,'(a,f10.3)') 'stretch_factor: ', stretch_fac
  CALL message(routine, message_text)

  WRITE(message_text,'(a,i6)') 'nflat_gradp: ', nflat_gradp(jg)
  CALL message(routine, message_text)

  WRITE(message_text,'(a,i6)') 'num_levels: ', p_patch(jg)%nlev
  CALL message(routine, message_text)

  mask_prog_halo_c_cb = LOGICAL(p_nh_state(jg)%metrics%mask_prog_halo_c(:, 1), c_bool)
  !$ACC ENTER DATA COPYIN(c_owner_mask_cb, mask_prog_halo_c_cb)

  CALL solve_nh_init( c_lin_e                       = p_int_state(jg)%c_lin_e(:,:,1), &
                      c_intp                        = p_int_state(jg)%cells_aw_verts(:,:,1), &
                      e_flx_avg                     = p_int_state(jg)%e_flx_avg(:,:,1), &
                      geofac_grdiv                  = p_int_state(jg)%geofac_grdiv(:,:,1), &
                      geofac_rot                    = p_int_state(jg)%geofac_rot(:,:,1), &
                      pos_on_tplane_e_1             = p_int_state(jg)%pos_on_tplane_e(:,:,1,1), &
                      pos_on_tplane_e_2             = p_int_state(jg)%pos_on_tplane_e(:,:,2,1), &
                      rbf_vec_coeff_e               = p_int_state(jg)%rbf_vec_coeff_e(:,:,1), &
                      e_bln_c_s                     = p_int_state(jg)%e_bln_c_s(:,:,1), &
                      rbf_vec_coeff_v               = p_int_state(jg)%rbf_vec_coeff_v(:,:,:,1), &
                      geofac_div                    = p_int_state(jg)%geofac_div(:,:,1), &
                      geofac_n2s                    = p_int_state(jg)%geofac_n2s(:,:,1), &
                      geofac_grg_x                  = p_int_state(jg)%geofac_grg(:,:,1,1), &
                      geofac_grg_y                  = p_int_state(jg)%geofac_grg(:,:,1,2), &
                      nudgecoeff_e                  = p_int_state(jg)%nudgecoeff_e(:,1), &
                      mask_prog_halo_c              = mask_prog_halo_c_cb, &
                      rayleigh_w                    = p_nh_state(jg)%metrics%rayleigh_w, &
                      exner_exfac                   = p_nh_state(jg)%metrics%exner_exfac(:,:,1), &
                      exner_ref_mc                  = p_nh_state(jg)%metrics%exner_ref_mc(:,:,1), &
                      wgtfac_c                      = p_nh_state(jg)%metrics%wgtfac_c(:,:,1), &
                      wgtfacq_c                     = p_nh_state(jg)%metrics%wgtfacq_c(:,:,1), &
                      inv_ddqz_z_full               = p_nh_state(jg)%metrics%inv_ddqz_z_full(:,:,1), &
                      rho_ref_mc                    = p_nh_state(jg)%metrics%rho_ref_mc(:,:,1), &
                      theta_ref_mc                  = p_nh_state(jg)%metrics%theta_ref_mc(:,:,1), &
                      vwind_expl_wgt                = p_nh_state(jg)%metrics%vwind_expl_wgt(:,1), &
                      d_exner_dz_ref_ic             = p_nh_state(jg)%metrics%d_exner_dz_ref_ic(:,:,1), &
                      ddqz_z_half                   = p_nh_state(jg)%metrics%ddqz_z_half(:,:,1), &
                      theta_ref_ic                  = p_nh_state(jg)%metrics%theta_ref_ic(:,:,1), &
                      d2dexdz2_fac1_mc              = p_nh_state(jg)%metrics%d2dexdz2_fac1_mc(:,:,1), &
                      d2dexdz2_fac2_mc              = p_nh_state(jg)%metrics%d2dexdz2_fac2_mc(:,:,1), &
                      rho_ref_me                    = p_nh_state(jg)%metrics%rho_ref_me(:,:,1), &
                      theta_ref_me                  = p_nh_state(jg)%metrics%theta_ref_me(:,:,1), &
                      ddxn_z_full                   = p_nh_state(jg)%metrics%ddxn_z_full(:,:,1), &
                      zdiff_gradp                   = p_nh_state(jg)%metrics%zdiff_gradp(:,:,:,1), &
                      vertidx_gradp                 = p_nh_state(jg)%metrics%vertidx_gradp(:,:,:,1), &
                      pg_edgeidx                    = p_nh_state(jg)%metrics%pg_edgeidx, &
                      pg_vertidx                    = p_nh_state(jg)%metrics%pg_vertidx, &
                      pg_exdist                     = p_nh_state(jg)%metrics%pg_exdist, &
                      ddqz_z_full_e                 = p_nh_state(jg)%metrics%ddqz_z_full_e(:,:,1), &
                      ddxt_z_full                   = p_nh_state(jg)%metrics%ddxt_z_full(:,:,1), &
                      wgtfac_e                      = p_nh_state(jg)%metrics%wgtfac_e(:,:,1), &
                      wgtfacq_e                     = p_nh_state(jg)%metrics%wgtfacq_e(:,:,1), &
                      vwind_impl_wgt                = p_nh_state(jg)%metrics%vwind_impl_wgt(:,1), &
                      hmask_dd3d                    = p_nh_state(jg)%metrics%hmask_dd3d(:,1), &
                      scalfac_dd3d                  = p_nh_state(jg)%metrics%scalfac_dd3d(:), &
                      coeff1_dwdz                   = p_nh_state(jg)%metrics%coeff1_dwdz(:,:,1), &
                      coeff2_dwdz                   = p_nh_state(jg)%metrics%coeff2_dwdz(:,:,1), &
                      coeff_gradekin                = p_nh_state(jg)%metrics%coeff_gradekin(:,:,1), &
                      c_owner_mask                  = c_owner_mask_cb, &
                      itime_scheme                  = itime_scheme, &
                      iadv_rhotheta                 = iadv_rhotheta, &
                      igradp_method                 = igradp_method, &
                      rayleigh_type                 = rayleigh_type, &
                      divdamp_order                 = divdamp_order, &
                      divdamp_type                  = divdamp_type, &
                      l_vert_nested                 = LOGICAL(lvert_nest, c_bool), &
                      ldeepatmo                     = LOGICAL(ldeepatmo, c_bool), &
                      iau_init                      = LOGICAL(init_mode == MODE_IAU, c_bool), &
                      extra_diffu                   = LOGICAL(lextra_diffu, c_bool), &
                      rhotheta_offctr               = rhotheta_offctr, &
                      veladv_offctr                 = veladv_offctr, &
                      nudge_max_coeff               = nudge_max_coeff, &
                      divdamp_fac                   = divdamp_fac, &
                      divdamp_fac2                  = divdamp_fac2, &
                      divdamp_fac3                  = divdamp_fac3, &
                      divdamp_fac4                  = divdamp_fac4, &
                      divdamp_z                     = divdamp_z, &
                      divdamp_z2                    = divdamp_z2, &
                      divdamp_z3                    = divdamp_z3, &
                      divdamp_z4                    = divdamp_z4, &
                      nflat_gradp                   = nflat_gradp(jg), &
                      backend                       = icon4py_backend, &
                      rc                            = rc )

  CALL icon4py_check_rc(rc, routine, 'solve_nh_init')

  END SUBROUTINE icon4py_init

  SUBROUTINE icon4py_special_fields_delete

  IF (icon4py_mode == ICON4PY_MODE_OFF) RETURN

  !$ACC EXIT DATA DELETE(cell_center_lat, cell_center_lon, edges_center_lat, edges_center_lon, primal_normal_x, primal_normal_y) &
  !$ACC   DELETE(primal_normal_vert_x, primal_normal_vert_y, dual_normal_vert_x, dual_normal_vert_y) &
  !$ACC   DELETE(primal_normal_cell_x, primal_normal_cell_y, dual_normal_cell_x, dual_normal_cell_y, vct_b) &
  !$ACC   DELETE(c2e_cc, e2c_cc, c2e2c_cc, e2c2e_cc, e2v_cc, v2e_cc, v2c_cc, e2c2v_cc, c2v_cc) &
  !$ACC   DELETE(c_owner_mask_cb, mask_prog_halo_c_cb)

  IF (ALLOCATED(cell_center_lat)) DEALLOCATE(cell_center_lat)
  IF (ALLOCATED(cell_center_lon)) DEALLOCATE(cell_center_lon)
  IF (ALLOCATED(edges_center_lat)) DEALLOCATE(edges_center_lat)
  IF (ALLOCATED(edges_center_lon)) DEALLOCATE(edges_center_lon)
  IF (ALLOCATED(primal_normal_x)) DEALLOCATE(primal_normal_x)
  IF (ALLOCATED(primal_normal_y)) DEALLOCATE(primal_normal_y)
  IF (ALLOCATED(primal_normal_vert_x)) DEALLOCATE(primal_normal_vert_x)
  IF (ALLOCATED(primal_normal_vert_y)) DEALLOCATE(primal_normal_vert_y)
  IF (ALLOCATED(dual_normal_vert_x)) DEALLOCATE(dual_normal_vert_x)
  IF (ALLOCATED(dual_normal_vert_y)) DEALLOCATE(dual_normal_vert_y)
  IF (ALLOCATED(primal_normal_cell_x)) DEALLOCATE(primal_normal_cell_x)
  IF (ALLOCATED(primal_normal_cell_y)) DEALLOCATE(primal_normal_cell_y)
  IF (ALLOCATED(dual_normal_cell_x)) DEALLOCATE(dual_normal_cell_x)
  IF (ALLOCATED(dual_normal_cell_y)) DEALLOCATE(dual_normal_cell_y)

  IF (ALLOCATED(c2e_cc))             DEALLOCATE(c2e_cc)
  IF (ALLOCATED(e2c_cc))             DEALLOCATE(e2c_cc)
  IF (ALLOCATED(c2e2c_cc))           DEALLOCATE(c2e2c_cc)
  IF (ALLOCATED(e2c2e_cc))           DEALLOCATE(e2c2e_cc)
  IF (ALLOCATED(e2v_cc))             DEALLOCATE(e2v_cc)
  IF (ALLOCATED(v2e_cc))             DEALLOCATE(v2e_cc)
  IF (ALLOCATED(v2c_cc))             DEALLOCATE(v2c_cc)
  IF (ALLOCATED(e2c2v_cc))           DEALLOCATE(e2c2v_cc)
  IF (ALLOCATED(c2v_cc))             DEALLOCATE(c2v_cc)
  IF (ALLOCATED(c_owner_mask_cb))    DEALLOCATE(c_owner_mask_cb)
  IF (ALLOCATED(e_owner_mask_cb))    DEALLOCATE(e_owner_mask_cb)
  IF (ALLOCATED(v_owner_mask_cb))    DEALLOCATE(v_owner_mask_cb)
  IF (ALLOCATED(mask_prog_halo_c_cb)) DEALLOCATE(mask_prog_halo_c_cb)

  END SUBROUTINE icon4py_special_fields_delete
  !-------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  !>
  !! Dispatcher selecting the implementation according to the icon4py_mode
  !! runtime namelist option.
  !!
  SUBROUTINE icon4py_solve_nh (p_nh, p_patch, p_int, prep_adv, nnow, nnew, l_init, l_recompute, lsave_mflx, &
                               lprep_adv, lclean_mflx, idyn_timestep, jstep, dtime, lis_iau_active, iau_wgt_dyn, lacc)

    TYPE(t_nh_state),    TARGET, INTENT(INOUT) :: p_nh
    TYPE(t_int_state),   TARGET, INTENT(IN)    :: p_int
    TYPE(t_patch),       TARGET, INTENT(INOUT) :: p_patch
    TYPE(t_prepare_adv), TARGET, INTENT(INOUT) :: prep_adv
    LOGICAL,                   INTENT(IN)    :: l_init, l_recompute, lsave_mflx, lprep_adv, lclean_mflx
    INTEGER,                   INTENT(IN)    :: idyn_timestep, jstep, nnow, nnew
    REAL(wp),                  INTENT(IN)    :: dtime
    LOGICAL,                   INTENT(IN)    :: lis_iau_active
    REAL(wp),                  INTENT(IN)    :: iau_wgt_dyn
    LOGICAL, INTENT(IN), OPTIONAL            :: lacc

    SELECT CASE (icon4py_mode)
    CASE (ICON4PY_MODE_OFF)
      CALL solve_nh(p_nh, p_patch, p_int, prep_adv, nnow, nnew, l_init, l_recompute, &
        &           lsave_mflx, lprep_adv, lclean_mflx, idyn_timestep, jstep, dtime, lacc=lacc)
    CASE (ICON4PY_MODE_SUBSTITUTE)
      CALL icon4py_solve_nh_substitute(p_nh, p_patch, p_int, prep_adv, nnow, nnew, l_init, l_recompute, lsave_mflx, &
                                       lprep_adv, lclean_mflx, idyn_timestep, jstep, dtime, lis_iau_active, iau_wgt_dyn, lacc)
    CASE (ICON4PY_MODE_VERIFY)
      CALL icon4py_solve_nh_verify(p_nh, p_patch, p_int, prep_adv, nnow, nnew, l_init, l_recompute, lsave_mflx, &
                                   lprep_adv, lclean_mflx, idyn_timestep, jstep, dtime, lis_iau_active, iau_wgt_dyn, lacc)
    END SELECT

  END SUBROUTINE icon4py_solve_nh


  SUBROUTINE icon4py_solve_nh_substitute (p_nh, p_patch, p_int, prep_adv, nnow, nnew, l_init, l_recompute, lsave_mflx, &
                               lprep_adv, lclean_mflx, idyn_timestep, jstep, dtime, lis_iau_active, iau_wgt_dyn, lacc)

    TYPE(t_nh_state),    TARGET, INTENT(INOUT) :: p_nh
    TYPE(t_int_state),   TARGET, INTENT(IN)    :: p_int
    TYPE(t_patch),       TARGET, INTENT(INOUT) :: p_patch
    TYPE(t_prepare_adv), TARGET, INTENT(INOUT) :: prep_adv
    REAL(vp), DIMENSION(:,:), POINTER          :: vn_incr_ptr => NULL()
    REAL(vp), DIMENSION(:,:), POINTER          :: rho_incr_ptr => NULL()
    REAL(vp), DIMENSION(:,:), POINTER          :: exner_incr_ptr => NULL()
    REAL(vp), DIMENSION(1)                     :: max_vcfl_size1_array

    ! Initialization switch that has to be .TRUE. at the initial time step only (not for restart)
    LOGICAL,                   INTENT(IN)    :: l_init
    ! Switch to recompute velocity tendencies after a physics call irrespective of the time scheme option
    LOGICAL,                   INTENT(IN)    :: l_recompute
    ! Switch if mass flux needs to be saved for nest boundary interpolation tendency computation
    LOGICAL,                   INTENT(IN)    :: lsave_mflx
    ! Switch if preparations for tracer advection shall be computed
    LOGICAL,                   INTENT(IN)    :: lprep_adv
    ! Switch if mass fluxes computed for tracer advection need to be reinitialized
    LOGICAL,                   INTENT(IN)    :: lclean_mflx
    ! Counter of dynamics time step within a large time step (ranges from 1 to ndyn_substeps)
    INTEGER,                   INTENT(IN)    :: idyn_timestep
    ! Time step count since last boundary interpolation (ranges from 0 to 2*ndyn_substeps-1)
    INTEGER,                   INTENT(IN)    :: jstep
    ! Time levels
    INTEGER,                   INTENT(IN)    :: nnow, nnew
    ! Dynamics time step
    REAL(wp),                  INTENT(IN)    :: dtime
    ! Switch if incremental analysis update is ative or not
    LOGICAL,                   INTENT(IN)    :: lis_iau_active
    ! Incremental analysis update weight for dynamics
    REAL(wp),                  INTENT(IN)    :: iau_wgt_dyn

    LOGICAL, INTENT(IN), OPTIONAL            :: lacc ! If true, use openacc

    INTEGER(c_int)           :: rc

      CALL assert_acc_device_only("icon4py_solve_nh_substitute", lacc)

      WRITE(message_text,'(a)') 'dycore_run ...'
      CALL message('icon4py_solve_nh_substitute', message_text)

      IF (ASSOCIATED(p_nh%diag%vn_incr)) vn_incr_ptr => p_nh%diag%vn_incr(:,:,1)
      IF (ASSOCIATED(p_nh%diag%rho_incr)) rho_incr_ptr => p_nh%diag%rho_incr(:,:,1)
      IF (ASSOCIATED(p_nh%diag%exner_incr)) exner_incr_ptr => p_nh%diag%exner_incr(:,:,1)
      !$ACC WAIT
      max_vcfl_size1_array(1) = p_nh%diag%max_vcfl_dyn ! max_vcfl is passed as a size 1 (host) array
      CALL solve_nh_run( rho_now             = p_nh%prog(nnow)%rho(:,:,1), &
                         rho_new             = p_nh%prog(nnew)%rho(:,:,1), &
                         exner_now           = p_nh%prog(nnow)%exner(:,:,1), &
                         exner_new           = p_nh%prog(nnew)%exner(:,:,1), &
                         w_now               = p_nh%prog(nnow)%w(:,:,1), &
                         w_new               = p_nh%prog(nnew)%w(:,:,1), &
                         theta_v_now         = p_nh%prog(nnow)%theta_v(:,:,1), &
                         theta_v_new         = p_nh%prog(nnew)%theta_v(:,:,1), &
                         vn_now              = p_nh%prog(nnow)%vn(:,:,1), &
                         vn_new              = p_nh%prog(nnew)%vn(:,:,1), &
                         w_concorr_c         = p_nh%diag%w_concorr_c(:,:,1), &
                         ddt_vn_apc_ntl1     = p_nh%diag%ddt_vn_apc_pc(:,:,1,nnow), &
                         ddt_vn_apc_ntl2     = p_nh%diag%ddt_vn_apc_pc(:,:,1,nnew), &
                         ddt_w_adv_ntl1      = p_nh%diag%ddt_w_adv_pc(:,:,1,nnow), &
                         ddt_w_adv_ntl2      = p_nh%diag%ddt_w_adv_pc(:,:,1,nnew), &
                         theta_v_ic          = p_nh%diag%theta_v_ic(:,:,1), &
                         rho_ic              = p_nh%diag%rho_ic(:,:,1), &
                         exner_pr            = p_nh%diag%exner_pr(:,:,1), &
                         exner_dyn_incr      = p_nh%diag%exner_dyn_incr(:,:,1), &
                         ddt_exner_phy       = p_nh%diag%ddt_exner_phy(:,:,1), &
                         grf_tend_rho        = p_nh%diag%grf_tend_rho(:,:,1), &
                         grf_tend_thv        = p_nh%diag%grf_tend_thv(:,:,1), &
                         grf_tend_w          = p_nh%diag%grf_tend_w(:,:,1), &
                         mass_fl_e           = p_nh%diag%mass_fl_e(:,:,1), &
                         ddt_vn_phy          = p_nh%diag%ddt_vn_phy(:,:,1), &
                         grf_tend_vn         = p_nh%diag%grf_tend_vn(:,:,1), &
                         vn_ie               = p_nh%diag%vn_ie(:,:,1), &
                         vt                  = p_nh%diag%vt(:,:,1), &
                         vn_incr             = vn_incr_ptr, &
                         rho_incr            = rho_incr_ptr, &
                         exner_incr          = exner_incr_ptr, &
                         mass_flx_me         = prep_adv%mass_flx_me(:,:,1), &
                         mass_flx_ic         = prep_adv%mass_flx_ic(:,:,1), &
                         vol_flx_ic          = prep_adv%vol_flx_ic(:,:,1), &
                         vn_traj             = prep_adv%vn_traj(:,:,1), &
                         dtime               = dtime, &
                         max_vcfl_size1_array= max_vcfl_size1_array, &
                         lprep_adv           = LOGICAL(lprep_adv, c_bool), &
                         at_initial_timestep = LOGICAL(l_init, c_bool), &
                         divdamp_fac_o2      = REAL(divdamp_fac_o2, c_double), &
                         ndyn_substeps_var   = ndyn_substeps_var(1), &
                         idyn_timestep       = idyn_timestep, &
                         is_iau_active       = LOGICAL(lis_iau_active, c_bool), &
                         iau_wgt_dyn         = iau_wgt_dyn, &
                         rc                  = rc )
      !$ACC WAIT
      CALL icon4py_check_rc(rc, 'icon4py_solve_nh_substitute', 'solve_nh_run')
      p_nh%diag%max_vcfl_dyn = max_vcfl_size1_array(1) ! retrieve max_vcfl from the array

  END SUBROUTINE icon4py_solve_nh_substitute


  SUBROUTINE icon4py_solve_nh_verify (p_nh, p_patch, p_int, prep_adv, nnow, nnew, l_init, l_recompute, lsave_mflx, &
                       lprep_adv, lclean_mflx, idyn_timestep, jstep, dtime, lis_iau_active, iau_wgt_dyn, lacc)

    TYPE(t_nh_state),    TARGET, INTENT(INOUT) :: p_nh
    TYPE(t_int_state),   TARGET, INTENT(IN)    :: p_int
    TYPE(t_patch),       TARGET, INTENT(INOUT) :: p_patch
    TYPE(t_prepare_adv), TARGET, INTENT(INOUT) :: prep_adv
    REAL(vp), DIMENSION(:,:), POINTER          :: vn_incr_ptr => NULL()
    REAL(vp), DIMENSION(:,:), POINTER          :: rho_incr_ptr => NULL()
    REAL(vp), DIMENSION(:,:), POINTER          :: exner_incr_ptr => NULL()
    REAL(vp), DIMENSION(1)                     :: max_vcfl_size1_array

    LOGICAL,                   INTENT(IN)    :: l_init
    LOGICAL,                   INTENT(IN)    :: l_recompute
    LOGICAL,                   INTENT(IN)    :: lsave_mflx
    LOGICAL,                   INTENT(IN)    :: lprep_adv
    LOGICAL,                   INTENT(IN)    :: lclean_mflx
    INTEGER,                   INTENT(IN)    :: idyn_timestep
    INTEGER,                   INTENT(IN)    :: jstep
    INTEGER,                   INTENT(IN)    :: nnow, nnew
    REAL(wp),                  INTENT(IN)    :: dtime
    LOGICAL,                   INTENT(IN)    :: lis_iau_active
    REAL(wp),                  INTENT(IN)    :: iau_wgt_dyn

    LOGICAL, INTENT(IN), OPTIONAL            :: lacc

    INTEGER(c_int)           :: rc

    REAL(wp), DIMENSION(nproma, p_patch%nlev) :: rho_now_before
    REAL(wp), DIMENSION(nproma, p_patch%nlev) :: rho_new_before
    REAL(wp), DIMENSION(nproma, p_patch%nlev) :: exner_now_before
    REAL(wp), DIMENSION(nproma, p_patch%nlev) :: exner_new_before
    REAL(wp), DIMENSION(nproma, p_patch%nlevp1) :: w_now_before
    REAL(wp), DIMENSION(nproma, p_patch%nlevp1) :: w_new_before
    REAL(wp), DIMENSION(nproma, p_patch%nlev) :: theta_v_now_before
    REAL(wp), DIMENSION(nproma, p_patch%nlev) :: theta_v_new_before
    REAL(wp), DIMENSION(nproma, p_patch%nlev) :: vn_now_before
    REAL(wp), DIMENSION(nproma, p_patch%nlev) :: vn_new_before
    REAL(wp), DIMENSION(nproma, p_patch%nlevp1) :: w_concorr_c_before
    REAL(wp), DIMENSION(nproma, p_patch%nlev) :: ddt_vn_apc_nnow_before
    REAL(wp), DIMENSION(nproma, p_patch%nlev) :: ddt_vn_apc_nnew_before
    REAL(wp), DIMENSION(nproma, p_patch%nlevp1) :: ddt_w_adv_nnow_before
    REAL(wp), DIMENSION(nproma, p_patch%nlevp1) :: ddt_w_adv_nnew_before
    REAL(wp), DIMENSION(nproma, p_patch%nlevp1) :: theta_v_ic_before
    REAL(wp), DIMENSION(nproma, p_patch%nlevp1) :: rho_ic_before
    REAL(wp), DIMENSION(nproma, p_patch%nlev) :: exner_pr_before
    REAL(wp), DIMENSION(nproma, p_patch%nlev) :: exner_dyn_incr_before
    REAL(wp), DIMENSION(nproma, p_patch%nlev) :: ddt_exner_phy_before
    REAL(wp), DIMENSION(nproma, p_patch%nlev) :: grf_tend_rho_before
    REAL(wp), DIMENSION(nproma, p_patch%nlev) :: grf_tend_thv_before
    REAL(wp), DIMENSION(nproma, p_patch%nlevp1) :: grf_tend_w_before
    REAL(wp), DIMENSION(nproma, p_patch%nlev) :: mass_fl_e_before
    REAL(wp), DIMENSION(nproma, p_patch%nlev) :: ddt_vn_phy_before
    REAL(wp), DIMENSION(nproma, p_patch%nlev) :: grf_tend_vn_before
    REAL(wp), DIMENSION(nproma, p_patch%nlevp1) :: vn_ie_before
    REAL(wp), DIMENSION(nproma, p_patch%nlev) :: vt_before
    REAL(wp), DIMENSION(:,:), ALLOCATABLE, TARGET :: vn_incr_before
    REAL(wp), DIMENSION(:,:), ALLOCATABLE, TARGET :: rho_incr_before
    REAL(wp), DIMENSION(:,:), ALLOCATABLE, TARGET :: exner_incr_before
    REAL(wp), DIMENSION(nproma, p_patch%nlev) :: mass_flx_me_before
    REAL(wp), DIMENSION(nproma, p_patch%nlevp1) :: mass_flx_ic_before
    REAL(wp), DIMENSION(nproma, p_patch%nlevp1) :: vol_flx_ic_before
    REAL(wp), DIMENSION(nproma, p_patch%nlev) :: vn_traj_before

    INTEGER :: i_startblk, i_endblk
    INTEGER :: ic_start, ic_end
    INTEGER :: ie_start, ie_end
    INTEGER :: ic_start_without_halo, ic_end_without_halo
    INTEGER :: rl_start, rl_end

    IF (ASSOCIATED(p_nh%diag%vn_incr)) ALLOCATE(vn_incr_before(nproma, p_patch%nlev))
    IF (ASSOCIATED(p_nh%diag%rho_incr)) ALLOCATE(rho_incr_before(nproma, p_patch%nlev))
    IF (ASSOCIATED(p_nh%diag%exner_incr)) ALLOCATE(exner_incr_before(nproma, p_patch%nlev))

      CALL assert_acc_device_only("icon4py_solve_nh_verify", lacc)

      WRITE(message_text,'(a)') 'dycore_run ...'
      CALL message('icon4py_solve_nh_verify', message_text)

      CALL message('Icon4Py', 'Dycore granule verifying')

      ! Array arguments
      WRITE(message_text,'(a,i6,a,i6,a,f10.3,a,f10.3)') 'rho_now shape: ', SIZE(p_nh%prog(nnow)%rho, 1), ' x ', SIZE(p_nh%prog(nnow)%rho, 2), ' min: ', MINVAL(p_nh%prog(nnow)%rho(:,:,1)), ' max: ', MAXVAL(p_nh%prog(nnow)%rho(:,:,1))
      CALL message('', message_text)

      WRITE(message_text,'(a,i6,a,i6,a,f10.3,a,f10.3)') 'rho_new shape: ', SIZE(p_nh%prog(nnew)%rho, 1), ' x ', SIZE(p_nh%prog(nnew)%rho, 2), ' min: ', MINVAL(p_nh%prog(nnew)%rho(:,:,1)), ' max: ', MAXVAL(p_nh%prog(nnew)%rho(:,:,1))
      CALL message('', message_text)

      WRITE(message_text,'(a,i6,a,i6,a,f10.3,a,f10.3)') 'exner_now shape: ', SIZE(p_nh%prog(nnow)%exner, 1), ' x ', SIZE(p_nh%prog(nnow)%exner, 2), ' min: ', MINVAL(p_nh%prog(nnow)%exner(:,:,1)), ' max: ', MAXVAL(p_nh%prog(nnow)%exner(:,:,1))
      CALL message('', message_text)

      WRITE(message_text,'(a,i6,a,i6,a,f10.3,a,f10.3)') 'exner_new shape: ', SIZE(p_nh%prog(nnew)%exner, 1), ' x ', SIZE(p_nh%prog(nnew)%exner, 2), ' min: ', MINVAL(p_nh%prog(nnew)%exner(:,:,1)), ' max: ', MAXVAL(p_nh%prog(nnew)%exner(:,:,1))
      CALL message('', message_text)

      WRITE(message_text,'(a,i6,a,i6,a,f10.3,a,f10.3)') 'w_now shape: ', SIZE(p_nh%prog(nnow)%w, 1), ' x ', SIZE(p_nh%prog(nnow)%w, 2), ' min: ', MINVAL(p_nh%prog(nnow)%w(:,:,1)), ' max: ', MAXVAL(p_nh%prog(nnow)%w(:,:,1))
      CALL message('', message_text)

      WRITE(message_text,'(a,i6,a,i6,a,f10.3,a,f10.3)') 'w_new shape: ', SIZE(p_nh%prog(nnew)%w, 1), ' x ', SIZE(p_nh%prog(nnew)%w, 2), ' min: ', MINVAL(p_nh%prog(nnew)%w(:,:,1)), ' max: ', MAXVAL(p_nh%prog(nnew)%w(:,:,1))
      CALL message('', message_text)

      WRITE(message_text,'(a,i6,a,i6,a,f10.3,a,f10.3)') 'theta_v_now shape: ', SIZE(p_nh%prog(nnow)%theta_v, 1), ' x ', SIZE(p_nh%prog(nnow)%theta_v, 2), ' min: ', MINVAL(p_nh%prog(nnow)%theta_v(:,:,1)), ' max: ', MAXVAL(p_nh%prog(nnow)%theta_v(:,:,1))
      CALL message('', message_text)

      WRITE(message_text,'(a,i6,a,i6,a,f10.3,a,f10.3)') 'theta_v_new shape: ', SIZE(p_nh%prog(nnew)%theta_v, 1), ' x ', SIZE(p_nh%prog(nnew)%theta_v, 2), ' min: ', MINVAL(p_nh%prog(nnew)%theta_v(:,:,1)), ' max: ', MAXVAL(p_nh%prog(nnew)%theta_v(:,:,1))
      CALL message('', message_text)

      WRITE(message_text,'(a,i6,a,i6,a,f10.3,a,f10.3)') 'vn_now shape: ', SIZE(p_nh%prog(nnow)%vn, 1), ' x ', SIZE(p_nh%prog(nnow)%vn, 2), ' min: ', MINVAL(p_nh%prog(nnow)%vn(:,:,1)), ' max: ', MAXVAL(p_nh%prog(nnow)%vn(:,:,1))
      CALL message('', message_text)

      WRITE(message_text,'(a,i6,a,i6,a,f10.3,a,f10.3)') 'vn_new shape: ', SIZE(p_nh%prog(nnew)%vn, 1), ' x ', SIZE(p_nh%prog(nnew)%vn, 2), ' min: ', MINVAL(p_nh%prog(nnew)%vn(:,:,1)), ' max: ', MAXVAL(p_nh%prog(nnew)%vn(:,:,1))
      CALL message('', message_text)

      WRITE(message_text,'(a,i6,a,i6,a,f10.3,a,f10.3)') 'w_concorr_c shape: ', SIZE(p_nh%diag%w_concorr_c, 1), ' x ', SIZE(p_nh%diag%w_concorr_c, 2), ' min: ', MINVAL(p_nh%diag%w_concorr_c(:,:,1)), ' max: ', MAXVAL(p_nh%diag%w_concorr_c(:,:,1))
      CALL message('', message_text)

      WRITE(message_text,'(a,i6,a,i6,a,f10.3,a,f10.3)') 'ddt_vn_apc_nnow shape: ', SIZE(p_nh%diag%ddt_vn_apc_pc, 1), ' x ', SIZE(p_nh%diag%ddt_vn_apc_pc, 2), ' min: ', MINVAL(p_nh%diag%ddt_vn_apc_pc(:,:,1,nnow)), ' max: ', MAXVAL(p_nh%diag%ddt_vn_apc_pc(:,:,1,nnow))
      CALL message('', message_text)

      WRITE(message_text,'(a,i6,a,i6,a,f10.3,a,f10.3)') 'ddt_vn_apc_nnew shape: ', SIZE(p_nh%diag%ddt_vn_apc_pc, 1), ' x ', SIZE(p_nh%diag%ddt_vn_apc_pc, 2), ' min: ', MINVAL(p_nh%diag%ddt_vn_apc_pc(:,:,1,nnew)), ' max: ', MAXVAL(p_nh%diag%ddt_vn_apc_pc(:,:,1,nnew))
      CALL message('', message_text)

      WRITE(message_text,'(a,i6,a,i6,a,f10.3,a,f10.3)') 'ddt_w_adv_nnow shape: ', SIZE(p_nh%diag%ddt_w_adv_pc, 1), ' x ', SIZE(p_nh%diag%ddt_w_adv_pc, 2), ' min: ', MINVAL(p_nh%diag%ddt_w_adv_pc(:,:,1,nnow)), ' max: ', MAXVAL(p_nh%diag%ddt_w_adv_pc(:,:,1,nnow))
      CALL message('', message_text)

      WRITE(message_text,'(a,i6,a,i6,a,f10.3,a,f10.3)') 'ddt_w_adv_nnew shape: ', SIZE(p_nh%diag%ddt_w_adv_pc, 1), ' x ', SIZE(p_nh%diag%ddt_w_adv_pc, 2), ' min: ', MINVAL(p_nh%diag%ddt_w_adv_pc(:,:,1,nnew)), ' max: ', MAXVAL(p_nh%diag%ddt_w_adv_pc(:,:,1,nnew))
      CALL message('', message_text)

      WRITE(message_text,'(a,i6,a,i6,a,f10.3,a,f10.3)') 'theta_v_ic shape: ', SIZE(p_nh%diag%theta_v_ic, 1), ' x ', SIZE(p_nh%diag%theta_v_ic, 2), ' min: ', MINVAL(p_nh%diag%theta_v_ic(:,:,1)), ' max: ', MAXVAL(p_nh%diag%theta_v_ic(:,:,1))
      CALL message('', message_text)

      WRITE(message_text,'(a,i6,a,i6,a,f10.3,a,f10.3)') 'rho_ic shape: ', SIZE(p_nh%diag%rho_ic, 1), ' x ', SIZE(p_nh%diag%rho_ic, 2), ' min: ', MINVAL(p_nh%diag%rho_ic(:,:,1)), ' max: ', MAXVAL(p_nh%diag%rho_ic(:,:,1))
      CALL message('', message_text)

      WRITE(message_text,'(a,i6,a,i6,a,f10.3,a,f10.3)') 'exner_pr shape: ', SIZE(p_nh%diag%exner_pr, 1), ' x ', SIZE(p_nh%diag%exner_pr, 2), ' min: ', MINVAL(p_nh%diag%exner_pr(:,:,1)), ' max: ', MAXVAL(p_nh%diag%exner_pr(:,:,1))
      CALL message('', message_text)

      WRITE(message_text,'(a,i6,a,i6,a,f10.3,a,f10.3)') 'exner_dyn_incr shape: ', SIZE(p_nh%diag%exner_dyn_incr, 1), ' x ', SIZE(p_nh%diag%exner_dyn_incr, 2), ' min: ', MINVAL(p_nh%diag%exner_dyn_incr(:,:,1)), ' max: ', MAXVAL(p_nh%diag%exner_dyn_incr(:,:,1))
      CALL message('', message_text)

      WRITE(message_text,'(a,i6,a,i6,a,f10.3,a,f10.3)') 'ddt_exner_phy shape: ', SIZE(p_nh%diag%ddt_exner_phy, 1), ' x ', SIZE(p_nh%diag%ddt_exner_phy, 2), ' min: ', MINVAL(p_nh%diag%ddt_exner_phy(:,:,1)), ' max: ', MAXVAL(p_nh%diag%ddt_exner_phy(:,:,1))
      CALL message('', message_text)

      WRITE(message_text,'(a,i6,a,i6,a,f10.3,a,f10.3)') 'grf_tend_rho shape: ', SIZE(p_nh%diag%grf_tend_rho, 1), ' x ', SIZE(p_nh%diag%grf_tend_rho, 2), ' min: ', MINVAL(p_nh%diag%grf_tend_rho(:,:,1)), ' max: ', MAXVAL(p_nh%diag%grf_tend_rho(:,:,1))
      CALL message('', message_text)

      WRITE(message_text,'(a,i6,a,i6,a,f10.3,a,f10.3)') 'grf_tend_thv shape: ', SIZE(p_nh%diag%grf_tend_thv, 1), ' x ', SIZE(p_nh%diag%grf_tend_thv, 2), ' min: ', MINVAL(p_nh%diag%grf_tend_thv(:,:,1)), ' max: ', MAXVAL(p_nh%diag%grf_tend_thv(:,:,1))
      CALL message('', message_text)

      WRITE(message_text,'(a,i6,a,i6,a,f10.3,a,f10.3)') 'grf_tend_w shape: ', SIZE(p_nh%diag%grf_tend_w, 1), ' x ', SIZE(p_nh%diag%grf_tend_w, 2), ' min: ', MINVAL(p_nh%diag%grf_tend_w(:,:,1)), ' max: ', MAXVAL(p_nh%diag%grf_tend_w(:,:,1))
      CALL message('', message_text)

      WRITE(message_text,'(a,i6,a,i6,a,f10.3,a,f10.3)') 'mass_fl_e shape: ', SIZE(p_nh%diag%mass_fl_e, 1), ' x ', SIZE(p_nh%diag%mass_fl_e, 2), ' min: ', MINVAL(p_nh%diag%mass_fl_e(:,:,1)), ' max: ', MAXVAL(p_nh%diag%mass_fl_e(:,:,1))
      CALL message('', message_text)

      WRITE(message_text,'(a,i6,a,i6,a,f10.3,a,f10.3)') 'ddt_vn_phy shape: ', SIZE(p_nh%diag%ddt_vn_phy, 1), ' x ', SIZE(p_nh%diag%ddt_vn_phy, 2), ' min: ', MINVAL(p_nh%diag%ddt_vn_phy(:,:,1)), ' max: ', MAXVAL(p_nh%diag%ddt_vn_phy(:,:,1))
      CALL message('', message_text)

      WRITE(message_text,'(a,i6,a,i6,a,f10.3,a,f10.3)') 'grf_tend_vn shape: ', SIZE(p_nh%diag%grf_tend_vn, 1), ' x ', SIZE(p_nh%diag%grf_tend_vn, 2), ' min: ', MINVAL(p_nh%diag%grf_tend_vn(:,:,1)), ' max: ', MAXVAL(p_nh%diag%grf_tend_vn(:,:,1))
      CALL message('', message_text)

      WRITE(message_text,'(a,i6,a,i6,a,f10.3,a,f10.3)') 'vn_ie shape: ', SIZE(p_nh%diag%vn_ie, 1), ' x ', SIZE(p_nh%diag%vn_ie, 2), ' min: ', MINVAL(p_nh%diag%vn_ie(:,:,1)), ' max: ', MAXVAL(p_nh%diag%vn_ie(:,:,1))
      CALL message('', message_text)

      WRITE(message_text,'(a,i6,a,i6,a,f10.3,a,f10.3)') 'vt shape: ', SIZE(p_nh%diag%vt, 1), ' x ', SIZE(p_nh%diag%vt, 2), ' min: ', MINVAL(p_nh%diag%vt(:,:,1)), ' max: ', MAXVAL(p_nh%diag%vt(:,:,1))
      CALL message('', message_text)

      WRITE(message_text,'(a,i6,a,i6,a,f10.3,a,f10.3)') 'vn_incr shape: ', SIZE(p_nh%diag%vn_incr, 1), ' x ', SIZE(p_nh%diag%vn_incr, 2), ' min: ', MINVAL(p_nh%diag%vn_incr(:,:,1)), ' max: ', MAXVAL(p_nh%diag%vn_incr(:,:,1))
      CALL message('', message_text)

      WRITE(message_text,'(a,i6,a,i6,a,f10.3,a,f10.3)') 'rho_incr shape: ', SIZE(p_nh%diag%rho_incr, 1), ' x ', SIZE(p_nh%diag%rho_incr, 2), ' min: ', MINVAL(p_nh%diag%rho_incr(:,:,1)), ' max: ', MAXVAL(p_nh%diag%rho_incr(:,:,1))
      CALL message('', message_text)

      WRITE(message_text,'(a,i6,a,i6,a,f10.3,a,f10.3)') 'exner_incr shape: ', SIZE(p_nh%diag%exner_incr, 1), ' x ', SIZE(p_nh%diag%exner_incr, 2), ' min: ', MINVAL(p_nh%diag%exner_incr(:,:,1)), ' max: ', MAXVAL(p_nh%diag%exner_incr(:,:,1))
      CALL message('', message_text)

      WRITE(message_text,'(a,i6,a,i6,a,f10.3,a,f10.3)') 'mass_flx_me shape: ', SIZE(prep_adv%mass_flx_me, 1), ' x ', SIZE(prep_adv%mass_flx_me, 2), ' min: ', MINVAL(prep_adv%mass_flx_me(:,:,1)), ' max: ', MAXVAL(prep_adv%mass_flx_me(:,:,1))
      CALL message('', message_text)

      WRITE(message_text,'(a,i6,a,i6,a,f10.3,a,f10.3)') 'mass_flx_ic shape: ', SIZE(prep_adv%mass_flx_ic, 1), ' x ', SIZE(prep_adv%mass_flx_ic, 2), ' min: ', MINVAL(prep_adv%mass_flx_ic(:,:,1)), ' max: ', MAXVAL(prep_adv%mass_flx_ic(:,:,1))
      CALL message('', message_text)

      WRITE(message_text,'(a,i6,a,i6,a,f10.3,a,f10.3)') 'vol_flx_ic shape: ', SIZE(prep_adv%vol_flx_ic, 1), ' x ', SIZE(prep_adv%vol_flx_ic, 2), ' min: ', MINVAL(prep_adv%vol_flx_ic(:,:,1)), ' max: ', MAXVAL(prep_adv%vol_flx_ic(:,:,1))
      CALL message('', message_text)

      WRITE(message_text,'(a,i6,a,i6,a,f10.3,a,f10.3)') 'vn_traj shape: ', SIZE(prep_adv%vn_traj, 1), ' x ', SIZE(prep_adv%vn_traj, 2), ' min: ', MINVAL(prep_adv%vn_traj(:,:,1)), ' max: ', MAXVAL(prep_adv%vn_traj(:,:,1))
      CALL message('', message_text)

      ! Scalar arguments
      WRITE(message_text,'(a,f10.3)') 'dtime: ', dtime
      CALL message('', message_text)

      WRITE(message_text,'(a,l1)') 'lprep_adv: ', LOGICAL(lprep_adv, c_int)
      CALL message('', message_text)

      WRITE(message_text,'(a,l1)') 'clean_mflx: ', lclean_mflx
      CALL message('', message_text)

      WRITE(message_text,'(a,l1)') 'recompute: ', l_recompute
      CALL message('', message_text)

      WRITE(message_text,'(a,l1)') 'linit: ', l_init
      CALL message('', message_text)

      WRITE(message_text,'(a,f10.3)') 'divdamp_fac_o2: ', REAL(divdamp_fac_o2, c_double)
      CALL message('', message_text)

      WRITE(message_text,'(a,f10.3)') 'ndyn_substeps_var: ', REAL(ndyn_substeps_var(1), c_double)
      CALL message('', message_text)

      WRITE(message_text,'(a,i6)') 'idyn_timestep: ', idyn_timestep
      CALL message('', message_text)

      WRITE(message_text,'(a,i6)') 'nnow: ', nnow
      CALL message('', message_text)

      WRITE(message_text,'(a,i6)') 'nnew: ', nnew
      CALL message('', message_text)

      !$ACC DATA &
      !$ACC   CREATE(rho_now_before) &
      !$ACC   CREATE(rho_new_before) &
      !$ACC   CREATE(exner_now_before) &
      !$ACC   CREATE(exner_new_before) &
      !$ACC   CREATE(w_now_before) &
      !$ACC   CREATE(w_new_before) &
      !$ACC   CREATE(theta_v_now_before) &
      !$ACC   CREATE(theta_v_new_before) &
      !$ACC   CREATE(vn_now_before) &
      !$ACC   CREATE(vn_new_before) &
      !$ACC   CREATE(w_concorr_c_before) &
      !$ACC   CREATE(ddt_vn_apc_nnow_before) &
      !$ACC   CREATE(ddt_vn_apc_nnew_before) &
      !$ACC   CREATE(ddt_w_adv_nnow_before) &
      !$ACC   CREATE(ddt_w_adv_nnew_before) &
      !$ACC   CREATE(theta_v_ic_before) &
      !$ACC   CREATE(rho_ic_before) &
      !$ACC   CREATE(exner_pr_before) &
      !$ACC   CREATE(exner_dyn_incr_before) &
      !$ACC   CREATE(ddt_exner_phy_before) &
      !$ACC   CREATE(grf_tend_rho_before) &
      !$ACC   CREATE(grf_tend_thv_before) &
      !$ACC   CREATE(grf_tend_w_before) &
      !$ACC   CREATE(mass_fl_e_before) &
      !$ACC   CREATE(ddt_vn_phy_before) &
      !$ACC   CREATE(grf_tend_vn_before) &
      !$ACC   CREATE(vn_ie_before) &
      !$ACC   CREATE(vt_before) &
      !$ACC   CREATE(mass_flx_me_before) &
      !$ACC   CREATE(mass_flx_ic_before) &
      !$ACC   CREATE(vol_flx_ic_before) &
      !$ACC   CREATE(vn_traj_before)
      !$ACC DATA CREATE(vn_incr_before) IF(ALLOCATED(vn_incr_before))
      !$ACC DATA CREATE(rho_incr_before) IF(ALLOCATED(rho_incr_before))
      !$ACC DATA CREATE(exner_incr_before) IF(ALLOCATED(exner_incr_before))

      CALL copy(p_nh%prog(nnow)%rho(:,:,1),          rho_now_before, lacc=.TRUE.)
      CALL copy(p_nh%prog(nnew)%rho(:,:,1),          rho_new_before, lacc=.TRUE.)
      CALL copy(p_nh%prog(nnow)%exner(:,:,1),        exner_now_before, lacc=.TRUE.)
      CALL copy(p_nh%prog(nnew)%exner(:,:,1),        exner_new_before, lacc=.TRUE.)
      CALL copy(p_nh%prog(nnow)%w(:,:,1),            w_now_before, lacc=.TRUE.)
      CALL copy(p_nh%prog(nnew)%w(:,:,1),            w_new_before, lacc=.TRUE.)
      CALL copy(p_nh%prog(nnow)%theta_v(:,:,1),      theta_v_now_before, lacc=.TRUE.)
      CALL copy(p_nh%prog(nnew)%theta_v(:,:,1),      theta_v_new_before, lacc=.TRUE.)
      CALL copy(p_nh%prog(nnow)%vn(:,:,1),           vn_now_before, lacc=.TRUE.)
      CALL copy(p_nh%prog(nnew)%vn(:,:,1),           vn_new_before, lacc=.TRUE.)
      CALL copy(p_nh%diag%w_concorr_c(:,:,1),        w_concorr_c_before, lacc=.TRUE.)
      CALL copy(p_nh%diag%ddt_vn_apc_pc(:,:,1,nnow), ddt_vn_apc_nnow_before, lacc=.TRUE.)
      CALL copy(p_nh%diag%ddt_vn_apc_pc(:,:,1,nnew), ddt_vn_apc_nnew_before, lacc=.TRUE.)
      CALL copy(p_nh%diag%ddt_w_adv_pc(:,:,1,nnow),  ddt_w_adv_nnow_before, lacc=.TRUE.)
      CALL copy(p_nh%diag%ddt_w_adv_pc(:,:,1,nnew),  ddt_w_adv_nnew_before, lacc=.TRUE.)
      CALL copy(p_nh%diag%theta_v_ic(:,:,1),         theta_v_ic_before, lacc=.TRUE.)
      CALL copy(p_nh%diag%rho_ic(:,:,1),             rho_ic_before, lacc=.TRUE.)
      CALL copy(p_nh%diag%exner_pr(:,:,1),           exner_pr_before, lacc=.TRUE.)
      CALL copy(p_nh%diag%exner_dyn_incr(:,:,1),     exner_dyn_incr_before, lacc=.TRUE.)
      CALL copy(p_nh%diag%ddt_exner_phy(:,:,1),      ddt_exner_phy_before, lacc=.TRUE.)
      CALL copy(p_nh%diag%grf_tend_rho(:,:,1),       grf_tend_rho_before, lacc=.TRUE.)
      CALL copy(p_nh%diag%grf_tend_thv(:,:,1),       grf_tend_thv_before, lacc=.TRUE.)
      CALL copy(p_nh%diag%grf_tend_w(:,:,1),         grf_tend_w_before, lacc=.TRUE.)
      CALL copy(p_nh%diag%mass_fl_e(:,:,1),          mass_fl_e_before, lacc=.TRUE.)
      CALL copy(p_nh%diag%ddt_vn_phy(:,:,1),         ddt_vn_phy_before, lacc=.TRUE.)
      CALL copy(p_nh%diag%grf_tend_vn(:,:,1),        grf_tend_vn_before, lacc=.TRUE.)
      CALL copy(p_nh%diag%vn_ie(:,:,1),              vn_ie_before, lacc=.TRUE.)
      CALL copy(p_nh%diag%vt(:,:,1),                 vt_before, lacc=.TRUE.)
      if (ALLOCATED(vn_incr_before)) then
          CALL copy(p_nh%diag%vn_incr(:,:,1),            vn_incr_before, lacc=.TRUE.)
      endif
      if (ALLOCATED(rho_incr_before)) then
          CALL copy(p_nh%diag%rho_incr(:,:,1),           rho_incr_before, lacc=.TRUE.)
      endif
      if (ALLOCATED(exner_incr_before)) then
          CALL copy(p_nh%diag%exner_incr(:,:,1),         exner_incr_before, lacc=.TRUE.)
      endif
      CALL copy(prep_adv%mass_flx_me(:,:,1),         mass_flx_me_before, lacc=.TRUE.)
      CALL copy(prep_adv%mass_flx_ic(:,:,1),         mass_flx_ic_before, lacc=.TRUE.)
      CALL copy(prep_adv%vol_flx_ic(:,:,1),          vol_flx_ic_before, lacc=.TRUE.)
      CALL copy(prep_adv%vn_traj(:,:,1),             vn_traj_before, lacc=.TRUE.)

      CALL solve_nh(p_nh, p_patch, p_int, prep_adv,                 &
        &           nnow, nnew, l_init, l_recompute,                &
        &           lsave_mflx, lprep_adv, lclean_mflx,             &
        &           idyn_timestep, jstep, dtime, lacc)

      IF(ALLOCATED(vn_incr_before)) vn_incr_ptr => vn_incr_before
      IF(ALLOCATED(rho_incr_before)) rho_incr_ptr => rho_incr_before
      IF(ALLOCATED(exner_incr_before)) exner_incr_ptr => exner_incr_before
      !$ACC WAIT
      max_vcfl_size1_array(1) = p_nh%diag%max_vcfl_dyn ! max_vcfl is passed as a size 1 (host) array
      CALL solve_nh_run( rho_now               = rho_now_before,                        &
                         rho_new               = rho_new_before,                        &
                         exner_now             = exner_now_before,                      &
                         exner_new             = exner_new_before,                      &
                         w_now                 = w_now_before,                          &
                         w_new                 = w_new_before,                          &
                         theta_v_now           = theta_v_now_before,                    &
                         theta_v_new           = theta_v_new_before,                    &
                         vn_now                = vn_now_before,                         &
                         vn_new                = vn_new_before,                         &
                         w_concorr_c           = w_concorr_c_before,                    &
                         ddt_vn_apc_ntl1       = ddt_vn_apc_nnow_before,                &
                         ddt_vn_apc_ntl2       = ddt_vn_apc_nnew_before,                &
                         ddt_w_adv_ntl1        = ddt_w_adv_nnow_before,                 &
                         ddt_w_adv_ntl2        = ddt_w_adv_nnew_before,                 &
                         theta_v_ic            = theta_v_ic_before,                     &
                         rho_ic                = rho_ic_before,                         &
                         exner_pr              = exner_pr_before,                       &
                         exner_dyn_incr        = exner_dyn_incr_before,                 &
                         ddt_exner_phy         = ddt_exner_phy_before,                  &
                         grf_tend_rho          = grf_tend_rho_before,                   &
                         grf_tend_thv          = grf_tend_thv_before,                   &
                         grf_tend_w            = grf_tend_w_before,                     &
                         mass_fl_e             = mass_fl_e_before,                      &
                         ddt_vn_phy            = ddt_vn_phy_before,                     &
                         grf_tend_vn           = grf_tend_vn_before,                    &
                         vn_ie                 = vn_ie_before,                          &
                         vt                    = vt_before,                             &
                         vn_incr               = vn_incr_ptr,                           &
                         rho_incr              = rho_incr_ptr,                          &
                         exner_incr            = exner_incr_ptr,                        &
                         mass_flx_me           = mass_flx_me_before,                    &
                         mass_flx_ic           = mass_flx_ic_before,                    &
                         vol_flx_ic            = vol_flx_ic_before,                     &
                         vn_traj               = vn_traj_before,                        &
                         dtime                 = dtime,                                 &
                         max_vcfl_size1_array  = max_vcfl_size1_array,                  &
                         lprep_adv             = LOGICAL(lprep_adv, c_bool),            &
                         at_initial_timestep   = LOGICAL(l_init, c_bool),               &
                         divdamp_fac_o2        = REAL(divdamp_fac_o2, c_double),        &
                         ndyn_substeps_var     = ndyn_substeps_var(1),                  &
                         idyn_timestep         = idyn_timestep,                         &
                         is_iau_active         = LOGICAL(lis_iau_active, c_bool),        &
                         iau_wgt_dyn           = iau_wgt_dyn,                           &
                         rc                    = rc )
      !$ACC WAIT
      CALL icon4py_check_rc(rc, 'icon4py_solve_nh_verify', 'solve_nh_run')
      p_nh%diag%max_vcfl_dyn = max_vcfl_size1_array(1) ! retrieve max_vcfl from the array

      rl_start = 1
      rl_end   = min_rlcell

      i_startblk = p_patch%cells%start_block(rl_start)
      i_endblk   = p_patch%cells%end_block(rl_end)

      CALL get_indices_c(p_patch, 1, i_startblk, i_endblk, ic_start, ic_end, rl_start, rl_end)

      rl_start = 1
      rl_end   = min_rledge

      i_startblk = p_patch%edges%start_block(rl_start)
      i_endblk   = p_patch%edges%end_block(rl_end)

      CALL get_indices_e(p_patch, 1, i_startblk, i_endblk, ie_start, ie_end, rl_start, rl_end)

      rl_start = grf_bdywidth_c+1
      rl_end   = min_rlcell_int

      i_startblk = p_patch%cells%start_block(rl_start)
      i_endblk   = p_patch%cells%end_block(rl_end)

      CALL get_indices_c(p_patch, 1, i_startblk, i_endblk, ic_start_without_halo, ic_end_without_halo, rl_start, rl_end)

      CALL print_verify_header('domain')

      CALL verify_field_wp( &
         field_name='rho_now', &
         field_ref=p_nh%prog(nnow)%rho(:,:,1), &
         field=rho_now_before, &
         i_start=ic_start, &
         i_end=ic_end)
      CALL verify_field_wp( &
         field_name='rho_new', &
         field_ref=p_nh%prog(nnew)%rho(:,:,1), &
         field=rho_new_before, &
         i_start=ic_start, &
         i_end=ic_end)
      CALL verify_field_wp( &
         field_name='exner_now', &
         field_ref=p_nh%prog(nnow)%exner(:,:,1), &
         field=exner_now_before, &
         i_start=ic_start, &
         i_end=ic_end)
      CALL verify_field_wp( &
         field_name='exner_new', &
         field_ref=p_nh%prog(nnew)%exner(:,:,1), &
         field=exner_new_before, &
         i_start=ic_start, &
         i_end=ic_end)
      CALL verify_field_wp( &
         field_name='w_now', &
         field_ref=p_nh%prog(nnow)%w(:,:,1), &
         field=w_now_before, &
         i_start=ic_start, &
         i_end=ic_end)
      CALL verify_field_wp( &
         field_name='w_new', &
         field_ref=p_nh%prog(nnew)%w(:,:,1), &
         field=w_new_before, &
         i_start=ic_start, &
         i_end=ic_end)
      CALL verify_field_wp( &
         field_name='theta_v_now', &
         field_ref=p_nh%prog(nnow)%theta_v(:,:,1), &
         field=theta_v_now_before, &
         i_start=ic_start, &
         i_end=ic_end)
      CALL verify_field_wp( &
         field_name='theta_v_new', &
         field_ref=p_nh%prog(nnew)%theta_v(:,:,1), &
         field=theta_v_new_before, &
         i_start=ic_start, &
         i_end=ic_end)
      CALL verify_field_wp( &
         field_name='vn_now', &
         field_ref=p_nh%prog(nnow)%vn(:,:,1), &
         field=vn_now_before, &
         i_start=ie_start, &
         i_end=ie_end)
      CALL verify_field_wp( &
         field_name='vn_new', &
         field_ref=p_nh%prog(nnew)%vn(:,:,1), &
         field=vn_new_before, &
         i_start=ie_start, &
         i_end=ie_end)
      CALL verify_field_wp( &
         field_name='w_concorr_c', &
         field_ref=p_nh%diag%w_concorr_c(:,:,1), &
         field=w_concorr_c_before, &
         i_start=ic_start, &
         i_end=ic_end)
      CALL verify_field_wp( &
         field_name='ddt_vn_apc_nnow', &
         field_ref=p_nh%diag%ddt_vn_apc_pc(:,:,1,nnow), &
         field=ddt_vn_apc_nnow_before, &
         i_start=ie_start, &
         i_end=ie_end)
      CALL verify_field_wp( &
         field_name='ddt_vn_apc_nnew', &
         field_ref=p_nh%diag%ddt_vn_apc_pc(:,:,1,nnew), &
         field=ddt_vn_apc_nnew_before, &
         i_start=ie_start, &
         i_end=ie_end)
      CALL verify_field_wp( &
         field_name='ddt_w_adv_nnow', &
         field_ref=p_nh%diag%ddt_w_adv_pc(:,:,1,nnow), &
          field=ddt_w_adv_nnow_before, &
         i_start=ic_start_without_halo, &
         i_end=ic_end_without_halo)
      CALL verify_field_wp( &
         field_name='ddt_w_adv_nnew', &
         field_ref=p_nh%diag%ddt_w_adv_pc(:,:,1,nnew), &
         field=ddt_w_adv_nnew_before, &
         i_start=ic_start_without_halo, &
         i_end=ic_end_without_halo)
      CALL verify_field_wp( &
         field_name='theta_v_ic', &
         field_ref=p_nh%diag%theta_v_ic(:,:,1), &
         field=theta_v_ic_before, &
         i_start=ic_start, &
         i_end=ic_end)
      CALL verify_field_wp( &
         field_name='rho_ic', &
         field_ref=p_nh%diag%rho_ic(:,:,1), &
         field=rho_ic_before, &
         i_start=ic_start, &
         i_end=ic_end_without_halo)
      CALL verify_field_wp( &
         field_name='exner_pr', &
         field_ref=p_nh%diag%exner_pr(:,:,1), &
         field=exner_pr_before, &
         i_start=ic_start, &
         i_end=ic_end)
      CALL verify_field_wp( &
         field_name='exner_dyn_incr', &
         field_ref=p_nh%diag%exner_dyn_incr(:,:,1), &
         field=exner_dyn_incr_before, &
         i_start=ic_start, &
         i_end=ic_end)
      CALL verify_field_wp( &
         field_name='ddt_exner_phy', &
         field_ref=p_nh%diag%ddt_exner_phy(:,:,1), &
         field=ddt_exner_phy_before, &
         i_start=ic_start, &
         i_end=ic_end)
      CALL verify_field_wp( &
         field_name='grf_tend_rho', &
         field_ref=p_nh%diag%grf_tend_rho(:,:,1), &
         field=grf_tend_rho_before, &
         i_start=ic_start, &
         i_end=ic_end)
      CALL verify_field_wp( &
         field_name='grf_tend_thv', &
         field_ref=p_nh%diag%grf_tend_thv(:,:,1), &
         field=grf_tend_thv_before, &
         i_start=ic_start, &
         i_end=ic_end)
      CALL verify_field_wp( &
         field_name='grf_tend_w', &
         field_ref=p_nh%diag%grf_tend_w(:,:,1), &
         field=grf_tend_w_before, &
         i_start=ic_start, &
         i_end=ic_end)
      CALL verify_field_wp( &
         field_name='mass_fl_e', &
         field_ref=p_nh%diag%mass_fl_e(:,:,1), &
         field=mass_fl_e_before, &
         i_start=ie_start, &
         i_end=ie_end)
      CALL verify_field_wp( &
         field_name='ddt_vn_phy', &
         field_ref=p_nh%diag%ddt_vn_phy(:,:,1), &
         field=ddt_vn_phy_before, &
         i_start=ie_start, &
         i_end=ie_end)
      CALL verify_field_wp( &
         field_name='grf_tend_vn', &
         field_ref=p_nh%diag%grf_tend_vn(:,:,1), &
         field=grf_tend_vn_before, &
         i_start=ie_start, &
         i_end=ie_end)
      CALL verify_field_wp( &
         field_name='vn_ie', &
         field_ref=p_nh%diag%vn_ie(:,:,1), &
         field=vn_ie_before, &
         i_start=ie_start, &
         i_end=ie_end)
      CALL verify_field_wp( &
         field_name='vt', &
         field_ref=p_nh%diag%vt(:,:,1), &
         field=vt_before, &
         i_start=ie_start, &
         i_end=ie_end)
      IF( ALLOCATED(vn_incr_before) ) CALL verify_field_wp( &
         field_name='vn_incr', &
         field_ref=p_nh%diag%vn_incr(:,:,1), &
         field=vn_incr_before, &
         i_start=ie_start, &
         i_end=ie_end)
      IF( ALLOCATED(rho_incr_before) ) CALL verify_field_wp( &
         field_name='rho_incr', &
         field_ref=p_nh%diag%rho_incr(:,:,1), &
         field=rho_incr_before, &
         i_start=ic_start, &
         i_end=ic_end)
      IF( ALLOCATED(exner_incr_before) ) CALL verify_field_wp( &
         field_name='exner_incr', &
         field_ref=p_nh%diag%exner_incr(:,:,1), &
         field=exner_incr_before, &
         i_start=ic_start, &
         i_end=ic_end)
      CALL verify_field_wp( &
         field_name='mass_flx_me', &
         field_ref=prep_adv%mass_flx_me(:,:,1), &
         field=mass_flx_me_before, &
         i_start=ie_start, &
         i_end=ie_end)
      CALL verify_field_wp( &
         field_name='mass_flx_ic', &
         field_ref=prep_adv%mass_flx_ic(:,:,1), &
         field=mass_flx_ic_before, &
         i_start=ic_start, &
         i_end=ic_end)
      CALL verify_field_wp( &
         field_name='vol_flx_ic', &
         field_ref=prep_adv%vol_flx_ic(:,:,1), &
         field=vol_flx_ic_before, &
         i_start=ic_start, &
         i_end=ic_end)
      CALL verify_field_wp( &
         field_name='vn_traj', &
         field_ref=prep_adv%vn_traj(:,:,1), &
         field=vn_traj_before, &
         i_start=ie_start, &
         i_end=ie_end)

      CALL print_verify_footer()

      CALL print_verify_header('padding')

      CALL verify_field_wp( &
         field_name='rho_now', &
         field_ref=p_nh%prog(nnow)%rho(:,:,1), &
         field=rho_now_before, &
         i_start=ic_end+1, &
         i_end=nproma)
      CALL verify_field_wp( &
         field_name='rho_new', &
         field_ref=p_nh%prog(nnew)%rho(:,:,1), &
         field=rho_new_before, &
         i_start=ic_end+1, &
         i_end=nproma)
      CALL verify_field_wp( &
         field_name='exner_now', &
         field_ref=p_nh%prog(nnow)%exner(:,:,1), &
         field=exner_now_before, &
         i_start=ic_end+1, &
         i_end=nproma)
      CALL verify_field_wp( &
         field_name='exner_new', &
         field_ref=p_nh%prog(nnew)%exner(:,:,1), &
         field=exner_new_before, &
         i_start=ic_end+1, &
         i_end=nproma)
      CALL verify_field_wp( &
         field_name='w_now', &
         field_ref=p_nh%prog(nnow)%w(:,:,1), &
         field=w_now_before, &
         i_start=ic_end+1, &
         i_end=nproma)
      CALL verify_field_wp( &
         field_name='w_new', &
         field_ref=p_nh%prog(nnew)%w(:,:,1), &
         field=w_new_before, &
         i_start=ic_end+1, &
         i_end=nproma)
      CALL verify_field_wp( &
         field_name='theta_v_now', &
         field_ref=p_nh%prog(nnow)%theta_v(:,:,1), &
         field=theta_v_now_before, &
         i_start=ic_end+1, &
         i_end=nproma)
      CALL verify_field_wp( &
         field_name='theta_v_new', &
         field_ref=p_nh%prog(nnew)%theta_v(:,:,1), &
         field=theta_v_new_before, &
         i_start=ic_end+1, &
         i_end=nproma)
      CALL verify_field_wp( &
         field_name='vn_now', &
         field_ref=p_nh%prog(nnow)%vn(:,:,1), &
         field=vn_now_before, &
         i_start=ie_end+1, &
         i_end=nproma)
      CALL verify_field_wp( &
         field_name='vn_new', &
         field_ref=p_nh%prog(nnew)%vn(:,:,1), &
         field=vn_new_before, &
         i_start=ie_end+1, &
         i_end=nproma)
      CALL verify_field_wp( &
         field_name='w_concorr_c', &
         field_ref=p_nh%diag%w_concorr_c(:,:,1), &
         field=w_concorr_c_before, &
         i_start=ic_end+1, &
         i_end=nproma)
      CALL verify_field_wp( &
         field_name='ddt_vn_apc_nnow', &
         field_ref=p_nh%diag%ddt_vn_apc_pc(:,:,1,nnow), &
         field=ddt_vn_apc_nnow_before, &
         i_start=ie_end+1, &
         i_end=nproma)
      CALL verify_field_wp( &
         field_name='ddt_vn_apc_nnew', &
         field_ref=p_nh%diag%ddt_vn_apc_pc(:,:,1,nnew), &
         field=ddt_vn_apc_nnew_before, &
         i_start=ie_end+1, &
         i_end=nproma)
      CALL verify_field_wp( &
         field_name='ddt_w_adv_nnow', &
         field_ref=p_nh%diag%ddt_w_adv_pc(:,:,1,nnow), &
         field=ddt_w_adv_nnow_before, &
         i_start=ic_end+1, &
         i_end=nproma)
      CALL verify_field_wp( &
         field_name='ddt_w_adv_nnew', &
         field_ref=p_nh%diag%ddt_w_adv_pc(:,:,1,nnew), &
         field=ddt_w_adv_nnew_before, &
         i_start=ic_end+1, &
         i_end=nproma)
      CALL verify_field_wp( &
         field_name='theta_v_ic', &
         field_ref=p_nh%diag%theta_v_ic(:,:,1), &
         field=theta_v_ic_before, &
         i_start=ic_end+1, &
         i_end=nproma)
      CALL verify_field_wp( &
         field_name='rho_ic', &
         field_ref=p_nh%diag%rho_ic(:,:,1), &
         field=rho_ic_before, &
         i_start=ic_end+1, &
         i_end=nproma)
      CALL verify_field_wp( &
         field_name='exner_pr', &
         field_ref=p_nh%diag%exner_pr(:,:,1), &
         field=exner_pr_before, &
         i_start=ic_end+1, &
         i_end=nproma)
      CALL verify_field_wp( &
         field_name='exner_dyn_incr', &
         field_ref=p_nh%diag%exner_dyn_incr(:,:,1), &
         field=exner_dyn_incr_before, &
         i_start=ic_end+1, &
         i_end=nproma)
      CALL verify_field_wp( &
         field_name='ddt_exner_phy', &
         field_ref=p_nh%diag%ddt_exner_phy(:,:,1), &
         field=ddt_exner_phy_before, &
         i_start=ic_end+1, &
         i_end=nproma)
      CALL verify_field_wp( &
         field_name='grf_tend_rho', &
         field_ref=p_nh%diag%grf_tend_rho(:,:,1), &
         field=grf_tend_rho_before, &
         i_start=ic_end+1, &
         i_end=nproma)
      CALL verify_field_wp( &
         field_name='grf_tend_thv', &
         field_ref=p_nh%diag%grf_tend_thv(:,:,1), &
         field=grf_tend_thv_before, &
         i_start=ic_end+1, &
         i_end=nproma)
      CALL verify_field_wp( &
         field_name='grf_tend_w', &
         field_ref=p_nh%diag%grf_tend_w(:,:,1), &
         field=grf_tend_w_before, &
         i_start=ic_end+1, &
         i_end=nproma)
      CALL verify_field_wp( &
         field_name='mass_fl_e', &
         field_ref=p_nh%diag%mass_fl_e(:,:,1), &
         field=mass_fl_e_before, &
         i_start=ie_end+1, &
         i_end=nproma)
      CALL verify_field_wp( &
         field_name='ddt_vn_phy', &
         field_ref=p_nh%diag%ddt_vn_phy(:,:,1), &
         field=ddt_vn_phy_before, &
         i_start=ie_end+1, &
         i_end=nproma)
      CALL verify_field_wp( &
         field_name='grf_tend_vn', &
         field_ref=p_nh%diag%grf_tend_vn(:,:,1), &
         field=grf_tend_vn_before, &
         i_start=ie_end+1, &
         i_end=nproma)
      CALL verify_field_wp( &
         field_name='vn_ie', &
         field_ref=p_nh%diag%vn_ie(:,:,1), &
         field=vn_ie_before, &
         i_start=ie_end+1, &
         i_end=nproma)
      CALL verify_field_wp( &
         field_name='vt', &
         field_ref=p_nh%diag%vt(:,:,1), &
         field=vt_before, &
         i_start=ie_end+1, &
         i_end=nproma)
      IF( ALLOCATED(vn_incr_before) ) CALL verify_field_wp( &
         field_name='vn_incr', &
         field_ref=p_nh%diag%vn_incr(:,:,1), &
         field=vn_incr_before, &
         i_start=ie_end+1, &
         i_end=nproma)
      IF( ALLOCATED(rho_incr_before) ) CALL verify_field_wp( &
         field_name='rho_incr', &
         field_ref=p_nh%diag%rho_incr(:,:,1), &
         field=rho_incr_before, &
         i_start=ic_end+1, &
         i_end=nproma)
      IF( ALLOCATED(exner_incr_before) ) CALL verify_field_wp( &
         field_name='exner_incr', &
         field_ref=p_nh%diag%exner_incr(:,:,1), &
         field=exner_incr_before, &
         i_start=ic_end+1, &
         i_end=nproma)
      CALL verify_field_wp( &
         field_name='mass_flx_me', &
         field_ref=prep_adv%mass_flx_me(:,:,1), &
         field=mass_flx_me_before, &
         i_start=ie_end+1, &
         i_end=nproma)
      CALL verify_field_wp( &
         field_name='mass_flx_ic', &
         field_ref=prep_adv%mass_flx_ic(:,:,1), &
         field=mass_flx_ic_before, &
         i_start=ic_end+1, &
         i_end=nproma)
      CALL verify_field_wp( &
         field_name='vol_flx_ic', &
         field_ref=prep_adv%vol_flx_ic(:,:,1), &
         field=vol_flx_ic_before, &
         i_start=ic_end+1, &
         i_end=nproma)
      CALL verify_field_wp( &
         field_name='vn_traj', &
         field_ref=prep_adv%vn_traj(:,:,1), &
         field=vn_traj_before, &
         i_start=ie_end+1, &
         i_end=nproma)

      CALL print_verify_footer()

      !$ACC END DATA
      !$ACC END DATA
      !$ACC END DATA
      !$ACC END DATA

  END SUBROUTINE icon4py_solve_nh_verify

  !-------------------------------------------------------------------------
  !>
  !! Dispatcher selecting the implementation according to the icon4py_mode
  !! runtime namelist option.
  !!
  SUBROUTINE icon4py_diffusion(p_nh_prog, p_nh_diag, p_nh_metrics, p_patch, p_int, &
      &                        dtime, linit, lacc)

    TYPE(t_patch), TARGET, INTENT(inout) :: p_patch
    TYPE(t_int_state),INTENT(in),TARGET  :: p_int
    TYPE(t_nh_prog), INTENT(inout)       :: p_nh_prog
    TYPE(t_nh_diag), INTENT(inout)       :: p_nh_diag
    TYPE(t_nh_metrics),INTENT(in),TARGET :: p_nh_metrics
    REAL(wp), INTENT(IN)                 :: dtime
    LOGICAL,  INTENT(IN)                 :: linit
    LOGICAL,  INTENT(IN)                 :: lacc

    SELECT CASE (icon4py_mode)
    CASE (ICON4PY_MODE_OFF)
      CALL diffusion(p_nh_prog, p_nh_diag, p_nh_metrics, p_patch, p_int, dtime, linit, lacc=lacc)
    CASE (ICON4PY_MODE_SUBSTITUTE)
      CALL icon4py_diffusion_substitute(p_nh_prog, p_nh_diag, p_nh_metrics, p_patch, p_int, dtime, linit, lacc)
    CASE (ICON4PY_MODE_VERIFY)
      CALL icon4py_diffusion_verify(p_nh_prog, p_nh_diag, p_nh_metrics, p_patch, p_int, dtime, linit, lacc)
    END SELECT

  END SUBROUTINE icon4py_diffusion


  SUBROUTINE  icon4py_diffusion_substitute(p_nh_prog, p_nh_diag, p_nh_metrics, p_patch, p_int, &
      &                         dtime, linit, lacc)

    TYPE(t_patch), TARGET, INTENT(inout) :: p_patch      !< single patch
    TYPE(t_int_state),INTENT(in),TARGET  :: p_int        !< single interpolation state
    TYPE(t_nh_prog), INTENT(inout)       :: p_nh_prog    !< single nh prognostic state
    TYPE(t_nh_diag), INTENT(inout)       :: p_nh_diag    !< single nh diagnostic state
    TYPE(t_nh_metrics),INTENT(in),TARGET :: p_nh_metrics !< single nh metric state
    REAL(vp), DIMENSION(:,:), POINTER    :: hdef_ic_ptr => NULL()
    REAL(vp), DIMENSION(:,:), POINTER    :: div_ic_ptr => NULL()
    REAL(vp), DIMENSION(:,:), POINTER    :: dwdx_ptr => NULL()
    REAL(vp), DIMENSION(:,:), POINTER    :: dwdy_ptr => NULL()
    REAL(wp), INTENT(IN)                 :: dtime        !< time step
    LOGICAL,  INTENT(IN)                 :: linit        !< initial call or runtime call
    LOGICAL,  INTENT(IN)                 :: lacc         ! If true, use openacc (if _OPENACC is enabled)

    INTEGER(c_int)                       :: rc

    CALL assert_acc_device_only("icon4py_diffusion_substitute", lacc)

    !$ACC DATA PRESENT(p_nh_diag%hdef_ic) IF(ASSOCIATED(p_nh_diag%hdef_ic))
    !$ACC DATA PRESENT(p_nh_diag%div_ic) IF(ASSOCIATED(p_nh_diag%div_ic))
    !$ACC DATA PRESENT(p_nh_diag%dwdx) IF(ASSOCIATED(p_nh_diag%dwdx))
    !$ACC DATA PRESENT(p_nh_diag%dwdy) IF(ASSOCIATED(p_nh_diag%dwdy))

    !$ACC DATA PRESENT(p_nh_prog%w, p_nh_prog%exner, p_nh_prog%vn) &
    !$ACC   PRESENT(p_nh_prog%theta_v, p_nh_prog%rho)

    IF (ASSOCIATED(p_nh_diag%hdef_ic)) hdef_ic_ptr => p_nh_diag%hdef_ic(:,:,1)
    IF (ASSOCIATED(p_nh_diag%div_ic)) div_ic_ptr => p_nh_diag%div_ic(:,:,1)
    IF (ASSOCIATED(p_nh_diag%dwdx)) dwdx_ptr => p_nh_diag%dwdx(:,:,1)
    IF (ASSOCIATED(p_nh_diag%dwdy)) dwdy_ptr => p_nh_diag%dwdy(:,:,1)
    !$ACC WAIT
    CALL diffusion_run( &
          w = p_nh_prog%w(:,:,1), &
          vn = p_nh_prog%vn(:,:,1), &
          exner = p_nh_prog%exner(:,:,1), &
          theta_v = p_nh_prog%theta_v(:,:,1), &
          rho = p_nh_prog%rho(:,:,1), &
          hdef_ic = hdef_ic_ptr, &
          div_ic = div_ic_ptr, &
          dwdx = dwdx_ptr, &
          dwdy = dwdy_ptr, &
          dtime = dtime, &
          linit = LOGICAL(linit, c_bool), &
          rc = rc)
    !$ACC WAIT
    CALL icon4py_check_rc(rc, 'icon4py_diffusion_substitute', 'diffusion_run')

    !$ACC END DATA
    !$ACC END DATA
    !$ACC END DATA
    !$ACC END DATA
    !$ACC END DATA

  END SUBROUTINE icon4py_diffusion_substitute


  SUBROUTINE  icon4py_diffusion_verify(p_nh_prog, p_nh_diag, p_nh_metrics, p_patch, p_int, &
      &                                dtime, linit, lacc)

    TYPE(t_patch), TARGET, INTENT(inout) :: p_patch
    TYPE(t_int_state),INTENT(in),TARGET  :: p_int
    TYPE(t_nh_prog), INTENT(inout)       :: p_nh_prog
    TYPE(t_nh_diag), INTENT(inout)       :: p_nh_diag
    TYPE(t_nh_metrics),INTENT(in),TARGET :: p_nh_metrics
    REAL(vp), DIMENSION(:,:), POINTER    :: hdef_ic_ptr => NULL()
    REAL(vp), DIMENSION(:,:), POINTER    :: div_ic_ptr => NULL()
    REAL(vp), DIMENSION(:,:), POINTER    :: dwdx_ptr => NULL()
    REAL(vp), DIMENSION(:,:), POINTER    :: dwdy_ptr => NULL()
    REAL(wp), INTENT(IN)                 :: dtime
    LOGICAL,  INTENT(IN)                 :: linit
    LOGICAL,  INTENT(IN)                 :: lacc

    INTEGER(c_int)                       :: rc

    REAL(wp), DIMENSION(nproma, p_patch%nlevp1, p_patch%nblks_c), TARGET :: w_before
    REAL(wp), DIMENSION(nproma, p_patch%nlev, p_patch%nblks_e), TARGET   :: vn_before
    REAL(wp), DIMENSION(nproma, p_patch%nlev, p_patch%nblks_c), TARGET   :: exner_before
    REAL(wp), DIMENSION(nproma, p_patch%nlev, p_patch%nblks_c), TARGET   :: theta_v_before
    REAL(wp), DIMENSION(nproma, p_patch%nlev, p_patch%nblks_c), TARGET   :: rho_before
    REAL(vp), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: hdef_ic_before
    REAL(vp), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: div_ic_before
    REAL(vp), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: dwdx_before
    REAL(vp), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: dwdy_before

    INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx
    INTEGER :: rl_start, rl_end

    IF (ASSOCIATED(p_nh_diag%hdef_ic)) ALLOCATE(hdef_ic_before(nproma, p_patch%nlevp1, p_patch%nblks_c))
    IF (ASSOCIATED(p_nh_diag%div_ic)) ALLOCATE(div_ic_before(nproma, p_patch%nlevp1, p_patch%nblks_c))
    IF (ASSOCIATED(p_nh_diag%dwdx)) ALLOCATE(dwdx_before(nproma, p_patch%nlevp1, p_patch%nblks_c))
    IF (ASSOCIATED(p_nh_diag%dwdy)) ALLOCATE(dwdy_before(nproma, p_patch%nlevp1, p_patch%nblks_c))

    CALL assert_acc_device_only("icon4py_diffusion_verify", lacc)

    !$ACC DATA PRESENT(p_nh_diag%hdef_ic) IF(ASSOCIATED(p_nh_diag%hdef_ic))
    !$ACC DATA PRESENT(p_nh_diag%div_ic) IF(ASSOCIATED(p_nh_diag%div_ic))
    !$ACC DATA PRESENT(p_nh_diag%dwdx) IF(ASSOCIATED(p_nh_diag%dwdx))
    !$ACC DATA PRESENT(p_nh_diag%dwdy) IF(ASSOCIATED(p_nh_diag%dwdy))

    !$ACC DATA PRESENT(p_nh_prog%w, p_nh_prog%exner, p_nh_prog%vn) &
    !$ACC   PRESENT(p_nh_prog%theta_v, p_nh_prog%rho)

    CALL message('Icon4py', 'Diffusion granule verifying', .TRUE.)

    !Before field copy
    !$ACC DATA &
    !$ACC   CREATE(w_before) &
    !$ACC   CREATE(vn_before) &
    !$ACC   CREATE(theta_v_before) &
    !$ACC   CREATE(exner_before) &
    !$ACC   CREATE(rho_before)
    !$ACC DATA CREATE(hdef_ic_before) IF(ALLOCATED(hdef_ic_before))
    !$ACC DATA CREATE(div_ic_before) IF(ALLOCATED(div_ic_before))
    !$ACC DATA CREATE(dwdx_before) IF(ALLOCATED(dwdx_before))
    !$ACC DATA CREATE(dwdy_before) IF(ALLOCATED(dwdy_before))

    CALL copy(p_nh_prog%w(: , : , :), w_before(:, :, :), lacc=.TRUE.)
    CALL copy(p_nh_prog%vn(: , : , :), vn_before(:, :, :), lacc=.TRUE.)
    CALL copy(p_nh_prog%exner(: , : , :), exner_before(:, :, :), lacc=.TRUE.)
    CALL copy(p_nh_prog%theta_v(: , : , :), theta_v_before(:, :, :), lacc=.TRUE.)
    CALL copy(p_nh_prog%rho(: , : , :), rho_before(:, :, :), lacc=.TRUE.)
    IF (ALLOCATED(hdef_ic_before)) then
        CALL copy(p_nh_diag%hdef_ic(: , : , :), hdef_ic_before(:, :, :), lacc=.TRUE.)
    ENDIF
    IF (ALLOCATED(div_ic_before)) then
        CALL copy(p_nh_diag%div_ic(: , : , :), div_ic_before(:, :, :), lacc=.TRUE.)
    ENDIF
    IF (ALLOCATED(dwdx_before)) then
        CALL copy(p_nh_diag%dwdx(: , : , :), dwdx_before(:, :, :), lacc=.TRUE.)
    ENDIF
    IF (ALLOCATED(dwdy_before)) then
        CALL copy(p_nh_diag%dwdy(: , : , :), dwdy_before(:, :, :), lacc=.TRUE.)
    ENDIF

    CALL diffusion(p_nh_prog, p_nh_diag, p_nh_metrics, p_patch, p_int, dtime, linit, .TRUE.)

    IF( ALLOCATED( hdef_ic_before ) ) hdef_ic_ptr => hdef_ic_before(:,:,1)
    IF( ALLOCATED( div_ic_before ) ) div_ic_ptr => div_ic_before(:,:,1)
    IF( ALLOCATED( dwdx_before ) ) dwdx_ptr => dwdx_before(:,:,1)
    IF( ALLOCATED( dwdy_before ) ) dwdy_ptr => dwdy_before(:,:,1)
    !$ACC WAIT
    CALL diffusion_run( &
          w = w_before(:,:,1), &
          vn = vn_before(:,:,1), &
          exner = exner_before(:,:,1), &
          theta_v = theta_v_before(:,:,1), &
          rho = rho_before(:,:,1), &
          hdef_ic = hdef_ic_ptr, &
          div_ic = div_ic_ptr, &
          dwdx = dwdx_ptr, &
          dwdy = dwdy_ptr, &
          dtime = dtime, &
          linit = LOGICAL(linit, c_bool), &
          rc = rc)
    !$ACC WAIT
    CALL icon4py_check_rc(rc, 'icon4py_diffusion_verify', 'diffusion_run')

    CALL print_verify_header('domain')

    rl_start = 1
    rl_end   = min_rledge

    i_startblk = p_patch%edges%start_block(rl_start)
    i_endblk   = p_patch%edges%end_block(rl_end)

    CALL get_indices_e(p_patch, 1, i_startblk, i_endblk, i_startidx, i_endidx, rl_start, rl_end)

    CALL verify_field_wp( &
        &   field_name='vn', &
        &   field_ref=p_nh_prog%vn(:, :, 1), &
        &   field=vn_before(:, :, 1), &
        &   i_start=i_startidx, &
        &   i_end=i_endidx)

    rl_start = 1
    rl_end   = min_rlcell

    i_startblk = p_patch%cells%start_block(rl_start)
    i_endblk   = p_patch%cells%end_block(rl_end)

    CALL get_indices_c(p_patch, 1, i_startblk, i_endblk, i_startidx, i_endidx, rl_start, rl_end)

    CALL verify_field_wp( &
        &   field_name='w', &
        &   field_ref=p_nh_prog%w(:, :, 1), &
        &   field=w_before(:, :, 1), &
        &   i_start=i_startidx, &
        &   i_end=i_endidx)
    CALL verify_field_wp( &
        &   field_name='exner', &
        &   field_ref=p_nh_prog%exner(:, :, 1), &
        &   field=exner_before(:, :, 1), &
        &   i_start=i_startidx, &
        &   i_end=i_endidx)
    CALL verify_field_wp( &
        &   field_name='theta_v', &
        &   field_ref=p_nh_prog%theta_v(:, :, 1), &
        &   field=theta_v_before(:, :, 1), &
        &   i_start=i_startidx, &
        &   i_end=i_endidx)
    CALL verify_field_wp(field_name='rho', &
        &   field_ref=p_nh_prog%rho(:, :, 1), &
        &   field=rho_before(:, :, 1), &
        &   i_start=i_startidx, &
        &   i_end=i_endidx)
    IF( ALLOCATED(hdef_ic_before) ) CALL verify_field_wp(&
        &   field_name='hdef_ic', &
        &   field_ref=p_nh_diag%hdef_ic(:, :, 1), &
        &   field=hdef_ic_before(:, :, 1), &
        &   i_start=i_startidx, &
        &   i_end=i_endidx)
    IF( ALLOCATED(div_ic_before) ) CALL verify_field_wp( &
        &   field_name='div_ic', &
        &   field_ref=p_nh_diag%div_ic(:, :, 1), &
        &   field=div_ic_before(:, :, 1), &
        &   i_start=i_startidx, &
        &   i_end=i_endidx)
    IF( ALLOCATED(dwdx_before) ) CALL verify_field_wp( &
        &   field_name='dwdx', &
        &   field_ref=p_nh_diag%dwdx(:, :, 1), &
        &   field=dwdx_before(:, :, 1), &
        &   i_start=i_startidx, &
        &   i_end=i_endidx)
    IF( ALLOCATED(dwdy_before) ) CALL verify_field_wp( &
        &   field_name='dwdy', &
        &   field_ref=p_nh_diag%dwdy(:, :, 1), &
        &   field=dwdy_before(:, :, 1), &
        &   i_start=i_startidx, &
        &   i_end=i_endidx)

    CALL print_verify_footer()

    CALL print_verify_header('padding')

    rl_start = 1
    rl_end   = min_rledge

    i_startblk = p_patch%edges%start_block(rl_start)
    i_endblk   = p_patch%edges%end_block(rl_end)

    CALL get_indices_e(p_patch, 1, i_startblk, i_endblk, i_startidx, i_endidx, rl_start, rl_end)

    CALL verify_field_wp( &
        &   field_name='vn', &
        &   field_ref=p_nh_prog%vn(:, :, 1), &
        &   field=vn_before(:, :, 1), &
        &   i_start=i_endidx+1, &
        &   i_end=nproma)

    rl_start = 1
    rl_end   = min_rlcell

    i_startblk = p_patch%cells%start_block(rl_start)
    i_endblk   = p_patch%cells%end_block(rl_end)

    CALL get_indices_c(p_patch, 1, i_startblk, i_endblk, i_startidx, i_endidx, rl_start, rl_end)

    CALL verify_field_wp( &
        &   field_name='w', &
        &   field_ref=p_nh_prog%w(:, :, 1), &
        &   field=w_before(:, :, 1), &
        &   i_start=i_endidx+1, &
        &   i_end=nproma)
    CALL verify_field_wp( &
        &   field_name='exner', &
        &   field_ref=p_nh_prog%exner(:, :, 1), &
        &   field=exner_before(:, :, 1), &
        &   i_start=i_endidx+1, &
        &   i_end=nproma)
    CALL verify_field_wp( &
        &   field_name='theta_v', &
        &   field_ref=p_nh_prog%theta_v(:, :, 1), &
        &   field=theta_v_before(:, :, 1), &
        &   i_start=i_endidx+1, &
        &   i_end=nproma)
    CALL verify_field_wp(field_name='rho', &
        &   field_ref=p_nh_prog%rho(:, :, 1), &
        &   field=rho_before(:, :, 1), &
        &   i_start=i_endidx+1, &
        &   i_end=nproma)
    IF( ALLOCATED(hdef_ic_before) ) CALL verify_field_wp(&
        &   field_name='hdef_ic', &
        &   field_ref=p_nh_diag%hdef_ic(:, :, 1), &
        &   field=hdef_ic_before(:, :, 1), &
        &   i_start=i_endidx+1, &
        &   i_end=nproma)
    IF( ALLOCATED(div_ic_before) ) CALL verify_field_wp( &
        &   field_name='div_ic', &
        &   field_ref=p_nh_diag%div_ic(:, :, 1), &
        &   field=div_ic_before(:, :, 1), &
        &   i_start=i_endidx+1, &
        &   i_end=nproma)
    IF( ALLOCATED(dwdx_before) ) CALL verify_field_wp( &
        &   field_name='dwdx', &
        &   field_ref=p_nh_diag%dwdx(:, :, 1), &
        &   field=dwdx_before(:, :, 1), &
        &   i_start=i_endidx+1, &
        &   i_end=nproma)
    IF( ALLOCATED(dwdy_before) ) CALL verify_field_wp( &
        &   field_name='dwdy', &
        &   field_ref=p_nh_diag%dwdy(:, :, 1), &
        &   field=dwdy_before(:, :, 1), &
        &   i_start=i_endidx+1, &
        &   i_end=nproma)

    CALL print_verify_footer()

    !$ACC END DATA
    !$ACC END DATA
    !$ACC END DATA
    !$ACC END DATA
    !$ACC END DATA

    !$ACC END DATA
    !$ACC END DATA
    !$ACC END DATA
    !$ACC END DATA
    !$ACC END DATA

  END SUBROUTINE icon4py_diffusion_verify

  SUBROUTINE icon4py_check_rc(rc, routine, bindings_routine)
    INTEGER,          INTENT(IN) :: rc
    CHARACTER(LEN=*), INTENT(IN) :: routine, bindings_routine
    CHARACTER(LEN=200)           :: msg

    ! rc convention from icon4py's py2fgen codegen:
    !   0 -> cffi embedding init failed (cffi forces 0 here, see
    !        _cffi_start_and_call_python in _embedding.h of cffi).
    !   1 -> success.
    !   2 -> Python wrapper raised an exception (traceback logged Python-side).
    SELECT CASE (rc)
    CASE (1)
      ! ok
    CASE (0)
      CALL finish(routine, TRIM(bindings_routine)// &
        ' failed: cffi embedding init failed (Python environment with ICON4Py active?). '// &
        'See stderr for "Failed to initialize the Python-CFFI embedding logic".')
    CASE (2)
      CALL finish(routine, TRIM(bindings_routine)// &
        ' failed: Python wrapper raised an exception. See log for traceback.')
    CASE DEFAULT
      WRITE(msg,'(a,i0)') TRIM(bindings_routine)//' failed: unexpected rc = ', rc
      CALL finish(routine, msg)
    END SELECT
  END SUBROUTINE icon4py_check_rc
#endif

END MODULE mo_icon4py_interfaces
