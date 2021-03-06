
! This file is automatically created by parse_castro_params.py.  To update
! or add runtime parameters, please edit _cpp_parameters and then run
! mk_params.sh

! This module stores the runtime parameters and integer names for 
! indexing arrays.
!
! The Fortran-specific parameters are initialized in set_method_params(),
! and the ones that we are mirroring from C++ and obtaining through the
! ParmParse module are initialized in ca_set_castro_method_params().

module meth_params_module

  use amrex_error_module
  use amrex_fort_module, only: rt => amrex_real
  use state_sizes_module, only : nadv, NQAUX, NVAR, NGDNV, NQ, QVAR

  implicit none

  ! number of ghost cells for the hyperbolic solver
  integer, parameter     :: NHYP    = 4

  ! conservative variables
  integer, allocatable, save :: URHO, UMX, UMY, UMZ, UMR, UML, UMP, UEDEN, UEINT, UTEMP, UFA, UFS, UFX
  integer, allocatable, save :: USHK

  ! primitive variables
  integer, allocatable, save :: QRHO, QU, QV, QW, QPRES, QREINT, QTEMP, QGAME
  integer, allocatable, save :: QGAMC, QC, QDPDR, QDPDE
#ifdef RADIATION
  integer, allocatable, save :: QGAMCG, QCG, QLAMS
#endif
  integer, allocatable, save :: QFA, QFS, QFX

#ifdef RADIATION
  integer, save :: QRAD, QRADHI, QPTOT, QREITOT
  integer, save :: fspace_type
  logical, save :: do_inelastic_scattering
  logical, save :: comoving

  real(rt)        , save :: flatten_pp_threshold = -1.e0_rt
#endif

  integer, save, allocatable :: npassive
  integer, save, allocatable :: qpass_map(:), upass_map(:)

  ! These are used for the Godunov state
  ! Note that the velocity indices here are picked to be the same value
  ! as in the primitive variable array
  integer, save, allocatable :: GDRHO, GDU, GDV, GDW, GDPRES, GDGAME
#ifdef RADIATION
  integer, save, allocatable :: GDLAMS, GDERADS
#endif

  integer         , save :: numpts_1d

  real(rt)        , save, allocatable :: outflow_data_old(:,:)
  real(rt)        , save, allocatable :: outflow_data_new(:,:)
  real(rt)        , save :: outflow_data_old_time
  real(rt)        , save :: outflow_data_new_time
  logical         , save :: outflow_data_allocated
  real(rt)        , save :: max_dist

  character(len=:), allocatable :: gravity_type

  ! these flags are for interpreting the EXT_DIR BCs
  integer, parameter :: EXT_UNDEFINED = -1
  integer, parameter :: EXT_HSE = 1
  integer, parameter :: EXT_INTERP = 2 
  
  integer, allocatable, save :: xl_ext, yl_ext, zl_ext, xr_ext, yr_ext, zr_ext

  ! Create versions of these variables on the GPU
  ! the device update is then done in Castro_nd.f90

#ifdef AMREX_USE_CUDA
  attributes(managed) :: URHO, UMX, UMY, UMZ, UMR, UML, UMP, UEDEN, UEINT, UTEMP, UFA, UFS, UFX
  attributes(managed) :: USHK
  attributes(managed) :: QRHO, QU, QV, QW, QPRES, QREINT, QTEMP, QGAME
  attributes(managed) :: QGAMC, QC, QDPDR, QDPDE
#ifdef RADIATION
  attributes(managed) :: QGAMCG, QCG, QLAMS
#endif
  attributes(managed) :: QFA, QFS, QFX
  attributes(managed) :: npassive
  attributes(managed) :: qpass_map, upass_map
  attributes(managed) :: GDRHO, GDU, GDV, GDW, GDPRES, GDGAME
#ifdef RADIATION
  attributes(managed) :: GDLAMS, GDERADS
#endif
  attributes(managed) :: xl_ext, yl_ext, zl_ext, xr_ext, yr_ext, zr_ext
#endif

  !$acc declare &
  !$acc create(URHO, UMX, UMY, UMZ, UMR, UML, UMP, UEDEN, UEINT, UTEMP, UFA, UFS,UFX) &
  !$acc create(USHK) &
  !$acc create(QRHO, QU, QV, QW, QPRES, QREINT, QTEMP) &
  !$acc create(QC, QDPDR, QDPDE, QGAMC, QGAME) &
#ifdef RADIATION
  !$acc create(QGAMCG, QCG, QLAMS) &
  !$acc create(QRAD, QRADHI, QPTOT, QREITOT) &
  !$acc create(fspace_type, do_inelastic_scattering, comoving) &
#endif
  !$acc create(QFA, QFS, QFX) &
  !$acc create(xl_ext, yl_ext, zl_ext, xr_ext, yr_ext, zr_ext)

  ! Begin the declarations of the ParmParse parameters

  real(rt), allocatable, save :: difmag
  real(rt), allocatable, save :: small_dens
  real(rt), allocatable, save :: small_temp
  real(rt), allocatable, save :: small_pres
  real(rt), allocatable, save :: small_ener
  integer,  allocatable, save :: do_hydro
  integer,  allocatable, save :: do_ctu
  integer,  allocatable, save :: fourth_order
  integer,  allocatable, save :: hybrid_hydro
  integer,  allocatable, save :: ppm_type
  integer,  allocatable, save :: ppm_temp_fix
  integer,  allocatable, save :: ppm_predict_gammae
  integer,  allocatable, save :: ppm_reference_eigenvectors
  integer,  allocatable, save :: plm_iorder
  integer,  allocatable, save :: hybrid_riemann
  integer,  allocatable, save :: riemann_solver
  integer,  allocatable, save :: cg_maxiter
  real(rt), allocatable, save :: cg_tol
  integer,  allocatable, save :: cg_blend
  integer,  allocatable, save :: use_eos_in_riemann
  integer,  allocatable, save :: use_flattening
  integer,  allocatable, save :: transverse_use_eos
  integer,  allocatable, save :: transverse_reset_density
  integer,  allocatable, save :: transverse_reset_rhoe
  real(rt), allocatable, save :: dual_energy_eta1
  real(rt), allocatable, save :: dual_energy_eta2
  integer,  allocatable, save :: use_pslope
  integer,  allocatable, save :: fix_mass_flux
  integer,  allocatable, save :: limit_fluxes_on_small_dens
  integer,  allocatable, save :: density_reset_method
  integer,  allocatable, save :: allow_small_energy
  integer,  allocatable, save :: do_sponge
  integer,  allocatable, save :: sponge_implicit
  integer,  allocatable, save :: first_order_hydro
  character (len=:), allocatable, save :: xl_ext_bc_type
  character (len=:), allocatable, save :: xr_ext_bc_type
  character (len=:), allocatable, save :: yl_ext_bc_type
  character (len=:), allocatable, save :: yr_ext_bc_type
  character (len=:), allocatable, save :: zl_ext_bc_type
  character (len=:), allocatable, save :: zr_ext_bc_type
  integer,  allocatable, save :: hse_zero_vels
  integer,  allocatable, save :: hse_interp_temp
  integer,  allocatable, save :: hse_reflect_vels
  integer,  allocatable, save :: mol_order
  real(rt), allocatable, save :: cfl
  real(rt), allocatable, save :: dtnuc_e
  real(rt), allocatable, save :: dtnuc_X
  real(rt), allocatable, save :: dtnuc_X_threshold
  real(rt), allocatable, save :: dxnuc
  real(rt), allocatable, save :: dxnuc_max
  integer,  allocatable, save :: max_dxnuc_lev
  integer,  allocatable, save :: do_react
  real(rt), allocatable, save :: react_T_min
  real(rt), allocatable, save :: react_T_max
  real(rt), allocatable, save :: react_rho_min
  real(rt), allocatable, save :: react_rho_max
  integer,  allocatable, save :: disable_shock_burning
  real(rt), allocatable, save :: diffuse_cutoff_density
  real(rt), allocatable, save :: diffuse_cond_scale_fac
  integer,  allocatable, save :: do_grav
  integer,  allocatable, save :: grav_source_type
  integer,  allocatable, save :: do_rotation
  real(rt), allocatable, save :: rot_period
  real(rt), allocatable, save :: rot_period_dot
  integer,  allocatable, save :: rotation_include_centrifugal
  integer,  allocatable, save :: rotation_include_coriolis
  integer,  allocatable, save :: rotation_include_domegadt
  integer,  allocatable, save :: state_in_rotating_frame
  integer,  allocatable, save :: rot_source_type
  integer,  allocatable, save :: implicit_rotation_update
  integer,  allocatable, save :: rot_axis
  integer,  allocatable, save :: use_point_mass
  real(rt), allocatable, save :: point_mass
  integer,  allocatable, save :: point_mass_fix_solution
  integer,  allocatable, save :: do_acc
  integer,  allocatable, save :: grown_factor
  integer,  allocatable, save :: track_grid_losses
  real(rt), allocatable, save :: const_grav
  integer,  allocatable, save :: get_g_from_phi

#ifdef AMREX_USE_CUDA
attributes(managed) :: difmag
attributes(managed) :: small_dens
attributes(managed) :: small_temp
attributes(managed) :: small_pres
attributes(managed) :: small_ener
attributes(managed) :: do_hydro
attributes(managed) :: do_ctu
attributes(managed) :: fourth_order
attributes(managed) :: hybrid_hydro
attributes(managed) :: ppm_type
attributes(managed) :: ppm_temp_fix
attributes(managed) :: ppm_predict_gammae
attributes(managed) :: ppm_reference_eigenvectors
attributes(managed) :: plm_iorder
attributes(managed) :: hybrid_riemann
attributes(managed) :: riemann_solver
attributes(managed) :: cg_maxiter
attributes(managed) :: cg_tol
attributes(managed) :: cg_blend
attributes(managed) :: use_eos_in_riemann
attributes(managed) :: use_flattening
attributes(managed) :: transverse_use_eos
attributes(managed) :: transverse_reset_density
attributes(managed) :: transverse_reset_rhoe
attributes(managed) :: dual_energy_eta1
attributes(managed) :: dual_energy_eta2
attributes(managed) :: use_pslope
attributes(managed) :: fix_mass_flux
attributes(managed) :: limit_fluxes_on_small_dens
attributes(managed) :: density_reset_method
attributes(managed) :: allow_small_energy
attributes(managed) :: do_sponge
attributes(managed) :: sponge_implicit
attributes(managed) :: first_order_hydro






attributes(managed) :: hse_zero_vels
attributes(managed) :: hse_interp_temp
attributes(managed) :: hse_reflect_vels
attributes(managed) :: mol_order
attributes(managed) :: cfl
attributes(managed) :: dtnuc_e
attributes(managed) :: dtnuc_X
attributes(managed) :: dtnuc_X_threshold
attributes(managed) :: dxnuc
attributes(managed) :: dxnuc_max
attributes(managed) :: max_dxnuc_lev
attributes(managed) :: do_react
attributes(managed) :: react_T_min
attributes(managed) :: react_T_max
attributes(managed) :: react_rho_min
attributes(managed) :: react_rho_max
attributes(managed) :: disable_shock_burning
#ifdef DIFFUSION
attributes(managed) :: diffuse_cutoff_density
#endif
#ifdef DIFFUSION
attributes(managed) :: diffuse_cond_scale_fac
#endif
attributes(managed) :: do_grav
attributes(managed) :: grav_source_type
attributes(managed) :: do_rotation
#ifdef ROTATION
attributes(managed) :: rot_period
#endif
#ifdef ROTATION
attributes(managed) :: rot_period_dot
#endif
#ifdef ROTATION
attributes(managed) :: rotation_include_centrifugal
#endif
#ifdef ROTATION
attributes(managed) :: rotation_include_coriolis
#endif
#ifdef ROTATION
attributes(managed) :: rotation_include_domegadt
#endif
#ifdef ROTATION
attributes(managed) :: state_in_rotating_frame
#endif
#ifdef ROTATION
attributes(managed) :: rot_source_type
#endif
#ifdef ROTATION
attributes(managed) :: implicit_rotation_update
#endif
#ifdef ROTATION
attributes(managed) :: rot_axis
#endif
#ifdef POINTMASS
attributes(managed) :: use_point_mass
#endif
#ifdef POINTMASS
attributes(managed) :: point_mass
#endif
#ifdef POINTMASS
attributes(managed) :: point_mass_fix_solution
#endif
attributes(managed) :: do_acc
attributes(managed) :: grown_factor
attributes(managed) :: track_grid_losses
attributes(managed) :: const_grav
attributes(managed) :: get_g_from_phi
#endif

  !$acc declare &
  !$acc create(difmag) &
  !$acc create(small_dens) &
  !$acc create(small_temp) &
  !$acc create(small_pres) &
  !$acc create(small_ener) &
  !$acc create(do_hydro) &
  !$acc create(do_ctu) &
  !$acc create(fourth_order) &
  !$acc create(hybrid_hydro) &
  !$acc create(ppm_type) &
  !$acc create(ppm_temp_fix) &
  !$acc create(ppm_predict_gammae) &
  !$acc create(ppm_reference_eigenvectors) &
  !$acc create(plm_iorder) &
  !$acc create(hybrid_riemann) &
  !$acc create(riemann_solver) &
  !$acc create(cg_maxiter) &
  !$acc create(cg_tol) &
  !$acc create(cg_blend) &
  !$acc create(use_eos_in_riemann) &
  !$acc create(use_flattening) &
  !$acc create(transverse_use_eos) &
  !$acc create(transverse_reset_density) &
  !$acc create(transverse_reset_rhoe) &
  !$acc create(dual_energy_eta1) &
  !$acc create(dual_energy_eta2) &
  !$acc create(use_pslope) &
  !$acc create(fix_mass_flux) &
  !$acc create(limit_fluxes_on_small_dens) &
  !$acc create(density_reset_method) &
  !$acc create(allow_small_energy) &
  !$acc create(do_sponge) &
  !$acc create(sponge_implicit) &
  !$acc create(first_order_hydro) &
  !$acc create(hse_zero_vels) &
  !$acc create(hse_interp_temp) &
  !$acc create(hse_reflect_vels) &
  !$acc create(mol_order) &
  !$acc create(cfl) &
  !$acc create(dtnuc_e) &
  !$acc create(dtnuc_X) &
  !$acc create(dtnuc_X_threshold) &
  !$acc create(dxnuc) &
  !$acc create(dxnuc_max) &
  !$acc create(max_dxnuc_lev) &
  !$acc create(do_react) &
  !$acc create(react_T_min) &
  !$acc create(react_T_max) &
  !$acc create(react_rho_min) &
  !$acc create(react_rho_max) &
  !$acc create(disable_shock_burning) &
#ifdef DIFFUSION
  !$acc create(diffuse_cutoff_density) &
#endif
#ifdef DIFFUSION
  !$acc create(diffuse_cond_scale_fac) &
#endif
  !$acc create(do_grav) &
  !$acc create(grav_source_type) &
  !$acc create(do_rotation) &
#ifdef ROTATION
  !$acc create(rot_period) &
#endif
#ifdef ROTATION
  !$acc create(rot_period_dot) &
#endif
#ifdef ROTATION
  !$acc create(rotation_include_centrifugal) &
#endif
#ifdef ROTATION
  !$acc create(rotation_include_coriolis) &
#endif
#ifdef ROTATION
  !$acc create(rotation_include_domegadt) &
#endif
#ifdef ROTATION
  !$acc create(state_in_rotating_frame) &
#endif
#ifdef ROTATION
  !$acc create(rot_source_type) &
#endif
#ifdef ROTATION
  !$acc create(implicit_rotation_update) &
#endif
#ifdef ROTATION
  !$acc create(rot_axis) &
#endif
#ifdef POINTMASS
  !$acc create(use_point_mass) &
#endif
#ifdef POINTMASS
  !$acc create(point_mass) &
#endif
#ifdef POINTMASS
  !$acc create(point_mass_fix_solution) &
#endif
  !$acc create(do_acc) &
  !$acc create(grown_factor) &
  !$acc create(track_grid_losses) &
  !$acc create(const_grav) &
  !$acc create(get_g_from_phi)

  ! End the declarations of the ParmParse parameters

  real(rt)        , save :: rot_vec(3)

contains

  subroutine ca_set_castro_method_params() bind(C, name="ca_set_castro_method_params")

    use amrex_parmparse_module, only: amrex_parmparse_build, amrex_parmparse_destroy, amrex_parmparse

    use amrex_fort_module, only : rt => amrex_real
    implicit none

    type (amrex_parmparse) :: pp


    allocate(URHO, UMX, UMY, UMZ, UMR, UML, UMP, UEDEN, UEINT, UTEMP, UFA, UFS, UFX)
    allocate(USHK)
    allocate(QRHO, QU, QV, QW, QPRES, QREINT, QTEMP, QGAME)
    allocate(QGAMC, QC, QDPDR, QDPDE)
#ifdef RADIATION
    allocate(QGAMCG, QCG, QLAMS)
#endif
    allocate(QFA, QFS, QFX)
    allocate(npassive)
    allocate(GDRHO, GDU, GDV, GDW, GDPRES, GDGAME)
#ifdef RADIATION
    allocate(GDLAMS, GDERADS)
#endif
    allocate(xl_ext, yl_ext, zl_ext, xr_ext, yr_ext, zr_ext)

    allocate(const_grav)
    const_grav = 0.0d0;
    allocate(get_g_from_phi)
    get_g_from_phi = 0;

    call amrex_parmparse_build(pp, "gravity")
    call pp%query("const_grav", const_grav)
    call pp%query("get_g_from_phi", get_g_from_phi)
    call amrex_parmparse_destroy(pp)


#ifdef DIFFUSION
    allocate(diffuse_cutoff_density)
    diffuse_cutoff_density = -1.d200;
    allocate(diffuse_cond_scale_fac)
    diffuse_cond_scale_fac = 1.0d0;
#endif
#ifdef ROTATION
    allocate(rot_period)
    rot_period = -1.d200;
    allocate(rot_period_dot)
    rot_period_dot = 0.0d0;
    allocate(rotation_include_centrifugal)
    rotation_include_centrifugal = 1;
    allocate(rotation_include_coriolis)
    rotation_include_coriolis = 1;
    allocate(rotation_include_domegadt)
    rotation_include_domegadt = 1;
    allocate(state_in_rotating_frame)
    state_in_rotating_frame = 1;
    allocate(rot_source_type)
    rot_source_type = 4;
    allocate(implicit_rotation_update)
    implicit_rotation_update = 1;
    allocate(rot_axis)
    rot_axis = 3;
#endif
#ifdef POINTMASS
    allocate(use_point_mass)
    use_point_mass = 1;
    allocate(point_mass)
    point_mass = 0.0d0;
    allocate(point_mass_fix_solution)
    point_mass_fix_solution = 0;
#endif
    allocate(difmag)
    difmag = 0.1d0;
    allocate(small_dens)
    small_dens = -1.d200;
    allocate(small_temp)
    small_temp = -1.d200;
    allocate(small_pres)
    small_pres = -1.d200;
    allocate(small_ener)
    small_ener = -1.d200;
    allocate(do_hydro)
    do_hydro = -1;
    allocate(do_ctu)
    do_ctu = 1;
    allocate(fourth_order)
    fourth_order = 0;
    allocate(hybrid_hydro)
    hybrid_hydro = 0;
    allocate(ppm_type)
    ppm_type = 1;
    allocate(ppm_temp_fix)
    ppm_temp_fix = 0;
    allocate(ppm_predict_gammae)
    ppm_predict_gammae = 0;
    allocate(ppm_reference_eigenvectors)
    ppm_reference_eigenvectors = 0;
    allocate(plm_iorder)
    plm_iorder = 2;
    allocate(hybrid_riemann)
    hybrid_riemann = 0;
    allocate(riemann_solver)
    riemann_solver = 0;
    allocate(cg_maxiter)
    cg_maxiter = 12;
    allocate(cg_tol)
    cg_tol = 1.0d-5;
    allocate(cg_blend)
    cg_blend = 2;
    allocate(use_eos_in_riemann)
    use_eos_in_riemann = 0;
    allocate(use_flattening)
    use_flattening = 1;
    allocate(transverse_use_eos)
    transverse_use_eos = 0;
    allocate(transverse_reset_density)
    transverse_reset_density = 1;
    allocate(transverse_reset_rhoe)
    transverse_reset_rhoe = 0;
    allocate(dual_energy_eta1)
    dual_energy_eta1 = 1.0d0;
    allocate(dual_energy_eta2)
    dual_energy_eta2 = 1.0d-4;
    allocate(use_pslope)
    use_pslope = 1;
    allocate(fix_mass_flux)
    fix_mass_flux = 0;
    allocate(limit_fluxes_on_small_dens)
    limit_fluxes_on_small_dens = 0;
    allocate(density_reset_method)
    density_reset_method = 1;
    allocate(allow_small_energy)
    allow_small_energy = 1;
    allocate(do_sponge)
    do_sponge = 0;
    allocate(sponge_implicit)
    sponge_implicit = 1;
    allocate(first_order_hydro)
    first_order_hydro = 0;
    allocate(character(len=1)::xl_ext_bc_type)
    xl_ext_bc_type = "";
    allocate(character(len=1)::xr_ext_bc_type)
    xr_ext_bc_type = "";
    allocate(character(len=1)::yl_ext_bc_type)
    yl_ext_bc_type = "";
    allocate(character(len=1)::yr_ext_bc_type)
    yr_ext_bc_type = "";
    allocate(character(len=1)::zl_ext_bc_type)
    zl_ext_bc_type = "";
    allocate(character(len=1)::zr_ext_bc_type)
    zr_ext_bc_type = "";
    allocate(hse_zero_vels)
    hse_zero_vels = 0;
    allocate(hse_interp_temp)
    hse_interp_temp = 0;
    allocate(hse_reflect_vels)
    hse_reflect_vels = 0;
    allocate(mol_order)
    mol_order = 2;
    allocate(cfl)
    cfl = 0.8d0;
    allocate(dtnuc_e)
    dtnuc_e = 1.d200;
    allocate(dtnuc_X)
    dtnuc_X = 1.d200;
    allocate(dtnuc_X_threshold)
    dtnuc_X_threshold = 1.d-3;
    allocate(dxnuc)
    dxnuc = 1.d200;
    allocate(dxnuc_max)
    dxnuc_max = 1.d200;
    allocate(max_dxnuc_lev)
    max_dxnuc_lev = -1;
    allocate(do_react)
    do_react = -1;
    allocate(react_T_min)
    react_T_min = 0.0d0;
    allocate(react_T_max)
    react_T_max = 1.d200;
    allocate(react_rho_min)
    react_rho_min = 0.0d0;
    allocate(react_rho_max)
    react_rho_max = 1.d200;
    allocate(disable_shock_burning)
    disable_shock_burning = 0;
    allocate(do_grav)
    do_grav = -1;
    allocate(grav_source_type)
    grav_source_type = 4;
    allocate(do_rotation)
    do_rotation = -1;
    allocate(do_acc)
    do_acc = -1;
    allocate(grown_factor)
    grown_factor = 1;
    allocate(track_grid_losses)
    track_grid_losses = 0;

    call amrex_parmparse_build(pp, "castro")
#ifdef DIFFUSION
    call pp%query("diffuse_cutoff_density", diffuse_cutoff_density)
    call pp%query("diffuse_cond_scale_fac", diffuse_cond_scale_fac)
#endif
#ifdef ROTATION
    call pp%query("rotational_period", rot_period)
    call pp%query("rotational_dPdt", rot_period_dot)
    call pp%query("rotation_include_centrifugal", rotation_include_centrifugal)
    call pp%query("rotation_include_coriolis", rotation_include_coriolis)
    call pp%query("rotation_include_domegadt", rotation_include_domegadt)
    call pp%query("state_in_rotating_frame", state_in_rotating_frame)
    call pp%query("rot_source_type", rot_source_type)
    call pp%query("implicit_rotation_update", implicit_rotation_update)
    call pp%query("rot_axis", rot_axis)
#endif
#ifdef POINTMASS
    call pp%query("use_point_mass", use_point_mass)
    call pp%query("point_mass", point_mass)
    call pp%query("point_mass_fix_solution", point_mass_fix_solution)
#endif
    call pp%query("difmag", difmag)
    call pp%query("small_dens", small_dens)
    call pp%query("small_temp", small_temp)
    call pp%query("small_pres", small_pres)
    call pp%query("small_ener", small_ener)
    call pp%query("do_hydro", do_hydro)
    call pp%query("do_ctu", do_ctu)
    call pp%query("fourth_order", fourth_order)
    call pp%query("hybrid_hydro", hybrid_hydro)
    call pp%query("ppm_type", ppm_type)
    call pp%query("ppm_temp_fix", ppm_temp_fix)
    call pp%query("ppm_predict_gammae", ppm_predict_gammae)
    call pp%query("ppm_reference_eigenvectors", ppm_reference_eigenvectors)
    call pp%query("plm_iorder", plm_iorder)
    call pp%query("hybrid_riemann", hybrid_riemann)
    call pp%query("riemann_solver", riemann_solver)
    call pp%query("cg_maxiter", cg_maxiter)
    call pp%query("cg_tol", cg_tol)
    call pp%query("cg_blend", cg_blend)
    call pp%query("use_eos_in_riemann", use_eos_in_riemann)
    call pp%query("use_flattening", use_flattening)
    call pp%query("transverse_use_eos", transverse_use_eos)
    call pp%query("transverse_reset_density", transverse_reset_density)
    call pp%query("transverse_reset_rhoe", transverse_reset_rhoe)
    call pp%query("dual_energy_eta1", dual_energy_eta1)
    call pp%query("dual_energy_eta2", dual_energy_eta2)
    call pp%query("use_pslope", use_pslope)
    call pp%query("fix_mass_flux", fix_mass_flux)
    call pp%query("limit_fluxes_on_small_dens", limit_fluxes_on_small_dens)
    call pp%query("density_reset_method", density_reset_method)
    call pp%query("allow_small_energy", allow_small_energy)
    call pp%query("do_sponge", do_sponge)
    call pp%query("sponge_implicit", sponge_implicit)
    call pp%query("first_order_hydro", first_order_hydro)
    call pp%query("xl_ext_bc_type", xl_ext_bc_type)
    call pp%query("xr_ext_bc_type", xr_ext_bc_type)
    call pp%query("yl_ext_bc_type", yl_ext_bc_type)
    call pp%query("yr_ext_bc_type", yr_ext_bc_type)
    call pp%query("zl_ext_bc_type", zl_ext_bc_type)
    call pp%query("zr_ext_bc_type", zr_ext_bc_type)
    call pp%query("hse_zero_vels", hse_zero_vels)
    call pp%query("hse_interp_temp", hse_interp_temp)
    call pp%query("hse_reflect_vels", hse_reflect_vels)
    call pp%query("mol_order", mol_order)
    call pp%query("cfl", cfl)
    call pp%query("dtnuc_e", dtnuc_e)
    call pp%query("dtnuc_X", dtnuc_X)
    call pp%query("dtnuc_X_threshold", dtnuc_X_threshold)
    call pp%query("dxnuc", dxnuc)
    call pp%query("dxnuc_max", dxnuc_max)
    call pp%query("max_dxnuc_lev", max_dxnuc_lev)
    call pp%query("do_react", do_react)
    call pp%query("react_T_min", react_T_min)
    call pp%query("react_T_max", react_T_max)
    call pp%query("react_rho_min", react_rho_min)
    call pp%query("react_rho_max", react_rho_max)
    call pp%query("disable_shock_burning", disable_shock_burning)
    call pp%query("do_grav", do_grav)
    call pp%query("grav_source_type", grav_source_type)
    call pp%query("do_rotation", do_rotation)
    call pp%query("do_acc", do_acc)
    call pp%query("grown_factor", grown_factor)
    call pp%query("track_grid_losses", track_grid_losses)
    call amrex_parmparse_destroy(pp)



    !$acc update &
    !$acc device(difmag, small_dens, small_temp) &
    !$acc device(small_pres, small_ener, do_hydro) &
    !$acc device(do_ctu, fourth_order, hybrid_hydro) &
    !$acc device(ppm_type, ppm_temp_fix, ppm_predict_gammae) &
    !$acc device(ppm_reference_eigenvectors, plm_iorder, hybrid_riemann) &
    !$acc device(riemann_solver, cg_maxiter, cg_tol) &
    !$acc device(cg_blend, use_eos_in_riemann, use_flattening) &
    !$acc device(transverse_use_eos, transverse_reset_density, transverse_reset_rhoe) &
    !$acc device(dual_energy_eta1, dual_energy_eta2, use_pslope) &
    !$acc device(fix_mass_flux, limit_fluxes_on_small_dens, density_reset_method) &
    !$acc device(allow_small_energy, do_sponge, sponge_implicit) &
    !$acc device(first_order_hydro, hse_zero_vels, hse_interp_temp) &
    !$acc device(hse_reflect_vels, mol_order, cfl) &
    !$acc device(dtnuc_e, dtnuc_X, dtnuc_X_threshold) &
    !$acc device(dxnuc, dxnuc_max, max_dxnuc_lev) &
    !$acc device(do_react, react_T_min, react_T_max) &
    !$acc device(react_rho_min, react_rho_max, disable_shock_burning) &
    !$acc device(diffuse_cutoff_density, diffuse_cond_scale_fac, do_grav) &
    !$acc device(grav_source_type, do_rotation, rot_period) &
    !$acc device(rot_period_dot, rotation_include_centrifugal, rotation_include_coriolis) &
    !$acc device(rotation_include_domegadt, state_in_rotating_frame, rot_source_type) &
    !$acc device(implicit_rotation_update, rot_axis, use_point_mass) &
    !$acc device(point_mass, point_mass_fix_solution, do_acc) &
    !$acc device(grown_factor, track_grid_losses, const_grav) &
    !$acc device(get_g_from_phi)


    ! now set the external BC flags
    select case (xl_ext_bc_type)
    case ("hse", "HSE")
       xl_ext = EXT_HSE
    case ("interp", "INTERP")       
       xl_ext = EXT_INTERP
    case default
       xl_ext = EXT_UNDEFINED
    end select

    select case (yl_ext_bc_type)
    case ("hse", "HSE")
       yl_ext = EXT_HSE
    case ("interp", "INTERP")       
       yl_ext = EXT_INTERP
    case default
       yl_ext = EXT_UNDEFINED
    end select

    select case (zl_ext_bc_type)
    case ("hse", "HSE")
       zl_ext = EXT_HSE
    case ("interp", "INTERP")       
       zl_ext = EXT_INTERP
    case default
       zl_ext = EXT_UNDEFINED
    end select

    select case (xr_ext_bc_type)
    case ("hse", "HSE")
       xr_ext = EXT_HSE
    case ("interp", "INTERP")       
       xr_ext = EXT_INTERP
    case default
       xr_ext = EXT_UNDEFINED
    end select

    select case (yr_ext_bc_type)
    case ("hse", "HSE")
       yr_ext = EXT_HSE
    case ("interp", "INTERP")       
       yr_ext = EXT_INTERP
    case default
       yr_ext = EXT_UNDEFINED
    end select

    select case (zr_ext_bc_type)
    case ("hse", "HSE")
       zr_ext = EXT_HSE
    case ("interp", "INTERP")       
       zr_ext = EXT_INTERP
    case default
       zr_ext = EXT_UNDEFINED
    end select

    !$acc update device(xl_ext, yl_ext, zl_ext, xr_ext, yr_ext, zr_ext)


  end subroutine ca_set_castro_method_params


  subroutine ca_finalize_meth_params() bind(C, name="ca_finalize_meth_params")
    implicit none

    deallocate(URHO, UMX, UMY, UMZ, UMR, UML, UMP, UEDEN, UEINT, UTEMP, UFA, UFS, UFX)
    deallocate(USHK)
    deallocate(QRHO, QU, QV, QW, QPRES, QREINT, QTEMP, QGAME)
    deallocate(QGAMC, QC, QDPDR, QDPDE)
#ifdef RADIATION
    deallocate(QGAMCG, QCG, QLAMS)
#endif
    deallocate(QFA, QFS, QFX)
    deallocate(npassive)
    deallocate(GDRHO, GDU, GDV, GDW, GDPRES, GDGAME)
#ifdef RADIATION
    deallocate(GDLAMS, GDERADS)
    deallocate(xl_ext, yl_ext, zl_ext, xr_ext, yr_ext, zr_ext)
#endif

    if (allocated(difmag)) then
        deallocate(difmag)
    end if
    if (allocated(small_dens)) then
        deallocate(small_dens)
    end if
    if (allocated(small_temp)) then
        deallocate(small_temp)
    end if
    if (allocated(small_pres)) then
        deallocate(small_pres)
    end if
    if (allocated(small_ener)) then
        deallocate(small_ener)
    end if
    if (allocated(do_hydro)) then
        deallocate(do_hydro)
    end if
    if (allocated(do_ctu)) then
        deallocate(do_ctu)
    end if
    if (allocated(fourth_order)) then
        deallocate(fourth_order)
    end if
    if (allocated(hybrid_hydro)) then
        deallocate(hybrid_hydro)
    end if
    if (allocated(ppm_type)) then
        deallocate(ppm_type)
    end if
    if (allocated(ppm_temp_fix)) then
        deallocate(ppm_temp_fix)
    end if
    if (allocated(ppm_predict_gammae)) then
        deallocate(ppm_predict_gammae)
    end if
    if (allocated(ppm_reference_eigenvectors)) then
        deallocate(ppm_reference_eigenvectors)
    end if
    if (allocated(plm_iorder)) then
        deallocate(plm_iorder)
    end if
    if (allocated(hybrid_riemann)) then
        deallocate(hybrid_riemann)
    end if
    if (allocated(riemann_solver)) then
        deallocate(riemann_solver)
    end if
    if (allocated(cg_maxiter)) then
        deallocate(cg_maxiter)
    end if
    if (allocated(cg_tol)) then
        deallocate(cg_tol)
    end if
    if (allocated(cg_blend)) then
        deallocate(cg_blend)
    end if
    if (allocated(use_eos_in_riemann)) then
        deallocate(use_eos_in_riemann)
    end if
    if (allocated(use_flattening)) then
        deallocate(use_flattening)
    end if
    if (allocated(transverse_use_eos)) then
        deallocate(transverse_use_eos)
    end if
    if (allocated(transverse_reset_density)) then
        deallocate(transverse_reset_density)
    end if
    if (allocated(transverse_reset_rhoe)) then
        deallocate(transverse_reset_rhoe)
    end if
    if (allocated(dual_energy_eta1)) then
        deallocate(dual_energy_eta1)
    end if
    if (allocated(dual_energy_eta2)) then
        deallocate(dual_energy_eta2)
    end if
    if (allocated(use_pslope)) then
        deallocate(use_pslope)
    end if
    if (allocated(fix_mass_flux)) then
        deallocate(fix_mass_flux)
    end if
    if (allocated(limit_fluxes_on_small_dens)) then
        deallocate(limit_fluxes_on_small_dens)
    end if
    if (allocated(density_reset_method)) then
        deallocate(density_reset_method)
    end if
    if (allocated(allow_small_energy)) then
        deallocate(allow_small_energy)
    end if
    if (allocated(do_sponge)) then
        deallocate(do_sponge)
    end if
    if (allocated(sponge_implicit)) then
        deallocate(sponge_implicit)
    end if
    if (allocated(first_order_hydro)) then
        deallocate(first_order_hydro)
    end if
    if (allocated(xl_ext_bc_type)) then
        deallocate(xl_ext_bc_type)
    end if
    if (allocated(xr_ext_bc_type)) then
        deallocate(xr_ext_bc_type)
    end if
    if (allocated(yl_ext_bc_type)) then
        deallocate(yl_ext_bc_type)
    end if
    if (allocated(yr_ext_bc_type)) then
        deallocate(yr_ext_bc_type)
    end if
    if (allocated(zl_ext_bc_type)) then
        deallocate(zl_ext_bc_type)
    end if
    if (allocated(zr_ext_bc_type)) then
        deallocate(zr_ext_bc_type)
    end if
    if (allocated(hse_zero_vels)) then
        deallocate(hse_zero_vels)
    end if
    if (allocated(hse_interp_temp)) then
        deallocate(hse_interp_temp)
    end if
    if (allocated(hse_reflect_vels)) then
        deallocate(hse_reflect_vels)
    end if
    if (allocated(mol_order)) then
        deallocate(mol_order)
    end if
    if (allocated(cfl)) then
        deallocate(cfl)
    end if
    if (allocated(dtnuc_e)) then
        deallocate(dtnuc_e)
    end if
    if (allocated(dtnuc_X)) then
        deallocate(dtnuc_X)
    end if
    if (allocated(dtnuc_X_threshold)) then
        deallocate(dtnuc_X_threshold)
    end if
    if (allocated(dxnuc)) then
        deallocate(dxnuc)
    end if
    if (allocated(dxnuc_max)) then
        deallocate(dxnuc_max)
    end if
    if (allocated(max_dxnuc_lev)) then
        deallocate(max_dxnuc_lev)
    end if
    if (allocated(do_react)) then
        deallocate(do_react)
    end if
    if (allocated(react_T_min)) then
        deallocate(react_T_min)
    end if
    if (allocated(react_T_max)) then
        deallocate(react_T_max)
    end if
    if (allocated(react_rho_min)) then
        deallocate(react_rho_min)
    end if
    if (allocated(react_rho_max)) then
        deallocate(react_rho_max)
    end if
    if (allocated(disable_shock_burning)) then
        deallocate(disable_shock_burning)
    end if
    if (allocated(diffuse_cutoff_density)) then
        deallocate(diffuse_cutoff_density)
    end if
    if (allocated(diffuse_cond_scale_fac)) then
        deallocate(diffuse_cond_scale_fac)
    end if
    if (allocated(do_grav)) then
        deallocate(do_grav)
    end if
    if (allocated(grav_source_type)) then
        deallocate(grav_source_type)
    end if
    if (allocated(do_rotation)) then
        deallocate(do_rotation)
    end if
    if (allocated(rot_period)) then
        deallocate(rot_period)
    end if
    if (allocated(rot_period_dot)) then
        deallocate(rot_period_dot)
    end if
    if (allocated(rotation_include_centrifugal)) then
        deallocate(rotation_include_centrifugal)
    end if
    if (allocated(rotation_include_coriolis)) then
        deallocate(rotation_include_coriolis)
    end if
    if (allocated(rotation_include_domegadt)) then
        deallocate(rotation_include_domegadt)
    end if
    if (allocated(state_in_rotating_frame)) then
        deallocate(state_in_rotating_frame)
    end if
    if (allocated(rot_source_type)) then
        deallocate(rot_source_type)
    end if
    if (allocated(implicit_rotation_update)) then
        deallocate(implicit_rotation_update)
    end if
    if (allocated(rot_axis)) then
        deallocate(rot_axis)
    end if
    if (allocated(use_point_mass)) then
        deallocate(use_point_mass)
    end if
    if (allocated(point_mass)) then
        deallocate(point_mass)
    end if
    if (allocated(point_mass_fix_solution)) then
        deallocate(point_mass_fix_solution)
    end if
    if (allocated(do_acc)) then
        deallocate(do_acc)
    end if
    if (allocated(grown_factor)) then
        deallocate(grown_factor)
    end if
    if (allocated(track_grid_losses)) then
        deallocate(track_grid_losses)
    end if
    if (allocated(const_grav)) then
        deallocate(const_grav)
    end if
    if (allocated(get_g_from_phi)) then
        deallocate(get_g_from_phi)
    end if


    
  end subroutine ca_finalize_meth_params


#ifdef RADIATION
  subroutine ca_init_radhydro_pars(fsp_type_in, do_is_in, com_in,fppt) &
       bind(C, name="ca_init_radhydro_pars")

    use rad_params_module, only : ngroups

    use amrex_fort_module, only : rt => amrex_real

    implicit none

    integer, intent(in) :: fsp_type_in, do_is_in, com_in
    real(rt)        , intent(in) :: fppt


    if (ngroups .eq. 1) then
       fspace_type = 1
    else
       fspace_type = fsp_type_in
    end if

#ifndef AMREX_USE_GPU
    if (fsp_type_in .ne. 1 .and. fsp_type_in .ne. 2) then
       call amrex_error("Unknown fspace_type", fspace_type)
    end if
#endif
    
    do_inelastic_scattering = (do_is_in .ne. 0)
    
    if (com_in .eq. 1) then
       comoving = .true.
    else if (com_in .eq. 0) then
       comoving = .false.
    else
#ifndef AMREX_USE_GPU
       call amrex_error("Wrong value for comoving", fspace_type)
#endif
    end if
    
    flatten_pp_threshold = fppt
    
    !$acc update &
    !$acc device(QRAD, QRADHI, QPTOT, QREITOT) &
    !$acc device(fspace_type) &
    !$acc device(do_inelastic_scattering) &
    !$acc device(comoving)
    !$acc device(flatten_pp_threshold = -1.e0_rt)

  end subroutine ca_init_radhydro_pars
#endif

end module meth_params_module
