! This module stores the runtime parameters and integer names for 
! indexing arrays.
!
! These parameter are initialized in set_method_params() 

module meth_params_module

  implicit none

  integer         , save :: iorder        ! used only in uslope 

  ! number of ghost cells for the hyperbolic solver
  integer, parameter     :: NHYP    = 4

  ! NTHERM: number of thermodynamic variables
  integer         , save :: NTHERM, NVAR
  integer         , save :: URHO, UMX, UMY, UMZ, UEDEN, UEINT, UTEMP, UFA, UFS, UFX

  ! QTHERM: number of primitive variables
  integer         , save :: QTHERM, QVAR
  integer         , save :: QRHO, QU, QV, QW, QPRES, QREINT, QTEMP
  integer         , save :: QGAMC, QGAME
  integer         , save :: QFA, QFS, QFX

  ! These are only used when we use the SGS model.
  integer         , save :: UESGS,QESGS

  integer         , save :: nadv

  integer, save :: npassive
  integer, save, allocatable :: qpass_map(:), upass_map(:)

  ! These are used for the Godunov state
  ! Note that the velocity indices here are picked to be the same value
  ! as in the primitive variable array
  integer, save :: ngdnv, GDRHO, GDU, GDV, GDW, GDPRES, GDGAME, GDLAMS, GDERADS

  integer         , save :: numpts_1d

  double precision, save, allocatable :: outflow_data_old(:,:)
  double precision, save, allocatable :: outflow_data_new(:,:)
  double precision, save :: outflow_data_old_time
  double precision, save :: outflow_data_new_time
  logical         , save :: outflow_data_allocated
  double precision, save :: max_dist

  double precision, save :: diffuse_cutoff_density

  double precision, save :: const_grav

  logical, save :: get_g_from_phi
  
  character(len=:), allocatable :: gravity_type
  

  double precision, save :: difmag
  double precision, save :: small_dens
  double precision, save :: small_temp
  double precision, save :: small_pres
  double precision, save :: small_ener
  integer         , save :: hybrid_hydro
  integer         , save :: ppm_type
  integer         , save :: ppm_reference
  integer         , save :: ppm_trace_sources
  integer         , save :: ppm_temp_fix
  integer         , save :: ppm_tau_in_tracing
  integer         , save :: ppm_predict_gammae
  integer         , save :: ppm_reference_edge_limit
  integer         , save :: ppm_reference_eigenvectors
  integer         , save :: hybrid_riemann
  integer         , save :: use_colglaz
  integer         , save :: riemann_solver
  integer         , save :: cg_maxiter
  double precision, save :: cg_tol
  integer         , save :: use_flattening
  integer         , save :: ppm_flatten_before_integrals
  integer         , save :: transverse_use_eos
  integer         , save :: transverse_reset_density
  integer         , save :: transverse_reset_rhoe
  logical         , save :: dual_energy_update_E_from_e
  double precision, save :: dual_energy_eta1
  double precision, save :: dual_energy_eta2
  double precision, save :: dual_energy_eta3
  integer         , save :: use_pslope
  integer         , save :: normalize_species
  integer         , save :: fix_mass_flux
  integer         , save :: allow_negative_energy
  integer         , save :: do_sponge
  double precision, save :: dtnuc
  double precision, save :: dxnuc
  double precision, save :: react_T_min
  double precision, save :: react_T_max
  integer         , save :: do_grav
  integer         , save :: grav_source_type
  integer         , save :: do_rotation
  double precision, save :: rot_period
  double precision, save :: rot_period_dot
  integer         , save :: rot_source_type
  integer         , save :: rot_axis
  double precision, save :: point_mass
  integer         , save :: do_acc

  double precision, save :: rot_vec(3)

  !$acc declare &
  !$acc create(NTHERM, NVAR, URHO, UMX, UMY, UMZ, UEDEN, UEINT, UTEMP, UFA, UFS, UFX) &
  !$acc create(QTHERM, QVAR, QRHO, QU, QV, QW, QPRES, QREINT, QTEMP, QGAMC, QGAME, QFA, QFS, QFX) &
  !$acc create(UESGS, QESGS, nadv, small_dens, small_temp, small_pres, small_ener) &
  !$acc create(ngdnv, GDRHO, GDU, GDV, GDW, GDPRES, GDGAME, GDLAMS, GDERADS) &
  !$acc create(allow_negative_energy, do_acc, ppm_type, ppm_reference, ppm_trace_sources) &
  !$acc create(ppm_temp_fix, ppm_tau_in_tracing, ppm_predict_gammae, ppm_reference_edge_limit) &
  !$acc create(ppm_flatten_before_integrals, ppm_reference_eigenvectors, hybrid_riemann, use_colglaz) &
  !$acc create(use_flattening, transverse_use_eos, transverse_reset_density, transverse_reset_rhoe) &
  !$acc create(burning_timestep_factor, cg_maxiter, cg_tol, use_pslope, do_grav, grav_source_type) &
  !$acc create(do_sponge, normalize_species, fix_mass_flux, numpts_1d) &
  !$acc create(dual_energy_eta1, dual_energy_eta2, dual_energy_eta3, dual_energy_update_E_from_e) &
  !$acc create(outflow_data_old, outflow_data_new, outflow_data_old_time, outflow_data_new_time, outflow_data_allocated) &
  !$acc create(max_dist, do_rotation, rot_period, rot_period_dot, rot_source_type, rot_axis, const_grav) &
  !$acc create(npassive, qpass_map, upass_map, deterministic)

end module meth_params_module
