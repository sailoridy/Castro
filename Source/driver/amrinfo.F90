
! This module contains information about the state of the Amr class in C++.

module amrinfo_module

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer, allocatable :: amr_level
  integer :: amr_iteration = 0
  integer :: amr_ncycle = 0

  real(rt)         :: amr_time = 0.0
  real(rt)         :: amr_dt = 0.0

#ifdef CUDA
  attributes(managed) :: amr_level
#endif

contains
  
  subroutine amrinfo_init()
    if (.not. allocated(amr_level)) then
       allocate(amr_level)
       amr_level = 0
    endif
  end subroutine amrinfo_init

end module amrinfo_module
