subroutine amrex_probinit (init,name,namlen,problo,probhi) bind(c)

  use probdata_module
  use model_parser_module
  use amrex_error_module
  use prob_params_module, only : center

  use amrex_fort_module, only : rt => amrex_real
  implicit none
  integer init, namlen
  integer name(namlen)
  real(rt)         problo(2), probhi(2)

  integer untin,i,j,k,dir

  namelist /fortin/ &
       model_name

  !
  !     Build "probin" filename -- the name of file containing fortin namelist.
  !     
  integer, parameter :: maxlen = 127
  character probin*(maxlen)
  character model*(maxlen)
  integer ipp, ierr, ipp1

  if (namlen .gt. maxlen) call amrex_error("probin file name too long")

  do i = 1, namlen
     probin(i:i) = char(name(i))
  end do

  ! set namelist defaults

  ! Read namelists
  untin = 9 
  open(untin,file=probin(1:namlen),form='formatted',status='old')
  read(untin,fortin)
  close(unit=untin)


  ! read initial model
  call read_model_file(model_name)

  ! assume axisymmetric
  center(1) = 0.e0_rt
  center(2) = 0.5e0_rt*(problo(2)+probhi(2))

end subroutine amrex_probinit

! ::: -----------------------------------------------------------
! ::: This routine is called at problem setup time and is used
! ::: to initialize data on each grid.  
! ::: 
! ::: NOTE:  all arrays have one cell of ghost zones surrounding
! :::        the grid interior.  Values in these cells need not
! :::        be set here.
! ::: 
! ::: INPUTS/OUTPUTS:
! ::: 
! ::: level     => amr level of grid
! ::: time      => time at which to init data             
! ::: lo,hi     => index limits of grid interior (cell centered)
! ::: nstate    => number of state components.  You should know
! :::		   this already!
! ::: state     <=  Scalar array
! ::: delta     => cell size
! ::: xlo,xhi   => physical locations of lower left and upper
! :::              right hand corner of grid.  (does not include
! :::		   ghost region).
! ::: -----------------------------------------------------------
subroutine ca_initdata(level,time,lo,hi,nscal, &
     state,state_l1,state_l2,state_h1,state_h2, &
     delta,xlo,xhi)

  use probdata_module
  use interpolate_module
  use eos_module
  use meth_params_module, only : NVAR, URHO, UMX, UMY, UTEMP,&
       UEDEN, UEINT, UFS
  use network, only : nspec
  use model_parser_module
  use prob_params_module, only : center
  use eos_type_module
  use eos_module

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer level, nscal
  integer lo(2), hi(2)
  integer state_l1,state_l2,state_h1,state_h2
  real(rt)         xlo(2), xhi(2), time, delta(2)
  real(rt)         state(state_l1:state_h1,state_l2:state_h2,NVAR)

  real(rt)         xcen,ycen,dist,pres
  integer i,j,n

  type(eos_t) :: eos_state

  do j = lo(2), hi(2)
     ycen = xlo(2) + delta(2)*(float(j-lo(2)) + 0.5e0_rt) - center(2)

     do i = lo(1), hi(1)   
        xcen = xlo(1) + delta(1)*(float(i-lo(1)) + 0.5e0_rt) - center(1)

        dist = sqrt(xcen**2 + ycen**2)

        state(i,j,URHO)  = interpolate(dist,npts_model,model_r,model_state(:,idens_model))
        state(i,j,UTEMP) = interpolate(dist,npts_model,model_r,model_state(:,itemp_model))

        do n = 1, nspec
           state(i,j,UFS-1+n) = interpolate(dist,npts_model,model_r,model_state(:,ispec_model-1+n))
        enddo

     enddo
  enddo

  do j = lo(2), hi(2)     
     do i = lo(1), hi(1)   
        eos_state%rho = state(i,j,URHO)
        eos_state%T = state(i,j,UTEMP)
        eos_state%xn(:) = state(i,j,UFS:UFS-1+nspec)

        call eos(eos_input_rt, eos_state)

        state(i,j,UEDEN) = eos_state%e

        state(i,j,UEINT) = state(i,j,URHO) * state(i,j,UEDEN)
        state(i,j,UEDEN) = state(i,j,URHO) * state(i,j,UEDEN)

        do n = 1,nspec
           state(i,j,UFS+n-1) = state(i,j,URHO) * state(i,j,UFS+n-1)
        end do

     enddo
  enddo

  ! Initial velocities = 0
  state(:,:,UMX:UMY) = 0.e0_rt

end subroutine ca_initdata
