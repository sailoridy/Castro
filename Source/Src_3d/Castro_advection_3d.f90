module advection_module

  implicit none

  private

  public umeth3d, ctoprim, divu, consup, enforce_minimum_density, normalize_new_species, &
       normalize_species_fluxes
  
contains

! ::: ---------------------------------------------------------------
! ::: :: UMETH3D     Compute hyperbolic fluxes using unsplit second
! ::: ::               order Godunov integrator.
! ::: :: 
! ::: :: inputs/outputs
! ::: :: q           => (const)  input state, primitives
! ::: :: c           => (const)  sound speed
! ::: :: gamc        => (const)  cound speed gamma
! ::: :: csml        => (const)  local small c val
! ::: :: flatn       => (const)  flattening parameter
! ::: :: src         => (const)  source
! ::: :: nx          => (const)  number of cells in X direction
! ::: :: ny          => (const)  number of cells in Y direction
! ::: :: nz          => (const)  number of cells in Z direction
! ::: :: dx          => (const)  grid spacing in X direction
! ::: :: dy          => (const)  grid spacing in Y direction
! ::: :: dz          => (const)  grid spacing in Z direction
! ::: :: dt          => (const)  time stepsize
! ::: :: flux1      <=  (modify) flux in X direction on X edges
! ::: :: flux2      <=  (modify) flux in Y direction on Y edges
! ::: :: flux3      <=  (modify) flux in Z direction on Z edges
! ::: ----------------------------------------------------------------

  subroutine umeth3d(lo, hi, dx, dy, dz, dt, domlo, domhi, &
                     q, c, gamc, csml, flatn, qlo, qhi, &
                     srcQ, slo, shi, &
                     grav, gvlo, gvhi, &
                     rot, rlo, rhi, &
                     flux1, f1lo, f1hi, &
                     flux2, f2lo, f2hi, &
                     flux3, f3lo, f3hi, &
                     ugdx, ugxlo, ugxhi, &
                     ugdy, ugylo, ugyhi, &
                     ugdz, ugzlo, ugzhi, &
                     pdivu, plo, phi )

    use meth_params_module, only : QVAR, NVAR, QPRES, QRHO, QU, QW, QFS, QFX, QTEMP, QREINT, ppm_type, &
                                   use_pslope, ppm_trace_grav, ppm_trace_rot, ppm_temp_fix, &
                                   do_grav, do_rotation, hybrid_riemann
    use trace_ppm_module, only : trace_ppm
    use trace_module, only : tracexy, tracez
    use transverse_module
    use ppm_module, only : ppm
    use slope_module, only : uslope, pslope
    use network
    use eos_module
    use eos_type_module
    use riemann_module, only: cmpflx, shock
    use bl_constants_module

    implicit none

    integer, intent(in) :: lo(3), hi(3), domlo(3), domhi(3), qlo(3), qhi(3), slo(3), shi(3), rlo(3), rhi(3)
    integer, intent(in) :: gvlo(3), gvhi(3), f1lo(3), f1hi(3), f2lo(3), f2hi(3), f3lo(3), f3hi(3)
    integer, intent(in) :: ugxlo(3), ugxhi(3), ugylo(3), ugyhi(3), ugzlo(3), ugzhi(3), plo(3), phi(3)
    double precision, intent(in) :: dx, dy, dz, dt
    double precision,intent(in   ):: q    (qlo(1):qhi(1),qlo(2):qhi(2),qlo(3):qhi(3),QVAR)
    double precision,intent(in   ):: c    (qlo(1):qhi(1),qlo(2):qhi(2),qlo(3):qhi(3))
    double precision,intent(in   ):: gamc (qlo(1):qhi(1),qlo(2):qhi(2),qlo(3):qhi(3))
    double precision,intent(in   ):: csml (qlo(1):qhi(1),qlo(2):qhi(2),qlo(3):qhi(3))
    double precision,intent(in   ):: flatn(qlo(1):qhi(1),qlo(2):qhi(2),qlo(3):qhi(3))
    double precision,intent(in   ):: srcQ (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),QVAR)
    double precision,intent(in   ):: rot  (rlo(1):rhi(1),rlo(2):rhi(2),rlo(3):rhi(3),3)
    double precision,intent(in   ):: grav (gvlo(1):gvhi(1),gvlo(2):gvhi(2),gvlo(3):gvhi(3),3)
    double precision,intent(inout):: flux1(f1lo(1):f1hi(1),f1lo(2):f1hi(2),f1lo(3):f1hi(3),NVAR)
    double precision,intent(inout):: flux2(f2lo(1):f2hi(1),f2lo(2):f2hi(2),f2lo(3):f2hi(3),NVAR)
    double precision,intent(inout):: flux3(f3lo(1):f3hi(1),f3lo(2):f3hi(2),f3lo(3):f3hi(3),NVAR)
    double precision,intent(inout):: ugdx(ugxlo(1):ugxhi(1),ugxlo(2):ugxhi(2),ugxlo(3):ugxhi(3))
    double precision,intent(inout):: ugdy(ugylo(1):ugyhi(1),ugylo(2):ugyhi(2),ugylo(3):ugyhi(3))
    double precision,intent(inout):: ugdz(ugzlo(1):ugzhi(1),ugzlo(2):ugzhi(2),ugzlo(3):ugzhi(3))
    double precision,intent(  out):: pdivu(plo(1):phi(1),plo(2):phi(2),plo(3):phi(3))

    integer :: i,j,k,n,iwave,idim
    double precision :: dtdx, dtdy, dtdz, hdt
    double precision :: cdtdx, cdtdy, cdtdz
    double precision :: hdtdx, hdtdy, hdtdz

    ! Left and right state arrays (edge centered, cell centered)
    double precision, allocatable:: dqx(:,:,:,:), dqy(:,:,:,:), dqz(:,:,:,:)
    double precision, allocatable::qxm(:,:,:,:),qym(:,:,:,:), qzm(:,:,:,:)
    double precision, allocatable::qxp(:,:,:,:),qyp(:,:,:,:), qzp(:,:,:,:)
    
    double precision, allocatable::qmxy(:,:,:,:),qpxy(:,:,:,:)
    double precision, allocatable::qmxz(:,:,:,:),qpxz(:,:,:,:)
    
    double precision, allocatable::qmyx(:,:,:,:),qpyx(:,:,:,:)
    double precision, allocatable::qmyz(:,:,:,:),qpyz(:,:,:,:)
    
    double precision, allocatable::qmzx(:,:,:,:),qpzx(:,:,:,:)
    double precision, allocatable::qmzy(:,:,:,:),qpzy(:,:,:,:)
    
    double precision, allocatable::qxl(:,:,:,:),qxr(:,:,:,:)
    double precision, allocatable::qyl(:,:,:,:),qyr(:,:,:,:)
    double precision, allocatable::qzl(:,:,:,:),qzr(:,:,:,:)
    
    ! Work arrays to hold 3 planes of riemann state and conservative fluxes
    double precision, allocatable::   fx(:,:,:,:),  fy(:,:,:,:), fz(:,:,:,:)
    
    double precision, allocatable::fxy(:,:,:,:),fxz(:,:,:,:)
    double precision, allocatable::fyx(:,:,:,:),fyz(:,:,:,:)
    double precision, allocatable::fzx(:,:,:,:),fzy(:,:,:,:)
    
    double precision, allocatable:: pgdnvx(:,:,:), ugdnvx(:,:,:), gegdnvx(:,:,:)
    double precision, allocatable:: pgdnvxf(:,:,:), ugdnvxf(:,:,:), gegdnvxf(:,:,:)
    double precision, allocatable:: pgdnvtmpx(:,:,:), ugdnvtmpx(:,:,:), gegdnvtmpx(:,:,:)
    
    double precision, allocatable:: pgdnvy(:,:,:), ugdnvy(:,:,:), gegdnvy(:,:,:)
    double precision, allocatable:: pgdnvyf(:,:,:), ugdnvyf(:,:,:), gegdnvyf(:,:,:)
    double precision, allocatable:: pgdnvtmpy(:,:,:), ugdnvtmpy(:,:,:), gegdnvtmpy(:,:,:)
    
    double precision, allocatable:: pgdnvz(:,:,:), ugdnvz(:,:,:), gegdnvz(:,:,:)
    double precision, allocatable:: pgdnvzf(:,:,:), ugdnvzf(:,:,:), gegdnvzf(:,:,:)
    double precision, allocatable:: pgdnvtmpz1(:,:,:), ugdnvtmpz1(:,:,:), gegdnvtmpz1(:,:,:)
    double precision, allocatable:: pgdnvtmpz2(:,:,:), ugdnvtmpz2(:,:,:), gegdnvtmpz2(:,:,:)
    
    double precision, allocatable:: Ip(:,:,:,:,:,:), Im(:,:,:,:,:,:)
    double precision, allocatable:: Ip_g(:,:,:,:,:,:), Im_g(:,:,:,:,:,:)
    double precision, allocatable:: Ip_r(:,:,:,:,:,:), Im_r(:,:,:,:,:,:)
    double precision, allocatable:: Ip_gc(:,:,:,:,:,:), Im_gc(:,:,:,:,:,:)

    double precision, allocatable :: shk(:,:,:)
    
    type (eos_t) :: eos_state

    integer :: fglo(3), fghi(3), glo(3), ghi(3)

    fglo = lo-1  ! face + one ghost
    fghi = hi+2  

    glo = lo-1  ! one ghost,  this can be used for face-based arrays too
    ghi = hi+1  
    
    allocate ( pgdnvx(fglo(1):fghi(1),fglo(2):fghi(2),fglo(3):fghi(3)))
    allocate ( ugdnvx(fglo(1):fghi(1),fglo(2):fghi(2),fglo(3):fghi(3)))
    allocate (gegdnvx(fglo(1):fghi(1),fglo(2):fghi(2),fglo(3):fghi(3)))

    allocate ( pgdnvxf(fglo(1):fghi(1),fglo(2):fghi(2),fglo(3):fghi(3)))
    allocate ( ugdnvxf(fglo(1):fghi(1),fglo(2):fghi(2),fglo(3):fghi(3)))
    allocate (gegdnvxf(fglo(1):fghi(1),fglo(2):fghi(2),fglo(3):fghi(3)))

    allocate ( pgdnvtmpx(fglo(1):fghi(1),fglo(2):fghi(2),fglo(3):fghi(3)))
    allocate ( ugdnvtmpx(fglo(1):fghi(1),fglo(2):fghi(2),fglo(3):fghi(3)))
    allocate (gegdnvtmpx(fglo(1):fghi(1),fglo(2):fghi(2),fglo(3):fghi(3)))
    
    allocate ( pgdnvy(fglo(1):fghi(1),fglo(2):fghi(2),fglo(3):fghi(3)))
    allocate ( ugdnvy(fglo(1):fghi(1),fglo(2):fghi(2),fglo(3):fghi(3)))
    allocate (gegdnvy(fglo(1):fghi(1),fglo(2):fghi(2),fglo(3):fghi(3)))

    allocate ( pgdnvyf(fglo(1):fghi(1),fglo(2):fghi(2),fglo(3):fghi(3)))
    allocate ( ugdnvyf(fglo(1):fghi(1),fglo(2):fghi(2),fglo(3):fghi(3)))
    allocate (gegdnvyf(fglo(1):fghi(1),fglo(2):fghi(2),fglo(3):fghi(3)))

    allocate ( pgdnvtmpy(fglo(1):fghi(1),fglo(2):fghi(2),fglo(3):fghi(3)))
    allocate ( ugdnvtmpy(fglo(1):fghi(1),fglo(2):fghi(2),fglo(3):fghi(3)))
    allocate (gegdnvtmpy(fglo(1):fghi(1),fglo(2):fghi(2),fglo(3):fghi(3)))

    allocate ( pgdnvz(fglo(1):fghi(1),fglo(2):fghi(2),fglo(3):fghi(3)))
    allocate ( ugdnvz(fglo(1):fghi(1),fglo(2):fghi(2),fglo(3):fghi(3)))
    allocate (gegdnvz(fglo(1):fghi(1),fglo(2):fghi(2),fglo(3):fghi(3)))

    allocate ( pgdnvzf(fglo(1):fghi(1),fglo(2):fghi(2),fglo(3):fghi(3)))
    allocate ( ugdnvzf(fglo(1):fghi(1),fglo(2):fghi(2),fglo(3):fghi(3)))
    allocate (gegdnvzf(fglo(1):fghi(1),fglo(2):fghi(2),fglo(3):fghi(3)))

    allocate ( pgdnvtmpz1(fglo(1):fghi(1),fglo(2):fghi(2),fglo(3):fghi(3)))
    allocate ( ugdnvtmpz1(fglo(1):fghi(1),fglo(2):fghi(2),fglo(3):fghi(3)))
    allocate (gegdnvtmpz1(fglo(1):fghi(1),fglo(2):fghi(2),fglo(3):fghi(3)))

    allocate ( pgdnvtmpz2(fglo(1):fghi(1),fglo(2):fghi(2),fglo(3):fghi(3)))
    allocate ( ugdnvtmpz2(fglo(1):fghi(1),fglo(2):fghi(2),fglo(3):fghi(3)))
    allocate (gegdnvtmpz2(fglo(1):fghi(1),fglo(2):fghi(2),fglo(3):fghi(3)))
    
    allocate ( dqx(fglo(1):fghi(1),fglo(2):fghi(2),fglo(3):fghi(3),QVAR))
    allocate ( dqy(fglo(1):fghi(1),fglo(2):fghi(2),fglo(3):fghi(3),QVAR))
    allocate ( dqz(fglo(1):fghi(1),fglo(2):fghi(2),fglo(3):fghi(3),QVAR))
    
    allocate ( qxm(fglo(1):fghi(1),fglo(2):fghi(2),fglo(3):fghi(3),QVAR))
    allocate ( qxp(fglo(1):fghi(1),fglo(2):fghi(2),fglo(3):fghi(3),QVAR))

    allocate ( qmxy(fglo(1):fghi(1),fglo(2):fghi(2),fglo(3):fghi(3),QVAR))
    allocate ( qpxy(fglo(1):fghi(1),fglo(2):fghi(2),fglo(3):fghi(3),QVAR))

    allocate ( qmxz(fglo(1):fghi(1),fglo(2):fghi(2),fglo(3):fghi(3),QVAR))
    allocate ( qpxz(fglo(1):fghi(1),fglo(2):fghi(2),fglo(3):fghi(3),QVAR))

    allocate ( qym(fglo(1):fghi(1),fglo(2):fghi(2),fglo(3):fghi(3),QVAR))
    allocate ( qyp(fglo(1):fghi(1),fglo(2):fghi(2),fglo(3):fghi(3),QVAR))

    allocate ( qmyx(fglo(1):fghi(1),fglo(2):fghi(2),fglo(3):fghi(3),QVAR))
    allocate ( qpyx(fglo(1):fghi(1),fglo(2):fghi(2),fglo(3):fghi(3),QVAR))

    allocate ( qmyz(fglo(1):fghi(1),fglo(2):fghi(2),fglo(3):fghi(3),QVAR))
    allocate ( qpyz(fglo(1):fghi(1),fglo(2):fghi(2),fglo(3):fghi(3),QVAR))

    allocate ( qzm(fglo(1):fghi(1),fglo(2):fghi(2),fglo(3):fghi(3),QVAR))
    allocate ( qzp(fglo(1):fghi(1),fglo(2):fghi(2),fglo(3):fghi(3),QVAR))

    allocate ( qxl(fglo(1):fghi(1),fglo(2):fghi(2),fglo(3):fghi(3),QVAR))
    allocate ( qxr(fglo(1):fghi(1),fglo(2):fghi(2),fglo(3):fghi(3),QVAR))
    allocate ( qyl(fglo(1):fghi(1),fglo(2):fghi(2),fglo(3):fghi(3),QVAR))
    allocate ( qyr(fglo(1):fghi(1),fglo(2):fghi(2),fglo(3):fghi(3),QVAR))
    allocate ( qzl(fglo(1):fghi(1),fglo(2):fghi(2),fglo(3):fghi(3),QVAR))
    allocate ( qzr(fglo(1):fghi(1),fglo(2):fghi(2),fglo(3):fghi(3),QVAR))

    allocate ( qmzx(fglo(1):fghi(1),fglo(2):fghi(2),fglo(3):fghi(3),QVAR))
    allocate ( qpzx(fglo(1):fghi(1),fglo(2):fghi(2),fglo(3):fghi(3),QVAR))

    allocate ( qmzy(fglo(1):fghi(1),fglo(2):fghi(2),fglo(3):fghi(3),QVAR))
    allocate ( qpzy(fglo(1):fghi(1),fglo(2):fghi(2),fglo(3):fghi(3),QVAR))

    allocate ( fx(glo(1):ghi(1),glo(2):ghi(2),glo(3):ghi(3),NVAR))
    allocate ( fy(glo(1):ghi(1),glo(2):ghi(2),glo(3):ghi(3),NVAR))
    allocate ( fz(glo(1):ghi(1),glo(2):ghi(2),glo(3):ghi(3),NVAR))

    allocate ( fxy(glo(1):ghi(1),glo(2):ghi(2),glo(3):ghi(3),NVAR))
    allocate ( fxz(glo(1):ghi(1),glo(2):ghi(2),glo(3):ghi(3),NVAR))

    allocate ( fyx(glo(1):ghi(1),glo(2):ghi(2),glo(3):ghi(3),NVAR))
    allocate ( fyz(glo(1):ghi(1),glo(2):ghi(2),glo(3):ghi(3),NVAR))

    allocate ( fzx(glo(1):ghi(1),glo(2):ghi(2),glo(3):ghi(3),NVAR))
    allocate ( fzy(glo(1):ghi(1),glo(2):ghi(2),glo(3):ghi(3),NVAR))

    ! x-index, y-index, z-index, dim, characteristics, variables
    allocate ( Ip(glo(1):ghi(1),glo(2):ghi(2),glo(3):ghi(3),3,3,QVAR))
    allocate ( Im(glo(1):ghi(1),glo(2):ghi(2),glo(3):ghi(3),3,3,QVAR))
    
    ! for gravity (last index is x,y,z component)
    allocate ( Ip_g(glo(1):ghi(1),glo(2):ghi(2),glo(3):ghi(3),3,3,3))
    allocate ( Im_g(glo(1):ghi(1),glo(2):ghi(2),glo(3):ghi(3),3,3,3))

    ! for rotation (last index is x,y,z component)
    allocate ( Ip_r(glo(1):ghi(1),glo(2):ghi(2),glo(3):ghi(3),3,3,3))
    allocate ( Im_r(glo(1):ghi(1),glo(2):ghi(2),glo(3):ghi(3),3,3,3))

    ! for gamc -- needed for the reference state in eigenvectors
    allocate ( Ip_gc(glo(1):ghi(1),glo(2):ghi(2),glo(3):ghi(3),3,3,1))
    allocate ( Im_gc(glo(1):ghi(1),glo(2):ghi(2),glo(3):ghi(3),3,3,1))

    ! for the hybrid Riemann solver
    allocate(shk(glo(1):ghi(1),glo(2):ghi(2),glo(3):ghi(3)))
    
    ! Local constants
    dtdx = dt/dx
    dtdy = dt/dy
    dtdz = dt/dz
    hdt = HALF*dt
    hdtdx = HALF*dtdx
    hdtdy = HALF*dtdy
    hdtdz = HALF*dtdz
    cdtdx = dtdx*THIRD
    cdtdy = dtdy*THIRD
    cdtdz = dtdz*THIRD

    ! Initialize pdivu to zero
    pdivu(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)) = ZERO

    ! multidimensional shock detection -- this will be used to do the
    ! hybrid Riemann solver
    if (hybrid_riemann == 1) then
       call shock(q,qlo(1),qlo(2),qlo(3),qhi(1),qhi(2),qhi(3), &
            shk,glo(1),glo(2),glo(3),ghi(1),ghi(2),ghi(3), &
            lo(1),lo(2),lo(3),hi(1),hi(2),hi(3),dx,dy,dz)
    else
       shk(:,:,:) = ZERO
    endif
    
    if (ppm_type .gt. 0) then

       do n=1,QVAR
          call ppm(q(:,:,:,n), qlo, qhi, &
                   q(:,:,:,QU:QW), c, qlo, qhi, &
                   flatn, qlo, qhi, &
                   Ip(:,:,:,:,:,n),Im(:,:,:,:,:,n), glo, ghi, &
                   lo,hi,dx,dy,dz,dt)
       end do
       
       if (do_grav .eq. 1 .and. ppm_trace_grav .eq. 1) then
          do n=1,3
             call ppm(grav(:,:,:,n), gvlo, gvhi, &
                      q(:,:,:,QU:),c, qlo, qhi, &
                      flatn, qlo, qhi, &
                      Ip_g(:,:,:,:,:,n),Im_g(:,:,:,:,:,n), glo, ghi, &
                      lo,hi,dx,dy,dz,dt)
          enddo
       endif

       if (do_rotation .eq. 1 .and. ppm_trace_rot .eq. 1) then
          do n=1,3
             call ppm(rot(:,:,:,n), slo, shi, &
                  q(:,:,:,QU:),c, qlo, qhi, &
                  flatn, qlo, qhi, &
                  Ip_r(:,:,:,:,:,n),Im_r(:,:,:,:,:,n), glo, ghi, &
                  lo,hi,dx,dy,dz,dt)
          enddo
       endif


       if (ppm_temp_fix /= 1) then
          call ppm(gamc(:,:,:), qlo, qhi, &
               q(:,:,:,QU:),c, qlo, qhi, &
               flatn, qlo, qhi, &
               Ip_gc(:,:,:,:,:,1),Im_gc(:,:,:,:,:,1), glo, qhi, &
               lo,hi,dx,dy,dz,dt)
       endif
       
       ! temperature-based PPM
       if (ppm_temp_fix == 1) then
          do k = glo(3), ghi(3)
          do j = glo(2), ghi(2)
          do i = glo(1), ghi(1)
             do idim = 1, 3
                do iwave = 1, 3
                   eos_state % rho = Ip(i,j,k,idim,iwave,QRHO)
                   eos_state % T   = Ip(i,j,k,idim,iwave,QTEMP)
                   eos_state % xn  = Ip(i,j,k,idim,iwave,QFS:QFS-1+nspec)
                   eos_state % aux = Ip(i,j,k,idim,iwave,QFX:QFX-1+naux)
                   
                   call eos(eos_input_rt, eos_state, .false.)
                   
                   Ip(i,j,k,idim,iwave,QPRES) = eos_state%p
                   Ip(i,j,k,idim,iwave,QREINT) = Ip(i,j,k,idim,iwave,QRHO)*eos_state%e
                   Ip_gc(i,j,k,idim,iwave,1) = eos_state%gam1
                   
                   
                   eos_state % rho = Im(i,j,k,idim,iwave,QRHO)
                   eos_state % T   = Im(i,j,k,idim,iwave,QTEMP)
                   eos_state % xn  = Im(i,j,k,idim,iwave,QFS:QFS-1+nspec)
                   eos_state % aux = Im(i,j,k,idim,iwave,QFX:QFX-1+naux)
                   
                   call eos(eos_input_rt, eos_state, .false.)
                   
                   Im(i,j,k,idim,iwave,QPRES) = eos_state%p
                   Im(i,j,k,idim,iwave,QREINT) = Im(i,j,k,idim,iwave,QRHO)*eos_state%e
                   Im_gc(i,j,k,idim,iwave,1) = eos_state%gam1
                   
                enddo
             enddo
          enddo
          enddo
          enddo
          
       endif

       ! Compute U_x and U_y
       ! Inputs: q, c, gamc, flatn             : lo-4:hi+4
       !         Ip,Im,Ip_g,Im_g,Ip_r,Im_r,... : lo-1:hi+1
       ! Outputs: qxm, qxp                     : xface, +-1 at y & z
       !          qym, qyp                     : yface, +-1 at x & z
       !          qzm, qzp                     : zface, +-1 at x & y
       call trace_ppm(q,c,gamc,flatn, qlo,qhi, &
                      Ip,Im,Ip_g,Im_g,Ip_r,Im_r,Ip_gc,Im_gc, glo, qhi, &
                      qxm,qxp,qym,qyp,qzm,qzp, fglo, fghi, &
                      lo,hi,dt)
       
    else
       
    !    ! Compute all slopes at kc (k3d)
    !    call uslope(q,flatn,qlo(1),qlo(2),qlo(3),qhi(1),qhi(2),qhi(3), &
    !         dqx,dqy,dqz,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
    !         ilo1,ilo2,ihi1,ihi2,kc,k3d,QVAR)
       
    !    if (use_pslope .eq. 1) &
    !         call pslope(q(:,:,:,QPRES),q(:,:,:,QRHO), &
    !         flatn,qlo(1),qlo(2),qlo(3),qhi(1),qhi(2),qhi(3), &
    !         dqx(:,:,:,QPRES),dqy(:,:,:,QPRES),dqz(:,:,:,QPRES), &
    !         ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
    !         grav,gv_l1,gv_l2,gv_l3,gv_h1,gv_h2,gv_h3, &
    !         ilo1,ilo2,ihi1,ihi2,kc,k3d,dx,dy,dz)
       
    !    ! Compute U_x and U_y at kc (k3d)
    !    call tracexy(q,c,qlo(1),qlo(2),qlo(3),qhi(1),qhi(2),qhi(3), &
    !         dqx,dqy,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
    !         qxm,qxp,qym,qyp,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
    !         ilo1,ilo2,ihi1,ihi2,dx,dy,dt,kc,k3d)
       
    end if
    
    ! Inputs: qxm, qxp                     : xface, +-1 at y & z
    !         gamc, csml, c                : +-4
    !         shk                          : +-1
    ! Outputs: fx, ugdnvx, pgdnvx, gegdnvx : xface, +-1 at y & z
    call cmpflx(qxm,qxp, fglo, fghi, &
                fx, glo, ghi, &
                ugdnvx,pgdnvx,gegdnvx, fglo, fghi, &
                gamc,csml,c, qlo, qhi, &
                shk, glo, ghi, &
                1, (/lo(1),lo(2)-1,lo(3)-1/), hi+1, domlo, domhi)

    ! Inputs: qym, qyp                     : yface, +-1 at x & z
    !         qzm, qzp                     : zface, +-1 at x & y
    !         fx, ugdnvx, pgdnvx, gegdnvx  : xface, +-1 at y & z
    !         gamc                         : +-4
    ! Outputs: qmyx, qpyx                  : yface, +-1 at z
    !          qmzx, qpzx                  : zface, +-1 at y
    call transx(qym,qmyx,qyp,qpyx, &
                qzm,qmzx,qzp,qpzx, fglo, fghi, &
                fx, glo, ghi, &
                ugdnvx,pgdnvx,gegdnvx, fglo, fghi, &
                gamc, qlo, qhi, &
                cdtdx, lo, hi)

    ! Inputs: qym, qyp                     : yface, +-1 at x & z
    !         gamc, csml, c                : +-4
    !         shk                          : +-1
    ! Outputs: fy, ugdnvy, pgdnvy, gegdnvy : yface, +-1 at x & z
    call cmpflx(qym,qyp, fglo, fghi, &
                fy, glo, ghi, &
                ugdnvy,pgdnvy,gegdnvy, fglo, fghi, &
                gamc,csml,c, qlo, qhi, &
                shk, glo, ghi, &
                2, (/lo(1)-1,lo(2),lo(3)-1/), hi+1, domlo, domhi)

    ! Inputs: qxm, qxp                     : xface, +-1 at y & z
    !         qzm, qzp                     : zface, +-1 at x & y
    !         fy, ugdnvy, pgdnvy, gegdnvy  : yface, +-1 at x & z
    !         gamc                         : +-4
    ! Outputs: qmxy, qpxy                  : xface, +-1 at z
    !          qmzy, qpzy                  : zface, +-1 at x
    call transy(qxm,qmxy,qxp,qpxy, &
                qzm,qmzy,qzp,qpzy, fglo, fghi, &
                fy, glo, ghi, &
                ugdnvy,pgdnvy,gegdnvy, fglo, fghi, &
                gamc, qlo, qhi, &
                cdtdy, lo, hi)

    ! Inputs: qzm, qzp                     : zface, +-1 at x & y
    !         gamc, csml, c                : +-4
    !         shk                          : +-1
    ! Outputs: fz, ugdnvz, pgdnvz, gegdnvz : zface, +-1 at x & y
    call cmpflx(qzm,qzp, fglo, fghi, &
                fz, glo, ghi, &
                ugdnvz,pgdnvz,gegdnvz, fglo, fghi, &
                gamc,csml,c, qlo, qhi, &
                shk, glo, ghi, &
                3, (/lo(1)-1,lo(2)-1,lo(3)/), hi+1, domlo, domhi)


       
    !    ! Compute U'^y_x at kc (k3d)
    !    call transy1(qxm,qmxy,qxp,qpxy,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
    !                 fy,ilo1-1,ilo2,1,ihi1+1,ihi2+1,2, &
    !                 ugdnvy,pgdnvy,gegdnvy,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
    !                 gamc,qlo(1),qlo(2),qlo(3),qhi(1),qhi(2),qhi(3), &
    !                 cdtdy,ilo1-1,ihi1+1,ilo2,ihi2,kc,k3d)

    !    ! Compute U'^x_y at kc (k3d)
    !    call transx1(qym,qmyx,qyp,qpyx,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
    !                 fx,ilo1,ilo2-1,1,ihi1+1,ihi2+1,2, &
    !                 ugdnvx,pgdnvx,gegdnvx,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
    !                 gamc,qlo(1),qlo(2),qlo(3),qhi(1),qhi(2),qhi(3), &
    !                 cdtdx,ilo1,ihi1,ilo2-1,ihi2+1,kc,k3d)

    !    ! Compute F^{x|y} at kc (k3d)
    !    call cmpflx(qmxy,qpxy,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
    !                fxy,ilo1,ilo2-1,1,ihi1+1,ihi2+1,2, &
    !                ugdnvtmpx,pgdnvtmpx,gegdnvtmpx,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
    !                gamc,csml,c,qlo(1),qlo(2),qlo(3),qhi(1),qhi(2),qhi(3), &
    !                shk,ilo1-1,ilo2-1,ilo3-1,ihi1+1,ihi2+1,ihi3+1, &
    !                1,ilo1,ihi1+1,ilo2,ihi2,kc,kc,k3d,domlo,domhi)

    !    ! Compute F^{y|x} at kc (k3d)
    !    call cmpflx(qmyx,qpyx,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
    !                fyx,ilo1-1,ilo2,1,ihi1+1,ihi2+1,2, &
    !                ugdnvtmpy,pgdnvtmpy,gegdnvtmpy,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
    !                gamc,csml,c,qlo(1),qlo(2),qlo(3),qhi(1),qhi(2),qhi(3), &
    !                shk,ilo1-1,ilo2-1,ilo3-1,ihi1+1,ihi2+1,ihi3+1, &
    !                2,ilo1,ihi1,ilo2,ihi2+1,kc,kc,k3d,domlo,domhi)

    !    if (k3d.ge.ilo3) then
          
    !       ! Compute U'^y_z at kc (k3d)
    !       call transy2(qzm,qmzy,qzp,qpzy,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
    !                    fy,ilo1-1,ilo2,1,ihi1+1,ihi2+1,2, &
    !                    ugdnvy,pgdnvy,gegdnvy,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
    !                    gamc,qlo(1),qlo(2),qlo(3),qhi(1),qhi(2),qhi(3), &
    !                    cdtdy,ilo1-1,ihi1+1,ilo2,ihi2,kc,km,k3d)

    !       ! Compute U'^x_z at kc (k3d)
    !       call transx2(qzm,qmzx,qzp,qpzx,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
    !                    fx,ilo1,ilo2-1,1,ihi1+1,ihi2+1,2, &
    !                    ugdnvx,pgdnvx,gegdnvx,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
    !                    gamc,qlo(1),qlo(2),qlo(3),qhi(1),qhi(2),qhi(3), &
    !                    cdtdx,ilo1,ihi1,ilo2-1,ihi2+1,kc,km,k3d)

    !       ! Compute F^{z|x} at kc (k3d)
    !       call cmpflx(qmzx,qpzx,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
    !                   fzx,ilo1,ilo2-1,1,ihi1,ihi2+1,2, &
    !                   ugdnvtmpz1,pgdnvtmpz1,gegdnvtmpz1,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
    !                   gamc,csml,c,qlo(1),qlo(2),qlo(3),qhi(1),qhi(2),qhi(3), &
    !                   shk,ilo1-1,ilo2-1,ilo3-1,ihi1+1,ihi2+1,ihi3+1, &
    !                   3,ilo1,ihi1,ilo2-1,ihi2+1,kc,kc,k3d,domlo,domhi)

    !       ! Compute F^{z|y} at kc (k3d)
    !       call cmpflx(qmzy,qpzy,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
    !                   fzy,ilo1-1,ilo2,1,ihi1+1,ihi2,2, &
    !                   ugdnvtmpz2,pgdnvtmpz2,gegdnvtmpz2,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
    !                   gamc,csml,c,qlo(1),qlo(2),qlo(3),qhi(1),qhi(2),qhi(3), &
    !                   shk,ilo1-1,ilo2-1,ilo3-1,ihi1+1,ihi2+1,ihi3+1, &                       
    !                   3,ilo1-1,ihi1+1,ilo2,ihi2,kc,kc,k3d,domlo,domhi)
          
    !       ! Compute U''_z at kc (k3d)
    !       call transxy(qzm,qzl,qzp,qzr,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
    !                    fxy,ilo1,ilo2-1,1,ihi1+1,ihi2+1,2, &
    !                    fyx,ilo1-1,ilo2,1,ihi1+1,ihi2+1,2, &
    !                    ugdnvtmpx,pgdnvtmpx,gegdnvtmpx,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
    !                    ugdnvtmpy,pgdnvtmpy,gegdnvtmpy,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
    !                    gamc,qlo(1),qlo(2),qlo(3),qhi(1),qhi(2),qhi(3), &
    !                    srcQ,src_l1,src_l2,src_l3,src_h1,src_h2,src_h3, &
    !                    grav,gv_l1,gv_l2,gv_l3,gv_h1,gv_h2,gv_h3,&
    !                    rot,rt_l1,rt_l2,rt_l3,rt_h1,rt_h2,rt_h3,&
    !                    hdt,hdtdx,hdtdy,ilo1,ihi1,ilo2,ihi2,kc,km,k3d)

    !       ! Compute F^z at kc (k3d) -- note that flux3 is indexed by k3d, not kc
    !       call cmpflx(qzl,qzr,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
    !                   flux3,fd3_l1,fd3_l2,fd3_l3,fd3_h1,fd3_h2,fd3_h3, &
    !                   ugdnvzf,pgdnvzf,gegdnvzf,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
    !                   gamc,csml,c,qlo(1),qlo(2),qlo(3),qhi(1),qhi(2),qhi(3), &
    !                   shk,ilo1-1,ilo2-1,ilo3-1,ihi1+1,ihi2+1,ihi3+1, &
    !                   3,ilo1,ihi1,ilo2,ihi2,kc,k3d,k3d,domlo,domhi)

    !       do j=ilo2-1,ihi2+1
    !          do i=ilo1-1,ihi1+1
    !             ugdz(i,j,k3d) = ugdnvzf(i,j,kc)
    !          end do
    !       end do

    !       if (k3d .ge. ilo3+1 .and. k3d .le. ihi3+1) then
    !          do j = ilo2,ihi2
    !             do i = ilo1,ihi1
    !                pdivu(i,j,k3d-1) = pdivu(i,j,k3d-1) +  &
    !                     HALF*(pgdnvzf(i,j,kc)+pgdnvzf(i,j,km)) * &
    !                           (ugdnvzf(i,j,kc)-ugdnvzf(i,j,km))/dz
    !             end do
    !          end do
    !       end if
          
    !       if (k3d.gt.ilo3) then

    !          ! Compute U'^z_x and U'^z_y at km (k3d-1) -- note flux3 has physical index
    !          call transz(qxm,qmxz,qxp,qpxz, &
    !                      qym,qmyz,qyp,qpyz,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
    !                      fz,ilo1-1,ilo2-1,1,ihi1+1,ihi2+1,2, &
    !                      ugdnvz,pgdnvz,gegdnvz,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
    !                      gamc,qlo(1),qlo(2),qlo(3),qhi(1),qhi(2),qhi(3), &
    !                      cdtdz,ilo1-1,ihi1+1,ilo2-1,ihi2+1,km,kc,k3d)
         
    !          ! Compute F^{x|z} at km (k3d-1)
    !          call cmpflx(qmxz,qpxz,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
    !                      fxz,ilo1,ilo2-1,1,ihi1+1,ihi2+1,2, &
    !                      ugdnvx,pgdnvx,gegdnvx,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
    !                      gamc,csml,c,qlo(1),qlo(2),qlo(3),qhi(1),qhi(2),qhi(3), &
    !                      shk,ilo1-1,ilo2-1,ilo3-1,ihi1+1,ihi2+1,ihi3+1, &
    !                      1,ilo1,ihi1+1,ilo2-1,ihi2+1,km,km,k3d-1,domlo,domhi)

    !          ! Compute F^{y|z} at km (k3d-1)
    !          call cmpflx(qmyz,qpyz,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
    !                      fyz,ilo1-1,ilo2,1,ihi1+1,ihi2+1,2, &
    !                      ugdnvy,pgdnvy,gegdnvy,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
    !                      gamc,csml,c,qlo(1),qlo(2),qlo(3),qhi(1),qhi(2),qhi(3), &
    !                      shk,ilo1-1,ilo2-1,ilo3-1,ihi1+1,ihi2+1,ihi3+1, &
    !                      2,ilo1-1,ihi1+1,ilo2,ihi2+1,km,km,k3d-1,domlo,domhi)

    !          ! Compute U''_x at km (k3d-1)
    !          call transyz(qxm,qxl,qxp,qxr,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
    !                       fyz,ilo1-1,ilo2,1,ihi1+1,ihi2+1,2, &
    !                       fzy,ilo1-1,ilo2,1,ihi1+1,ihi2,2, &
    !                       ugdnvy,pgdnvy,gegdnvy,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
    !                       ugdnvtmpz2,pgdnvtmpz2,gegdnvtmpz2,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
    !                       gamc,qlo(1),qlo(2),qlo(3),qhi(1),qhi(2),qhi(3), &
    !                       srcQ,src_l1,src_l2,src_l3,src_h1,src_h2,src_h3, &
    !                       grav,gv_l1,gv_l2,gv_l3,gv_h1,gv_h2,gv_h3,&
    !                       rot,rt_l1,rt_l2,rt_l3,rt_h1,rt_h2,rt_h3,&
    !                       hdt,hdtdy,hdtdz,ilo1-1,ihi1+1,ilo2,ihi2,km,kc,k3d-1)

    !          ! Compute U''_y at km (k3d-1)
    !          call transxz(qym,qyl,qyp,qyr,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
    !                       fxz,ilo1,ilo2-1,1,ihi1+1,ihi2+1,2, &
    !                       fzx,ilo1,ilo2-1,1,ihi1,ihi2+1,2, &
    !                       ugdnvx,pgdnvx,gegdnvx,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
    !                       ugdnvtmpz1,pgdnvtmpz1,gegdnvtmpz1,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
    !                       gamc,qlo(1),qlo(2),qlo(3),qhi(1),qhi(2),qhi(3), &
    !                       srcQ,src_l1,src_l2,src_l3,src_h1,src_h2,src_h3, &
    !                       grav,gv_l1,gv_l2,gv_l3,gv_h1,gv_h2,gv_h3,&
    !                       rot,rt_l1,rt_l2,rt_l3,rt_h1,rt_h2,rt_h3,&
    !                       hdt,hdtdx,hdtdz,ilo1,ihi1,ilo2-1,ihi2+1,km,kc,k3d-1)

    !          ! Compute F^x at km (k3d-1)
    !          call cmpflx(qxl,qxr,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
    !                      flux1,fd1_l1,fd1_l2,fd1_l3,fd1_h1,fd1_h2,fd1_h3, &
    !                      ugdnvxf,pgdnvxf,gegdnvxf,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
    !                      gamc,csml,c,qlo(1),qlo(2),qlo(3),qhi(1),qhi(2),qhi(3), &
    !                      shk,ilo1-1,ilo2-1,ilo3-1,ihi1+1,ihi2+1,ihi3+1, &
    !                      1,ilo1,ihi1+1,ilo2,ihi2,km,k3d-1,k3d-1,domlo,domhi)
             
    !          do j=ilo2-1,ihi2+1
    !             do i=ilo1-1,ihi1+2
    !                ugdx(i,j,k3d-1) = ugdnvxf(i,j,km)
    !             end do
    !          end do
             
    !          ! Compute F^y at km (k3d-1)
    !          call cmpflx(qyl,qyr,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
    !                      flux2,fd2_l1,fd2_l2,fd2_l3,fd2_h1,fd2_h2,fd2_h3, &
    !                      ugdnvyf,pgdnvyf,gegdnvyf,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
    !                      gamc,csml,c,qlo(1),qlo(2),qlo(3),qhi(1),qhi(2),qhi(3), &
    !                      shk,ilo1-1,ilo2-1,ilo3-1,ihi1+1,ihi2+1,ihi3+1, &
    !                      2,ilo1,ihi1,ilo2,ihi2+1,km,k3d-1,k3d-1,domlo,domhi)

    !          do j=ilo2-1,ihi2+2
    !             do i=ilo1-1,ihi1+1
    !                ugdy(i,j,k3d-1) = ugdnvyf(i,j,km)
    !             end do
    !          end do

    !          do j = ilo2,ihi2
    !             do i = ilo1,ihi1
    !                pdivu(i,j,k3d-1) = pdivu(i,j,k3d-1) +  &
    !                     HALF*(pgdnvxf(i+1,j,km) + pgdnvxf(i,j,km)) *  &
    !                           (ugdnvxf(i+1,j,km)-ugdnvxf(i,j,km))/dx + &
    !                     HALF*(pgdnvyf(i,j+1,km) + pgdnvyf(i,j,km)) *  &
    !                           (ugdnvyf(i,j+1,km)-ugdnvyf(i,j,km))/dy
    !             end do
    !          end do
               
    !       end if
    !    end if
    ! enddo

    ! ! Deallocate arrays
    ! deallocate(pgdnvx,ugdnvx,gegdnvx)
    ! deallocate(pgdnvxf,ugdnvxf,gegdnvxf)
    ! deallocate(pgdnvtmpx,ugdnvtmpx,gegdnvtmpx)
    ! deallocate(pgdnvy,ugdnvy,gegdnvy)
    ! deallocate(pgdnvyf,ugdnvyf,gegdnvyf)
    ! deallocate(pgdnvtmpy,ugdnvtmpy,gegdnvtmpy)
    ! deallocate(pgdnvz,ugdnvz,gegdnvz)
    ! deallocate(pgdnvtmpz1,ugdnvtmpz1,gegdnvtmpz1)
    ! deallocate(pgdnvtmpz2,ugdnvtmpz2,gegdnvtmpz2)
    ! deallocate(pgdnvzf,ugdnvzf,gegdnvzf)
    ! deallocate(dqx,dqy,dqz)
    ! deallocate(qxm,qxp)
    ! deallocate(qmxy,qpxy)
    ! deallocate(qmxz,qpxz)
    ! deallocate(qym,qyp)
    ! deallocate(qmyx,qpyx)
    ! deallocate(qmyz,qpyz)
    ! deallocate(qzm,qzp)
    ! deallocate(qxl,qxr,qyl,qyr,qzl,qzr)
    ! deallocate(qmzx,qpzx)
    ! deallocate(qmzy,qpzy)
    ! deallocate(fx,fy,fz)
    ! deallocate(fxy,fxz)
    ! deallocate(fyx,fyz)
    ! deallocate(fzx,fzy)
    ! deallocate(Ip,Im)
    ! deallocate(Ip_g,Im_g)
    ! deallocate(Ip_gc,Im_gc)
      
  end subroutine umeth3d

! ::: 
! ::: ------------------------------------------------------------------
! ::: 

  subroutine ctoprim(lo,hi, &
                     uin, ulo, uhi, &
                     q,c,gamc,csml,flatn, qlo, qhi, &
                     src, srclo, srchi, &
                     srcQ, slo, shi, &
                     courno,dx,dy,dz,dt,ngp,ngf)
    !
    !     Will give primitive variables on lo-ngp:hi+ngp, and flatn on lo-ngf:hi+ngf
    !     if use_flattening=1.  Declared dimensions of q,c,gamc,csml,flatn are given
    !     by DIMS(q).  This declared region is assumed to encompass lo-ngp:hi+ngp.
    !     Also, uflaten call assumes ngp>=ngf+3 (ie, primitve data is used by the
    !     routine that computes flatn).  
    !
    use network, only : nspec, naux
    use eos_module
    use eos_type_module
    use meth_params_module, only : NVAR, URHO, UMX, UMY, UMZ, &
                                   UEDEN, UEINT, UESGS, UTEMP, UFA, UFS, UFX, &
                                   QVAR, QRHO, QU, QV, QW, &
                                   QREINT, QESGS, QPRES, QTEMP, QGAME, QFA, QFS, QFX, &
                                   nadv, allow_negative_energy, small_temp, use_flattening, &
                                   npassive, upass_map, qpass_map, dual_energy_eta1
    
    use flatten_module
    use bl_constants_module

    implicit none

    double precision, parameter:: small = 1.d-8

    integer, intent(in) :: lo(3), hi(3), ulo(3), uhi(3), qlo(3), qhi(3), srclo(3), srchi(3), slo(3), shi(3)
    double precision, intent(in ) :: uin  (ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3),NVAR)
    double precision, intent(out) :: q    (qlo(1):qhi(1),qlo(2):qhi(2),qlo(3):qhi(3),QVAR)
    double precision, intent(out) :: c    (qlo(1):qhi(1),qlo(2):qhi(2),qlo(3):qhi(3))
    double precision, intent(out) :: gamc (qlo(1):qhi(1),qlo(2):qhi(2),qlo(3):qhi(3))
    double precision, intent(out) :: csml (qlo(1):qhi(1),qlo(2):qhi(2),qlo(3):qhi(3))
    double precision, intent(out) :: flatn(qlo(1):qhi(1),qlo(2):qhi(2),qlo(3):qhi(3))
    double precision, intent(in ) :: src  (srclo(1):srchi(1),srclo(2):srchi(2),srclo(3):srchi(3),NVAR)
    double precision, intent(out) :: srcQ (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),QVAR)
    double precision, intent(in   ) :: dx, dy, dz, dt
    double precision, intent(inout) :: courno
    integer         , intent(in   ) :: ngp, ngf

    double precision, allocatable:: dpdrho(:,:,:)
    double precision, allocatable:: dpde(:,:,:)
!    double precision, allocatable:: dpdX_er(:,:,:,:)

    integer          :: i, j, k
    integer          :: pt_index(3)
    integer          :: n, nq
    integer          :: iadv, ispec, iaux
    double precision :: courx, coury, courz, courmx, courmy, courmz
    double precision :: kineng

    integer :: ipassive

    type (eos_t) :: eos_state

    allocate( dpdrho(qlo(1):qhi(1),qlo(2):qhi(2),qlo(3):qhi(3)))
    allocate(   dpde(qlo(1):qhi(1),qlo(2):qhi(2),qlo(3):qhi(3)))
!    allocate(dpdX_er(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3,nspec))

    !
    ! Make q (all but p), except put e in slot for rho.e, fix after eos call.
    ! The temperature is used as an initial guess for the eos call and will be overwritten.
    !
    do       k = qlo(3), qhi(3)
       do    j = qlo(2), qhi(2)
          do i = qlo(1), qhi(1)
             
             if (uin(i,j,k,URHO) .le. ZERO) then
                print *,'   '
                print *,'>>> Error: Castro_advection_3d::ctoprim ',i,j,k
                print *,'>>> ... negative density ',uin(i,j,k,URHO)
                call bl_error("Error:: Castro_advection_3d.f90 :: ctoprim")
             end if

             q(i,j,k,QRHO) = uin(i,j,k,URHO)
             q(i,j,k,QU) = uin(i,j,k,UMX)/uin(i,j,k,URHO)
             q(i,j,k,QV) = uin(i,j,k,UMY)/uin(i,j,k,URHO)
             q(i,j,k,QW) = uin(i,j,k,UMZ)/uin(i,j,k,URHO)

             ! Get the internal energy, which we'll use for determining the pressure.
             ! We use a dual energy formalism. If (E - K) < eta1 and eta1 is suitably small, 
             ! then we risk serious numerical truncation error in the internal energy.
             ! Therefore we'll use the result of the separately updated internal energy equation.
             ! Otherwise, we'll set e = E - K.

             kineng = HALF * q(i,j,k,QRHO) * (q(i,j,k,QU)**2 + q(i,j,k,QV)**2 + q(i,j,k,QW)**2)

             if ( (uin(i,j,k,UEDEN) - kineng) / uin(i,j,k,UEDEN) .lt. dual_energy_eta1) then
                q(i,j,k,QREINT) = (uin(i,j,k,UEDEN) - kineng) / q(i,j,k,QRHO)
             else
                q(i,j,k,QREINT) = uin(i,j,k,UEINT) / q(i,j,k,QRHO)
             endif

             q(i,j,k,QTEMP) = uin(i,j,k,UTEMP)
             
             ! convert "rho K" to "K"
             if (QESGS .gt. -1) &
                  q(i,j,k,QESGS) = uin(i,j,k,UESGS)/q(i,j,k,QRHO)

          enddo
       enddo
    enddo

    ! Load passively-advected quatities, c, into q, assuming they 
    ! arrived in uin as rho.c
    do ipassive = 1, npassive
       n = upass_map(ipassive)
       nq = qpass_map(ipassive)
       do       k = qlo(3), qhi(3)
          do    j = qlo(2), qhi(2)
             do i = qlo(1), qhi(1)
                q(i,j,k,nq) = uin(i,j,k,n)/q(i,j,k,QRHO)
             enddo
          enddo
       enddo
    enddo

    ! Get gamc, p, T, c, csml using q state
    do       k = qlo(3), qhi(3)
       do    j = qlo(2), qhi(2)
          do i = qlo(1), qhi(1)
             
             pt_index(:) = (/i, j, k/)

             eos_state % T   = q(i,j,k,QTEMP)
             eos_state % rho = q(i,j,k,QRHO)
             eos_state % xn  = q(i,j,k,QFS:QFS+nspec-1)
             eos_state % aux = q(i,j,k,QFX:QFX+naux-1)             

             ! If necessary, reset the energy using small_temp
             if ((allow_negative_energy .eq. 0) .and. (q(i,j,k,QREINT) .lt. ZERO)) then
                q(i,j,k,QTEMP) = small_temp
                eos_state % T =  q(i,j,k,QTEMP)

                call eos(eos_input_rt, eos_state, .false., pt_index = pt_index)
                q(i,j,k,QREINT) = eos_state % e

                if (q(i,j,k,QREINT) .lt. ZERO) then
                   print *,'   '
                   print *,'>>> Error: Castro_advection_3d::ctoprim ',i,j,k
                   print *,'>>> ... new e from eos (input_rt) call is negative ' &
                        ,q(i,j,k,QREINT)
                   print *,'    '
                   call bl_error("Error:: Castro_advection_3d.f90 :: ctoprim")
                end if
             end if

             eos_state % e = q(i,j,k,QREINT)

             call eos(eos_input_re, eos_state, .false., pt_index = pt_index)

             q(i,j,k,QTEMP)  = eos_state % T
             q(i,j,k,QREINT) = eos_state % e
             q(i,j,k,QPRES)  = eos_state % p

             dpdrho(i,j,k) = eos_state % dpdr_e
             dpde(i,j,k)   = eos_state % dpde
             c(i,j,k)      = eos_state % cs
             gamc(i,j,k)   = eos_state % gam1

             csml(i,j,k) = max(small, small * c(i,j,k))

             ! convert "e" back to "rho e"
             q(i,j,k,QREINT) = q(i,j,k,QREINT)*q(i,j,k,QRHO)

             q(i,j,k,QGAME) = q(i,j,k,QPRES)/q(i,j,k,QREINT) + ONE

          end do
       end do
    end do

    ! compute srcQ terms
    do       k = slo(3), shi(3)
       do    j = slo(2), shi(2)
          do i = slo(1), shi(1)
             
             srcQ(i,j,k,QRHO  ) = src(i,j,k,URHO)
             srcQ(i,j,k,QU    ) = (src(i,j,k,UMX) - q(i,j,k,QU) * srcQ(i,j,k,QRHO)) / q(i,j,k,QRHO)
             srcQ(i,j,k,QV    ) = (src(i,j,k,UMY) - q(i,j,k,QV) * srcQ(i,j,k,QRHO)) / q(i,j,k,QRHO)
             srcQ(i,j,k,QW    ) = (src(i,j,k,UMZ) - q(i,j,k,QW) * srcQ(i,j,k,QRHO)) / q(i,j,k,QRHO)
             srcQ(i,j,k,QREINT) = src(i,j,k,UEDEN) - q(i,j,k,QU)*src(i,j,k,UMX) &
                                                   - q(i,j,k,QV)*src(i,j,k,UMY) &
                                                   - q(i,j,k,QW)*src(i,j,k,UMZ) &
                                    + HALF * (q(i,j,k,QU)**2 + q(i,j,k,QV)**2 + q(i,j,k,QW)**2) * srcQ(i,j,k,QRHO)

             srcQ(i,j,k,QPRES ) = dpde(i,j,k)*(srcQ(i,j,k,QREINT) - &
                  q(i,j,k,QREINT)*srcQ(i,j,k,QRHO)/q(i,j,k,QRHO)) /q(i,j,k,QRHO) + &
                  dpdrho(i,j,k)*srcQ(i,j,k,QRHO)! + &
!                                    sum(dpdX_er(i,j,k,:)*(src(i,j,k,UFS:UFS+nspec-1) - &
!                                                          q(i,j,k,QFS:QFS+nspec-1)*srcQ(i,j,k,QRHO))) &
!                                    /q(i,j,k,QRHO)

             if (QESGS .gt. -1) &
                  srcQ(i,j,k,QESGS) = src(i,j,k,UESGS)/q(i,j,k,QRHO) - q(i,j,k,QESGS) * srcQ(i,j,k,QRHO)

          enddo
       enddo
    enddo

    do ipassive = 1, npassive
       n = upass_map(ipassive)
       nq = qpass_map(ipassive)

       do       k = slo(3), shi(3)
          do    j = slo(2), shi(2)
             do i = slo(1), shi(1)
                srcQ(i,j,k,nq) = ( src(i,j,k,n) - q(i,j,k,nq) * srcQ(i,j,k,QRHO) ) / &
                     q(i,j,k,QRHO)
             enddo
          enddo
       enddo

    enddo

    ! Compute running max of Courant number over grids
    courmx = courno
    courmy = courno
    courmz = courno
    do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)
             
             courx = ( c(i,j,k)+abs(q(i,j,k,QU)) ) * dt/dx
             coury = ( c(i,j,k)+abs(q(i,j,k,QV)) ) * dt/dy
             courz = ( c(i,j,k)+abs(q(i,j,k,QW)) ) * dt/dz
             
             courmx = max( courmx, courx )
             courmy = max( courmy, coury )
             courmz = max( courmz, courz )
             
             if (courx .gt. ONE) then
                print *,'   '
                call bl_warning("Warning:: Castro_advection_3d.f90 :: CFL violation in ctoprim")
                print *,'>>> ... (u+c) * dt / dx > 1 ', courx
                print *,'>>> ... at cell (i,j,k)   : ',i,j,k
                print *,'>>> ... u, c                ',q(i,j,k,QU), c(i,j,k)
                print *,'>>> ... density             ',q(i,j,k,QRHO)
             end if
             
             if (coury .gt. ONE) then
                print *,'   '
                call bl_warning("Warning:: Castro_advection_3d.f90 :: CFL violation in ctoprim")
                print *,'>>> ... (v+c) * dt / dx > 1 ', coury
                print *,'>>> ... at cell (i,j,k)   : ',i,j,k
                print *,'>>> ... v, c                ',q(i,j,k,QV), c(i,j,k)
                print *,'>>> ... density             ',q(i,j,k,QRHO)
             end if
             
             if (courz .gt. ONE) then
                print *,'   '
                call bl_warning("Warning:: Castro_advection_3d.f90 :: CFL violation in ctoprim")
                print *,'>>> ... (w+c) * dt / dx > 1 ', courz
                print *,'>>> ... at cell (i,j,k)   : ',i,j,k
                print *,'>>> ... w, c                ',q(i,j,k,QW), c(i,j,k)
                print *,'>>> ... density             ',q(i,j,k,QRHO)
             end if
             
          enddo
       enddo
    enddo

    courno = max( courmx, courmy, courmz )

    ! Compute flattening coef for slope calculations
    if (use_flattening == 1) then
       call uflaten(lo-ngf,hi+ngf, &
                    q(:,:,:,QPRES), &
                    q(:,:,:,QU), &
                    q(:,:,:,QV), &
                    q(:,:,:,QW), &
                    flatn,qlo,qhi)
    else
       flatn = ONE
    endif

    deallocate(dpdrho,dpde)
    
  end subroutine ctoprim

! ::: 
! ::: ------------------------------------------------------------------
! ::: 

  subroutine consup(lo,hi, &
                    uin,uin_l1,uin_l2,uin_l3,uin_h1,uin_h2,uin_h3, &
                    uout,uout_l1,uout_l2,uout_l3,uout_h1,uout_h2,uout_h3, &
                    src ,src_l1,src_l2,src_l3,src_h1,src_h2,src_h3, &
                    flux1,flux1_l1,flux1_l2,flux1_l3,flux1_h1,flux1_h2,flux1_h3, &
                    flux2,flux2_l1,flux2_l2,flux2_l3,flux2_h1,flux2_h2,flux2_h3, &
                    flux3,flux3_l1,flux3_l2,flux3_l3,flux3_h1,flux3_h2,flux3_h3, &
                    area1,area1_l1,area1_l2,area1_l3,area1_h1,area1_h2,area1_h3, &
                    area2,area2_l1,area2_l2,area2_l3,area2_h1,area2_h2,area2_h3, &
                    area3,area3_l1,area3_l2,area3_l3,area3_h1,area3_h2,area3_h3, &
                    vol,vol_l1,vol_l2,vol_l3,vol_h1,vol_h2,vol_h3, &
                    div, dlo, dhi, &
                    pdivu, plo, phi, &
                    dx,dy,dz,dt,E_added_flux,&
                    xmom_added_flux,ymom_added_flux,zmom_added_flux)

    use network, only : nspec, naux
    use eos_module
    use meth_params_module, only : difmag, NVAR, URHO, UMX, UMY, UMZ, &
         UEDEN, UEINT, UTEMP, normalize_species
    use bl_constants_module

    implicit none

    integer, intent(in) :: lo(3), hi(3), plo(3), phi(3), dlo(3), dhi(3)
    integer, intent(in) :: uin_l1,uin_l2,uin_l3,uin_h1,uin_h2,uin_h3
    integer, intent(in) ::  uout_l1, uout_l2, uout_l3, uout_h1, uout_h2, uout_h3
    integer, intent(in) ::   src_l1,  src_l2,  src_l3,  src_h1,  src_h2,  src_h3 
    integer, intent(in) :: flux1_l1,flux1_l2,flux1_l3,flux1_h1,flux1_h2,flux1_h3
    integer, intent(in) :: flux2_l1,flux2_l2,flux2_l3,flux2_h1,flux2_h2,flux2_h3
    integer, intent(in) :: flux3_l1,flux3_l2,flux3_l3,flux3_h1,flux3_h2,flux3_h3
    integer, intent(in) :: area1_l1,area1_l2,area1_l3,area1_h1,area1_h2,area1_h3
    integer, intent(in) :: area2_l1,area2_l2,area2_l3,area2_h1,area2_h2,area2_h3
    integer, intent(in) :: area3_l1,area3_l2,area3_l3,area3_h1,area3_h2,area3_h3
    integer, intent(in) :: vol_l1,vol_l2,vol_l3,vol_h1,vol_h2,vol_h3

    double precision,intent(in   ):: uin(uin_l1:uin_h1,uin_l2:uin_h2,uin_l3:uin_h3,NVAR)
    double precision,intent(inout):: uout(uout_l1:uout_h1,uout_l2:uout_h2,uout_l3:uout_h3,NVAR)
    double precision,intent(in   )::   src(src_l1:src_h1,src_l2:src_h2,src_l3:src_h3,NVAR)
    double precision,intent(inout):: flux1(flux1_l1:flux1_h1,flux1_l2:flux1_h2,flux1_l3:flux1_h3,NVAR)
    double precision,intent(inout):: flux2(flux2_l1:flux2_h1,flux2_l2:flux2_h2,flux2_l3:flux2_h3,NVAR)
    double precision,intent(inout):: flux3(flux3_l1:flux3_h1,flux3_l2:flux3_h2,flux3_l3:flux3_h3,NVAR)
    double precision,intent(in   ):: area1(area1_l1:area1_h1,area1_l2:area1_h2,area1_l3:area1_h3)
    double precision,intent(in   ):: area2(area2_l1:area2_h1,area2_l2:area2_h2,area2_l3:area2_h3)
    double precision,intent(in   ):: area3(area3_l1:area3_h1,area3_l2:area3_h2,area3_l3:area3_h3)
    double precision,intent(in   ):: vol(vol_l1:vol_h1,vol_l2:vol_h2,vol_l3:vol_h3)
    double precision,intent(in   ):: div(dlo(1):dhi(1),dlo(2):dhi(2),dlo(3):dhi(3))
    double precision,intent(in   ):: pdivu(plo(1):phi(1),plo(2):phi(2),plo(3):phi(3))
    double precision, intent(in) :: dx, dy, dz, dt
    double precision, intent(inout) :: E_added_flux, xmom_added_flux, ymom_added_flux, zmom_added_flux

    double precision :: div1
    integer          :: i, j, k, n

    do n = 1, NVAR
         
       if ( n.eq.UTEMP ) then
          
          flux1(:,:,:,n) = ZERO
          flux2(:,:,:,n) = ZERO
          flux3(:,:,:,n) = ZERO
          
       else

          do k = lo(3),hi(3)
             do j = lo(2),hi(2)
                do i = lo(1),hi(1)+1
                   div1 = FOURTH*(div(i,j,k) + div(i,j+1,k) + div(i,j,k+1) + div(i,j+1,k+1))
                   div1 = difmag*min(ZERO,div1)
                   flux1(i,j,k,n) = flux1(i,j,k,n) + dx*div1*(uin(i,j,k,n)-uin(i-1,j,k,n))
                   flux1(i,j,k,n) = flux1(i,j,k,n) * area1(i,j,k) * dt
                enddo
             enddo
          enddo

          do k = lo(3),hi(3)
             do j = lo(2),hi(2)+1
                do i = lo(1),hi(1)
                   div1 = FOURTH*(div(i,j,k) + div(i+1,j,k) + div(i,j,k+1) + div(i+1,j,k+1))
                   div1 = difmag*min(ZERO,div1)
                   flux2(i,j,k,n) = flux2(i,j,k,n) + dy*div1*(uin(i,j,k,n)-uin(i,j-1,k,n))
                   flux2(i,j,k,n) = flux2(i,j,k,n) * area2(i,j,k) * dt
                enddo
             enddo
          enddo

          do k = lo(3),hi(3)+1
             do j = lo(2),hi(2)
                do i = lo(1),hi(1)
                   div1 = FOURTH*(div(i,j,k) + div(i+1,j,k) + div(i,j+1,k) + div(i+1,j+1,k))
                   div1 = difmag*min(ZERO,div1)
                   flux3(i,j,k,n) = flux3(i,j,k,n) + dz*div1*(uin(i,j,k,n)-uin(i,j,k-1,n))
                   flux3(i,j,k,n) = flux3(i,j,k,n) * area3(i,j,k) * dt
                enddo
             enddo
          enddo
          
       endif

    enddo

    if (normalize_species .eq. 1) &
         call normalize_species_fluxes( &
                  flux1,flux1_l1,flux1_l2,flux1_l3,flux1_h1,flux1_h2,flux1_h3, &
                  flux2,flux2_l1,flux2_l2,flux2_l3,flux2_h1,flux2_h2,flux2_h3, &
                  flux3,flux3_l1,flux3_l2,flux3_l3,flux3_h1,flux3_h2,flux3_h3, &
                  lo,hi)

    do n = 1, NVAR

       ! pass temperature through
       if (n .eq. UTEMP) then
          do k = lo(3),hi(3)
             do j = lo(2),hi(2)
                do i = lo(1),hi(1)
                   uout(i,j,k,n) = uin(i,j,k,n)
                enddo
             enddo
          enddo
       else 
          ! update everything else with fluxes and source terms
          do k = lo(3),hi(3)
             do j = lo(2),hi(2)
                do i = lo(1),hi(1)
                   uout(i,j,k,n) = uin(i,j,k,n) &
                          + ( flux1(i,j,k,n) - flux1(i+1,j,k,n) &
                          +   flux2(i,j,k,n) - flux2(i,j+1,k,n) &
                          +   flux3(i,j,k,n) - flux3(i,j,k+1,n)) / vol(i,j,k) &
                          +   dt * src(i,j,k,n)
                   !
                   ! Add the source term to (rho e)
                   !
                   if (n .eq. UEINT) then
                      uout(i,j,k,n) = uout(i,j,k,n) - dt * pdivu(i,j,k)
                   else if (n .eq. UMX) then
                      xmom_added_flux = xmom_added_flux + &
                           ( flux1(i,j,k,n) - flux1(i+1,j,k,n) &
                         +   flux2(i,j,k,n) - flux2(i,j+1,k,n) &
                         +   flux3(i,j,k,n) - flux3(i,j,k+1,n)) / vol(i,j,k)
                   else if (n .eq. UMY) then
                      ymom_added_flux = ymom_added_flux + &
                           ( flux1(i,j,k,n) - flux1(i+1,j,k,n) &
                         +   flux2(i,j,k,n) - flux2(i,j+1,k,n) &
                         +   flux3(i,j,k,n) - flux3(i,j,k+1,n)) / vol(i,j,k)
                   else if (n .eq. UMZ) then
                      zmom_added_flux = zmom_added_flux + &
                           ( flux1(i,j,k,n) - flux1(i+1,j,k,n) &
                         +   flux2(i,j,k,n) - flux2(i,j+1,k,n) &
                         +   flux3(i,j,k,n) - flux3(i,j,k+1,n)) / vol(i,j,k)
                   else if (n .eq. UEDEN) then
                      E_added_flux = E_added_flux + &
                           ( flux1(i,j,k,n) - flux1(i+1,j,k,n) &
                         +   flux2(i,j,k,n) - flux2(i,j+1,k,n) &
                         +   flux3(i,j,k,n) - flux3(i,j,k+1,n)) / vol(i,j,k) 
                   endif
                enddo
             enddo
          enddo
       endif
         
    enddo

  end subroutine consup

! ::: 
! ::: ------------------------------------------------------------------
! ::: 

  subroutine divu(lo,hi,q,qlo,qhi,dx,dy,dz,div,dlo,dhi)
    
    use meth_params_module, only : QU, QV, QW, QVAR
    use bl_constants_module
    
    implicit none

    integer, intent(in) :: lo(3),hi(3),qlo(3),qhi(3),dlo(3),dhi(3)
    double precision, intent(in ) :: dx, dy, dz
    double precision, intent(out) :: div(dlo(1):dhi(1),dlo(2):dhi(2),dlo(3):dhi(3))
    double precision, intent(in ) :: q  (qlo(1):qhi(1),qlo(2):qhi(2),qlo(3):qhi(3),QVAR)

    integer          :: i, j, k
    double precision :: ux, vy, wz

    do k=lo(3),hi(3)+1
       do j=lo(2),hi(2)+1
          do i=lo(1),hi(1)+1
             
             ux = FOURTH*( &
                    + q(i  ,j  ,k  ,QU) - q(i-1,j  ,k  ,QU) &
                    + q(i  ,j  ,k-1,QU) - q(i-1,j  ,k-1,QU) &
                    + q(i  ,j-1,k  ,QU) - q(i-1,j-1,k  ,QU) &
                    + q(i  ,j-1,k-1,QU) - q(i-1,j-1,k-1,QU) )/ dx

             vy = FOURTH*( &
                    + q(i  ,j  ,k  ,QV) - q(i  ,j-1,k  ,QV) &
                    + q(i  ,j  ,k-1,QV) - q(i  ,j-1,k-1,QV) &
                    + q(i-1,j  ,k  ,QV) - q(i-1,j-1,k  ,QV) &
                    + q(i-1,j  ,k-1,QV) - q(i-1,j-1,k-1,QV) )/ dy

             wz = FOURTH*( &
                    + q(i  ,j  ,k  ,QW) - q(i  ,j  ,k-1,QW) &
                    + q(i  ,j-1,k  ,QW) - q(i  ,j-1,k-1,QW) &
                    + q(i-1,j  ,k  ,QW) - q(i-1,j  ,k-1,QW) &
                    + q(i-1,j-1,k  ,QW) - q(i-1,j-1,k-1,QW) )/ dz

             div(i,j,k) = ux + vy + wz

          enddo
       enddo
    enddo
    
  end subroutine divu

! ::
! :: ----------------------------------------------------------
! ::

  subroutine normalize_species_fluxes(flux1,flux1_l1,flux1_l2,flux1_l3, &
                                      flux1_h1,flux1_h2,flux1_h3, &
                                      flux2,flux2_l1,flux2_l2,flux2_l3, &
                                      flux2_h1,flux2_h2,flux2_h3, &
                                      flux3,flux3_l1,flux3_l2,flux3_l3, &
                                      flux3_h1,flux3_h2,flux3_h3, &
                                      lo,hi)
    
    use network, only : nspec
    use meth_params_module, only : NVAR, URHO, UFS
    use bl_constants_module

    implicit none

    integer          :: lo(3),hi(3)
    integer          :: flux1_l1,flux1_l2,flux1_l3,flux1_h1,flux1_h2,flux1_h3
    integer          :: flux2_l1,flux2_l2,flux2_l3,flux2_h1,flux2_h2,flux2_h3
    integer          :: flux3_l1,flux3_l2,flux3_l3,flux3_h1,flux3_h2,flux3_h3
    double precision :: flux1(flux1_l1:flux1_h1,flux1_l2:flux1_h2,flux1_l3:flux1_h3,NVAR)
    double precision :: flux2(flux2_l1:flux2_h1,flux2_l2:flux2_h2,flux2_l3:flux2_h3,NVAR)
    double precision :: flux3(flux3_l1:flux3_h1,flux3_l2:flux3_h2,flux3_l3:flux3_h3,NVAR)
    
    ! Local variables
    integer          :: i,j,k,n
    double precision :: sum,fac
    
    do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)+1
             sum = ZERO
             do n = UFS, UFS+nspec-1
                sum = sum + flux1(i,j,k,n)
             end do
             if (sum .ne. ZERO) then
                fac = flux1(i,j,k,URHO) / sum
             else
                fac = ONE
             end if
             do n = UFS, UFS+nspec-1
                flux1(i,j,k,n) = flux1(i,j,k,n) * fac
             end do
          end do
       end do
    end do

    do k = lo(3),hi(3)
       do j = lo(2),hi(2)+1
          do i = lo(1),hi(1)
             sum = ZERO
             do n = UFS, UFS+nspec-1
                sum = sum + flux2(i,j,k,n)
             end do
             if (sum .ne. ZERO) then
                fac = flux2(i,j,k,URHO) / sum
             else
                fac = ONE
             end if
             do n = UFS, UFS+nspec-1
                flux2(i,j,k,n) = flux2(i,j,k,n) * fac
             end do
          end do
       end do
    end do

    do k = lo(3),hi(3)+1
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)
             sum = ZERO
             do n = UFS, UFS+nspec-1
                sum = sum + flux3(i,j,k,n)
             end do
             if (sum .ne. ZERO) then
                fac = flux3(i,j,k,URHO) / sum
             else
                fac = ONE
             end if
             do n = UFS, UFS+nspec-1
                flux3(i,j,k,n) = flux3(i,j,k,n) * fac
             end do
          end do
       end do
    end do

  end subroutine normalize_species_fluxes

! ::
! :: ----------------------------------------------------------
! ::

  subroutine enforce_minimum_density(uin,uin_l1,uin_l2,uin_l3,uin_h1,uin_h2,uin_h3, &
                                     uout,uout_l1,uout_l2,uout_l3, &
                                     uout_h1,uout_h2,uout_h3, &
                                     lo,hi,mass_added,eint_added,eden_added,verbose)
    
    use network, only : nspec, naux
    use meth_params_module, only : NVAR, URHO, UMX, UMY, UMZ, UTEMP, UEDEN, UEINT, UFS, UFX, &
                                     UFA, small_dens, small_temp, nadv
    use bl_constants_module
    use eos_type_module, only : eos_t
    use eos_module, only : eos
    use eos_data_module, only : eos_input_rt

    implicit none

    integer, intent(in) :: lo(3), hi(3), verbose
    integer, intent(in) ::  uin_l1,  uin_l2,  uin_l3,  uin_h1,  uin_h2,  uin_h3
    integer, intent(in) :: uout_l1, uout_l2, uout_l3, uout_h1, uout_h2, uout_h3
    double precision,intent(in   )::  uin( uin_l1: uin_h1, uin_l2: uin_h2, uin_l3: uin_h3,NVAR)
    double precision,intent(inout):: uout(uout_l1:uout_h1,uout_l2:uout_h2,uout_l3:uout_h3,NVAR)
    double precision,intent(inout):: mass_added, eint_added, eden_added
    
    ! Local variables
    integer          :: i,ii,j,jj,k,kk,n
    integer          :: i_set, j_set, k_set
    double precision :: max_dens
    
    double precision :: initial_mass, final_mass
    double precision :: initial_eint, final_eint
    double precision :: initial_eden, final_eden

    type (eos_t) :: eos_state
    
    initial_mass = ZERO
      final_mass = ZERO

    initial_eint = ZERO
      final_eint = ZERO

    initial_eden = ZERO
      final_eden = ZERO

    max_dens = ZERO

    do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)
             
             initial_mass = initial_mass + uout(i,j,k,URHO )
             initial_eint = initial_eint + uout(i,j,k,UEINT)
             initial_eden = initial_eden + uout(i,j,k,UEDEN)
             
             if (uout(i,j,k,URHO) .eq. ZERO) then
                
                print *,'DENSITY EXACTLY ZERO AT CELL ',i,j,k
                print *,'  in grid ',lo(1),lo(2),lo(3),hi(1),hi(2),hi(3)
                call bl_error("Error:: Castro_3d.f90 :: enforce_minimum_density")
                
             else if (uout(i,j,k,URHO) < small_dens) then
                
                max_dens = uout(i,j,k,URHO)
                i_set = i
                j_set = j
                k_set = k
                do kk = -1,1
                   do jj = -1,1
                      do ii = -1,1
                         if (i+ii.ge.lo(1) .and. j+jj.ge.lo(2) .and. k+kk.ge.lo(3) .and. &
                             i+ii.le.hi(1) .and. j+jj.le.hi(2) .and. k+kk.le.hi(3)) then
                              if (uout(i+ii,j+jj,k+kk,URHO) .gt. max_dens) then
                                  i_set = i+ii
                                  j_set = j+jj
                                  k_set = k+kk
                                  max_dens = uout(i_set,j_set,k_set,URHO)
                              endif
                         endif
                      end do
                   end do
                end do

                ! If no neighboring zones are above small_dens, our only recourse 
                ! is to set the density equal to small_dens, and the temperature 
                ! equal to small_temp. We set the velocities to zero, 
                ! though any choice here would be arbitrary.

                if (max_dens < small_dens) then

                   do n = UFS, UFS+nspec-1
                      uout(i,j,k,n) = uout(i_set,j_set,k_set,n) * (small_dens / uout(i,j,k,URHO))
                   end do
                   do n = UFX, UFX+naux-1
                      uout(i,j,k,n) = uout(i_set,j_set,k_set,n) * (small_dens / uout(i,j,k,URHO))
                   end do
                   do n = UFA, UFA+nadv-1
                      uout(i,j,k,n) = uout(i_set,j_set,k_set,n) * (small_dens / uout(i,j,k,URHO))
                   end do

                   eos_state % rho = small_dens
                   eos_state % T   = small_temp
                   eos_state % xn  = uout(i,j,k,UFS:UFS+nspec-1) / uout(i,j,k,URHO)

                   call eos(eos_input_rt, eos_state)

                   uout(i,j,k,URHO ) = eos_state % rho
                   uout(i,j,k,UTEMP) = eos_state % T

                   uout(i,j,k,UMX  ) = ZERO
                   uout(i,j,k,UMY  ) = ZERO
                   uout(i,j,k,UMZ  ) = ZERO

                   uout(i,j,k,UEINT) = eos_state % rho * eos_state % e
                   uout(i,j,k,UEDEN) = uout(i,j,k,UEINT)

                endif
                
                if (verbose .gt. 0) then
                   if (uout(i,j,k,URHO) < ZERO) then
                      print *,'   '
                      print *,'>>> RESETTING NEG.  DENSITY AT ',i,j,k
                      print *,'>>> FROM ',uout(i,j,k,URHO),' TO ',uout(i_set,j_set,k_set,URHO)
                      print *,'   '
                   else
                      print *,'   '
                      print *,'>>> RESETTING SMALL DENSITY AT ',i,j,k
                      print *,'>>> FROM ',uout(i,j,k,URHO),' TO ',uout(i_set,j_set,k_set,URHO)
                      print *,'   '
                   end if
                end if
                
                uout(i,j,k,URHO ) = uout(i_set,j_set,k_set,URHO )
                uout(i,j,k,UTEMP) = uout(i_set,j_set,k_set,UTEMP)
                uout(i,j,k,UEINT) = uout(i_set,j_set,k_set,UEINT)
                uout(i,j,k,UEDEN) = uout(i_set,j_set,k_set,UEDEN)
                uout(i,j,k,UMX  ) = uout(i_set,j_set,k_set,UMX  )
                uout(i,j,k,UMY  ) = uout(i_set,j_set,k_set,UMY  )
                uout(i,j,k,UMZ  ) = uout(i_set,j_set,k_set,UMZ  )
   
                do n = UFS, UFS+nspec-1
                   uout(i,j,k,n) = uout(i_set,j_set,k_set,n)
                end do
                do n = UFX, UFX+naux-1
                   uout(i,j,k,n) = uout(i_set,j_set,k_set,n)
                end do
                do n = UFA, UFA+nadv-1
                   uout(i,j,k,n) = uout(i_set,j_set,k_set,n)
                end do
                
             end if

             final_mass = final_mass + uout(i,j,k,URHO )
             final_eint = final_eint + uout(i,j,k,UEINT)
             final_eden = final_eden + uout(i,j,k,UEDEN)                
             
          enddo
       enddo
    enddo

    if ( max_dens /= ZERO ) then
       mass_added = mass_added + final_mass - initial_mass
       eint_added = eint_added + final_eint - initial_eint
       eden_added = eden_added + final_eden - initial_eden
    endif
    
  end subroutine enforce_minimum_density

! :::
! ::: ------------------------------------------------------------------
! :::

  subroutine normalize_new_species(u,u_l1,u_l2,u_l3,u_h1,u_h2,u_h3,lo,hi)

    use network, only : nspec
    use meth_params_module, only : NVAR, URHO, UFS
    use bl_constants_module

    implicit none

    integer          :: lo(3), hi(3)
    integer          :: u_l1,u_l2,u_l3,u_h1,u_h2,u_h3
    double precision :: u(u_l1:u_h1,u_l2:u_h2,u_l3:u_h3,NVAR)
    
    ! Local variables
    integer          :: i,j,k,n
    double precision :: fac,sum
    
    do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)
             sum = ZERO
             do n = UFS, UFS+nspec-1
                sum = sum + u(i,j,k,n)
             end do
             if (sum .ne. ZERO) then
                fac = u(i,j,k,URHO) / sum
             else
                fac = ONE
             end if
             do n = UFS, UFS+nspec-1
                u(i,j,k,n) = u(i,j,k,n) * fac
             end do
          end do
       end do
    end do
    
  end subroutine normalize_new_species

end module advection_module
