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
    double precision :: dxinv, dyinv, dzinv
    double precision :: dtdx, dtdy, dtdz, hdt
    double precision :: cdtdx, cdtdy, cdtdz
    double precision :: hdtdx, hdtdy, hdtdz

    double precision, allocatable::dqx(:,:,:,:), dqy(:,:,:,:), dqz(:,:,:,:)

    double precision, allocatable:: Ip(:,:,:,:,:,:), Im(:,:,:,:,:,:)
    double precision, allocatable:: Ip_g(:,:,:,:,:,:), Im_g(:,:,:,:,:,:)
    double precision, allocatable:: Ip_r(:,:,:,:,:,:), Im_r(:,:,:,:,:,:)
    double precision, allocatable:: Ip_gc(:,:,:,:,:,:), Im_gc(:,:,:,:,:,:)

    double precision, allocatable :: shk(:,:,:)

    ! Left and right state arrays
    double precision, dimension(:,:,:,:), pointer :: &
         qxm,qym,qzm,qxp,qyp,qzp, ql, qr, &
         qmxy,qpxy,qmxz,qpxz,qmyx,qpyx,qmyz,qpyz,qmzx,qpzx,qmzy,qpzy, &
         qxl,qxr,qyl,qyr,qzl,qzr
    
    double precision, dimension(:,:,:,:), pointer:: &
         fx,fy,fz,fxy,fxz,fyx,fyz,fzx,fzy

    double precision, dimension(:,:,:), pointer:: &
         pgdnvx,ugdnvx,gegdnvx, &
         pgdnvy,ugdnvy,gegdnvy, &
         pgdnvz,ugdnvz,gegdnvz, &
         pgdnvxy,ugdnvxy,gegdnvxy, &
         pgdnvxz,ugdnvxz,gegdnvxz, &
         pgdnvyx,ugdnvyx,gegdnvyx, &
         pgdnvyz,ugdnvyz,gegdnvyz, &
         pgdnvzx,ugdnvzx,gegdnvzx, &
         pgdnvzy,ugdnvzy,gegdnvzy
    
    double precision, dimension(:,:,:,:), pointer :: ftmp1,ftmp2
    double precision, dimension(:,:,:), pointer :: pgdnvtmp1,ugdnvtmp1,gegdnvtmp1, &
         &                                         pgdnvtmp2,ugdnvtmp2,gegdnvtmp2
    
    type (eos_t) :: eos_state

    integer :: fglo(3), fghi(3), glo(3), ghi(3)

    ! Local constants
    dxinv = ONE/dx
    dyinv = ONE/dy
    dzinv = ONE/dz
    dtdx = dt*dxinv
    dtdy = dt*dyinv
    dtdz = dt*dzinv
    hdt = HALF*dt
    hdtdx = HALF*dtdx
    hdtdy = HALF*dtdy
    hdtdz = HALF*dtdz
    cdtdx = dtdx*THIRD
    cdtdy = dtdy*THIRD
    cdtdz = dtdz*THIRD

    fglo = lo-1  ! face + one ghost
    fghi = hi+2  

    glo = lo-1  ! one ghost,  this can be used for face-based arrays too
    ghi = hi+1  
    
    allocate ( qxm(fglo(1):fghi(1),fglo(2):fghi(2),fglo(3):fghi(3),QVAR))
    allocate ( qxp(fglo(1):fghi(1),fglo(2):fghi(2),fglo(3):fghi(3),QVAR))

    allocate ( qym(fglo(1):fghi(1),fglo(2):fghi(2),fglo(3):fghi(3),QVAR))
    allocate ( qyp(fglo(1):fghi(1),fglo(2):fghi(2),fglo(3):fghi(3),QVAR))

    allocate ( qzm(fglo(1):fghi(1),fglo(2):fghi(2),fglo(3):fghi(3),QVAR))
    allocate ( qzp(fglo(1):fghi(1),fglo(2):fghi(2),fglo(3):fghi(3),QVAR))

    ! for the hybrid Riemann solver
    allocate(shk(glo(1):ghi(1),glo(2):ghi(2),glo(3):ghi(3)))
    
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
               Ip_gc(:,:,:,:,:,1),Im_gc(:,:,:,:,:,1), glo, ghi, &
               lo,hi,dx,dy,dz,dt)
       else
          ! temperature-based PPM
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
                      Ip,Im,Ip_g,Im_g,Ip_r,Im_r,Ip_gc,Im_gc, glo, ghi, &
                      qxm,qxp,qym,qyp,qzm,qzp, fglo, fghi, &
                      lo,hi,dt)
       
       deallocate(Ip,Im,Ip_g,Im_g,Ip_r,Im_r,Ip_gc,Im_gc)

    else

       allocate ( dqx(fglo(1):fghi(1),fglo(2):fghi(2),fglo(3):fghi(3),QVAR))
       allocate ( dqy(fglo(1):fghi(1),fglo(2):fghi(2),fglo(3):fghi(3),QVAR))
       allocate ( dqz(fglo(1):fghi(1),fglo(2):fghi(2),fglo(3):fghi(3),QVAR))
           
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
       
       deallocate(dqx,dqy,dqz)

    end if

    allocate ( qmxy(fglo(1):fghi(1),fglo(2):fghi(2),fglo(3):fghi(3),QVAR))
    allocate ( qpxy(fglo(1):fghi(1),fglo(2):fghi(2),fglo(3):fghi(3),QVAR))

    allocate ( qmxz(fglo(1):fghi(1),fglo(2):fghi(2),fglo(3):fghi(3),QVAR))
    allocate ( qpxz(fglo(1):fghi(1),fglo(2):fghi(2),fglo(3):fghi(3),QVAR))

    allocate ( qmyx(fglo(1):fghi(1),fglo(2):fghi(2),fglo(3):fghi(3),QVAR))
    allocate ( qpyx(fglo(1):fghi(1),fglo(2):fghi(2),fglo(3):fghi(3),QVAR))

    allocate ( qmyz(fglo(1):fghi(1),fglo(2):fghi(2),fglo(3):fghi(3),QVAR))
    allocate ( qpyz(fglo(1):fghi(1),fglo(2):fghi(2),fglo(3):fghi(3),QVAR))

    allocate ( qmzx(fglo(1):fghi(1),fglo(2):fghi(2),fglo(3):fghi(3),QVAR))
    allocate ( qpzx(fglo(1):fghi(1),fglo(2):fghi(2),fglo(3):fghi(3),QVAR))

    allocate ( qmzy(fglo(1):fghi(1),fglo(2):fghi(2),fglo(3):fghi(3),QVAR))
    allocate ( qpzy(fglo(1):fghi(1),fglo(2):fghi(2),fglo(3):fghi(3),QVAR))

    allocate ( ql(fglo(1):fghi(1),fglo(2):fghi(2),fglo(3):fghi(3),QVAR))
    allocate ( qr(fglo(1):fghi(1),fglo(2):fghi(2),fglo(3):fghi(3),QVAR))

    allocate ( ftmp1(glo(1):ghi(1),glo(2):ghi(2),glo(3):ghi(3),NVAR))
    allocate ( ftmp2(glo(1):ghi(1),glo(2):ghi(2),glo(3):ghi(3),NVAR))

    allocate ( pgdnvtmp1(fglo(1):fghi(1),fglo(2):fghi(2),fglo(3):fghi(3)))
    allocate ( ugdnvtmp1(fglo(1):fghi(1),fglo(2):fghi(2),fglo(3):fghi(3)))
    allocate (gegdnvtmp1(fglo(1):fghi(1),fglo(2):fghi(2),fglo(3):fghi(3)))

    allocate ( pgdnvtmp2(fglo(1):fghi(1),fglo(2):fghi(2),fglo(3):fghi(3)))
    allocate ( ugdnvtmp2(fglo(1):fghi(1),fglo(2):fghi(2),fglo(3):fghi(3)))
    allocate (gegdnvtmp2(fglo(1):fghi(1),fglo(2):fghi(2),fglo(3):fghi(3)))

    !-------------------------------------------------------------------------!
    ! Some notes on the work index (i.e., lo and hi arguments near the end    !
    !                               of the argument list).                    !
    ! * For cmpflx, we use face index in the flux direction and cell-centered !
    !   index for others.                                                     !
    ! * For trans*, we use cell-centered index of the valid region.           !
    !-------------------------------------------------------------------------!
    
    fx      =>      ftmp1
    ugdnvx  =>  ugdnvtmp1
    pgdnvx  =>  pgdnvtmp1
    gegdnvx => gegdnvtmp1

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
    ! Outputs: qmyx, qpyx                  : yface, +-0 at x, +-1 at z
    !          qmzx, qpzx                  : zface, +-0 at x, +-1 at y
    call transx(qym,qmyx,qyp,qpyx, &
                qzm,qmzx,qzp,qpzx, fglo, fghi, &
                fx, glo, ghi, &
                ugdnvx,pgdnvx,gegdnvx, fglo, fghi, &
                gamc, qlo, qhi, &
                cdtdx, lo, hi)

    nullify(fx,ugdnvx,pgdnvx,gegdnvx)

    fy      =>      ftmp1
    ugdnvy  =>  ugdnvtmp1
    pgdnvy  =>  pgdnvtmp1
    gegdnvy => gegdnvtmp1

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
    ! Outputs: qmxy, qpxy                  : xface, +-0 at y, +-1 at z
    !          qmzy, qpzy                  : zface, +-0 at y, +-1 at x
    call transy(qxm,qmxy,qxp,qpxy, &
                qzm,qmzy,qzp,qpzy, fglo, fghi, &
                fy, glo, ghi, &
                ugdnvy,pgdnvy,gegdnvy, fglo, fghi, &
                gamc, qlo, qhi, &
                cdtdy, lo, hi)

    nullify(fy,ugdnvy,pgdnvy,gegdnvy)

    fz      =>      ftmp1
    ugdnvz  =>  ugdnvtmp1
    pgdnvz  =>  pgdnvtmp1
    gegdnvz => gegdnvtmp1

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

    ! Inputs: qxm, qxp                     : xface, +-1 at y & z
    !         qym, qyp                     : yface, +-1 at x & z
    !         fz, ugdnvz, pgdnvz, gegdnvz  : zface, +-1 at x & y
    !         gamc                         : +-4
    ! Outputs: qmxz, qpxz                  : xface, +-0 at z, +-1 at y
    !          qmyz, qpyz                  : yface, +-0 at z, +-1 at x
    call transz(qxm,qmxz,qxp,qpxz, &
                qym,qmyz,qyp,qpyz, fglo, fghi, &
                fz, glo, ghi, &
                ugdnvz,pgdnvz,gegdnvz, fglo, fghi, &
                gamc, qlo, qhi, &
                cdtdz, lo, hi)

    nullify(fz,ugdnvz,pgdnvz,gegdnvz)

    ! We now have qx?, qy?, qz?
    !         and q?zx, q?yx, q?zy, q?xy, q?yz, q?xz

    !
    ! Use qx?, q?yz, q?zy to compute final x-flux
    !

    fyz      =>      ftmp1
    ugdnvyz  =>  ugdnvtmp1
    pgdnvyz  =>  pgdnvtmp1
    gegdnvyz => gegdnvtmp1

    ! Inputs: qmyz, qpyz                       : yface, +-1 at x, +-0 at z
    !         gamc, csml, c                    : +-4
    !         shk                              : +-1
    ! Outputs: fyz, ugdnvyz, pgdnvyz, gegdnvyz : yface, +-1 at x, +-0 at z
    call cmpflx(qmyz,qpyz, fglo, fghi, &
                fyz, glo, ghi, &
                ugdnvyz,pgdnvyz,gegdnvyz, fglo, fghi, &
                gamc,csml,c, qlo, qhi, &
                shk, glo, ghi, &
                2, (/lo(1)-1,lo(2),lo(3)/), (/hi(1)+1,hi(2)+1,hi(3)/), domlo, domhi)

    fzy      =>      ftmp2
    ugdnvzy  =>  ugdnvtmp2
    pgdnvzy  =>  pgdnvtmp2
    gegdnvzy => gegdnvtmp2

    ! Inputs: qmzy, qpzy                       : zface, +-1 at x, +-0 at y
    !         gamc, csml, c                    : +-4
    !         shk                              : +-1
    ! Outputs: fzy, ugdnvzy, pgdnvzy, gegdnvzy : zface, +-1 at x, +-0 at y
    call cmpflx(qmzy,qpzy, fglo, fghi, &
                fzy, glo, ghi, &
                ugdnvzy,pgdnvzy,gegdnvzy, fglo, fghi, &
                gamc,csml,c, qlo, qhi, &
                shk, glo, ghi, &
                3, (/lo(1)-1,lo(2),lo(3)/), (/hi(1)+1,hi(2),hi(3)+1/), domlo, domhi)

    qxl => ql
    qxr => qr
 
    ! Inputs: qxm, qxp                        : xface, +-1 at y & z
    !         fyz, ugdnvyz, pgdnvyz, gegdnvyz : yface, +-1 at x, +-0 at z
    !         fzy, ugdnvzy, pgdnvzy, gegdnvzy : zface, +-1 at x, +-0 at y
    !         gamc, grav, rot                 : +-4
    !         srcQ                            : +-1
    ! Outputs: qxl, qxr                       : xface, +-0 at y & z
    call transyz(qxm,qxl,qxp,qxr, fglo, fghi, &
                 fyz, fzy, glo, ghi, &
                 ugdnvyz,pgdnvyz,gegdnvyz,ugdnvzy,pgdnvzy,gegdnvzy, fglo, fghi, &
                 gamc, qlo, qhi, &
                 srcQ, slo, shi, &
                 grav, gvlo, gvhi, &
                 rot, rlo, rhi, &
                 hdt,hdtdy,hdtdz,lo,hi)

    nullify(fyz, ugdnvyz, pgdnvyz, gegdnvyz)
    nullify(fzy, ugdnvzy, pgdnvzy, gegdnvzy)

    ugdnvx  =>  ugdnvtmp1
    pgdnvx  =>  pgdnvtmp1
    gegdnvx => gegdnvtmp1

    ! Inputs: qxl, qxr                        : xface, +-0 at y & z
    !         gamc, csml, c                   : +-4
    !         shk                             : +-1
    ! Outputs: flux1, ugdnvx, pgdnvx, gegdnvx : xface, +-0 at y & z
    call cmpflx(qxl,qxr, fglo, fghi, &
                flux1, f1lo, f1hi, &
                ugdnvx,pgdnvx,gegdnvx, fglo, fghi, &
                gamc,csml,c, qlo, qhi, &
                shk, glo, ghi, &
                1, lo, (/hi(1)+1,hi(2),hi(3)/), domlo, domhi)

    do       k = lo(3), hi(3)
       do    j = lo(2), hi(2)
          do i = lo(1), hi(1)
             pdivu(i,j,k) = HALF*(pgdnvx(i,j,k)+pgdnvx(i+1,j,k)) * &
                                 (ugdnvx(i+1,j,k)-ugdnvx(i,j,k)) * dxinv
          end do
          do i = lo(1), hi(1)+1
             ugdx(i,j,k) = ugdx(i,j,k)
          end do
       end do
    end do

    nullify(ugdnvx,pgdnvx,gegdnvx)
    nullify(qxl,qxr)

    !
    ! Use qy?, q?zx, q?xz to compute final y-flux
    !

    fzx      =>      ftmp1
    ugdnvzx  =>  ugdnvtmp1
    pgdnvzx  =>  pgdnvtmp1
    gegdnvzx => gegdnvtmp1

    ! Inputs: qmzx, qpzx                       : zface, +-0 at x, +-1 at y
    !         gamc, csml, c                    : +-4
    !         shk                              : +-1
    ! Outputs: fzx, ugdnvzx, pgdnvzx, gegdnvzx : zface, +-0 at x, +-1 at y
    call cmpflx(qmzx,qpzx, fglo, fghi, &
                fzx, glo, ghi, &
                ugdnvzx,pgdnvzx,gegdnvzx, fglo, fghi, &
                gamc,csml,c, qlo, qhi, &
                shk, glo, ghi, &
                3, (/lo(1),lo(2)-1,lo(3)/), (/hi(1),hi(2)+1,hi(3)+1/), domlo, domhi)

    fxz      =>      ftmp2
    ugdnvxz  =>  ugdnvtmp2
    pgdnvxz  =>  pgdnvtmp2
    gegdnvxz => gegdnvtmp2

    ! Inputs: qmxz, qpxz                       : xface, +-1 at y, +-0 at z
    !         gamc, csml, c                    : +-4
    !         shk                              : +-1
    ! Outputs: fxz, ugdnvxz, pgdnvxz, gegdnvxz : xface, +-1 at y, +-0 at z
    call cmpflx(qmxz,qpxz, fglo, fghi, &
                fxz, glo, ghi, &
                ugdnvxz,pgdnvxz,gegdnvxz, fglo, fghi, &
                gamc,csml,c, qlo, qhi, &
                shk, glo, ghi, &
                1, (/lo(1),lo(2)-1,lo(3)/), (/hi(1)+1,hi(2)+1,hi(3)/), domlo, domhi)

    qyl => ql
    qyr => qr
 
    ! Inputs: qym, qyp                        : yface, +-1 at x & z
    !         fxz, ugdnvxz, pgdnvxz, gegdnvxz : xface, +-1 at y, +-0 at z
    !         fzx, ugdnvzx, pgdnvzx, gegdnvzx : zface, +-0 at x, +-1 at y
    !         gamc, grav, rot                 : +-4
    !         srcQ                            : +-1
    ! Outputs: qyl, qyr                       : yface, +-0 at x & z
    call transxz(qym,qyl,qyp,qyr, fglo, fghi, &
                 fxz, fzx, glo, ghi, &
                 ugdnvxz,pgdnvxz,gegdnvxz,ugdnvzx,pgdnvzx,gegdnvzx, fglo, fghi, &
                 gamc, qlo, qhi, &
                 srcQ, slo, shi, &
                 grav, gvlo, gvhi, &
                 rot, rlo, rhi, &
                 hdt,hdtdx,hdtdz,lo,hi)

    nullify(fzx, ugdnvzx, pgdnvzx, gegdnvzx)
    nullify(fxz, ugdnvxz, pgdnvxz, gegdnvxz)

    ugdnvy  =>  ugdnvtmp1
    pgdnvy  =>  pgdnvtmp1
    gegdnvy => gegdnvtmp1

    ! Inputs: qyl, qyr                        : yface, +-0 at x & y
    !         gamc, csml, c                   : +-4
    !         shk                             : +-1
    ! Outputs: flux2, ugdnvy, pgdnvy, gegdnvy : yface, +-0 at x & y
    call cmpflx(qyl,qyr, fglo, fghi, &
                flux2, f2lo, f2hi, &
                ugdnvy,pgdnvy,gegdnvy, fglo, fghi, &
                gamc,csml,c, qlo, qhi, &
                shk, glo, ghi, &
                2, lo, (/hi(1),hi(2)+1,hi(3)/), domlo, domhi)

    do       k = lo(3), hi(3)
       do    j = lo(2), hi(2)
          do i = lo(1), hi(1)
             pdivu(i,j,k) = pdivu(i,j,k) + &
                            HALF*(pgdnvy(i,j,k)+pgdnvy(i,j+1,k)) * &
                                 (ugdnvy(i,j+1,k)-ugdnvy(i,j,k)) * dyinv
          end do
       end do
       do    j = lo(2), hi(2)+1
          do i = lo(1), hi(1)
             ugdy(i,j,k) = ugdy(i,j,k)
          end do
       end do
    end do

    nullify(ugdnvy,pgdnvy,gegdnvy)
    nullify(qyl,qyr)

    !
    ! Use qz?, q?xy, q?yx to compute final z-flux
    !

    fxy      =>      ftmp1
    ugdnvxy  =>  ugdnvtmp1
    pgdnvxy  =>  pgdnvtmp1
    gegdnvxy => gegdnvtmp1

    ! Inputs: qmxy, qpxy                       : xface, +-0 at y, +-1 at z
    !         gamc, csml, c                    : +-4
    !         shk                              : +-1
    ! Outputs: fxy, ugdnvxy, pgdnvxy, gegdnvxy : xface, +-0 at y, +-1 at z
    call cmpflx(qmxy,qpxy, fglo, fghi, &
                fxy, glo, ghi, &
                ugdnvxy,pgdnvxy,gegdnvxy, fglo, fghi, &
                gamc,csml,c, qlo, qhi, &
                shk, glo, ghi, &
                1, (/lo(1),lo(2),lo(3)-1/), (/hi(1)+1,hi(2),hi(3)+1/), domlo, domhi)

    fyx      =>      ftmp2
    ugdnvyx  =>  ugdnvtmp2
    pgdnvyx  =>  pgdnvtmp2
    gegdnvyx => gegdnvtmp2

    ! Inputs: qmyx, qpyx                       : yface, +-0 at x, +-1 at z
    !         gamc, csml, c                    : +-4
    !         shk                              : +-1
    ! Outputs: fyx, ugdnvyx, pgdnvyx, gegdnvyx : yface, +-0 at x, +-1 at z
    call cmpflx(qmyx,qpyx, fglo, fghi, &
                fyx, glo, ghi, &
                ugdnvyx,pgdnvyx,gegdnvyx, fglo, fghi, &
                gamc,csml,c, qlo, qhi, &
                shk, glo, ghi, &
                2, (/lo(1),lo(2),lo(3)-1/), (/hi(1),hi(2)+1,hi(3)+1/), domlo, domhi)
    
    qzl => ql
    qzr => qr
 
    ! Inputs: qzm, qzp                        : zface, +-1 at x & y
    !         fxy, ugdnvxy, pgdnvxy, gegdnvxy : xface, +-0 at y, +-1 at z
    !         fyx, ugdnvyx, pgdnvyx, gegdnvyx : yface, +-0 at x, +-1 at z
    !         gamc, grav, rot                 : +-4
    !         srcQ                            : +-1
    ! Outputs: qzl, qzr                       : zface, +-0 at x & y
    call transxy(qzm,qzl,qzp,qzr, fglo, fghi, &
                 fxy, fyx, glo, ghi, &
                 ugdnvxy,pgdnvxy,gegdnvxy,ugdnvyx,pgdnvyx,gegdnvyx, fglo, fghi, &
                 gamc, qlo, qhi, &
                 srcQ, slo, shi, &
                 grav, gvlo, gvhi, &
                 rot, rlo, rhi, &
                 hdt,hdtdx,hdtdy,lo,hi)

    nullify(fxy, ugdnvxy, pgdnvxy, gegdnvxy)
    nullify(fyx, ugdnvyx, pgdnvyx, gegdnvyx)

    ugdnvz  =>  ugdnvtmp1
    pgdnvz  =>  pgdnvtmp1
    gegdnvz => gegdnvtmp1

    ! Inputs: qzl, qzr                        : zface, +-0 at x & y
    !         gamc, csml, c                   : +-4
    !         shk                             : +-1
    ! Outputs: flux3, ugdnvz, pgdnvz, gegdnvz : zface, +-0 at x & y
    call cmpflx(qzl,qzr, fglo, fghi, &
                flux3, f3lo, f3hi, &
                ugdnvz,pgdnvz,gegdnvz, fglo, fghi, &
                gamc,csml,c, qlo, qhi, &
                shk, glo, ghi, &
                3, lo, (/hi(1),hi(2),hi(3)+1/), domlo, domhi)

    do       k = lo(3), hi(3)
       do    j = lo(2), hi(2)
          do i = lo(1), hi(1)
             pdivu(i,j,k) = pdivu(i,j,k) + &
                            HALF*(pgdnvz(i,j,k)+pgdnvz(i,j,k+1)) * &
                                 (ugdnvz(i,j,k+1)-ugdnvz(i,j,k)) * dzinv
          end do
       end do
    end do
    ugdz(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)+1) = ugdnvz(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)+1)
    
    nullify(ugdnvz,pgdnvz,gegdnvz)
    nullify(qzl,qzr)

    deallocate(qxm,qym,qzm,qxp,qyp,qzp,ql,qr)
    deallocate(qmxy,qpxy,qmxz,qpxz,qmyx,qpyx,qmyz,qpyz,qmzx,qpzx,qmzy,qpzy)
    deallocate(ftmp1,ftmp2)
    deallocate(pgdnvtmp1,pgdnvtmp2)
    deallocate(ugdnvtmp1,ugdnvtmp2)
    deallocate(gegdnvtmp1,gegdnvtmp2)

    deallocate(shk)
      
  end subroutine umeth3d

! ::: 
! ::: ------------------------------------------------------------------
! ::: 

  subroutine ctoprim(glo,ghi,lo,hi,flo,fhi, &
                     uin, ulo, uhi, &
                     q,c,gamc,csml,flatn, qlo, qhi, &
                     src, srclo, srchi, &
                     srcQ, slo, shi, &
                     courno,dx,dt) bind(c,name='ctoprim')
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

    integer, intent(in) :: glo(3), ghi(3), lo(3), hi(3), flo(3), fhi(3)
    integer, intent(in) :: ulo(3), uhi(3), qlo(3), qhi(3), srclo(3), srchi(3), slo(3), shi(3)
    double precision, intent(in ) :: uin  (ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3),NVAR)
    double precision, intent(out) :: q    (qlo(1):qhi(1),qlo(2):qhi(2),qlo(3):qhi(3),QVAR)
    double precision, intent(out) :: c    (qlo(1):qhi(1),qlo(2):qhi(2),qlo(3):qhi(3))
    double precision, intent(out) :: gamc (qlo(1):qhi(1),qlo(2):qhi(2),qlo(3):qhi(3))
    double precision, intent(out) :: csml (qlo(1):qhi(1),qlo(2):qhi(2),qlo(3):qhi(3))
    double precision, intent(out) :: flatn(qlo(1):qhi(1),qlo(2):qhi(2),qlo(3):qhi(3))
    double precision, intent(in ) :: src  (srclo(1):srchi(1),srclo(2):srchi(2),srclo(3):srchi(3),NVAR)
    double precision, intent(out) :: srcQ (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),QVAR)
    double precision, intent(in   ) :: dx(3), dt
    double precision, intent(inout) :: courno

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
    do       k = glo(3), ghi(3)
       do    j = glo(2), ghi(2)
          do i = glo(1), ghi(1)
             
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
       do       k = glo(3), ghi(3)
          do    j = glo(2), ghi(2)
             do i = glo(1), ghi(1)
                q(i,j,k,nq) = uin(i,j,k,n)/q(i,j,k,QRHO)
             enddo
          enddo
       enddo
    enddo

    ! Get gamc, p, T, c, csml using q state
    do       k = glo(3), ghi(3)
       do    j = glo(2), ghi(2)
          do i = glo(1), ghi(1)
             
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
             
             courx = ( c(i,j,k)+abs(q(i,j,k,QU)) ) * dt/dx(1)
             coury = ( c(i,j,k)+abs(q(i,j,k,QV)) ) * dt/dx(2)
             courz = ( c(i,j,k)+abs(q(i,j,k,QW)) ) * dt/dx(3)
             
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
       call uflaten(flo,fhi, &
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
