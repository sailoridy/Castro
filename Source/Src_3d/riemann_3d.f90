module riemann_module

  use bl_constants_module

  implicit none

  private

  public cmpflx, shock

contains

! :::
! ::: ------------------------------------------------------------------
! :::

  subroutine cmpflx(qm,qp,qlo,qhi, &
                    flx,flo,fhi, &
                    ugdnv,pgdnv,gegdnv,gdlo,gdhi, &
                    gamc,csml,c,gclo,gchi, &
                    shk,slo,shi,&
                    idir,lo,hi,domlo,domhi)
    ! note that lo(idir) and hi(idir) are face index, whereas other lo and hi are cell-centerd

    use meth_params_module, only : QVAR, NVAR, use_colglaz, hybrid_riemann
    implicit none

    integer, intent(in) :: qlo(3), qhi(3), flo(3), fhi(3), gdlo(3), gdhi(3), gclo(3), gchi(3), &
         slo(3), shi(3), idir, lo(3), hi(3), domlo(3), domhi(3)
    double precision,intent(inout)::    qm( qlo(1): qhi(1), qlo(2): qhi(2), qlo(3): qhi(3),QVAR)
    double precision,intent(inout)::    qp( qlo(1): qhi(1), qlo(2): qhi(2), qlo(3): qhi(3),QVAR)
    double precision,intent(inout):: ugdnv(gdlo(1):gdhi(1),gdlo(2):gdhi(2),gdlo(3):gdhi(3))
    double precision,intent(inout):: pgdnv(gdlo(1):gdhi(1),gdlo(2):gdhi(2),gdlo(3):gdhi(3))
    double precision,intent(inout)::gegdnv(gdlo(1):gdhi(1),gdlo(2):gdhi(2),gdlo(3):gdhi(3))
    double precision,intent(in   )::  gamc(gclo(1):gchi(1),gclo(2):gchi(2),gclo(3):gchi(3))
    double precision,intent(in   )::  csml(gclo(1):gchi(1),gclo(2):gchi(2),gclo(3):gchi(3))
    double precision,intent(in   )::     c(gclo(1):gchi(1),gclo(2):gchi(2),gclo(3):gchi(3))
    double precision,intent(in   )::   shk( slo(1): shi(1), slo(2): shi(2), slo(3): shi(3))
    double precision,intent(inout)::   flx( flo(1): fhi(1), flo(2): fhi(2), flo(3): fhi(3),NVAR)

    integer :: i,j,k, is_shock
    double precision :: cl, cr

    if (use_colglaz == 1) then
       call riemanncg(qm,qp,qlo,qhi, &
                      flx,flo,fhi, &
                      ugdnv,pgdnv,gegdnv,gdlo,gdhi, &
                      gamc,csml,c,gclo,gchi, &
                      idir,lo,hi,domlo,domhi)
    else
       call riemannus(qm,qp,qlo,qhi, &
                      flx,flo,fhi, &
                      ugdnv,pgdnv,gegdnv,gdlo,gdhi, &
                      gamc,csml,c,gclo,gchi, &
                      idir,lo,hi,domlo,domhi)
    end if

    if (hybrid_riemann == 1) then
       ! correct the fluxes using an HLL scheme if we are in a shock
       ! and doing the hybrid approach
       do    k = lo(3), hi(3)
          do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             select case (idir)
             case (1)
                is_shock = shk(i-1,j,k) + shk(i,j,k)
             case (2)
                is_shock = shk(i,j-1,k) + shk(i,j,k)
             case (3)
                is_shock = shk(i,j,k-1) + shk(i,j,k)
             end select

             if (is_shock >= 1) then

                select case (idir)
                case (1)
                   cl = c(i-1,j,k)
                   cr = c(i,j,k)
                case (2)
                   cl = c(i,j-1,k)
                   cr = c(i,j,k)
                case (3)
                   cl = c(i,j,k-1)
                   cr = c(i,j,k)
                end select

                call HLL(qm(i,j,k,:), qp(i,j,k,:), cl, cr, &
                         idir, flx(i,j,k,:))

             endif

          enddo
          enddo
       enddo

    endif

  end subroutine cmpflx


  subroutine shock(q,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                   shk,s_l1,s_l2,s_l3,s_h1,s_h2,s_h3, &
                   ilo1,ilo2,ilo3,ihi1,ihi2,ihi3,dx,dy,dz)

    use prob_params_module, only : coord_type
    use meth_params_module, only : QU, QV, QW, QPRES, QVAR

    integer, intent(in) :: qd_l1, qd_l2, qd_l3, qd_h1, qd_h2, qd_h3
    integer, intent(in) :: s_l1, s_l2, s_l3, s_h1, s_h2, s_h3
    integer, intent(in) :: ilo1, ilo2, ilo3, ihi1, ihi2, ihi3
    double precision, intent(in) :: dx, dy, dz
    double precision, intent(in) :: q(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
    double precision, intent(inout) :: shk(s_l1:s_h1,s_l2:s_h2,s_l3:s_h3)

    integer :: i, j, k

    double precision :: divU
    double precision :: px_pre, px_post, py_pre, py_post, pz_pre, pz_post
    double precision :: e_x, e_y, e_z, d
    double precision :: p_pre, p_post, pjump

    double precision, parameter :: small = 1.d-10
    double precision, parameter :: eps = 0.33d0

    ! This is a basic multi-dimensional shock detection algorithm.
    ! This implementation follows Flash, which in turn follows
    ! AMRA and a Woodward (1995) (supposedly -- couldn't locate that).
    !
    ! The spirit of this follows the shock detection in Colella &
    ! Woodward (1984)

    if (coord_type /= 0) then
       call bl_error("ERROR: invalid geometry in shock()")
    endif

    do k = ilo3-1, ihi3+1
       do j = ilo2-1, ihi2+1
          do i = ilo1-1, ihi1+1

             ! construct div{U}
             divU = HALF*(q(i+1,j,k,QU) - q(i-1,j,k,QU))/dx + &
                    HALF*(q(i,j+1,k,QV) - q(i,j-1,k,QV))/dy + &
                    HALF*(q(i,j,k+1,QW) - q(i,j,k-1,QW))/dz

             ! find the pre- and post-shock pressures in each direction
             if (q(i+1,j,k,QPRES) - q(i-1,j,k,QPRES) < ZERO) then
                px_pre  = q(i+1,j,k,QPRES)
                px_post = q(i-1,j,k,QPRES)
             else
                px_pre  = q(i-1,j,k,QPRES)
                px_post = q(i+1,j,k,QPRES)
             endif

             if (q(i,j+1,k,QPRES) - q(i,j-1,k,QPRES) < ZERO) then
                py_pre  = q(i,j+1,k,QPRES)
                py_post = q(i,j-1,k,QPRES)
             else
                py_pre  = q(i,j-1,k,QPRES)
                py_post = q(i,j+1,k,QPRES)
             endif

             if (q(i,j,k+1,QPRES) - q(i,j,k-1,QPRES) < ZERO) then
                pz_pre  = q(i,j,k+1,QPRES)
                pz_post = q(i,j,k-1,QPRES)
             else
                pz_pre  = q(i,j,k-1,QPRES)
                pz_post = q(i,j,k+1,QPRES)
             endif

             ! use compression to create unit vectors for the shock direction
             e_x = (q(i+1,j,k,QU) - q(i-1,j,k,QU))**2
             e_y = (q(i,j+1,k,QV) - q(i,j-1,k,QV))**2
             e_z = (q(i,j,k+1,QW) - q(i,j,k-1,QW))**2
             d = ONE/(e_x + e_y + e_z + small)

             e_x = e_x*d
             e_y = e_y*d
             e_z = e_z*d

             ! project the pressures onto the shock direction
             p_pre  = e_x*px_pre + e_y*py_pre + e_z*pz_pre
             p_post = e_x*px_post + e_y*py_post + e_z*pz_post

             ! test for compression + pressure jump to flag a shock
             pjump = eps - (p_post - p_pre)/p_pre

             if (pjump < ZERO .and. divU < ZERO) then
                shk(i,j,k) = ONE
             else
                shk(i,j,k) = ZERO
             endif

          enddo
       enddo
    enddo

  end subroutine shock



! :::
! ::: ------------------------------------------------------------------
! :::

  subroutine riemanncg(ql,qr,qlo,qhi, &
                       uflx,flo,fhi, &
                       ugdnv,pgdnv,gegdnv,gdlo,gdhi, &
                       gamc,csml,c,gclo,gchi, &
                       idir,lo,hi,domlo,domhi)
    ! note that lo(idir) and hi(idir) are face index, whereas other lo and hi are cell-centerd

    ! this implements the approximate Riemann solver of Colella & Glaz (1985)

    use bl_error_module
    use network, only : nspec, naux
    use eos_type_module
    use eos_module
    use prob_params_module, only : physbc_lo, physbc_hi, Symmetry, SlipWall, NoSlipWall
    use meth_params_module, only : QVAR, NVAR, QRHO, QU, QV, QW, &
                                   QPRES, QREINT, QESGS, QFA, QFS, &
                                   QFX, URHO, UMX, UMY, UMZ, UEDEN, UEINT, &
                                   UESGS, UFA, UFS, UFX, &
                                   nadv, small_dens, small_pres, small_temp, &
                                   cg_maxiter, cg_tol, &
                                   npassive, upass_map, qpass_map, ppm_temp_fix, allow_negative_energy
    use bl_constants_module

    implicit none

    double precision, parameter:: small = 1.d-8

    integer, intent(in) :: qlo(3), qhi(3), flo(3), fhi(3), gdlo(3), gdhi(3), gclo(3), gchi(3), &
         idir, lo(3), hi(3), domlo(3), domhi(3)
    double precision,intent(inout)::    ql( qlo(1): qhi(1), qlo(2): qhi(2), qlo(3): qhi(3),QVAR)
    double precision,intent(inout)::    qr( qlo(1): qhi(1), qlo(2): qhi(2), qlo(3): qhi(3),QVAR)
    double precision,intent(inout):: ugdnv(gdlo(1):gdhi(1),gdlo(2):gdhi(2),gdlo(3):gdhi(3))
    double precision,intent(inout):: pgdnv(gdlo(1):gdhi(1),gdlo(2):gdhi(2),gdlo(3):gdhi(3))
    double precision,intent(inout)::gegdnv(gdlo(1):gdhi(1),gdlo(2):gdhi(2),gdlo(3):gdhi(3))
    double precision,intent(in   )::  gamc(gclo(1):gchi(1),gclo(2):gchi(2),gclo(3):gchi(3))
    double precision,intent(in   )::  csml(gclo(1):gchi(1),gclo(2):gchi(2),gclo(3):gchi(3))
    double precision,intent(in   )::     c(gclo(1):gchi(1),gclo(2):gchi(2),gclo(3):gchi(3))
    double precision,intent(inout)::  uflx( flo(1): fhi(1), flo(2): fhi(2), flo(3): fhi(3),NVAR)

    integer :: i,j,k,n,nq,ipassive,iadv, ispec, iaux

    double precision :: gamcl, gamcr, cav, smallc
    double precision :: rgdnv,v1gdnv,v2gdnv,ustar,gamgdnv
    double precision :: rl, ul, v1l, v2l, pl, rel
    double precision :: rr, ur, v1r, v2r, pr, rer
    double precision :: wl, wr, rhoetot, scr
    double precision :: rstar, cstar, pstar
    double precision :: ro, uo, po, co, gamco
    double precision :: sgnm, spin, spout, ushock, frac
    double precision :: wsmall, csmall,qavg
    double precision :: rho_K_contrib

    double precision :: gcl, gcr
    double precision :: clsq, clsql, clsqr, wlsq, wosq, wrsq, wo
    double precision :: zm, zp
    double precision :: denom, dpditer, dpjmp
    double precision :: gamc_bar, game_bar
    double precision :: gamel, gamer, gameo, gamstar, gmin, gmax, gdot

    integer :: iter, iter_max
    double precision :: tol
    double precision :: err

    logical :: converged

    double precision :: pstnm1
    double precision :: taul, taur, tauo
    double precision :: ustarm, ustarp, ustnm1, ustnp1

    double precision, parameter :: weakwv = 1.d-3

    double precision, allocatable :: pstar_hist(:)

    type (eos_t) :: eos_state

    tol = cg_tol
    iter_max = cg_maxiter

    allocate (pstar_hist(iter_max))

    do    k = lo(3), hi(3)
       do j = lo(2), hi(2)
       do i = lo(1), hi(1)

         ! pick left and right states based on direction
          if(idir.eq.1) then
             ul  = ql(i,j,k,QU)
             v1l = ql(i,j,k,QV)
             v2l = ql(i,j,k,QW)
             !
             ur  = qr(i,j,k,QU)
             v1r = qr(i,j,k,QV)
             v2r = qr(i,j,k,QW)
             !
             gamcl  = gamc(i-1,j,k)
             gamcr  = gamc(i  ,j,k)
             cav    = HALF*(c(i-1,j,k)+c(i,j,k))
             smallc = max(csml(i-1,j,k),csml(i,j,k)) 
          elseif(idir.eq.2) then
             ul  = ql(i,j,k,QV)
             v1l = ql(i,j,k,QU)
             v2l = ql(i,j,k,QW)
             !
             ur  = qr(i,j,k,QV)
             v1r = qr(i,j,k,QU)
             v2r = qr(i,j,k,QW)
             !
             gamcl  = gamc(i,j-1,k)
             gamcr  = gamc(i,j  ,k)
             cav    = HALF*(c(i,j-1,k)+c(i,j,k))
             smallc = max(csml(i,j-1,k),csml(i,j,k)) 
          else
             ul  = ql(i,j,k,QW)
             v1l = ql(i,j,k,QU)
             v2l = ql(i,j,k,QV)
             !
             ur  = qr(i,j,k,QW)
             v1r = qr(i,j,k,QU)
             v2r = qr(i,j,k,QV)
             !
             gamcl  = gamc(i,j,k-1)
             gamcr  = gamc(i,j,k  )
             cav    = HALF*(c(i,j,k-1)+c(i,j,k))
             smallc = max(csml(i,j,k-1),csml(i,j,k)) 
          endif

          if (ppm_temp_fix == 2) then
             ! recompute the thermodynamics on the interface to make it
             ! all consistent

             ! we want to take the edge states of rho, p, and X, and get
             ! new values for gamc and (rho e) on the edges that are
             ! thermodynamically consistent.

             ! this is an initial guess for iterations, since we
             ! can't be certain that temp is on interfaces
             eos_state%T = 10000.0d0

             ! minus state
             eos_state%rho = ql(i,j,k,QRHO)
             eos_state%p   = ql(i,j,k,QPRES)
             eos_state%e   = ql(i,j,k,QREINT)/ql(i,j,k,QRHO)
             eos_state%xn  = ql(i,j,k,QFS:QFS-1+nspec)
             eos_state%aux = ql(i,j,k,QFX:QFX-1+naux)

             ! Protect against negative energies

             if (allow_negative_energy .eq. 0 .and. eos_state % e < ZERO) then
                eos_state % T = small_temp
                call eos(eos_input_rt, eos_state, .false.)                
             else
                call eos(eos_input_re, eos_state, .false.)
             endif
             
             ql(i,j,k,QREINT) = ql(i,j,k,QRHO)*eos_state%e
             ql(i,j,k,QPRES) = eos_state%p
             gamcl = eos_state%gam1

             ! plus state
             eos_state%rho = qr(i,j,k,QRHO)
             eos_state%p   = qr(i,j,k,QPRES)
             eos_state%e   = qr(i,j,k,QREINT)/qr(i,j,k,QRHO)
             eos_state%xn  = qr(i,j,k,QFS:QFS-1+nspec)
             eos_state%aux = qr(i,j,k,QFX:QFX-1+naux)
             
             ! Protect against negative energies
             
             if (allow_negative_energy .eq. 0 .and. eos_state % e < ZERO) then
                eos_state % T = small_temp
                call eos(eos_input_rt, eos_state, .false.)                
             else              
                call eos(eos_input_re, eos_state, .false.)
             endif
             
             qr(i,j,k,QREINT) = qr(i,j,k,QRHO)*eos_state%e
             qr(i,j,k,QPRES) = eos_state%p
             gamcr = eos_state%gam1
          end if

          ! left state
          rl = max(ql(i,j,k,QRHO),small_dens)
          pl  = ql(i,j,k,QPRES)
          rel = ql(i,j,k,QREINT)
          gcl = gamcl

          ! sometime we come in here with negative energy or pressure
          ! note: reset both in either case, to remain thermo
          ! consistent
          if (rel <= ZERO .or. pl < small_pres) then
             print *, "WARNING: (rho e)_l < 0 or pl < small_pres in Riemann: ", rel, pl, small_pres

             eos_state%T   = small_temp
             eos_state%rho = rl
             eos_state%xn  = ql(i,j,k,QFS:QFS-1+nspec)
             eos_state%aux = ql(i,j,k,QFX:QFX-1+naux)

             call eos(eos_input_rt, eos_state, .false.)

             rel = rl*eos_state%e
             pl  = eos_state%p
             gcl = eos_state%gam1
          endif

          ! right state
          rr = max(qr(i,j,k,QRHO),small_dens)
          pr  = qr(i,j,k,QPRES)
          rer = qr(i,j,k,QREINT)
          gcr = gamcr

          if (rer <= ZERO .or. pr < small_pres) then
             print *, "WARNING: (rho e)_r < 0 or pr < small_pres in Riemann: ", rer, pr, small_pres

             eos_state % T   = small_temp
             eos_state % rho = rr
             eos_state % xn  = qr(i,j,k,QFS:QFS-1+nspec)
             eos_state % aux = qr(i,j,k,QFX:QFX-1+naux)

             call eos(eos_input_rt, eos_state, .false.)

             rer = rr*eos_state%e
             pr  = eos_state%p
             gcr = eos_state%gam1
          endif


          ! common quantities
          taul = ONE/rl
          taur = ONE/rr

          ! lagrangian sound speeds
          clsql = gcl*pl*rl
          clsqr = gcr*pr*rr


          ! Note: in the original Colella & Glaz paper, they predicted
          ! gamma_e to the interfaces using a special (non-hyperbolic)
          ! evolution equation.  In Castro, we instead bring (rho e)
          ! to the edges, so we construct the necessary gamma_e here from
          ! what we have on the interfaces.
          gamel = pl/rel + ONE
          gamer = pr/rer + ONE

          ! these should consider a wider average of the cell-centered
          ! gammas
          gmin = min(gamel, gamer, ONE, FOUR3RD)
          gmax = max(gamel, gamer, TWO, FIVE3RD)

          game_bar = HALF*(gamel + gamer)
          gamc_bar = HALF*(gcl + gcr)

          gdot = TWO*(ONE - game_bar/gamc_bar)*(game_bar - ONE)

          csmall = smallc
          wsmall = small_dens*csmall
          wl = max(wsmall,sqrt(abs(clsql)))
          wr = max(wsmall,sqrt(abs(clsqr)))

          ! make an initial guess for pstar -- this is a two-shock
          ! approximation
          !pstar = ((wr*pl + wl*pr) + wl*wr*(ul - ur))/(wl + wr)
          pstar = pl + ( (pr - pl) - wr*(ur - ul) )*wl/(wl+wr)
          pstar = max(pstar,small_pres)

          ! get the shock speeds -- this computes W_s from CG Eq. 34
          call wsqge(pl,taul,gamel,gdot,  &
                     gamstar,pstar,wlsq,clsql,gmin,gmax)

          call wsqge(pr,taur,gamer,gdot,  &
                     gamstar,pstar,wrsq,clsqr,gmin,gmax)

          pstnm1 = pstar

          wl = sqrt(wlsq)
          wr = sqrt(wrsq)

          ! R-H jump conditions give ustar across each wave -- these
          ! should be equal when we are done iterating.  Our notation
          ! here is a little funny, comparing to CG, ustarp = u*_L and
          ! ustarm = u*_R.
          ustarp = ul - (pstar-pl)/wl
          ustarm = ur + (pstar-pr)/wr

          ! revise our pstar guess
          !pstar = ((wr*pl + wl*pr) + wl*wr*(ul - ur))/(wl + wr)
          pstar = pl + ( (pr - pl) - wr*(ur - ul) )*wl/(wl+wr)
          pstar = max(pstar,small_pres)

          ! sectant iteration
          converged = .false.
          iter = 1
          do while ((iter <= iter_max .and. .not. converged) .or. iter <= 2)

             call wsqge(pl,taul,gamel,gdot,  &
                        gamstar,pstar,wlsq,clsql,gmin,gmax)

             call wsqge(pr,taur,gamer,gdot,  &
                        gamstar,pstar,wrsq,clsqr,gmin,gmax)

             ! NOTE: these are really the inverses of the wave speeds!
             wl = ONE / sqrt(wlsq)
             wr = ONE / sqrt(wrsq)

             ustnm1 = ustarm
             ustnp1 = ustarp

             ustarm = ur-(pr-pstar)*wr
             ustarp = ul+(pl-pstar)*wl

             dpditer=abs(pstnm1-pstar)

             ! Here we are going to do the Secant iteration version in
             ! CG.  Note that what we call zp and zm here are not
             ! actually the Z_p = |dp*/du*_p| defined in CG, by rather
             ! simply |du*_p| (or something that looks like dp/Z!).
             zp=abs(ustarp-ustnp1)
             if(zp-weakwv*cav <= ZERO)then
                zp = dpditer*wl
             endif

             zm=abs(ustarm-ustnm1)
             if(zm-weakwv*cav <= ZERO)then
                zm = dpditer*wr
             endif

             ! the new pstar is found via CG Eq. 18
             denom=dpditer/max(zp+zm,small*cav)
             pstnm1 = pstar
             pstar=pstar-denom*(ustarm-ustarp)
             pstar=max(pstar,small_pres)

             err = abs(pstar - pstnm1)
             if (err < tol*pstar) converged = .true.

             pstar_hist(iter) = pstar

             iter = iter + 1

          enddo

          if (.not. converged) then
             print *, 'pstar history: '
             do iter = 1, iter_max
                print *, iter, pstar_hist(iter)
             enddo

             print *, ' '
             print *, 'left state  (r,u,p,re,gc): ', rl, ul, pl, rel, gcl
             print *, 'right state (r,u,p,re,gc): ', rr, ur, pr, rer, gcr
             call bl_error("ERROR: non-convergence in the Riemann solver")
          endif


          ! we converged!  construct the single ustar for the region
          ! between the left and right waves, using the updated wave speeds
          ustarm = ur-(pr-pstar)*wr  ! careful -- here wl, wr are 1/W
          ustarp = ul+(pl-pstar)*wl

          ustar = HALF* ( ustarp + ustarm )


          ! sample the solution -- here we look first at the direction
          ! that the contact is moving.  This tells us if we need to
          ! worry about the L/L* states or the R*/R states.
          if (ustar .gt. ZERO) then
             ro = rl
             uo = ul
             po = pl
             tauo = taul
             gamco = gcl
             gameo = gamel

          else if (ustar .lt. ZERO) then
             ro = rr
             uo = ur
             po = pr
             tauo = taur
             gamco = gcr
             gameo = gamer
          else
             ro = HALF*(rl+rr)
             uo = HALF*(ul+ur)
             po = HALF*(pl+pr)
             tauo = HALF*(taul+taur)
             gamco = HALF*(gcl+gcr)
             gameo = HALF*(gamel + gamer)
          endif

          ! use tau = 1/rho as the independent variable here
          ro = max(small_dens,ONE/tauo)
          tauo = ONE/ro

          co = sqrt(abs(gamco*po/ro))
          co = max(csmall,co)
          clsq = (co*ro)**2

          ! now that we know which state (left or right) we need to worry
          ! about, get the value of gamstar and wosq across the wave we
          ! are dealing with.
          call wsqge(po,tauo,gameo,gdot,   &
                     gamstar,pstar,wosq,clsq,gmin,gmax)

          sgnm = sign(ONE,ustar)

          wo = sqrt(wosq)
          dpjmp = pstar - po

          ! is this max really necessary?
          !rstar=max(ONE-ro*dpjmp/wosq, (gameo-ONE)/(gameo+ONE))
          rstar=ONE-ro*dpjmp/wosq
          rstar=ro/rstar
          rstar = max(small_dens,rstar)

          cstar = sqrt(abs(gamco*pstar/rstar))
          cstar = max(cstar,csmall)


          spout = co - sgnm*uo
          spin = cstar - sgnm*ustar

          !ushock = HALF*(spin + spout)
          ushock = wo/ro - sgnm*uo

          if (pstar-po .ge. ZERO) then
             spin = ushock
             spout = ushock
          endif
          !if (spout-spin .eq. ZERO) then
          !   scr = small*cav
          !else
          !   scr = spout-spin
          !endif
          !frac = (ONE + (spout + spin)/scr)*HALF
          !frac = max(ZERO,min(ONE,frac))

          frac = HALF*(ONE + (spin + spout)/max(spout-spin,spin+spout, small*cav))

          ! the transverse velocity states only depend on the
          ! direction that the contact moves
          if (ustar .gt. ZERO) then
             v1gdnv = v1l
             v2gdnv = v2l
          else if (ustar .lt. ZERO) then
             v1gdnv = v1r
             v2gdnv = v2r
          else
             v1gdnv = HALF*(v1l+v1r)
             v2gdnv = HALF*(v2l+v2r)
          endif

          ! linearly interpolate between the star and normal state -- this covers the
          ! case where we are inside the rarefaction fan.
          rgdnv = frac*rstar + (ONE - frac)*ro
          ugdnv(i,j,k) = frac*ustar + (ONE - frac)*uo
          pgdnv(i,j,k) = frac*pstar + (ONE - frac)*po
          gamgdnv =  frac*gamstar + (ONE-frac)*gameo

          ! now handle the cases where instead we are fully in the
          ! star or fully in the original (l/r) state
          if (spout .lt. ZERO) then
             rgdnv = ro
             ugdnv(i,j,k) = uo
             pgdnv(i,j,k) = po
             gamgdnv = gameo
          endif
          if (spin .ge. ZERO) then
             rgdnv = rstar
             ugdnv(i,j,k) = ustar
             pgdnv(i,j,k) = pstar
             gamgdnv = gamstar
          endif

          gegdnv(i,j,k) = gamgdnv

          pgdnv(i,j,k) = max(pgdnv(i,j,k),small_pres)

          ! Enforce that fluxes through a symmetry plane or wall are hard zero.
          if (idir .eq. 1) then
             if (i.eq.domlo(1) .and. &
                 (physbc_lo(1) .eq. Symmetry .or.  physbc_lo(1) .eq. SlipWall .or. &
                  physbc_lo(1) .eq. NoSlipWall) ) &
                  ugdnv(i,j,k) = ZERO
             if (i.eq.domhi(1)+1 .and. &
                 (physbc_hi(1) .eq. Symmetry .or.  physbc_hi(1) .eq. SlipWall .or. &
                  physbc_hi(1) .eq. NoSlipWall) ) &
                  ugdnv(i,j,k) = ZERO
          end if
          if (idir .eq. 2) then
             if (j.eq.domlo(2) .and. &
                 (physbc_lo(2) .eq. Symmetry .or.  physbc_lo(2) .eq. SlipWall .or. &
                  physbc_lo(2) .eq. NoSlipWall) ) &
                  ugdnv(i,j,k) = ZERO
             if (j.eq.domhi(2)+1 .and. &
                 (physbc_hi(2) .eq. Symmetry .or.  physbc_hi(2) .eq. SlipWall .or. &
                  physbc_hi(2) .eq. NoSlipWall) ) &
                  ugdnv(i,j,k) = ZERO
          end if
          if (idir .eq. 3) then
             if (k.eq.domlo(3) .and. &
                 (physbc_lo(3) .eq. Symmetry .or.  physbc_lo(3) .eq. SlipWall .or. &
                  physbc_lo(3) .eq. NoSlipWall) ) &
                  ugdnv(i,j,k) = ZERO
             if (k.eq.domhi(3)+1 .and. &
                 (physbc_hi(3) .eq. Symmetry .or.  physbc_hi(3) .eq. SlipWall .or. &
                  physbc_hi(3) .eq. NoSlipWall) ) &
                  ugdnv(i,j,k) = ZERO
          end if

          ! Compute fluxes, order as conserved state (not q)
          uflx(i,j,k,URHO) = rgdnv*ugdnv(i,j,k)

          if(idir.eq.1) then
             uflx(i,j,k,UMX) = uflx(i,j,k,URHO)*ugdnv(i,j,k) + pgdnv(i,j,k)
             uflx(i,j,k,UMY) = uflx(i,j,k,URHO)*v1gdnv
             uflx(i,j,k,UMZ) = uflx(i,j,k,URHO)*v2gdnv
          elseif(idir.eq.2) then
             uflx(i,j,k,UMX) = uflx(i,j,k,URHO)*v1gdnv
             uflx(i,j,k,UMY) = uflx(i,j,k,URHO)*ugdnv(i,j,k) + pgdnv(i,j,k)
             uflx(i,j,k,UMZ) = uflx(i,j,k,URHO)*v2gdnv
          else
             uflx(i,j,k,UMX) = uflx(i,j,k,URHO)*v1gdnv
             uflx(i,j,k,UMY) = uflx(i,j,k,URHO)*v2gdnv
             uflx(i,j,k,UMZ) = uflx(i,j,k,URHO)*ugdnv(i,j,k) + pgdnv(i,j,k)
          endif

          ! compute the total energy from the internal, p/(gamma - 1), and the kinetic
          rhoetot = pgdnv(i,j,k)/(gamgdnv - ONE) + &
               HALF*rgdnv*(ugdnv(i,j,k)**2 + v1gdnv**2 + v2gdnv**2)

          uflx(i,j,k,UEDEN) = ugdnv(i,j,k)*(rhoetot + pgdnv(i,j,k))
          uflx(i,j,k,UEINT) = ugdnv(i,j,k)*pgdnv(i,j,k)/(gamgdnv - ONE)


          ! Treat K as a passively advected quantity but allow it to
          ! affect fluxes of (rho E) and momenta.
          if (UESGS .gt. -1) then
             n  = UESGS
             nq = QESGS
             if (ustar .gt. ZERO) then
                qavg = ql(i,j,k,nq)
             else if (ustar .lt. ZERO) then
                qavg = qr(i,j,k,nq)
             else
                qavg = HALF * (ql(i,j,k,nq) + qr(i,j,k,nq))
             endif

             uflx(i,j,k,n) = uflx(i,j,k,URHO)*qavg

             rho_K_contrib =  TWO3RD * rgdnv * qavg

             if(idir.eq.1) then
                uflx(i,j,k,UMX) = uflx(i,j,k,UMX) + rho_K_contrib
             elseif(idir.eq.2) then
                uflx(i,j,k,UMY) = uflx(i,j,k,UMY) + rho_K_contrib
             elseif(idir.eq.3) then
                uflx(i,j,k,UMZ) = uflx(i,j,k,UMZ) + rho_K_contrib
             endif

             uflx(i,j,k,UEDEN) = uflx(i,j,k,UEDEN) + ugdnv(i,j,k) * rho_K_contrib
          end if

          ! advected quantities -- only the contact matters
          do ipassive = 1, npassive
             n  = upass_map(ipassive)
             nq = qpass_map(ipassive)

             if (ustar .gt. ZERO) then
                uflx(i,j,k,n) = uflx(i,j,k,URHO)*ql(i,j,k,nq)
             else if (ustar .lt. ZERO) then
                uflx(i,j,k,n) = uflx(i,j,k,URHO)*qr(i,j,k,nq)
             else
                qavg = HALF * (ql(i,j,k,nq) + qr(i,j,k,nq))
                uflx(i,j,k,n) = uflx(i,j,k,URHO)*qavg
             endif
          enddo

       enddo
       enddo
    enddo

  end subroutine riemanncg

  subroutine wsqge(p,v,gam,gdot,gstar,pstar,wsq,csq,gmin,gmax)

    use bl_constants_module

    implicit none

    double precision p,v,gam,gdot,gstar,pstar,wsq,csq,gmin,gmax

    double precision, parameter :: smlp1 = 1.d-10
    double precision, parameter :: small = 1.d-7

    double precision :: alpha, beta

    ! First predict a value of game across the shock

    ! CG Eq. 31
    gstar=(pstar-p)*gdot/(pstar+p) + gam
    gstar=max(gmin,min(gmax,gstar))

    ! Now use that predicted value of game with the R-H jump conditions
    ! to compute the wave speed.

    ! CG Eq. 34
    ! wsq = (HALF*(gstar-ONE)*(pstar+p)+pstar)
    ! temp = ((gstar-gam)/(gam-ONE))

    ! if (pstar-p == ZERO) then
    !    divide=small
    ! else
    !    divide=pstar-p
    ! endif

    ! temp=temp/divide
    ! wsq = wsq/(v - temp*p*v)

    alpha = pstar-(gstar-ONE)*p/(gam-ONE)
    if (alpha == ZERO) alpha = smlp1*(pstar + p)

    beta = pstar + HALF*(gstar-ONE)*(pstar+p)

    wsq = (pstar-p)*beta/(v*alpha)

    if (abs(pstar - p) < smlp1*(pstar + p)) then
       wsq = csq
    endif
    wsq=max(wsq,(HALF*(gam-ONE)/gam)*csq)

    return
  end subroutine wsqge

! :::
! ::: ------------------------------------------------------------------
! :::

  subroutine riemannus(ql,qr,qlo,qhi, &
                       uflx,flo,fhi, &
                       ugdnv,pgdnv,gegdnv,gdlo,gdhi, &
                       gamc,csml,c,gclo,gchi, &
                       idir,lo,hi,domlo,domhi)
    ! note that lo(idir) and hi(idir) are face index, whereas other lo and hi are cell-centerd

    use network, only : nspec, naux
    use eos_type_module
    use eos_module
    use prob_params_module, only : physbc_lo, physbc_hi, Symmetry, SlipWall, NoSlipWall
    use meth_params_module, only : QVAR, NVAR, QRHO, QU, QV, QW, QPRES, QREINT, QESGS, QFA, QFS, &
         QFX, URHO, UMX, UMY, UMZ, UEDEN, UEINT, UESGS, UFA, UFS, UFX, &
         nadv, small_dens, small_pres, small_temp, npassive, upass_map, qpass_map, ppm_temp_fix, &
         allow_negative_energy
    use bl_constants_module
    implicit none

    double precision, parameter:: small = 1.d-8

    integer, intent(in) :: qlo(3), qhi(3), flo(3), fhi(3), gdlo(3), gdhi(3), gclo(3), gchi(3), &
         idir, lo(3), hi(3), domlo(3), domhi(3)
    double precision,intent(inout)::    ql( qlo(1): qhi(1), qlo(2): qhi(2), qlo(3): qhi(3),QVAR)
    double precision,intent(inout)::    qr( qlo(1): qhi(1), qlo(2): qhi(2), qlo(3): qhi(3),QVAR)
    double precision,intent(inout):: ugdnv(gdlo(1):gdhi(1),gdlo(2):gdhi(2),gdlo(3):gdhi(3))
    double precision,intent(inout):: pgdnv(gdlo(1):gdhi(1),gdlo(2):gdhi(2),gdlo(3):gdhi(3))
    double precision,intent(inout)::gegdnv(gdlo(1):gdhi(1),gdlo(2):gdhi(2),gdlo(3):gdhi(3))
    double precision,intent(in   )::  gamc(gclo(1):gchi(1),gclo(2):gchi(2),gclo(3):gchi(3))
    double precision,intent(in   )::  csml(gclo(1):gchi(1),gclo(2):gchi(2),gclo(3):gchi(3))
    double precision,intent(in   )::     c(gclo(1):gchi(1),gclo(2):gchi(2),gclo(3):gchi(3))
    double precision,intent(inout)::  uflx( flo(1): fhi(1), flo(2): fhi(2), flo(3): fhi(3),NVAR)

    integer :: i,j,k,n,nq,ipassive,iadv, ispec, iaux

    double precision :: gamcl, gamcr, cav, smallc
    double precision :: rgdnv,v1gdnv,v2gdnv,regdnv,ustar
    double precision :: rl, ul, v1l, v2l, pl, rel
    double precision :: rr, ur, v1r, v2r, pr, rer
    double precision :: wl, wr, rhoetot, scr
    double precision :: rstar, cstar, estar, pstar
    double precision :: ro, uo, po, reo, co, gamco, entho
    double precision :: sgnm, spin, spout, ushock, frac
    double precision :: wsmall, csmall,qavg
    double precision :: rho_K_contrib
    type(eos_t) :: eos_state

    do    k = lo(3), hi(3)
       do j = lo(2), hi(2)
       do i = lo(1), hi(1)

         ! pick left and right states based on direction
          if(idir.eq.1) then
             ul  = ql(i,j,k,QU)
             v1l = ql(i,j,k,QV)
             v2l = ql(i,j,k,QW)
             !
             ur  = qr(i,j,k,QU)
             v1r = qr(i,j,k,QV)
             v2r = qr(i,j,k,QW)
             !
             gamcl  = gamc(i-1,j,k)
             gamcr  = gamc(i  ,j,k)
             cav    = HALF*(c(i-1,j,k)+c(i,j,k))
             smallc = max(csml(i-1,j,k),csml(i,j,k)) 
          elseif(idir.eq.2) then
             ul  = ql(i,j,k,QV)
             v1l = ql(i,j,k,QU)
             v2l = ql(i,j,k,QW)
             !
             ur  = qr(i,j,k,QV)
             v1r = qr(i,j,k,QU)
             v2r = qr(i,j,k,QW)
             !
             gamcl  = gamc(i,j-1,k)
             gamcr  = gamc(i,j  ,k)
             cav    = HALF*(c(i,j-1,k)+c(i,j,k))
             smallc = max(csml(i,j-1,k),csml(i,j,k)) 
          else
             ul  = ql(i,j,k,QW)
             v1l = ql(i,j,k,QU)
             v2l = ql(i,j,k,QV)
             !
             ur  = qr(i,j,k,QW)
             v1r = qr(i,j,k,QU)
             v2r = qr(i,j,k,QV)
             !
             gamcl  = gamc(i,j,k-1)
             gamcr  = gamc(i,j,k  )
             cav    = HALF*(c(i,j,k-1)+c(i,j,k))
             smallc = max(csml(i,j,k-1),csml(i,j,k)) 
          endif

          if (ppm_temp_fix == 2) then
             ! recompute the thermodynamics on the interface to make it
             ! all consistent

             ! we want to take the edge states of rho, p, and X, and get
             ! new values for gamc and (rho e) on the edges that are
             ! thermodynamically consistent.

             ! this is an initial guess for iterations, since we
             ! can't be certain that temp is on interfaces
             eos_state%T = 10000.0d0

             ! minus state
             eos_state%rho = ql(i,j,k,QRHO)
             eos_state%p   = ql(i,j,k,QPRES)
             eos_state%e   = ql(i,j,k,QREINT)/ql(i,j,k,QRHO)
             eos_state%xn  = ql(i,j,k,QFS:QFS-1+nspec)
             eos_state%aux = ql(i,j,k,QFX:QFX-1+naux)

             ! Protect against negative energies

             if (allow_negative_energy .eq. 0 .and. eos_state % e < ZERO) then
                eos_state % T = small_temp
                call eos(eos_input_rt, eos_state, .false.)                
             else
                call eos(eos_input_re, eos_state, .false.)
             endif
             
             ql(i,j,k,QREINT) = ql(i,j,k,QRHO)*eos_state%e
             ql(i,j,k,QPRES) = eos_state%p
             gamcl = eos_state%gam1

             ! plus state
             eos_state%rho = qr(i,j,k,QRHO)
             eos_state%p   = qr(i,j,k,QPRES)
             eos_state%e   = qr(i,j,k,QREINT)/qr(i,j,k,QRHO)
             eos_state%xn  = qr(i,j,k,QFS:QFS-1+nspec)
             eos_state%aux = qr(i,j,k,QFX:QFX-1+naux)
             
             ! Protect against negative energies
             
             if (allow_negative_energy .eq. 0 .and. eos_state % e < ZERO) then
                eos_state % T = small_temp
                call eos(eos_input_rt, eos_state, .false.)                
             else              
                call eos(eos_input_re, eos_state, .false.)
             endif
             
             qr(i,j,k,QREINT) = qr(i,j,k,QRHO)*eos_state%e
             qr(i,j,k,QPRES) = eos_state%p
             gamcr = eos_state%gam1
          end if

          rl = max(ql(i,j,k,QRHO),small_dens)
          pl  = max(ql(i,j,k,QPRES ),small_pres)
          rel =     ql(i,j,k,QREINT)

          rr = max(qr(i,j,k,QRHO),small_dens)
          pr  = max(qr(i,j,k,QPRES),small_pres)
          rer =     qr(i,j,k,QREINT)

          csmall = smallc
          wsmall = small_dens*csmall
          wl = max(wsmall,sqrt(abs(gamcl*pl*rl)))
          wr = max(wsmall,sqrt(abs(gamcr*pr*rr)))

          pstar = ((wr*pl + wl*pr) + wl*wr*(ul - ur))/(wl + wr)
          ustar = ((wl*ul + wr*ur) + (pl - pr))/(wl + wr)
          pstar = max(pstar,small_pres)

          if (ustar .gt. ZERO) then
             ro = rl
             uo = ul
             po = pl
             reo = rel
             gamco = gamcl
          else if (ustar .lt. ZERO) then
             ro = rr
             uo = ur
             po = pr
             reo = rer
             gamco = gamcr
          else
             ro = HALF*(rl+rr)
             uo = HALF*(ul+ur)
             po = HALF*(pl+pr)
             reo = HALF*(rel+rer)
             gamco = HALF*(gamcl+gamcr)
          endif
          ro = max(small_dens,ro)

          co = sqrt(abs(gamco*po/ro))
          co = max(csmall,co)
          entho = (reo/ro + po/ro)/co**2
          rstar = ro + (pstar - po)/co**2
          rstar = max(small_dens,rstar)
          estar = reo + (pstar - po)*entho
          cstar = sqrt(abs(gamco*pstar/rstar))
          cstar = max(cstar,csmall)

          sgnm = sign(ONE,ustar)
          spout = co - sgnm*uo
          spin = cstar - sgnm*ustar
          ushock = HALF*(spin + spout)
          if (pstar-po .ge. ZERO) then
             spin = ushock
             spout = ushock
          endif
          if (spout-spin .eq. ZERO) then
             scr = small*cav
          else
             scr = spout-spin
          endif
          frac = (ONE + (spout + spin)/scr)*HALF
          frac = max(ZERO,min(ONE,frac))

          if (ustar .gt. ZERO) then
             v1gdnv = v1l
             v2gdnv = v2l
          else if (ustar .lt. ZERO) then
             v1gdnv = v1r
             v2gdnv = v2r
          else
             v1gdnv = HALF*(v1l+v1r)
             v2gdnv = HALF*(v2l+v2r)
          endif
          rgdnv = frac*rstar + (ONE - frac)*ro

          ugdnv(i,j,k) = frac*ustar + (ONE - frac)*uo
          pgdnv(i,j,k) = frac*pstar + (ONE - frac)*po

          regdnv = frac*estar + (ONE - frac)*reo
          if (spout .lt. ZERO) then
             rgdnv = ro
             ugdnv(i,j,k) = uo
             pgdnv(i,j,k) = po
             regdnv = reo
          endif
          if (spin .ge. ZERO) then
             rgdnv = rstar
             ugdnv(i,j,k) = ustar
             pgdnv(i,j,k) = pstar
             regdnv = estar
          endif

          gegdnv(i,j,k) = pgdnv(i,j,k)/regdnv + 1.0d0

          pgdnv(i,j,k) = max(pgdnv(i,j,k),small_pres)

          ! Enforce that fluxes through a symmetry plane or wall are hard zero.
          if (idir .eq. 1) then
             if (i.eq.domlo(1) .and. &
                  (physbc_lo(1) .eq. Symmetry .or.  physbc_lo(1) .eq. SlipWall .or. &
                  physbc_lo(1) .eq. NoSlipWall) ) &
                  ugdnv(i,j,k) = ZERO
             if (i.eq.domhi(1)+1 .and. &
                  (physbc_hi(1) .eq. Symmetry .or.  physbc_hi(1) .eq. SlipWall .or. &
                  physbc_hi(1) .eq. NoSlipWall) ) &
                  ugdnv(i,j,k) = ZERO
          end if
          if (idir .eq. 2) then
             if (j.eq.domlo(2) .and. &
                  (physbc_lo(2) .eq. Symmetry .or.  physbc_lo(2) .eq. SlipWall .or. &
                  physbc_lo(2) .eq. NoSlipWall) ) &
                  ugdnv(i,j,k) = ZERO
             if (j.eq.domhi(2)+1 .and. &
                  (physbc_hi(2) .eq. Symmetry .or.  physbc_hi(2) .eq. SlipWall .or. &
                  physbc_hi(2) .eq. NoSlipWall) ) &
                  ugdnv(i,j,k) = ZERO
          end if
          if (idir .eq. 3) then
             if (k.eq.domlo(3) .and. &
                  (physbc_lo(3) .eq. Symmetry .or.  physbc_lo(3) .eq. SlipWall .or. &
                  physbc_lo(3) .eq. NoSlipWall) ) &
                  ugdnv(i,j,k) = ZERO
             if (k.eq.domhi(3)+1 .and. &
                  (physbc_hi(3) .eq. Symmetry .or.  physbc_hi(3) .eq. SlipWall .or. &
                  physbc_hi(3) .eq. NoSlipWall) ) &
                  ugdnv(i,j,k) = ZERO
          end if

          ! Compute fluxes, order as conserved state (not q)
          uflx(i,j,k,URHO) = rgdnv*ugdnv(i,j,k)

          if(idir.eq.1) then
             uflx(i,j,k,UMX) = uflx(i,j,k,URHO)*ugdnv(i,j,k) + pgdnv(i,j,k)
             uflx(i,j,k,UMY) = uflx(i,j,k,URHO)*v1gdnv
             uflx(i,j,k,UMZ) = uflx(i,j,k,URHO)*v2gdnv
          elseif(idir.eq.2) then
             uflx(i,j,k,UMX) = uflx(i,j,k,URHO)*v1gdnv
             uflx(i,j,k,UMY) = uflx(i,j,k,URHO)*ugdnv(i,j,k) + pgdnv(i,j,k)
             uflx(i,j,k,UMZ) = uflx(i,j,k,URHO)*v2gdnv
          else
             uflx(i,j,k,UMX) = uflx(i,j,k,URHO)*v1gdnv
             uflx(i,j,k,UMY) = uflx(i,j,k,URHO)*v2gdnv
             uflx(i,j,k,UMZ) = uflx(i,j,k,URHO)*ugdnv(i,j,k) + pgdnv(i,j,k)
          endif

          rhoetot = regdnv + HALF*rgdnv*(ugdnv(i,j,k)**2 + v1gdnv**2 + v2gdnv**2)

          uflx(i,j,k,UEDEN) = ugdnv(i,j,k)*(rhoetot + pgdnv(i,j,k))
          uflx(i,j,k,UEINT) = ugdnv(i,j,k)*regdnv

          ! Treat K as a passively advected quantity but allow it to affect fluxes of (rho E) and momenta.
          if (UESGS .gt. -1) then
             n  = UESGS
             nq = QESGS
             if (ustar .gt. ZERO) then
                qavg = ql(i,j,k,nq)
             else if (ustar .lt. ZERO) then
                qavg = qr(i,j,k,nq)
             else
                qavg = HALF * (ql(i,j,k,nq) + qr(i,j,k,nq))
             endif

             uflx(i,j,k,n) = uflx(i,j,k,URHO)*qavg

             rho_K_contrib =  TWO3RD * rgdnv * qavg

             if(idir.eq.1) then
                uflx(i,j,k,UMX) = uflx(i,j,k,UMX) + rho_K_contrib
             elseif(idir.eq.2) then
                uflx(i,j,k,UMY) = uflx(i,j,k,UMY) + rho_K_contrib
             elseif(idir.eq.3) then
                uflx(i,j,k,UMZ) = uflx(i,j,k,UMZ) + rho_K_contrib
             endif

             uflx(i,j,k,UEDEN) = uflx(i,j,k,UEDEN) + ugdnv(i,j,k) * rho_K_contrib
          end if

          do ipassive = 1, npassive
             n  = upass_map(ipassive)
             nq = qpass_map(ipassive)

             if (ustar .gt. ZERO) then
                uflx(i,j,k,n) = uflx(i,j,k,URHO)*ql(i,j,k,nq)
             else if (ustar .lt. ZERO) then
                uflx(i,j,k,n) = uflx(i,j,k,URHO)*qr(i,j,k,nq)
             else
                qavg = HALF * (ql(i,j,k,nq) + qr(i,j,k,nq))
                uflx(i,j,k,n) = uflx(i,j,k,URHO)*qavg
             endif
          enddo

       enddo
       enddo
    enddo

  end subroutine riemannus


  subroutine HLL(ql, qr, cl, cr, idir, f)

    use meth_params_module, only : QVAR, NVAR, QRHO, QU, QV, QW, QPRES, QREINT, &
                                   URHO, UMX, UMY, UMZ, UEDEN, UEINT, &
                                   npassive, upass_map, qpass_map

    use network, only : nspec, naux

    double precision, intent(in) :: ql(QVAR), qr(QVAR), cl, cr
    double precision, intent(inout) :: f(NVAR)
    integer, intent(in) :: idir

    integer :: ivel, ivelt, iveltt, imom, imomt, imomtt
    double precision :: a1, a4, bd, bl, bm, bp, br
    double precision :: cavg, uavg
    double precision :: fl_tmp, fr_tmp
    double precision :: rhod, rhoEl, rhoEr, rhol_sqrt, rhor_sqrt
    integer :: n, nq

    integer :: ipassive

    double precision, parameter :: small = 1.d-10

    select case (idir)
    case (1)
       ivel = QU
       ivelt = QV
       iveltt = QW

       imom = UMX
       imomt = UMY
       imomtt = UMZ

    case (2)
       ivel = QV
       ivelt = QU
       iveltt = QW

       imom = UMY
       imomt = UMX
       imomtt = UMZ

    case (3)
       ivel = QW
       ivelt = QU
       iveltt = QV

       imom = UMZ
       imomt = UMX
       imomtt = UMY

    end select

    rhol_sqrt = sqrt(ql(QRHO))
    rhor_sqrt = sqrt(qr(QRHO))

    rhod = ONE/(rhol_sqrt + rhor_sqrt)



    ! compute the average sound speed. This uses an approximation from
    ! E88, eq. 5.6, 5.7 that assumes gamma falls between 1
    ! and 5/3
    cavg = sqrt( (rhol_sqrt*cl**2 + rhor_sqrt*cr**2)*rhod + &
         HALF*rhol_sqrt*rhor_sqrt*rhod**2*(qr(ivel) - ql(ivel))**2 )


    ! Roe eigenvalues (E91, eq. 5.3b)
    uavg = (rhol_sqrt*ql(ivel) + rhor_sqrt*qr(ivel))*rhod

    a1 = uavg - cavg
    a4 = uavg + cavg


    ! signal speeds (E91, eq. 4.5)
    bl = min(a1, ql(ivel) - cl)
    br = max(a4, qr(ivel) + cr)

    bm = min(ZERO, bl)
    bp = max(ZERO, br)

    bd = bp - bm

    if (abs(bd) < small*max(abs(bm),abs(bp))) return

    bd = ONE/bd


    ! compute the fluxes according to E91, eq. 4.4b -- note that the
    ! min/max above picks the correct flux if we are not in the star
    ! region

    ! density flux
    fl_tmp = ql(QRHO)*ql(ivel)
    fr_tmp = qr(QRHO)*qr(ivel)

    f(URHO) = (bp*fl_tmp - bm*fr_tmp)*bd + bp*bm*bd*(qr(QRHO) - ql(QRHO))


    ! normal momentum flux -- we handle that separately
    fl_tmp = ql(QRHO)*ql(ivel)**2 + ql(QPRES)
    fr_tmp = qr(QRHO)*qr(ivel)**2 + qr(QPRES)

    f(imom) = (bp*fl_tmp - bm*fr_tmp)*bd + bp*bm*bd*(qr(QRHO)*qr(ivel) - ql(QRHO)*ql(ivel))


    ! transverse momentum flux
    fl_tmp = ql(QRHO)*ql(ivel)*ql(ivelt)
    fr_tmp = qr(QRHO)*qr(ivel)*qr(ivelt)

    f(imomt) = (bp*fl_tmp - bm*fr_tmp)*bd + bp*bm*bd*(qr(QRHO)*qr(ivelt) - ql(QRHO)*ql(ivelt))


    fl_tmp = ql(QRHO)*ql(ivel)*ql(iveltt)
    fr_tmp = qr(QRHO)*qr(ivel)*qr(iveltt)

    f(imomtt) = (bp*fl_tmp - bm*fr_tmp)*bd + bp*bm*bd*(qr(QRHO)*qr(iveltt) - ql(QRHO)*ql(iveltt))


    ! total energy flux
    rhoEl = ql(QREINT) + HALF*ql(QRHO)*(ql(ivel)**2 + ql(ivelt)**2 + ql(iveltt)**2)
    fl_tmp = ql(ivel)*(rhoEl + ql(QPRES))

    rhoEr = qr(QREINT) + HALF*qr(QRHO)*(qr(ivel)**2 + qr(ivelt)**2 + qr(iveltt)**2)
    fr_tmp = qr(ivel)*(rhoEr + qr(QPRES))

    f(UEDEN) = (bp*fl_tmp - bm*fr_tmp)*bd + bp*bm*bd*(rhoEr - rhoEl)


    ! eint flux
    fl_tmp = ql(QREINT)*ql(ivel)
    fr_tmp = qr(QREINT)*qr(ivel)

    f(UEINT) = (bp*fl_tmp - bm*fr_tmp)*bd + bp*bm*bd*(qr(QREINT) - ql(QREINT))


    ! passively-advected scalar fluxes
    do ipassive = 1, npassive
       n  = upass_map(ipassive)
       nq = qpass_map(ipassive)

       fl_tmp = ql(QRHO)*ql(nq)*ql(ivel)
       fr_tmp = qr(QRHO)*qr(nq)*qr(ivel)

       f(n) = (bp*fl_tmp - bm*fr_tmp)*bd + bp*bm*bd*(qr(QRHO)*qr(nq) - ql(QRHO)*ql(nq))
    enddo

  end subroutine HLL

end module riemann_module
