module transverse_module

  use bl_constants_module

  implicit none

contains

  !
  ! lo and hi in this module are cell-centered index of valid region
  ! 

  !===========================================================================
  ! transx
  !===========================================================================
  subroutine transx(qym,qymo,qyp,qypo, &
                    qzm,qzmo,qzp,qzpo, qlo, qhi, &
                    fx, flo, fhi, &
                    ugdnvx,pgdnvx,gegdnvx, gdlo, gdhi, &
                    gamc, gclo, gchi, &
                    cdtdx, lo, hi)
    use network, only : nspec, naux
    use meth_params_module, only : QVAR, NVAR, QRHO, QU, QV, QW, &
                                   QPRES, QREINT, QGAME, QESGS, QFA, QFS, &
                                   URHO, UMX, UMY, UMZ, UEDEN, UEINT, UESGS, UFA, UFS, &
                                   nadv, small_pres, small_temp, &
                                   npassive, upass_map, qpass_map, &
                                   ppm_predict_gammae, &
                                   transverse_use_eos, transverse_reset_density, transverse_reset_rhoe
    use eos_module
    implicit none

    integer, intent(in) :: qlo(3), qhi(3), flo(3), fhi(3), gdlo(3), gdhi(3), &
         gclo(3), gchi(3), lo(3), hi(3)
    double precision, intent(in) :: cdtdx
    double precision,intent(in )::    qym( qlo(1): qhi(1), qlo(2): qhi(2), qlo(3): qhi(3),QVAR)
    double precision,intent(in )::    qyp( qlo(1): qhi(1), qlo(2): qhi(2), qlo(3): qhi(3),QVAR)
    double precision,intent(in )::    qzm( qlo(1): qhi(1), qlo(2): qhi(2), qlo(3): qhi(3),QVAR)
    double precision,intent(in )::    qzp( qlo(1): qhi(1), qlo(2): qhi(2), qlo(3): qhi(3),QVAR)
    double precision,intent(in )::     fx( flo(1): fhi(1), flo(2): fhi(2), flo(3): fhi(3),NVAR)
    double precision,intent(in ):: ugdnvx(gdlo(1):gdhi(1),gdlo(2):gdhi(2),gdlo(3):gdhi(3))
    double precision,intent(in ):: pgdnvx(gdlo(1):gdhi(1),gdlo(2):gdhi(2),gdlo(3):gdhi(3))
    double precision,intent(in )::gegdnvx(gdlo(1):gdhi(1),gdlo(2):gdhi(2),gdlo(3):gdhi(3))
    double precision,intent(in )::   gamc(gclo(1):gchi(1),gclo(2):gchi(2),gclo(3):gchi(3))
    double precision,intent(out)::   qymo( qlo(1): qhi(1), qlo(2): qhi(2), qlo(3): qhi(3),QVAR)
    double precision,intent(out)::   qypo( qlo(1): qhi(1), qlo(2): qhi(2), qlo(3): qhi(3),QVAR)
    double precision,intent(out)::   qzmo( qlo(1): qhi(1), qlo(2): qhi(2), qlo(3): qhi(3),QVAR)
    double precision,intent(out)::   qzpo( qlo(1): qhi(1), qlo(2): qhi(2), qlo(3): qhi(3),QVAR)

    integer :: i,j,k,n,nq,ipassive

    double precision rrnew, compu
    double precision rrry, rrly, rrrz, rrlz
    double precision rury, ruly, rurz, rulz
    double precision rvry, rvly, rvrz, rvlz
    double precision rwry, rwly, rwrz, rwlz
    double precision ekenry, ekenly, ekenrz, ekenlz
    double precision rery, rely, rerz, relz
    double precision rrnewry, rrnewly, rrnewrz, rrnewlz
    double precision runewry, runewly, runewrz, runewlz
    double precision rvnewry, rvnewly, rvnewrz, rvnewlz
    double precision rwnewry, rwnewly, rwnewrz, rwnewlz
    double precision renewry, renewly, renewrz, renewlz
    double precision pnewry, pnewly, pnewrz, pnewlz
    double precision rhoekenry, rhoekenly, rhoekenrz, rhoekenlz
    double precision pgp, pgm, ugp, ugm, gegp, gegm, dup, pav, du, dge, uav, geav
    
    type (eos_t) :: eos_state

    ! work on qy* first

    !-------------------------------------------------------------------------
    ! update all of the passively-advected quantities with the
    ! transerse term and convert back to the primitive quantity
    !-------------------------------------------------------------------------

    do ipassive = 1,npassive
       n  = upass_map(ipassive)
       nq = qpass_map(ipassive)

       do       k = lo(3)-1, hi(3)+1
          do    j = lo(2)  , hi(2)+1 
             do i = lo(1)  , hi(1)             
                rrnew = qyp(i,j,k,QRHO) - cdtdx*(fx(i+1,j,k,URHO) - fx(i,j,k,URHO))
                compu = qyp(i,j,k,QRHO)*qyp(i,j,k,nq) - cdtdx*(fx(i+1,j,k,n) - fx(i,j,k,n))
                qypo(i,j,k,nq) = compu/rrnew             
             enddo
          enddo
          do    j = lo(2)-1, hi(2)
             do i = lo(1)  , hi(1)
                rrnew = qym(i,j+1,k,QRHO) - cdtdx*(fx(i+1,j,k,URHO) - fx(i,j,k,URHO))
                compu = qym(i,j+1,k,QRHO)*qym(i,j+1,k,nq) - cdtdx*(fx(i+1,j,k,n) - fx(i,j,k,n))
                qymo(i,j+1,k,nq) = compu/rrnew
             enddo
          enddo
       enddo
    enddo

    !-------------------------------------------------------------------
    ! add the transverse flux difference in the x-direction to y-states
    ! for the fluid variables
    !-------------------------------------------------------------------

    do k = lo(3)-1, hi(3)+1

       !----------------------------------------------------------------
       ! qypo state
       !----------------------------------------------------------------

       do    j = lo(2), hi(2)+1 
          do i = lo(1), hi(1)             
             
             pgp  =  pgdnvx(i+1,j,k)
             pgm  =  pgdnvx(i  ,j,k)
             ugp  =  ugdnvx(i+1,j,k)
             ugm  =  ugdnvx(i  ,j,k)
             gegp = gegdnvx(i+1,j,k)
             gegm = gegdnvx(i  ,j,k)

             ! we need to augment our conserved system with either a p
             ! equation or gammae (if we have ppm_predict_gammae = 1) to
             ! be able to deal with the general EOS

             dup = pgp*ugp - pgm*ugm
             pav = HALF*(pgp+pgm)
             uav = HALF*(ugp+ugm)
             geav = HALF*(gegp+gegm)
             du = ugp-ugm
             dge = gegp-gegm
          
             ! Convert to conservation form
             rrry = qyp(i,j,k,QRHO)
             rury = rrry*qyp(i,j,k,QU)
             rvry = rrry*qyp(i,j,k,QV)
             rwry = rrry*qyp(i,j,k,QW)
             ekenry = HALF*rrry*(qyp(i,j,k,QU)**2 + qyp(i,j,k,QV)**2 + qyp(i,j,k,QW)**2)
             rery = qyp(i,j,k,QREINT) + ekenry
          
             ! Add transverse predictor
             rrnewry = rrry - cdtdx*(fx(i+1,j,k,URHO ) - fx(i,j,k,URHO))
             runewry = rury - cdtdx*(fx(i+1,j,k,UMX  ) - fx(i,j,k,UMX))          
             rvnewry = rvry - cdtdx*(fx(i+1,j,k,UMY  ) - fx(i,j,k,UMY))
             rwnewry = rwry - cdtdx*(fx(i+1,j,k,UMZ  ) - fx(i,j,k,UMZ))
             renewry = rery - cdtdx*(fx(i+1,j,k,UEDEN) - fx(i,j,k,UEDEN))
          
             ! Reset to original value if adding transverse terms made density negative
             if (transverse_reset_density == 1 .and. rrnewry .lt. ZERO) then
                rrnewry = rrry
                runewry = rury
                rvnewry = rvry
                rwnewry = rwry
                renewry = rery
             endif

             ! Convert back to primitive form
             qypo(i,j,k,QRHO) = rrnewry
             qypo(i,j,k,QU) = runewry/qypo(i,j,k,QRHO)
             qypo(i,j,k,QV) = rvnewry/qypo(i,j,k,QRHO)
             qypo(i,j,k,QW) = rwnewry/qypo(i,j,k,QRHO)

             ! note: we run the risk of (rho e) being negative here
             rhoekenry = HALF*(runewry**2 + rvnewry**2 + rwnewry**2)/qypo(i,j,k,QRHO)
             qypo(i,j,k,QREINT) = renewry - rhoekenry

             if (transverse_reset_rhoe == 1) then
                ! If it is negative, reset the internal energy by
                ! using the discretized expression for updating (rho e).

                if (qypo(i,j,k,QREINT) .le. ZERO) then
                   qypo(i,j,k,QREINT) = qyp(i,j,k,QREINT) - &
                        cdtdx*(fx(i+1,j,k,UEINT) - fx(i,j,k,UEINT) + pav*du)
                   
                   ! if we are still negative, then we need to reset
                   if (qypo(i,j,k,QREINT) < ZERO) then
                      eos_state % rho = qypo(i,j,k,QRHO)
                      eos_state % T = small_temp
                      eos_state % xn(:) = qypo(i,j,k,QFS:QFS-1+nspec)
                      
                      call eos(eos_input_rt, eos_state)
                      
                      qypo(i,j,k,QREINT) = qypo(i,j,k,QRHO)*eos_state % e
                      qypo(i,j,k,QPRES) = eos_state % p
                   endif
                endif
             endif

             if (ppm_predict_gammae == 0) then
                
                ! Optionally, use the EOS to calculate the pressure.

                if (transverse_use_eos .eq. 1) then
                   eos_state % rho = qypo(i,j,k,QRHO)
                   eos_state % e   = qypo(i,j,k,QREINT) / qypo(i,j,k,QRHO)
                   eos_state % T   = small_temp
                   eos_state % xn  = qypo(i,j,k,QFS:QFS+nspec-1)
                
                   call eos(eos_input_re, eos_state)
                
                   pnewry = eos_state % p
                   qypo(i,j,k,QREINT) = eos_state % e * eos_state % rho
                else
                   ! add the transverse term to the p evolution eq here
                   pnewry = qyp(i,j,k,QPRES) - cdtdx*(dup + pav*du*(gamc(i,j,k)-ONE))
                endif

                qypo(i,j,k,QPRES) = max(pnewry,small_pres)

             else

                ! Update gammae with its transverse terms
                qypo(i,j,k,QGAME) = qyp(i,j,k,QGAME) + &
                     cdtdx*( (geav-ONE)*(geav-gamc(i,j,k))*du - uav*dge )

                ! and compute the p edge state from this and (rho e)
                qypo(i,j,k,QPRES) = qypo(i,j,k,QREINT)*(qypo(i,j,k,QGAME)-ONE)
                qypo(i,j,k,QPRES) = max(qypo(i,j,k,QPRES),small_pres)
             endif

          end do
       end do

       !-------------------------------------------------------------------   
       ! qymo state
       !-------------------------------------------------------------------
       
       do    j = lo(2), hi(2)+1 
          do i = lo(1), hi(1)             
             
             pgp  =  pgdnvx(i+1,j,k)
             pgm  =  pgdnvx(i  ,j,k)
             ugp  =  ugdnvx(i+1,j,k)
             ugm  =  ugdnvx(i  ,j,k)
             gegp = gegdnvx(i+1,j,k)
             gegm = gegdnvx(i  ,j,k)

             ! we need to augment our conserved system with either a p
             ! equation or gammae (if we have ppm_predict_gammae = 1) to
             ! be able to deal with the general EOS

             dup = pgp*ugp - pgm*ugm
             pav = HALF*(pgp+pgm)
             uav = HALF*(ugp+ugm)
             geav = HALF*(gegp+gegm)
             du = ugp-ugm
             dge = gegp-gegm

             ! Convert to conservation form
             rrly = qym(i,j+1,k,QRHO)
             ruly = rrly*qym(i,j+1,k,QU)
             rvly = rrly*qym(i,j+1,k,QV)
             rwly = rrly*qym(i,j+1,k,QW)
             ekenly = HALF*rrly*(qym(i,j+1,k,QU)**2 + qym(i,j+1,k,QV)**2 + qym(i,j+1,k,QW)**2)
             rely = qym(i,j+1,k,QREINT) + ekenly
          
             ! Add transverse predictor
             rrnewly = rrly - cdtdx*(fx(i+1,j,k,URHO ) - fx(i,j,k,URHO))
             runewly = ruly - cdtdx*(fx(i+1,j,k,UMX  ) - fx(i,j,k,UMX))
             rvnewly = rvly - cdtdx*(fx(i+1,j,k,UMY  ) - fx(i,j,k,UMY))
             rwnewly = rwly - cdtdx*(fx(i+1,j,k,UMZ  ) - fx(i,j,k,UMZ))
             renewly = rely - cdtdx*(fx(i+1,j,k,UEDEN) - fx(i,j,k,UEDEN))

             ! Reset to original value if adding transverse terms made density negative
             if (transverse_reset_density == 1 .and. rrnewly .lt. ZERO) then
                rrnewly = rrly 
                runewly = ruly 
                rvnewly = rvly 
                rwnewly = rwly 
                renewly = rely 
             end if

             ! Convert back to primitive form
             qymo(i,j+1,k,QRHO) = rrnewly
             qymo(i,j+1,k,QU) = runewly/qymo(i,j+1,k,QRHO)
             qymo(i,j+1,k,QV) = rvnewly/qymo(i,j+1,k,QRHO)
             qymo(i,j+1,k,QW) = rwnewly/qymo(i,j+1,k,QRHO)

             ! note: we run the risk of (rho e) being negative here
             rhoekenly = HALF*(runewly**2 + rvnewly**2 + rwnewly**2)/qymo(i,j+1,k,QRHO)
             qymo(i,j+1,k,QREINT) = renewly - rhoekenly
             
             if (transverse_reset_rhoe == 1) then
                ! If it is negative, reset the internal energy by using the discretized
                ! expression for updating (rho e).

                if (qymo(i,j+1,k,QREINT) .le. ZERO) then
                   qymo(i,j+1,k,QREINT) = qym(i,j+1,k,QREINT) - &
                        cdtdx*(fx(i+1,j,k,UEINT) - fx(i,j,k,UEINT) + pav*du)
                   
                   ! if we are still negative, then we need to reset
                   if (qymo(i,j+1,k,QREINT) < ZERO) then
                      eos_state % rho = qymo(i,j+1,k,QRHO)
                      eos_state % T = small_temp
                      eos_state % xn(:) = qymo(i,j+1,k,QFS:QFS-1+nspec)
                      
                      call eos(eos_input_rt, eos_state)
                      
                      qymo(i,j+1,k,QREINT) = qymo(i,j+1,k,QRHO)*eos_state % e
                      qymo(i,j+1,k,QPRES) = eos_state % p
                   endif
                endif
             endif

             if (ppm_predict_gammae == 0) then

                ! Optionally, use the EOS to calculate the pressure.

                if (transverse_use_eos .eq. 1) then
                   eos_state % rho = qymo(i,j+1,k,QRHO)
                   eos_state % e   = qymo(i,j+1,k,QREINT) / qymo(i,j+1,k,QRHO)
                   eos_state % T   = small_temp
                   eos_state % xn  = qymo(i,j+1,k,QFS:QFS+nspec-1)

                   call eos(eos_input_re, eos_state)
  
                   pnewly = eos_state % p
                   qymo(i,j+1,k,QREINT) = eos_state % e * eos_state % rho
                else
                   pnewly = qym(i,j+1,k,QPRES) - cdtdx*(dup + pav*du*(gamc(i,j,k)-ONE))
                endif
                
                qymo(i,j+1,k,QPRES) = max(pnewly,small_pres)

             else
                
                ! Update gammae with its transverse terms
                qymo(i,j+1,k,QGAME) = qym(i,j+1,k,QGAME) + &
                     cdtdx*( (geav-ONE)*(geav-gamc(i,j,k))*du - uav*dge )

                ! and compute the p edge state from this and (rho e)
                qymo(i,j+1,k,QPRES) = qymo(i,j+1,k,QREINT)*(qymo(i,j+1,k,QGAME)-ONE)
                qymo(i,j+1,k,QPRES) = max(qymo(i,j+1,k,QPRES), small_pres)

             end if
          
          enddo
       enddo
    enddo

    ! work on qz*

    !-------------------------------------------------------------------------
    ! update all of the passively-advected quantities with the
    ! transerse term and convert back to the primitive quantity
    !-------------------------------------------------------------------------

    do ipassive = 1,npassive
       n  = upass_map(ipassive)
       nq = qpass_map(ipassive)

       do       k = lo(3)  , hi(3)+1
          do    j = lo(2)-1, hi(2)+1 
             do i = lo(1)  , hi(1)             
                rrnew = qzp(i,j,k,QRHO) - cdtdx*(fx(i+1,j,k,URHO) - fx(i,j,k,URHO))
                compu = qzp(i,j,k,QRHO)*qzp(i,j,k,nq) - cdtdx*(fx(i+1,j,k,n) - fx(i,j,k,n))
                qzpo(i,j,k,nq) = compu/rrnew             
             enddo
          enddo
       end do
       
       do       k = lo(3)-1, hi(3)
          do    j = lo(2)-1, hi(2)+1
             do i = lo(1)  , hi(1)
             rrnew = qzm(i,j,k+1,QRHO) - cdtdx*(fx(i+1,j,k,URHO) - fx(i,j,k,URHO))
             compu = qzm(i,j,k+1,QRHO)*qzm(i,j,k+1,nq) - cdtdx*(fx(i+1,j,k,n) - fx(i,j,k,n))
             qzmo(i,j,k+1,nq) = compu/rrnew
          enddo
          enddo
       enddo
    enddo

    !-------------------------------------------------------------------
    ! add the transverse flux difference in the x-direction to z-states
    ! for the fluid variables
    !-------------------------------------------------------------------

    !-------------------------------------------------------------------
    ! qzpo state
    !-------------------------------------------------------------------
    
    do       k = lo(3)  , hi(3)+1
       do    j = lo(2)-1, hi(2)+1 
          do i = lo(1)  , hi(1)             

             pgp = pgdnvx(i+1,j,k)
             pgm = pgdnvx(i,j,k)
             ugp = ugdnvx(i+1,j,k)
             ugm = ugdnvx(i,j,k)
             gegp = gegdnvx(i+1,j,k)
             gegm = gegdnvx(i,j,k)
             
             ! we need to augment our conserved system with either a p
             ! equation or gammae (if we have ppm_predict_gammae = 1) to
             ! be able to deal with the general EOS
             
             dup = pgp*ugp - pgm*ugm
             pav = HALF*(pgp+pgm)
             uav = HALF*(ugp+ugm)
             geav = HALF*(gegp+gegm)
             du = ugp-ugm
             dge = gegp-gegm
             
             ! Convert to conservation form
             rrrz = qzp(i,j,k,QRHO)
             rurz = rrrz*qzp(i,j,k,QU)
             rvrz = rrrz*qzp(i,j,k,QV)
             rwrz = rrrz*qzp(i,j,k,QW)
             ekenrz = HALF*rrrz*(qzp(i,j,k,QU)**2 + qzp(i,j,k,QV)**2 + qzp(i,j,k,QW)**2)
             rerz = qzp(i,j,k,QREINT) + ekenrz
             
             ! Add transverse predictor
             rrnewrz = rrrz - cdtdx*(fx(i+1,j,k,URHO) - fx(i,j,k,URHO))
             runewrz = rurz - cdtdx*(fx(i+1,j,k,UMX) - fx(i,j,k,UMX))
             rvnewrz = rvrz - cdtdx*(fx(i+1,j,k,UMY) - fx(i,j,k,UMY))
             rwnewrz = rwrz - cdtdx*(fx(i+1,j,k,UMZ) - fx(i,j,k,UMZ))
             renewrz = rerz - cdtdx*(fx(i+1,j,k,UEDEN) - fx(i,j,k,UEDEN))

             ! Reset to original value if adding transverse terms made density negative
             if (transverse_reset_density == 1 .and. rrnewrz .lt. ZERO) then
                rrnewrz = rrrz 
                runewrz = rurz 
                rvnewrz = rvrz 
                rwnewrz = rwrz 
                renewrz = rerz 
             endif
                 
             ! Convert back to primitive form
             qzpo(i,j,k,QRHO) = rrnewrz
             qzpo(i,j,k,QU) = runewrz/qzpo(i,j,k,QRHO)
             qzpo(i,j,k,QV) = rvnewrz/qzpo(i,j,k,QRHO)
             qzpo(i,j,k,QW) = rwnewrz/qzpo(i,j,k,QRHO)
             
             ! note: we run the risk of (rho e) being negative here
             rhoekenrz = HALF*(runewrz**2 + rvnewrz**2 + rwnewrz**2)/qzpo(i,j,k,QRHO)
             qzpo(i,j,k,QREINT) = renewrz - rhoekenrz
             
             if (transverse_reset_rhoe == 1) then
                ! If it is negative, reset the internal energy by using the discretized
                ! expression for updating (rho e).

                if (qzpo(i,j,k,QREINT) .le. ZERO) then
                   qzpo(i,j,k,QREINT) = qzp(i,j,k,QREINT) - &
                        cdtdx*(fx(i+1,j,k,UEINT) - fx(i,j,k,UEINT) + pav*du)
                   
                   ! if we are still negative, then we need to reset
                   if (qzpo(i,j,k,QREINT) < ZERO) then
                      eos_state % rho = qzpo(i,j,k,QRHO)
                      eos_state % T = small_temp
                      eos_state % xn(:) = qzpo(i,j,k,QFS:QFS-1+nspec)
                      
                      call eos(eos_input_rt, eos_state)
                      
                      qzpo(i,j,k,QREINT) = qzpo(i,j,k,QRHO)*eos_state % e
                      qzpo(i,j,k,QPRES) = eos_state % p
                   endif
                endif
             endif
             
             if (ppm_predict_gammae == 0) then
                
                ! Optionally, use the EOS to calculate the pressure.
                
                if (transverse_use_eos .eq. 1) then
                   eos_state % rho = qzpo(i,j,k,QRHO)
                   eos_state % e   = qzpo(i,j,k,QREINT) / qzpo(i,j,k,QRHO)
                   eos_state % T   = small_temp
                   eos_state % xn  = qzpo(i,j,k,QFS:QFS+nspec-1)
                   
                   call eos(eos_input_re, eos_state)
                   
                   pnewrz = eos_state % p
                   qzpo(i,j,k,QREINT) = eos_state % e * eos_state % rho
                else
                   ! add the transverse term to the p evolution eq here
                   pnewrz = qzp(i,j,k,QPRES) - cdtdx*(dup + pav*du*(gamc(i,j,k)-ONE))
                endif
                
                qzpo(i,j,k,QPRES) = max(pnewrz,small_pres)
                
             else
                
                ! Update gammae with its transverse terms
                qzpo(i,j,k,QGAME) = qzp(i,j,k,QGAME) + &
                     cdtdx*( (geav-ONE)*(geav-gamc(i,j,k))*du - uav*dge )
                
                ! and compute the p edge state from this and (rho e)
                qzpo(i,j,k,QPRES) = qzpo(i,j,k,QREINT)*(qzpo(i,j,k,QGAME)-ONE)
                qzpo(i,j,k,QPRES) = max(qzpo(i,j,k,QPRES), small_pres)
                
             endif

          end do
       end do
    end do

    !-------------------------------------------------------------------
    ! qzmo state
    !-------------------------------------------------------------------
    
    do       k = lo(3)-1, hi(3)
       do    j = lo(2)-1, hi(2)+1 
          do i = lo(1)  , hi(1)             
          
             pgp = pgdnvx(i+1,j,k)
             pgm = pgdnvx(i,j,k)
             ugp = ugdnvx(i+1,j,k)
             ugm = ugdnvx(i,j,k)
             gegp = gegdnvx(i+1,j,k)
             gegm = gegdnvx(i,j,k)

             ! we need to augment our conserved system with either a p
             ! equation or gammae (if we have ppm_predict_gammae = 1) to
             ! be able to deal with the general EOS
             
             dup = pgp*ugp - pgm*ugm
             pav = HALF*(pgp+pgm)
             uav = HALF*(ugp+ugm)
             geav = HALF*(gegp+gegm)
             du = ugp-ugm         
             dge = gegp-gegm
          
             ! Convert to conservation form
             rrlz = qzm(i,j,k+1,QRHO)
             rulz = rrlz*qzm(i,j,k+1,QU)
             rvlz = rrlz*qzm(i,j,k+1,QV)
             rwlz = rrlz*qzm(i,j,k+1,QW)
             ekenlz = HALF*rrlz*(qzm(i,j,k+1,QU)**2 + qzm(i,j,k+1,QV)**2 + qzm(i,j,k+1,QW)**2)
             relz = qzm(i,j,k+1,QREINT) + ekenlz
             
             ! Add transverse predictor
             rrnewlz = rrlz - cdtdx*(fx(i+1,j,k,URHO) - fx(i,j,k,URHO))
             runewlz = rulz - cdtdx*(fx(i+1,j,k,UMX) - fx(i,j,k,UMX))
             rvnewlz = rvlz - cdtdx*(fx(i+1,j,k,UMY) - fx(i,j,k,UMY))
             rwnewlz = rwlz - cdtdx*(fx(i+1,j,k,UMZ) - fx(i,j,k,UMZ))
             renewlz = relz - cdtdx*(fx(i+1,j,k,UEDEN) - fx(i,j,k,UEDEN))
             
             ! Reset to original value if adding transverse terms made density negative
             if (transverse_reset_density == 1 .and. rrnewlz .lt. ZERO) then
                rrnewlz = rrlz
                runewlz = rulz
                rvnewlz = rvlz
                rwnewlz = rwlz
                renewlz = relz
             endif

             ! Convert back to primitive form
             qzmo(i,j,k+1,QRHO) = rrnewlz
             qzmo(i,j,k+1,QU) = runewlz/qzmo(i,j,k+1,QRHO)
             qzmo(i,j,k+1,QV) = rvnewlz/qzmo(i,j,k+1,QRHO)
             qzmo(i,j,k+1,QW) = rwnewlz/qzmo(i,j,k+1,QRHO)
             
             ! note: we run the risk of (rho e) being negative here
             rhoekenlz = HALF*(runewlz**2 + rvnewlz**2 + rwnewlz**2)/qzmo(i,j,k+1,QRHO)
             qzmo(i,j,k+1,QREINT) = renewlz - rhoekenlz
             
             if (transverse_reset_rhoe == 1) then
                ! If it is negative, reset the internal energy by using the discretized
                ! expression for updating (rho e).

                if (qzmo(i,j,k+1,QREINT) .le. ZERO) then
                   qzmo(i,j,k+1,QREINT) = qzm(i,j,k+1,QREINT) - &
                        cdtdx*(fx(i+1,j,k,UEINT) - fx(i,j,k,UEINT) + pav*du)
                
                   ! if we are still negative, then we need to reset
                   if (qzmo(i,j,k+1,QREINT) < ZERO) then
                      eos_state % rho = qzmo(i,j,k+1,QRHO)
                      eos_state % T = small_temp
                      eos_state % xn(:) = qzmo(i,j,k+1,QFS:QFS-1+nspec)
                      
                      call eos(eos_input_rt, eos_state)
                      
                      qzmo(i,j,k+1,QREINT) = qzmo(i,j,k+1,QRHO)*eos_state % e
                      qzmo(i,j,k+1,QPRES) = eos_state % p
                   endif
                endif
             endif

             if (ppm_predict_gammae == 0) then
                
                ! Optionally, use the EOS to calculate the pressure.
                
                if (transverse_use_eos .eq. 1) then
                   eos_state % rho = qzmo(i,j,k+1,QRHO)
                   eos_state % e   = qzmo(i,j,k+1,QREINT) / qzmo(i,j,k+1,QRHO)
                   eos_state % T   = small_temp
                   eos_state % xn  = qzmo(i,j,k+1,QFS:QFS+nspec-1)
                   
                   call eos(eos_input_re, eos_state)
                   
                   pnewlz = eos_state % p
                   qzmo(i,j,k+1,QREINT) = eos_state % e * eos_state % rho
                else
                   ! add the transverse term to the p evolution eq here
                   pnewlz = qzm(i,j,k+1,QPRES) - cdtdx*(dup + pav*du*(gamc(i,j,k)-ONE))
                endif
                
                qzmo(i,j,k+1,QPRES) = max(pnewlz,small_pres)
                
             else
                
                ! Update gammae with its transverse terms
                qzmo(i,j,k+1,QGAME) = qzm(i,j,k+1,QGAME) + &
                     cdtdx*( (geav-ONE)*(geav-gamc(i,j,k))*du - uav*dge )
                
                ! and compute the p edge state from this and (rho e)
                qzmo(i,j,k+1,QPRES) = qzmo(i,j,k+1,QREINT)*(qzmo(i,j,k+1,QGAME)-ONE)
                qzmo(i,j,k+1,QPRES) = max(qzmo(i,j,k+1,QPRES), small_pres)
                
             endif

          end do
       end do
    end do

  end subroutine transx


  !===========================================================================
  ! transy
  !===========================================================================
  subroutine transy(qxm,qxmo,qxp,qxpo, &
                    qzm,qzmo,qzp,qzpo, qlo, qhi, &
                    fy, flo, fhi, &
                    ugdnvy,pgdnvy,gegdnvy, gdlo, gdhi, &
                    gamc, gclo, gchi, &
                    cdtdy, lo, hi)
    use network, only : nspec, naux
    use meth_params_module, only : QVAR, NVAR, QRHO, QU, QV, QW, &
                                   QPRES, QREINT, QGAME, QESGS, QFA, QFS, &
                                   URHO, UMX, UMY, UMZ, UEDEN, UEINT, UESGS, UFA, UFS, &
                                   nadv, small_pres, small_temp, &
                                   npassive, upass_map, qpass_map, &
                                   ppm_predict_gammae, &
                                   transverse_use_eos, transverse_reset_density, transverse_reset_rhoe
    use eos_module
    implicit none

    integer, intent(in) :: qlo(3), qhi(3), flo(3), fhi(3), gdlo(3), gdhi(3), &
         gclo(3), gchi(3), lo(3), hi(3)
    double precision, intent(in) :: cdtdy
    double precision,intent(in )::    qxm( qlo(1): qhi(1), qlo(2): qhi(2), qlo(3): qhi(3),QVAR)
    double precision,intent(in )::    qxp( qlo(1): qhi(1), qlo(2): qhi(2), qlo(3): qhi(3),QVAR)
    double precision,intent(in )::    qzm( qlo(1): qhi(1), qlo(2): qhi(2), qlo(3): qhi(3),QVAR)
    double precision,intent(in )::    qzp( qlo(1): qhi(1), qlo(2): qhi(2), qlo(3): qhi(3),QVAR)
    double precision,intent(in )::     fy( flo(1): fhi(1), flo(2): fhi(2), flo(3): fhi(3),NVAR)
    double precision,intent(in ):: ugdnvy(gdlo(1):gdhi(1),gdlo(2):gdhi(2),gdlo(3):gdhi(3))
    double precision,intent(in ):: pgdnvy(gdlo(1):gdhi(1),gdlo(2):gdhi(2),gdlo(3):gdhi(3))
    double precision,intent(in )::gegdnvy(gdlo(1):gdhi(1),gdlo(2):gdhi(2),gdlo(3):gdhi(3))
    double precision,intent(in )::   gamc(gclo(1):gchi(1),gclo(2):gchi(2),gclo(3):gchi(3))
    double precision,intent(out)::   qxmo( qlo(1): qhi(1), qlo(2): qhi(2), qlo(3): qhi(3),QVAR)
    double precision,intent(out)::   qxpo( qlo(1): qhi(1), qlo(2): qhi(2), qlo(3): qhi(3),QVAR)
    double precision,intent(out)::   qzmo( qlo(1): qhi(1), qlo(2): qhi(2), qlo(3): qhi(3),QVAR)
    double precision,intent(out)::   qzpo( qlo(1): qhi(1), qlo(2): qhi(2), qlo(3): qhi(3),QVAR)

    integer :: i,j,k,n,nq,ipassive

    double precision rrnew, compu
    double precision rrrx, rrlx, rrrz, rrlz
    double precision rurx, rulx, rurz, rulz
    double precision rvrx, rvlx, rvrz, rvlz
    double precision rwrx, rwlx, rwrz, rwlz
    double precision ekenrx, ekenlx, ekenrz, ekenlz
    double precision rerx, relx, rerz, relz
    double precision rrnewrx, rrnewlx, rrnewrz, rrnewlz    
    double precision runewrx, runewlx, runewrz, runewlz    
    double precision rvnewrx, rvnewlx, rvnewrz, rvnewlz    
    double precision rwnewrx, rwnewlx, rwnewrz, rwnewlz    
    double precision renewrx, renewlx, renewrz, renewlz    
    double precision pnewrx, pnewlx, pnewrz, pnewlz      
    double precision rhoekenrx, rhoekenlx, rhoekenrz, rhoekenlz
    double precision pgp, pgm, ugp, ugm, gegp, gegm, dup, pav, du, dge, uav, geav
    
    type (eos_t) :: eos_state

    ! work on qx* first

    !-------------------------------------------------------------------------    
    ! update all of the passively-advected quantities with the
    ! transerse term and convert back to the primitive quantity
    !-------------------------------------------------------------------------

    do ipassive = 1,npassive
       n  = upass_map(ipassive)
       nq = qpass_map(ipassive)

       do       k = lo(3)-1, hi(3)+1
          do    j = lo(2)  , hi(2) 
             do i = lo(1)  , hi(1)+1
                rrnew = qxp(i,j,k,QRHO) - cdtdy*(fy(i,j+1,k,URHO) - fy(i,j,k,URHO))
                compu = qxp(i,j,k,QRHO)*qxp(i,j,k,nq) - cdtdy*(fy(i,j+1,k,n) - fy(i,j,k,n))
                qxpo(i,j,k,nq) = compu/rrnew             
             enddo
             do i = lo(1)-1, hi(1)
                rrnew = qxm(i+1,j,k,QRHO) - cdtdy*(fy(i,j+1,k,URHO) - fy(i,j,k,URHO))
                compu = qxm(i+1,j,k,QRHO)*qxm(i+1,j,k,nq) - cdtdy*(fy(i,j+1,k,n) - fy(i,j,k,n))
                qxmo(i+1,j,k,nq) = compu/rrnew             
             end do
          end do
       end do
    enddo

    !-------------------------------------------------------------------
    ! add the transverse flux difference in the y-direction to x-states
    ! for the fluid variables
    !-------------------------------------------------------------------
          
    do    k = lo(3)-1, hi(3)+1
       do j = lo(2)  , hi(2) 

          !-------------------------------------------------------------------
          ! qxpo state
          !-------------------------------------------------------------------    
          
          do i = lo(1), hi(1)+1

             pgp = pgdnvy(i,j+1,k)
             pgm = pgdnvy(i,j,k)
             ugp = ugdnvy(i,j+1,k)
             ugm = ugdnvy(i,j,k)
             gegp = gegdnvy(i,j+1,k)
             gegm = gegdnvy(i,j,k)

             ! we need to augment our conserved system with either a p
             ! equation or gammae (if we have ppm_predict_gammae = 1) to
             ! be able to deal with the general EOS
             
             dup = pgp*ugp - pgm*ugm
             pav = HALF*(pgp+pgm)
             uav = HALF*(ugp+ugm)
             geav = HALF*(gegp+gegm)
             du = ugp-ugm
             dge = gegp-gegm
             
             ! Convert to conservation form
             rrrx = qxp(i,j,k,QRHO)
             rurx = rrrx*qxp(i,j,k,QU)
             rvrx = rrrx*qxp(i,j,k,QV)
             rwrx = rrrx*qxp(i,j,k,QW)
             ekenrx = HALF*rrrx*(qxp(i,j,k,QU)**2 + qxp(i,j,k,QV)**2 + qxp(i,j,k,QW)**2)
             rerx = qxp(i,j,k,QREINT) + ekenrx
                       
             ! Add transverse predictor
             rrnewrx = rrrx - cdtdy*(fy(i,j+1,k,URHO) - fy(i,j,k,URHO))
             runewrx = rurx - cdtdy*(fy(i,j+1,k,UMX) - fy(i,j,k,UMX))
             rvnewrx = rvrx - cdtdy*(fy(i,j+1,k,UMY) - fy(i,j,k,UMY))
             rwnewrx = rwrx - cdtdy*(fy(i,j+1,k,UMZ) - fy(i,j,k,UMZ))
             renewrx = rerx - cdtdy*(fy(i,j+1,k,UEDEN) - fy(i,j,k,UEDEN))          
             
             ! Reset to original value if adding transverse terms made density negative
             if (transverse_reset_density == 1 .and. rrnewrx .lt. ZERO) then
                rrnewrx = rrrx 
                runewrx = rurx 
                rvnewrx = rvrx 
                rwnewrx = rwrx 
                renewrx = rerx 
             endif

             ! Convert back to primitive form
             qxpo(i,j,k,QRHO) = rrnewrx
             qxpo(i,j,k,QU) = runewrx/qxpo(i,j,k,QRHO)
             qxpo(i,j,k,QV) = rvnewrx/qxpo(i,j,k,QRHO)
             qxpo(i,j,k,QW) = rwnewrx/qxpo(i,j,k,QRHO)
             
             ! note: we run the risk of (rho e) being negative here
             rhoekenrx = HALF*(runewrx**2 + rvnewrx**2 + rwnewrx**2)/qxpo(i,j,k,QRHO)
             qxpo(i,j,k,QREINT) = renewrx - rhoekenrx
             
             if (transverse_reset_rhoe == 1) then
                ! If it is negative, reset the internal energy by
                ! using the discretized expression for updating (rho e).
             
                if (qxpo(i,j,k,QREINT) .le. ZERO) then
                   qxpo(i,j,k,QREINT) = qxp(i,j,k,QREINT) - &
                        cdtdy*(fy(i,j+1,k,UEINT) - fy(i,j,k,UEINT) + pav*du)
                   
                   ! if we are still negative, then we need to reset
                   if (qxpo(i,j,k,QREINT) < ZERO) then
                      eos_state % rho = qxpo(i,j,k,QRHO)
                      eos_state % T = small_temp
                      eos_state % xn(:) = qxpo(i,j,k,QFS:QFS-1+nspec)
                      
                      call eos(eos_input_rt, eos_state)
                      
                      qxpo(i,j,k,QREINT) = qxpo(i,j,k,QRHO) * eos_state % e
                      qxpo(i,j,k,QPRES) = eos_state % p
                   endif
                endif
             endif

             if (ppm_predict_gammae == 0) then

                ! Optionally, use the EOS to calculate the pressure.             

                if (transverse_use_eos .eq. 1) then
                   eos_state % rho = qxpo(i,j,k,QRHO)
                   eos_state % e   = qxpo(i,j,k,QREINT) / qxpo(i,j,k,QRHO)
                   eos_state % T   = small_temp
                   eos_state % xn  = qxpo(i,j,k,QFS:QFS+nspec-1)

                   call eos(eos_input_re, eos_state)

                   pnewrx = eos_state % p
                   qxpo(i,j,k,QREINT) = eos_state % e * eos_state % rho
                else
                   ! add the transverse term to the p evolution eq here
                   pnewrx = qxp(i,j,k,QPRES) - cdtdy*(dup + pav*du*(gamc(i,j,k) - ONE))
                endif

                qxpo(i,j,k,QPRES) = max(pnewrx,small_pres)

             else

                ! Update gammae with its transverse terms
                qxpo(i,j,k,QGAME) = qxp(i,j,k,QGAME) + &
                     cdtdy*( (geav-ONE)*(geav-gamc(i,j,k))*du - uav*dge )

                ! and compute the p edge state from this and (rho e)
                qxpo(i,j,k,QPRES) = qxpo(i,j,k,QREINT)*(qxpo(i,j,k,QGAME)-ONE)
                qxpo(i,j,k,QPRES) = max(qxpo(i,j,k,QPRES), small_pres)
                
             endif

          end do

          !-------------------------------------------------------------------
          ! qxmo state
          !-------------------------------------------------------------------

          do i = lo(1)-1, hi(1)

             pgp = pgdnvy(i,j+1,k)
             pgm = pgdnvy(i,j,k)
             ugp = ugdnvy(i,j+1,k)
             ugm = ugdnvy(i,j,k)
             gegp = gegdnvy(i,j+1,k)
             gegm = gegdnvy(i,j,k)
             
             ! we need to augment our conserved system with either a p
             ! equation or gammae (if we have ppm_predict_gammae = 1) to
             ! be able to deal with the general EOS
             
             dup = pgp*ugp - pgm*ugm
             pav = HALF*(pgp+pgm)
             uav = HALF*(ugp+ugm)
             geav = HALF*(gegp+gegm)
             du = ugp-ugm
             dge = gegp-gegm

             ! Convert to conservation form
             rrlx = qxm(i+1,j,k,QRHO)
             rulx = rrlx*qxm(i+1,j,k,QU)
             rvlx = rrlx*qxm(i+1,j,k,QV)
             rwlx = rrlx*qxm(i+1,j,k,QW)
             ekenlx = HALF*rrlx*(qxm(i+1,j,k,QU)**2 + qxm(i+1,j,k,QV)**2 &
                  + qxm(i+1,j,k,QW)**2)
             relx = qxm(i+1,j,k,QREINT) + ekenlx
             
             ! Add transverse predictor
             rrnewlx = rrlx - cdtdy*(fy(i,j+1,k,URHO) - fy(i,j,k,URHO))
             runewlx = rulx - cdtdy*(fy(i,j+1,k,UMX) - fy(i,j,k,UMX))
             rvnewlx = rvlx - cdtdy*(fy(i,j+1,k,UMY) - fy(i,j,k,UMY))
             rwnewlx = rwlx - cdtdy*(fy(i,j+1,k,UMZ) - fy(i,j,k,UMZ))
             renewlx = relx - cdtdy*(fy(i,j+1,k,UEDEN)- fy(i,j,k,UEDEN))
             
             ! Reset to original value if adding transverse terms made density negative
             if (transverse_reset_density == 1 .and. rrnewlx .lt. ZERO) then
                rrnewlx = rrlx 
                runewlx = rulx 
                rvnewlx = rvlx 
                rwnewlx = rwlx 
                renewlx = relx 
             endif
          
             ! Convert back to primitive form
             qxmo(i+1,j,k,QRHO) = rrnewlx
             qxmo(i+1,j,k,QU) = runewlx/qxmo(i+1,j,k,QRHO)
             qxmo(i+1,j,k,QV) = rvnewlx/qxmo(i+1,j,k,QRHO)
             qxmo(i+1,j,k,QW) = rwnewlx/qxmo(i+1,j,k,QRHO)

             ! note: we run the risk of (rho e) being negative here
             rhoekenlx = HALF*(runewlx**2 + rvnewlx**2 + rwnewlx**2)/qxmo(i+1,j,k,QRHO)
             qxmo(i+1,j,k,QREINT) = renewlx - rhoekenlx

             if (transverse_reset_rhoe == 1) then
                ! If it is negative, reset the internal energy by using the discretized
                ! expression for updating (rho e).
                
                if (qxmo(i+1,j,k,QREINT) .le. ZERO) then
                   qxmo(i+1,j,k,QREINT) = qxm(i+1,j,k,QREINT) - &
                        cdtdy*(fy(i,j+1,k,UEINT) - fy(i,j,k,UEINT) + pav*du)
                   
                   ! if we are still negative, then we need to reset
                   if (qxmo(i+1,j,k,QREINT) < ZERO) then
                      eos_state % rho = qxmo(i+1,j,k,QRHO) 
                      eos_state % T = small_temp
                      eos_state % xn(:) = qxmo(i+1,j,k,QFS:QFS-1+nspec) 
                   
                      call eos(eos_input_rt, eos_state)
                      
                      qxmo(i+1,j,k,QREINT) = qxmo(i+1,j,k,QRHO)*eos_state % e
                      qxmo(i+1,j,k,QPRES) = eos_state % p
                   endif
                endif
             endif

             if (ppm_predict_gammae == 0) then

                ! Optionally, use the EOS to calculate the pressure.             

                if (transverse_use_eos .eq. 1) then
                   eos_state % rho = qxmo(i+1,j,k,QRHO)
                   eos_state % e   = qxmo(i+1,j,k,QREINT) / qxmo(i+1,j,k,QRHO)
                   eos_state % T   = small_temp
                   eos_state % xn  = qxmo(i+1,j,k,QFS:QFS+nspec-1)
                   
                   call eos(eos_input_re, eos_state)

                   pnewlx = eos_state % p
                   qxmo(i+1,j,k,QREINT) = eos_state % e * eos_state % rho
                else
                   ! add the transverse term to the p evolution eq here
                   pnewlx = qxm(i+1,j,k,QPRES) - cdtdy*(dup + pav*du*(gamc(i,j,k) - ONE))
                endif

                qxmo(i+1,j,k,QPRES) = max(pnewlx,small_pres)

             else

                ! Update gammae with its transverse terms
                qxmo(i+1,j,k,QGAME) = qxm(i+1,j,k,QGAME) + &
                     cdtdy*( (geav-ONE)*(geav-gamc(i,j,k))*du - uav*dge )

                ! and compute the p edge state from this and (rho e)
                qxmo(i+1,j,k,QPRES) = qxmo(i+1,j,k,QREINT)*(qxmo(i+1,j,k,QGAME)-ONE)
                qxmo(i+1,j,k,QPRES) = max(qxmo(i+1,j,k,QPRES), small_pres)

             endif

          enddo
       enddo
    enddo

    ! work on qz*

    !-------------------------------------------------------------------------
    ! update all of the passively-advected quantities with the
    ! transerse term and convert back to the primitive quantity
    !-------------------------------------------------------------------------

    do ipassive = 1,npassive
       n  = upass_map(ipassive)
       nq = qpass_map(ipassive)

       do    k = lo(3)  , hi(3)+1
          do j = lo(2)  , hi(2)
          do i = lo(1)-1, hi(1)+1             
             rrnew = qzp(i,j,k,QRHO) - cdtdy*(fy(i,j+1,k,URHO) - fy(i,j,k,URHO))
             compu = qzp(i,j,k,QRHO)*qzp(i,j,k,nq) - cdtdy*(fy(i,j+1,k,n) - fy(i,j,k,n))
             qzpo(i,j,k,nq) = compu/rrnew
          enddo
          enddo
       enddo
       do    k = lo(3)-1, hi(3)
          do j = lo(2)  , hi(2)
          do i = lo(1)-1, hi(1)+1             
             rrnew = qzm(i,j,k+1,QRHO) - cdtdy*(fy(i,j+1,k,URHO) - fy(i,j,k,URHO))
             compu = qzm(i,j,k+1,QRHO)*qzm(i,j,k+1,nq) - cdtdy*(fy(i,j+1,k,n) - fy(i,j,k,n))
             qzmo(i,j,k+1,nq) = compu/rrnew
          enddo
          enddo
       enddo
    enddo

    !-------------------------------------------------------------------
    ! add the transverse flux difference in the y-direction to z-states
    ! for the fluid variables
    !-------------------------------------------------------------------
    
    !-------------------------------------------------------------------
    ! qzpo states
    !-------------------------------------------------------------------          

    do    k = lo(3)  , hi(3)+1
       do j = lo(2)  , hi(2)
          do i = lo(1)-1, hi(1)+1             
          
             pgp = pgdnvy(i,j+1,k)
             pgm = pgdnvy(i,j,k)
             ugp = ugdnvy(i,j+1,k)
             ugm = ugdnvy(i,j,k)
             gegp = gegdnvy(i,j+1,k)
             gegm = gegdnvy(i,j,k)

             ! we need to augment our conserved system with either a p
             ! equation or gammae (if we have ppm_predict_gammae = 1) to
             ! be able to deal with the general EOS
             
             dup = pgp*ugp - pgm*ugm
             pav = HALF*(pgp+pgm)
             uav = HALF*(ugp+ugm)
             geav = HALF*(gegp+gegm)
             du = ugp-ugm
             dge = gegp-gegm
             
             ! Convert to conservation form
             rrrz = qzp(i,j,k,QRHO)
             rurz = rrrz*qzp(i,j,k,QU)
             rvrz = rrrz*qzp(i,j,k,QV)
             rwrz = rrrz*qzp(i,j,k,QW)
             ekenrz = HALF*rrrz*(qzp(i,j,k,QU)**2 + qzp(i,j,k,QV)**2 + qzp(i,j,k,QW)**2)
             rerz = qzp(i,j,k,QREINT) + ekenrz
             
             ! Add transverse predictor
             rrnewrz = rrrz - cdtdy*(fy(i,j+1,k,URHO) - fy(i,j,k,URHO))
             runewrz = rurz - cdtdy*(fy(i,j+1,k,UMX) - fy(i,j,k,UMX))
             rvnewrz = rvrz - cdtdy*(fy(i,j+1,k,UMY) - fy(i,j,k,UMY))
             rwnewrz = rwrz - cdtdy*(fy(i,j+1,k,UMZ) - fy(i,j,k,UMZ))
             renewrz = rerz - cdtdy*(fy(i,j+1,k,UEDEN) - fy(i,j,k,UEDEN))
             
             ! Reset to original value if adding transverse terms made density negative
             if (transverse_reset_density == 1 .and. rrnewrz .lt. ZERO) then
                rrnewrz = rrrz 
                runewrz = rurz 
                rvnewrz = rvrz 
                rwnewrz = rwrz 
                renewrz = rerz 
             endif
             
             ! Convert back to primitive form
             qzpo(i,j,k,QRHO) = rrnewrz
             qzpo(i,j,k,QU) = runewrz/qzpo(i,j,k,QRHO)
             qzpo(i,j,k,QV) = rvnewrz/qzpo(i,j,k,QRHO)
             qzpo(i,j,k,QW) = rwnewrz/qzpo(i,j,k,QRHO)
             
             ! note: we run the risk of (rho e) being negative here
             rhoekenrz = HALF*(runewrz**2 + rvnewrz**2 + rwnewrz**2)/qzpo(i,j,k,QRHO)
             qzpo(i,j,k,QREINT) = renewrz - rhoekenrz
             
             if (transverse_reset_rhoe == 1) then
                ! If it is negative, reset the internal energy by using the discretized
                ! expression for updating (rho e).
                
                if (qzpo(i,j,k,QREINT) .le. ZERO) then
                   qzpo(i,j,k,QREINT) = qzp(i,j,k,QREINT) - &
                        cdtdy*(fy(i,j+1,k,UEINT) - fy(i,j,k,UEINT) + pav*du)
                   
                   ! if we are still negative, then we need to reset
                   if (qzpo(i,j,k,QREINT) < ZERO) then
                      eos_state % rho = qzpo(i,j,k,QRHO)
                      eos_state % T = small_temp
                      eos_state % xn(:) = qzpo(i,j,k,QFS:QFS-1+nspec)
                      
                      call eos(eos_input_rt, eos_state)
                      
                      qzpo(i,j,k,QREINT) = qzpo(i,j,k,QRHO)*eos_state % e
                      qzpo(i,j,k,QPRES) = eos_state % p
                   endif
                endif
             endif

             if (ppm_predict_gammae == 0) then

                ! Optionally, use the EOS to calculate the pressure.
                
                if (transverse_use_eos .eq. 1) then
                   eos_state % rho = qzpo(i,j,k,QRHO)
                   eos_state % e   = qzpo(i,j,k,QREINT) / qzpo(i,j,k,QRHO)
                   eos_state % T   = small_temp
                   eos_state % xn  = qzpo(i,j,k,QFS:QFS+nspec-1)
                   
                   call eos(eos_input_re, eos_state)
                   
                   pnewrz = eos_state % p
                   qzpo(i,j,k,QREINT) = eos_state % e * eos_state % rho
                else
                   ! add the transverse term to the p evolution eq here
                   pnewrz = qzp(i,j,k,QPRES) - cdtdy*(dup + pav*du*(gamc(i,j,k) - ONE))
                endif
                
                qzpo(i,j,k,QPRES) = max(pnewrz,small_pres)
                
             else
                
                ! Update gammae with its transverse terms
                qzpo(i,j,k,QGAME) = qzp(i,j,k,QGAME) + &
                     cdtdy*( (geav-ONE)*(geav-gamc(i,j,k))*du - uav*dge )
                
                ! and compute the p edge state from this and (rho e)
                qzpo(i,j,k,QPRES) = qzpo(i,j,k,QREINT)*(qzpo(i,j,k,QGAME)-ONE)
                qzpo(i,j,k,QPRES) = max(qzpo(i,j,k,QPRES), small_pres)
                
             endif
             
          end do
       end do
    end do

    !-------------------------------------------------------------------
    ! qzmo states
    !-------------------------------------------------------------------
    
    do       k = lo(3)-1, hi(3)
       do    j = lo(2)  , hi(2)
          do i = lo(1)-1, hi(1)+1             

             pgp  =  pgdnvy(i,j+1,k)
             pgm  =  pgdnvy(i,j,k)
             ugp  =  ugdnvy(i,j+1,k)
             ugm  =  ugdnvy(i,j,k)
             gegp = gegdnvy(i,j+1,k)
             gegm = gegdnvy(i,j,k)

             ! we need to augment our conserved system with either a p
             ! equation or gammae (if we have ppm_predict_gammae = 1) to
             ! be able to deal with the general EOS
             
             dup = pgp*ugp - pgm*ugm
             pav = HALF*(pgp+pgm)
             uav = HALF*(ugp+ugm)
             geav = HALF*(gegp+gegm)
             du = ugp-ugm
             dge = gegp-gegm
                       
             ! Convert to conservation form
             rrlz = qzm(i,j,k+1,QRHO)
             rulz = rrlz*qzm(i,j,k+1,QU)
             rvlz = rrlz*qzm(i,j,k+1,QV)
             rwlz = rrlz*qzm(i,j,k+1,QW)
             ekenlz = HALF*rrlz*(qzm(i,j,k+1,QU)**2 + qzm(i,j,k+1,QV)**2 + qzm(i,j,k+1,QW)**2)
             relz = qzm(i,j,k+1,QREINT) + ekenlz
          
             ! Add transverse predictor
             rrnewlz = rrlz - cdtdy*(fy(i,j+1,k,URHO) - fy(i,j,k,URHO))
             runewlz = rulz - cdtdy*(fy(i,j+1,k,UMX) - fy(i,j,k,UMX))
             rvnewlz = rvlz - cdtdy*(fy(i,j+1,k,UMY) - fy(i,j,k,UMY))
             rwnewlz = rwlz - cdtdy*(fy(i,j+1,k,UMZ) - fy(i,j,k,UMZ))
             renewlz = relz - cdtdy*(fy(i,j+1,k,UEDEN)- fy(i,j,k,UEDEN))
             
             ! Reset to original value if adding transverse terms made density negative
             if (transverse_reset_density == 1 .and. rrnewlz .lt. ZERO) then
                rrnewlz = rrlz
                runewlz = rulz
                rvnewlz = rvlz
                rwnewlz = rwlz
                renewlz = relz
             endif
             
             ! Convert back to primitive form
             qzmo(i,j,k+1,QRHO) = rrnewlz
             qzmo(i,j,k+1,QU) = runewlz/qzmo(i,j,k+1,QRHO)
             qzmo(i,j,k+1,QV) = rvnewlz/qzmo(i,j,k+1,QRHO)
             qzmo(i,j,k+1,QW) = rwnewlz/qzmo(i,j,k+1,QRHO)
             
             ! note: we run the risk of (rho e) being negative here
             rhoekenlz = HALF*(runewlz**2 + rvnewlz**2 + rwnewlz**2)/qzmo(i,j,k+1,QRHO)
             qzmo(i,j,k+1,QREINT) = renewlz - rhoekenlz
             
             if (transverse_reset_rhoe == 1) then
                ! If it is negative, reset the internal energy by using the discretized
                ! expression for updating (rho e).
                
                if (qzmo(i,j,k+1,QREINT) .le. ZERO) then
                   qzmo(i,j,k+1,QREINT) = qzm(i,j,k+1,QREINT) - &
                        cdtdy*(fy(i,j+1,k,UEINT) - fy(i,j,k,UEINT) + pav*du)
                   
                   ! if we are still negative, then we need to reset
                   if (qzmo(i,j,k+1,QREINT) < ZERO) then
                      eos_state % rho = qzmo(i,j,k+1,QRHO)
                      eos_state % T = small_temp
                      eos_state % xn(:) = qzmo(i,j,k+1,QFS:QFS-1+nspec)
                      
                      call eos(eos_input_rt, eos_state)
                      
                      qzmo(i,j,k+1,QREINT) = qzmo(i,j,k+1,QRHO)*eos_state % e
                      qzmo(i,j,k+1,QPRES) = eos_state % p
                   endif
                endif
             endif
             
             if (ppm_predict_gammae == 0) then
                
                ! Optionally, use the EOS to calculate the pressure.
                
                if (transverse_use_eos .eq. 1) then
                   eos_state % rho = qzmo(i,j,k+1,QRHO)
                   eos_state % e   = qzmo(i,j,k+1,QREINT) / qzmo(i,j,k+1,QRHO)
                   eos_state % T   = small_temp
                   eos_state % xn  = qzmo(i,j,k+1,QFS:QFS+nspec-1)
                   
                   call eos(eos_input_re, eos_state)
                   
                   pnewlz = eos_state % p
                   qzmo(i,j,k+1,QREINT) = eos_state % e * eos_state % rho
                else
                   ! add the transverse term to the p evolution eq here
                   pnewlz = qzm(i,j,k+1,QPRES) - cdtdy*(dup + pav*du*(gamc(i,j,k) - ONE))
                endif
                
                qzmo(i,j,k+1,QPRES) = max(pnewlz,small_pres)
                
             else
                
                ! Update gammae with its transverse terms
                qzmo(i,j,k+1,QGAME) = qzm(i,j,k+1,QGAME) + &
                     cdtdy*( (geav-ONE)*(geav-gamc(i,j,k))*du - uav*dge )
                
                ! and compute the p edge state from this and (rho e)
                qzmo(i,j,k+1,QPRES) = qzmo(i,j,k+1,QREINT)*(qzmo(i,j,k+1,QGAME)-ONE)
                qzmo(i,j,k+1,QPRES) = max(qzmo(i,j,k+1,QPRES), small_pres)

             endif
             
          enddo
       enddo
    enddo

  end subroutine transy


  !===========================================================================
  ! transz
  !===========================================================================
  subroutine transz(qxm,qxmo,qxp,qxpo, &
                    qym,qymo,qyp,qypo, qlo, qhi, &
                    fz, flo, fhi, &
                    ugdnvz,pgdnvz,gegdnvz, gdlo, gdhi, &
                    gamc, gclo, gchi, &
                    cdtdz, lo, hi)
    use network, only : nspec, naux
    use meth_params_module, only : QVAR, NVAR, QRHO, QU, QV, QW, &
                                   QPRES, QREINT, QGAME, QESGS, QFA, QFS, &
                                   URHO, UMX, UMY, UMZ, UEDEN, UEINT, UESGS, UFA, UFS, &
                                   nadv, small_pres, small_temp, &
                                   npassive, upass_map, qpass_map, &
                                   ppm_predict_gammae, &
                                   transverse_use_eos, transverse_reset_density, transverse_reset_rhoe
    use eos_module
    implicit none

    integer, intent(in) :: qlo(3), qhi(3), flo(3), fhi(3), gdlo(3), gdhi(3), &
         gclo(3), gchi(3), lo(3), hi(3)
    double precision, intent(in) :: cdtdz
    double precision,intent(in )::    qxm( qlo(1): qhi(1), qlo(2): qhi(2), qlo(3): qhi(3),QVAR)
    double precision,intent(in )::    qxp( qlo(1): qhi(1), qlo(2): qhi(2), qlo(3): qhi(3),QVAR)
    double precision,intent(in )::    qym( qlo(1): qhi(1), qlo(2): qhi(2), qlo(3): qhi(3),QVAR)
    double precision,intent(in )::    qyp( qlo(1): qhi(1), qlo(2): qhi(2), qlo(3): qhi(3),QVAR)
    double precision,intent(in )::     fz( flo(1): fhi(1), flo(2): fhi(2), flo(3): fhi(3),NVAR)
    double precision,intent(in ):: ugdnvz(gdlo(1):gdhi(1),gdlo(2):gdhi(2),gdlo(3):gdhi(3))
    double precision,intent(in ):: pgdnvz(gdlo(1):gdhi(1),gdlo(2):gdhi(2),gdlo(3):gdhi(3))
    double precision,intent(in )::gegdnvz(gdlo(1):gdhi(1),gdlo(2):gdhi(2),gdlo(3):gdhi(3))
    double precision,intent(in )::   gamc(gclo(1):gchi(1),gclo(2):gchi(2),gclo(3):gchi(3))
    double precision,intent(out)::   qxmo( qlo(1): qhi(1), qlo(2): qhi(2), qlo(3): qhi(3),QVAR)
    double precision,intent(out)::   qxpo( qlo(1): qhi(1), qlo(2): qhi(2), qlo(3): qhi(3),QVAR)
    double precision,intent(out)::   qymo( qlo(1): qhi(1), qlo(2): qhi(2), qlo(3): qhi(3),QVAR)
    double precision,intent(out)::   qypo( qlo(1): qhi(1), qlo(2): qhi(2), qlo(3): qhi(3),QVAR)

    integer :: i,j,k,n,nq,ipassive
    
    double precision rrnew, compu
    double precision rrrx, rrry, rrlx, rrly
    double precision rurx, rury, rulx, ruly
    double precision rvrx, rvry, rvlx, rvly
    double precision rwrx, rwry, rwlx, rwly
    double precision ekenrx, ekenry, ekenlx, ekenly
    double precision rerx, rery, relx, rely
    double precision rrnewrx, rrnewry, rrnewlx, rrnewly
    double precision runewrx, runewry, runewlx, runewly
    double precision rvnewrx, rvnewry, rvnewlx, rvnewly
    double precision rwnewrx, rwnewry, rwnewlx, rwnewly
    double precision renewrx, renewry, renewlx, renewly
    double precision pnewrx, pnewry, pnewlx, pnewly
    double precision rhoekenrx, rhoekenry, rhoekenlx, rhoekenly
    double precision pgp, pgm, ugp, ugm, gegp, gegm, dup, pav, du, dge, uav, geav

    type (eos_t) :: eos_state

    ! work on qx* first

    !-------------------------------------------------------------------------    
    ! update all of the passively-advected quantities with the
    ! transerse term and convert back to the primitive quantity
    !-------------------------------------------------------------------------

    do ipassive = 1,npassive
       n  = upass_map(ipassive)
       nq = qpass_map(ipassive)

       do    k = lo(3)  , hi(3)
       do    j = lo(2)-1, hi(2)+1
          do i = lo(1)  , hi(1)+1
             rrnew = qxp(i,j,k,QRHO) - cdtdz*(fz(i,j,k+1,URHO) - fz(i,j,k,URHO))
             compu = qxp(i,j,k,QRHO)*qxp(i,j,k,nq) - cdtdz*(fz(i,j,k+1,n) - fz(i,j,k,n))
             qxpo(i,j,k,nq) = compu/rrnew
          end do
          do i = lo(1)-1, hi(1)
             rrnew = qxm(i+1,j,k,QRHO) - cdtdz*(fz(i,j,k+1,URHO) - fz(i,j,k,URHO))
             compu = qxm(i+1,j,k,QRHO)*qxm(i+1,j,k,nq) - cdtdz*(fz(i,j,k+1,n) - fz(i,j,k,n))
             qxmo(i+1,j,k,nq) = compu/rrnew
          enddo
       enddo
       enddo
    enddo

    !-------------------------------------------------------------------
    ! add transverse flux difference in the z-direction to x-states
    ! for the fluid variables
    !-------------------------------------------------------------------    

    do    k = lo(3)  , hi(3)
       do j = lo(2)-1, hi(2)+1
       do i = lo(1)-1, hi(1)+1

          pgp  =  pgdnvz(i,j,k+1)
          pgm  =  pgdnvz(i,j,k)
          ugp  =  ugdnvz(i,j,k+1)
          ugm  =  ugdnvz(i,j,k)
          gegp = gegdnvz(i,j,k+1)
          gegm = gegdnvz(i,j,k)
          
          ! Convert to conservation form
          rrrx = qxp(i,j,k,QRHO)
          rurx = rrrx*qxp(i,j,k,QU)
          rvrx = rrrx*qxp(i,j,k,QV)
          rwrx = rrrx*qxp(i,j,k,QW)
          ekenrx = HALF*rrrx*(qxp(i,j,k,QU)**2 + qxp(i,j,k,QV)**2 &
               + qxp(i,j,k,QW)**2)
          rerx = qxp(i,j,k,QREINT) + ekenrx
          
          rrlx = qxm(i+1,j,k,QRHO)
          rulx = rrlx*qxm(i+1,j,k,QU)
          rvlx = rrlx*qxm(i+1,j,k,QV)
          rwlx = rrlx*qxm(i+1,j,k,QW)
          ekenlx = HALF*rrlx*(qxm(i+1,j,k,QU)**2 + qxm(i+1,j,k,QV)**2 &
               + qxm(i+1,j,k,QW)**2)
          relx = qxm(i+1,j,k,QREINT) + ekenlx
          
          ! Add transverse predictor
          rrnewrx = rrrx - cdtdz*(fz(i,j,k+1,URHO) - fz(i,j,k,URHO))
          runewrx = rurx - cdtdz*(fz(i,j,k+1,UMX) - fz(i,j,k,UMX))
          rvnewrx = rvrx - cdtdz*(fz(i,j,k+1,UMY) - fz(i,j,k,UMY))
          rwnewrx = rwrx - cdtdz*(fz(i,j,k+1,UMZ) - fz(i,j,k,UMZ))
          renewrx = rerx - cdtdz*(fz(i,j,k+1,UEDEN) - fz(i,j,k,UEDEN))
          
          rrnewlx = rrlx - cdtdz*(fz(i,j,k+1,URHO) - fz(i,j,k,URHO))
          runewlx = rulx - cdtdz*(fz(i,j,k+1,UMX) - fz(i,j,k,UMX))
          rvnewlx = rvlx - cdtdz*(fz(i,j,k+1,UMY) - fz(i,j,k,UMY))
          rwnewlx = rwlx - cdtdz*(fz(i,j,k+1,UMZ) - fz(i,j,k,UMZ))
          renewlx = relx - cdtdz*(fz(i,j,k+1,UEDEN) - fz(i,j,k,UEDEN))

          ! Reset to original value if adding transverse terms made density negative
          if (transverse_reset_density == 1) then
             if (rrnewrx .lt. ZERO) then
                rrnewrx = rrrx 
                runewrx = rurx 
                rvnewrx = rvrx 
                rwnewrx = rwrx 
                renewrx = rerx 
             endif
             if (rrnewlx .lt. ZERO) then
                rrnewlx = rrlx 
                runewlx = rulx 
                rvnewlx = rvlx 
                rwnewlx = rwlx 
                renewlx = relx 
             endif
          endif

          ! we need to augment our conserved system with either a p
          ! equation or gammae (if we have ppm_predict_gammae = 1) to
          ! be able to deal with the general EOS

          dup = pgp*ugp - pgm*ugm
          pav = HALF*(pgp+pgm)
          uav = HALF*(ugp+ugm)
          geav = HALF*(gegp+gegm)
          du = ugp-ugm
          dge = gegp-gegm

          !-------------------------------------------------------------------
          ! qxpo state
          !-------------------------------------------------------------------
          
          ! Convert back to primitive form
          if (i.ge.lo(1)) then
             qxpo(i,j,k,QRHO) = rrnewrx
             qxpo(i,j,k,QU) = runewrx/qxpo(i,j,k,QRHO)
             qxpo(i,j,k,QV) = rvnewrx/qxpo(i,j,k,QRHO)
             qxpo(i,j,k,QW) = rwnewrx/qxpo(i,j,k,QRHO)

             ! note: we run the risk of (rho e) being negative here
             rhoekenrx = HALF*(runewrx**2 + rvnewrx**2 + rwnewrx**2)/qxpo(i,j,k,QRHO)
             qxpo(i,j,k,QREINT) = renewrx - rhoekenrx

             if (transverse_reset_rhoe == 1) then
                ! If it is negative, reset the internal energy by using the discretized
                ! expression for updating (rho e).
                
                if (qxpo(i,j,k,QREINT) .le. ZERO) then
                   qxpo(i,j,k,QREINT) = qxp(i,j,k,QREINT) - &
                        cdtdz*(fz(i,j,k+1,UEINT) - fz(i,j,k,UEINT) + pav*du)

                   ! if we are still negative, then we need to reset
                   if (qxpo(i,j,k,QREINT) < ZERO) then
                      eos_state % rho = qxpo(i,j,k,QRHO) 
                      eos_state % T = small_temp
                      eos_state % xn(:) = qxpo(i,j,k,QFS:QFS-1+nspec) 
                      
                      call eos(eos_input_rt, eos_state)
                      
                      qxpo(i,j,k,QREINT) = qxpo(i,j,k,QRHO)*eos_state % e
                      qxpo(i,j,k,QPRES) = eos_state % p
                   endif
                endif
             endif

             if (ppm_predict_gammae == 0) then
             
                ! Optionally, use the EOS to calculate the pressure.

                if (transverse_use_eos .eq. 1) then
                   eos_state % rho = qxpo(i,j,k,QRHO)
                   eos_state % e   = qxpo(i,j,k,QREINT) / qxpo(i,j,k,QRHO)
                   eos_state % T   = small_temp
                   eos_state % xn  = qxpo(i,j,k,QFS:QFS+nspec-1)

                   call eos(eos_input_re, eos_state)

                   pnewrx = eos_state % p
                   qxpo(i,j,k,QREINT) = eos_state % e * eos_state % rho
                else
                   ! add the transverse term to the p evolution eq here
                   pnewrx = qxp(i,j,k,QPRES) - cdtdz*(dup + pav*du*(gamc(i,j,k) - ONE))

                endif

                qxpo(i,j,k,QPRES) = max(pnewrx,small_pres)

             else

                ! Update gammae with its transverse terms
                qxpo(i,j,k,QGAME) = qxp(i,j,k,QGAME) + &
                     cdtdz*( (geav-ONE)*(geav-gamc(i,j,k))*du - uav*dge )

                ! and compute the p edge state from this and (rho e)
                qxpo(i,j,k,QPRES) = qxpo(i,j,k,QREINT)*(qxpo(i,j,k,QGAME)-ONE)
                qxpo(i,j,k,QPRES) = max(qxpo(i,j,k,QPRES), small_pres)

             endif

          end if

          !-------------------------------------------------------------------
          ! qxmo state
          !-------------------------------------------------------------------
          
          if (i.le.hi(1)) then
             qxmo(i+1,j,k,QRHO) = rrnewlx
             qxmo(i+1,j,k,QU) = runewlx/qxmo(i+1,j,k,QRHO)
             qxmo(i+1,j,k,QV) = rvnewlx/qxmo(i+1,j,k,QRHO)
             qxmo(i+1,j,k,QW) = rwnewlx/qxmo(i+1,j,k,QRHO)

             ! note: we run the risk of (rho e) being negative here
             rhoekenlx = HALF*(runewlx**2 + rvnewlx**2 + rwnewlx**2)/qxmo(i+1,j,k,QRHO)
             qxmo(i+1,j,k,QREINT) = renewlx - rhoekenlx

             if (transverse_reset_rhoe == 1) then
                ! If it is negative, reset the internal energy by using the discretized
                ! expression for updating (rho e).
                
                if (qxmo(i+1,j,k,QREINT) .le. ZERO) then
                   qxmo(i+1,j,k,QREINT) = qxm(i+1,j,k,QREINT) - &
                        cdtdz*(fz(i,j,k+1,UEINT) - fz(i,j,k,UEINT) + pav*du)
                   
                   ! if we are still negative, then we need to reset
                   if (qxmo(i+1,j,k,QREINT) < ZERO) then
                      eos_state % rho = qxmo(i+1,j,k,QRHO)
                      eos_state % T = small_temp
                      eos_state % xn(:) = qxmo(i+1,j,k,QFS:QFS-1+nspec)
                      
                      call eos(eos_input_rt, eos_state)
                      
                      qxmo(i+1,j,k,QREINT) = qxmo(i+1,j,k,QRHO)*eos_state % e 
                      qxmo(i+1,j,k,QPRES) = eos_state % p
                   endif
                endif
             endif

             if (ppm_predict_gammae == 0) then
             
                ! Optionally, use the EOS to calculate the pressure.

                if (transverse_use_eos .eq. 1) then
                   eos_state % rho = qxmo(i+1,j,k,QRHO)
                   eos_state % e   = qxmo(i+1,j,k,QREINT) / qxmo(i+1,j,k,QRHO)
                   eos_state % T   = small_temp
                   eos_state % xn  = qxmo(i+1,j,k,QFS:QFS+nspec-1)
                   
                   call eos(eos_input_re, eos_state)

                   pnewlx = eos_state % p
                   qxmo(i+1,j,k,QREINT) = eos_state % e * eos_state % rho
                else
                   ! add the transverse term to the p evolution eq here
                   pnewlx = qxm(i+1,j,k,QPRES) - cdtdz*(dup + pav*du*(gamc(i,j,k) - ONE))
                endif
                
                qxmo(i+1,j,k,QPRES) = max(pnewlx,small_pres)

             else
                
                ! Update gammae with its transverse terms             
                qxmo(i+1,j,k,QGAME) = qxm(i+1,j,k,QGAME) + &
                     cdtdz*( (geav-ONE)*(geav-gamc(i,j,k))*du - uav*dge )
                
                ! and compute the p edge state from this and (rho e)
                qxmo(i+1,j,k,QPRES) = qxmo(i+1,j,k,QREINT)*(qxmo(i+1,j,k,QGAME)-ONE)
                qxmo(i+1,j,k,QPRES) = max(qxmo(i+1,j,k,QPRES), small_pres)
                
             end if

          endif

       enddo
       enddo
    enddo

    ! work on qy*

    !-------------------------------------------------------------------------    
    ! update all of the passively-advected quantities with the
    ! transerse term and convert back to the primitive quantity
    !-------------------------------------------------------------------------

    do ipassive = 1,npassive
       n  = upass_map(ipassive)
       nq = qpass_map(ipassive)

       do    k = lo(3)  , hi(3)
          do j = lo(2)  , hi(2)+1
          do i = lo(1)-1, hi(1)+1
             rrnew = qyp(i,j,k,QRHO) - cdtdz*(fz(i,j,k+1,URHO) - fz(i,j,k,URHO))
             compu = qyp(i,j,k,QRHO)*qyp(i,j,k,nq) - cdtdz*(fz(i,j,k+1,n) - fz(i,j,k,n))
             qypo(i,j,k,nq) = compu/rrnew
          end do
          end do
          do j = lo(2)-1, hi(2)
          do i = lo(1)-1, hi(1)+1
             rrnew = qym(i,j+1,k,QRHO) - cdtdz*(fz(i,j,k+1,URHO) - fz(i,j,k,URHO))
             compu = qym(i,j+1,k,QRHO)*qym(i,j+1,k,nq) - cdtdz*(fz(i,j,k+1,n) - fz(i,j,k,n))
             qymo(i,j+1,k,nq) = compu/rrnew
          end do
          end do
       end do
    end do

    !-------------------------------------------------------------------
    ! add transverse flux difference in the z-direction to y-states
    ! for the fluid variables
    !-------------------------------------------------------------------    

    do    k = lo(3)  , hi(3)
       do j = lo(2)-1, hi(2)+1 
       do i = lo(1)-1, hi(1)+1             

          pgp  =  pgdnvz(i,j,k+1)
          pgm  =  pgdnvz(i,j,k)
          ugp  =  ugdnvz(i,j,k+1)
          ugm  =  ugdnvz(i,j,k)
          gegp = gegdnvz(i,j,k+1)
          gegm = gegdnvz(i,j,k)

          ! Convert to conservation form
          rrry = qyp(i,j,k,QRHO)
          rury = rrry*qyp(i,j,k,QU)
          rvry = rrry*qyp(i,j,k,QV)
          rwry = rrry*qyp(i,j,k,QW)
          ekenry = HALF*rrry*(qyp(i,j,k,QU)**2 + qyp(i,j,k,QV)**2 &
               + qyp(i,j,k,QW)**2)
          rery = qyp(i,j,k,QREINT) + ekenry
          
          rrly = qym(i,j+1,k,QRHO)
          ruly = rrly*qym(i,j+1,k,QU)
          rvly = rrly*qym(i,j+1,k,QV)
          rwly = rrly*qym(i,j+1,k,QW)
          ekenly = HALF*rrly*(qym(i,j+1,k,QU)**2 + qym(i,j+1,k,QV)**2 &
               + qym(i,j+1,k,QW)**2)
          rely = qym(i,j+1,k,QREINT) + ekenly
          
          ! Add transverse predictor
          rrnewry = rrry - cdtdz*(fz(i,j,k+1,URHO) - fz(i,j,k,URHO))
          runewry = rury - cdtdz*(fz(i,j,k+1,UMX) - fz(i,j,k,UMX))
          rvnewry = rvry - cdtdz*(fz(i,j,k+1,UMY) - fz(i,j,k,UMY))
          rwnewry = rwry - cdtdz*(fz(i,j,k+1,UMZ) - fz(i,j,k,UMZ))
          renewry = rery - cdtdz*(fz(i,j,k+1,UEDEN) - fz(i,j,k,UEDEN))
                    
          rrnewly = rrly - cdtdz*(fz(i,j,k+1,URHO) - fz(i,j,k,URHO))
          runewly = ruly - cdtdz*(fz(i,j,k+1,UMX) - fz(i,j,k,UMX))
          rvnewly = rvly - cdtdz*(fz(i,j,k+1,UMY) - fz(i,j,k,UMY))
          rwnewly = rwly - cdtdz*(fz(i,j,k+1,UMZ) - fz(i,j,k,UMZ))
          renewly = rely - cdtdz*(fz(i,j,k+1,UEDEN) - fz(i,j,k,UEDEN))

          ! Reset to original value if adding transverse terms made density negative
          if (transverse_reset_density == 1) then
             if (rrnewry .lt. ZERO) then
                rrnewry = rrry 
                runewry = rury 
                rvnewry = rvry 
                rwnewry = rwry 
                renewry = rery 
             endif
             if (rrnewly .lt. ZERO) then
                rrnewly = rrly 
                runewly = ruly 
                rvnewly = rvly 
                rwnewly = rwly 
                renewly = rely 
             endif
          endif

          ! we need to augment our conserved system with either a p
          ! equation or gammae (if we have ppm_predict_gammae = 1) to
          ! be able to deal with the general EOS

          dup = pgp*ugp - pgm*ugm
          pav = HALF*(pgp+pgm)
          uav = HALF*(ugp+ugm)
          geav = HALF*(gegp+gegm)
          du = ugp-ugm
          dge = gegp-gegm

          !-------------------------------------------------------------------
          ! qypo state
          !-------------------------------------------------------------------
          
          if (j.ge.lo(2)) then
             qypo(i,j,k,QRHO) = rrnewry
             qypo(i,j,k,QU) = runewry/qypo(i,j,k,QRHO)
             qypo(i,j,k,QV) = rvnewry/qypo(i,j,k,QRHO)
             qypo(i,j,k,QW) = rwnewry/qypo(i,j,k,QRHO)

             ! note: we run the risk of (rho e) being negative here
             rhoekenry = HALF*(runewry**2 + rvnewry**2 + rwnewry**2)/qypo(i,j,k,QRHO)
             qypo(i,j,k,QREINT) = renewry - rhoekenry

             if (transverse_reset_rhoe == 1) then
                ! If it is negative, reset the internal energy by using the discretized
                ! expression for updating (rho e).

                if (qypo(i,j,k,QREINT) .le. ZERO) then
                   qypo(i,j,k,QREINT) = qyp(i,j,k,QREINT) - &
                        cdtdz*(fz(i,j,k+1,UEINT) - fz(i,j,k,UEINT) + pav*du)
                   
                   ! if we are still negative, then we need to reset
                   if (qypo(i,j,k,QREINT) < ZERO) then
                      eos_state % rho = qypo(i,j,k,QRHO)
                      eos_state % T = small_temp
                      eos_state % xn(:) = qypo(i,j,k,QFS:QFS-1+nspec)
                      
                      call eos(eos_input_rt, eos_state)
                      
                      qypo(i,j,k,QREINT) = qypo(i,j,k,QRHO)*eos_state % e
                      qypo(i,j,k,QPRES) = eos_state % p
                   endif
                endif
             endif

             if (ppm_predict_gammae == 0) then

                ! Optionally, use the EOS to calculate the pressure.

                if (transverse_use_eos .eq. 1) then
                   eos_state % rho = qypo(i,j,k,QRHO)
                   eos_state % e   = qypo(i,j,k,QREINT) / qypo(i,j,k,QRHO)
                   eos_state % T   = small_temp
                   eos_state % xn  = qypo(i,j,k,QFS:QFS+nspec-1)

                   call eos(eos_input_re, eos_state)
                   
                   pnewry = eos_state % p
                   qypo(i,j,k,QREINT) = eos_state % e * eos_state % rho
                else
                   ! add the transverse term to the p evolution eq here
                   pnewry = qyp(i,j,k,QPRES) - cdtdz*(dup + pav*du*(gamc(i,j,k) - ONE))
                endif

                qypo(i,j,k,QPRES) = max(pnewry,small_pres)

             else

                ! Update gammae with its transverse terms
                qypo(i,j,k,QGAME) = qyp(i,j,k,QGAME) + &
                     cdtdz*( (geav-ONE)*(geav-gamc(i,j,k))*du - uav*dge )

                ! and compute the p edge state from this and (rho e)
                qypo(i,j,k,QPRES) = qypo(i,j,k,QREINT)*(qypo(i,j,k,QGAME)-ONE)
                qypo(i,j,k,QPRES) = max(qypo(i,j,k,QPRES), small_pres)

             endif

          end if

          !-------------------------------------------------------------------
          ! qymo state
          !-------------------------------------------------------------------
          
          if (j.le.hi(2)) then
             qymo(i,j+1,k,QRHO) = rrnewly
             qymo(i,j+1,k,QU) = runewly/qymo(i,j+1,k,QRHO)
             qymo(i,j+1,k,QV) = rvnewly/qymo(i,j+1,k,QRHO)
             qymo(i,j+1,k,QW) = rwnewly/qymo(i,j+1,k,QRHO)

             ! note: we run the risk of (rho e) being negative here
             rhoekenly = HALF*(runewly**2 + rvnewly**2 + rwnewly**2)/qymo(i,j+1,k,QRHO)
             qymo(i,j+1,k,QREINT) = renewly - rhoekenly

             if (transverse_reset_rhoe == 1) then
                ! If it is negative, reset the internal energy by using the discretized
                ! expression for updating (rho e).

                if (qymo(i,j+1,k,QREINT) .le. ZERO) then
                   qymo(i,j+1,k,QREINT) = qym(i,j+1,k,QREINT) - &
                        cdtdz*(fz(i,j,k+1,UEINT) - fz(i,j,k,UEINT) + pav*du)
                   
                   ! if we are still negative, then we need to reset
                   if (qymo(i,j+1,k,QREINT) < ZERO) then
                      eos_state % rho = qymo(i,j+1,k,QRHO)
                      eos_state % T = small_temp
                      eos_state % xn(:) = qymo(i,j+1,k,QFS:QFS-1+nspec)
                      
                      call eos(eos_input_rt, eos_state)
                      
                      qymo(i,j+1,k,QREINT) =  qymo(i,j+1,k,QRHO)*eos_state % e
                      qymo(i,j+1,k,QPRES) =  eos_state % p
                   endif
                endif
             endif

             if (ppm_predict_gammae == 0) then

                ! Optionally, use the EOS to calculate the pressure.             

                if (transverse_use_eos .eq. 1) then
                   eos_state % rho = qymo(i,j+1,k,QRHO)
                   eos_state % e   = qymo(i,j+1,k,QREINT) / qymo(i,j+1,k,QRHO)
                   eos_state % T   = small_temp
                   eos_state % xn  = qymo(i,j+1,k,QFS:QFS+nspec-1)
                   
                   call eos(eos_input_re, eos_state)
                
                   pnewly = eos_state % p
                   qymo(i,j+1,k,QREINT) = eos_state % e * eos_state % rho
                else
                   ! add the transverse term to the p evolution eq here
                   pnewly = qym(i,j+1,k,QPRES) - cdtdz*(dup + pav*du*(gamc(i,j,k) - ONE))
                endif

                qymo(i,j+1,k,QPRES) = max(pnewly,small_pres)

             else

                ! Update gammae with its transverse terms
                qymo(i,j+1,k,QGAME) = qym(i,j+1,k,QGAME) + &
                     cdtdz*( (geav-ONE)*(geav-gamc(i,j,k))*du - uav*dge )

                ! and compute the p edge state from this and (rho e)
                qymo(i,j+1,k,QPRES) = qymo(i,j+1,k,QREINT)*(qymo(i,j+1,k,QGAME)-ONE)
                qymo(i,j+1,k,QPRES) = max(qymo(i,j+1,k,QPRES), small_pres)
                
             endif

          endif

       enddo
       enddo
    enddo

  end subroutine transz


  !===========================================================================
  ! transxy
  !===========================================================================
  subroutine transxy(qm,qmo,qp,qpo, qlo, qhi, &
                     fxy,fyx, flo, fhi, &
                     ugdnvx,pgdnvx,gegdnvx,ugdnvy,pgdnvy,gegdnvy, gdlo, gdhi, &
                     gamc,gclo,gchi, &
                     srcQ,slo,shi, &
                     grav,gvlo,gvhi, &
                     rot,rlo,rhi, &
                     hdt,hdtdx,hdtdy,lo,hi)
    
    use network, only : nspec, naux
    use meth_params_module, only : QVAR, NVAR, QRHO, QU, QV, QW, &
                                   QPRES, QREINT, QGAME, QESGS, QFA, QFS, &
                                   URHO, UMX, UMY, UMZ, UEDEN, UEINT, UESGS, UFA, UFS, &
                                   nadv, small_pres, small_temp, &
                                   npassive, upass_map, qpass_map, &
                                   ppm_predict_gammae, &
                                   transverse_use_eos, transverse_reset_density, transverse_reset_rhoe, &
                                   ppm_type, ppm_trace_grav, rot_period, ppm_trace_rot, &
                                   do_grav, do_rotation
    use eos_module

    implicit none

    integer, intent(in) :: qlo(3),qhi(3),flo(3),fhi(3),gdlo(3),gdhi(3),gclo(3),gchi(3), &
         slo(3),shi(3),gvlo(3),gvhi(3),rlo(3),rhi(3),lo(3),hi(3)
    double precision, intent(in) :: hdt, hdtdx, hdtdy
    double precision,intent(in )::     qm( qlo(1): qhi(1), qlo(2): qhi(2), qlo(3): qhi(3),QVAR)
    double precision,intent(out)::    qmo( qlo(1): qhi(1), qlo(2): qhi(2), qlo(3): qhi(3),QVAR)
    double precision,intent(in )::     qp( qlo(1): qhi(1), qlo(2): qhi(2), qlo(3): qhi(3),QVAR)
    double precision,intent(out)::    qpo( qlo(1): qhi(1), qlo(2): qhi(2), qlo(3): qhi(3),QVAR)
    double precision,intent(in )::    fxy( flo(1): fhi(1), flo(2): fhi(2), flo(3): fhi(3),NVAR)
    double precision,intent(in )::    fyx( flo(1): fhi(1), flo(2): fhi(2), flo(3): fhi(3),NVAR)
    double precision,intent(in ):: ugdnvx(gdlo(1):gdhi(1),gdlo(2):gdhi(2),gdlo(3):gdhi(3))
    double precision,intent(in ):: pgdnvx(gdlo(1):gdhi(1),gdlo(2):gdhi(2),gdlo(3):gdhi(3))
    double precision,intent(in )::gegdnvx(gdlo(1):gdhi(1),gdlo(2):gdhi(2),gdlo(3):gdhi(3))
    double precision,intent(in ):: ugdnvy(gdlo(1):gdhi(1),gdlo(2):gdhi(2),gdlo(3):gdhi(3))
    double precision,intent(in ):: pgdnvy(gdlo(1):gdhi(1),gdlo(2):gdhi(2),gdlo(3):gdhi(3))
    double precision,intent(in )::gegdnvy(gdlo(1):gdhi(1),gdlo(2):gdhi(2),gdlo(3):gdhi(3))
    double precision,intent(in )::   gamc(gclo(1):gchi(1),gclo(2):gchi(2),gclo(3):gchi(3))
    double precision,intent(in )::   srcQ( slo(1): shi(1), slo(2): shi(2), slo(3): shi(3),QVAR)
    double precision,intent(in )::   grav(gvlo(1):gvhi(1),gvlo(2):gvhi(2),gvlo(3):gvhi(3),3)
    double precision,intent(in )::    rot( rlo(1): rhi(1), rlo(2): rhi(2), rlo(3): rhi(3),3)

    integer i, j, k, n, nq, ipassive    
    
    double precision rrr, rur, rvr, rwr, rer, ekenr, rhoekenr
    double precision rrl, rul, rvl, rwl, rel, ekenl, rhoekenl
    double precision rrnewr, runewr, rvnewr, rwnewr, renewr
    double precision rrnewl, runewl, rvnewl, rwnewl, renewl
    double precision pnewr, pnewl
    double precision pgxp, pgxm, ugxp, ugxm, gegxp, gegxm, duxp, pxav, dux, pxnew, gexnew
    double precision pgyp, pgym, ugyp, ugym, gegyp, gegym, duyp, pyav, duy, pynew, geynew
    double precision uxav, gexav, dgex, uyav, geyav, dgey
    double precision pgxpm, pgxmm, ugxpm, ugxmm, gegxpm, gegxmm, duxpm, pxavm, duxm, pxnewm, gexnewm
    double precision pgypm, pgymm, ugypm, ugymm, gegypm, gegymm, duypm, pyavm, duym, pynewm, geynewm
    double precision uxavm, gexavm, dgexm, uyavm, geyavm, dgeym
    double precision compr, compl, compnr, compnl
    
    type (eos_t) :: eos_state

    !-------------------------------------------------------------------------    
    ! update all of the passively-advected quantities with the
    ! transerse term and convert back to the primitive quantity
    !-------------------------------------------------------------------------

    do ipassive = 1,npassive
       n  = upass_map(ipassive)
       nq = qpass_map(ipassive)

       do    k = lo(3), hi(3)+1
          do j = lo(2), hi(2) 
          do i = lo(1), hi(1)
             rrr = qp(i,j,k,QRHO)
             compr = rrr*qp(i,j,k,nq)
             rrnewr = rrr - hdtdx*(fxy(i+1,j,k,URHO) - fxy(i,j,k,URHO)) &
                          - hdtdy*(fyx(i,j+1,k,URHO) - fyx(i,j,k,URHO))
             compnr = compr - hdtdx*(fxy(i+1,j,k,n) - fxy(i,j,k,n)) &
                            - hdtdy*(fyx(i,j+1,k,n) - fyx(i,j,k,n))
             qpo(i,j,k,nq) = compnr/rrnewr + hdt*srcQ(i,j,k,nq)
          end do
          end do
       end do

       do    k = lo(3)-1, hi(3)
          do j = lo(2)  , hi(2) 
          do i = lo(1)  , hi(1)
             rrl = qm(i,j,k+1,QRHO)
             compl = rrl*qm(i,j,k+1,nq)
             rrnewl = rrl - hdtdx*(fxy(i+1,j,k,URHO) - fxy(i,j,k,URHO)) &
                          - hdtdy*(fyx(i,j+1,k,URHO) - fyx(i,j,k,URHO))
             compnl = compl - hdtdx*(fxy(i+1,j,k,n) - fxy(i,j,k,n)) &
                            - hdtdy*(fyx(i,j+1,k,n) - fyx(i,j,k,n))
             qmo(i,j,k+1,nq) = compnl/rrnewl + hdt*srcQ(i,j,k,nq)
          enddo
          enddo
       enddo
    enddo

    do    k = lo(3)-1, hi(3)+1
       do j = lo(2)  , hi(2) 
       do i = lo(1)  , hi(1)

          !-------------------------------------------------------------------
          ! add the transverse xy and yx differences to the z-states for the
          ! fluid variables
          !-------------------------------------------------------------------          

          pgxp = pgdnvx(i+1,j,k)
          pgxm = pgdnvx(i,j,k)
          ugxp = ugdnvx(i+1,j,k)          
          ugxm = ugdnvx(i,j,k)
          gegxp = gegdnvx(i+1,j,k)          
          gegxm = gegdnvx(i,j,k)
          
          pgyp = pgdnvy(i,j+1,k)
          pgym = pgdnvy(i,j,k)
          ugyp = ugdnvy(i,j+1,k)
          ugym = ugdnvy(i,j,k)
          gegyp = gegdnvy(i,j+1,k)
          gegym = gegdnvy(i,j,k)
          
          pgxpm = pgdnvx(i+1,j,k)
          pgxmm = pgdnvx(i,j,k)
          ugxpm = ugdnvx(i+1,j,k)
          ugxmm = ugdnvx(i,j,k)
          gegxpm = gegdnvx(i+1,j,k)
          gegxmm = gegdnvx(i,j,k)
          
          pgypm = pgdnvy(i,j+1,k)
          pgymm = pgdnvy(i,j,k)
          ugypm = ugdnvy(i,j+1,k)
          ugymm = ugdnvy(i,j,k)
          gegypm = gegdnvy(i,j+1,k)
          gegymm = gegdnvy(i,j,k)
          
          ! Convert to conservation form
          rrr = qp(i,j,k,QRHO)
          rur = rrr*qp(i,j,k,QU)
          rvr = rrr*qp(i,j,k,QV)
          rwr = rrr*qp(i,j,k,QW)
          ekenr = HALF*rrr*(qp(i,j,k,QU)**2 + qp(i,j,k,QV)**2 + &
               qp(i,j,k,QW)**2)
          rer = qp(i,j,k,QREINT) + ekenr
          
          rrl = qm(i,j,k+1,QRHO)
          rul = rrl*qm(i,j,k+1,QU)
          rvl = rrl*qm(i,j,k+1,QV)
          rwl = rrl*qm(i,j,k+1,QW)
          ekenl = HALF*rrl*(qm(i,j,k+1,QU)**2 + qm(i,j,k+1,QV)**2 + &
               qm(i,j,k+1,QW)**2)
          rel = qm(i,j,k+1,QREINT) + ekenl
          
          ! Add transverse predictor
          rrnewr = rrr - hdtdx*(fxy(i+1,j,k,URHO) - fxy(i,j,k,URHO)) &
                       - hdtdy*(fyx(i,j+1,k,URHO) - fyx(i,j,k,URHO))
          runewr = rur - hdtdx*(fxy(i+1,j,k,UMX) - fxy(i,j,k,UMX)) &
                       - hdtdy*(fyx(i,j+1,k,UMX) - fyx(i,j,k,UMX))
          rvnewr = rvr - hdtdx*(fxy(i+1,j,k,UMY) - fxy(i,j,k,UMY)) &
                       - hdtdy*(fyx(i,j+1,k,UMY) - fyx(i,j,k,UMY))
          rwnewr = rwr - hdtdx*(fxy(i+1,j,k,UMZ) - fxy(i,j,k,UMZ)) &
                       - hdtdy*(fyx(i,j+1,k,UMZ) - fyx(i,j,k,UMZ))
          renewr = rer - hdtdx*(fxy(i+1,j,k,UEDEN) - fxy(i,j,k,UEDEN)) &
                       - hdtdy*(fyx(i,j+1,k,UEDEN) - fyx(i,j,k,UEDEN))


          rrnewl = rrl - hdtdx*(fxy(i+1,j,k,URHO) - fxy(i,j,k,URHO)) &
                       - hdtdy*(fyx(i,j+1,k,URHO) - fyx(i,j,k,URHO))
          runewl = rul - hdtdx*(fxy(i+1,j,k,UMX) - fxy(i,j,k,UMX)) &
                       - hdtdy*(fyx(i,j+1,k,UMX) - fyx(i,j,k,UMX))
          rvnewl = rvl - hdtdx*(fxy(i+1,j,k,UMY) - fxy(i,j,k,UMY)) &
                       - hdtdy*(fyx(i,j+1,k,UMY) - fyx(i,j,k,UMY))
          rwnewl = rwl - hdtdx*(fxy(i+1,j,k,UMZ) - fxy(i,j,k,UMZ)) &
                       - hdtdy*(fyx(i,j+1,k,UMZ) - fyx(i,j,k,UMZ))
          renewl = rel - hdtdx*(fxy(i+1,j,k,UEDEN) - fxy(i,j,k,UEDEN)) &
                       - hdtdy*(fyx(i,j+1,k,UEDEN) - fyx(i,j,k,UEDEN))

          ! Reset to original value if adding transverse terms made density negative
          if (transverse_reset_density == 1) then
             if (rrnewr .lt. ZERO) then
                rrnewr = rrr 
                runewr = rur 
                rvnewr = rvr 
                rwnewr = rwr 
                renewr = rer 
             endif
             if (rrnewl .lt. ZERO) then
                rrnewl = rrl 
                runewl = rul 
                rvnewl = rvl 
                rwnewl = rwl 
                renewl = rel 
             endif
          endif

          rhoekenr = HALF*(runewr**2 + rvnewr**2 + rwnewr**2)/rrnewr
          rhoekenl = HALF*(runewl**2 + rvnewl**2 + rwnewl**2)/rrnewl

          duxp = pgxp*ugxp - pgxm*ugxm
          pxav = HALF*(pgxp+pgxm)
          uxav = HALF*(ugxp+ugxm)
          gexav = HALF*(gegxp+gegxm)
          dux = ugxp-ugxm
          dgex = gegxp-gegxm
          pxnew = hdtdx*(duxp + pxav*dux*(gamc(i,j,k)-ONE))
          gexnew = hdtdx*( (gexav-ONE)*(gexav-gamc(i,j,k))*dux - uxav*dgex )

          duxpm = pgxpm*ugxpm - pgxmm*ugxmm
          pxavm = HALF*(pgxpm+pgxmm)
          uxavm = HALF*(ugxpm+ugxmm)
          gexavm = HALF*(gegxpm+gegxmm)
          duxm = ugxpm-ugxmm
          dgexm = gegxpm-gegxmm
          pxnewm = hdtdx*(duxpm + pxavm*duxm*(gamc(i,j,k)-ONE))
          gexnewm = hdtdx*( (gexavm-ONE)*(gexavm-gamc(i,j,k))*duxm - uxavm*dgexm )
          
          duyp = pgyp*ugyp - pgym*ugym
          pyav = HALF*(pgyp+pgym)
          uyav = HALF*(ugyp+ugym)
          geyav = HALF*(gegyp+gegym)
          duy = ugyp-ugym
          dgey = gegyp-gegym
          pynew = hdtdy*(duyp + pyav*duy*(gamc(i,j,k)-ONE))
          geynew = hdtdy*( (geyav-ONE)*(geyav-gamc(i,j,k))*duy - uyav*dgey )

          duypm = pgypm*ugypm - pgymm*ugymm
          pyavm = HALF*(pgypm+pgymm)
          uyavm = HALF*(ugypm+ugymm)
          geyavm = HALF*(gegypm+gegymm)
          duym = ugypm-ugymm
          dgeym = gegypm-gegymm
          pynewm = hdtdy*(duypm + pyavm*duym*(gamc(i,j,k)-ONE))
          geynewm = hdtdy*( (geyavm-ONE)*(geyavm-gamc(i,j,k))*duym - uyavm*dgeym )


          !-------------------------------------------------------------------
          ! qzpo state
          !-------------------------------------------------------------------          
          
          ! Convert back to primitive form
          qpo(i,j,k,QRHO  ) = rrnewr        + hdt*srcQ(i,j,k,QRHO)
          qpo(i,j,k,QU    ) = runewr/rrnewr + hdt*srcQ(i,j,k,QU) 
          qpo(i,j,k,QV    ) = rvnewr/rrnewr + hdt*srcQ(i,j,k,QV) 
          qpo(i,j,k,QW    ) = rwnewr/rrnewr + hdt*srcQ(i,j,k,QW) 

          ! note: we run the risk of (rho e) being negative here
          qpo(i,j,k,QREINT) = renewr - rhoekenr + hdt*srcQ(i,j,k,QREINT)

          if (transverse_reset_rhoe == 1) then
             ! If it is negative, reset the internal energy by using the discretized
             ! expression for updating (rho e).
             
             if (qpo(i,j,k,QREINT) .le. ZERO) then
                qpo(i,j,k,QREINT) = qp(i,j,k,QREINT) &
                     - hdtdx*(fxy(i+1,j,k,UEINT) - fxy(i,j,k,UEINT) + pxav*dux) &
                     - hdtdy*(fyx(i,j+1,k,UEINT) - fyx(i,j,k,UEINT) + pyav*duy) &
                     + hdt*srcQ(i,j,k,QREINT)
                
                ! if we are still negative, then we need to reset
                if (qpo(i,j,k,QREINT) < ZERO) then
                   eos_state % rho = qpo(i,j,k,QRHO)
                   eos_state % T = small_temp
                   eos_state % xn(:) = qpo(i,j,k,QFS:QFS-1+nspec)
                   
                   call eos(eos_input_rt, eos_state)
                   
                   qpo(i,j,k,QREINT) = qpo(i,j,k,QRHO)*eos_state % e
                   qpo(i,j,k,QPRES) = eos_state % p
                endif
             endif
          endif

          if (ppm_predict_gammae == 0) then
             ! Optionally, use the EOS to calculate the pressure.

             if (transverse_use_eos .eq. 1) then
                eos_state % rho = qpo(i,j,k,QRHO)
                eos_state % e   = qpo(i,j,k,QREINT) / qpo(i,j,k,QRHO)
                eos_state % T   = small_temp
                eos_state % xn  = qpo(i,j,k,QFS:QFS+nspec-1)
                
                call eos(eos_input_re, eos_state)

                pnewr = eos_state % p
                qpo(i,j,k,QPRES ) = pnewr
                qpo(i,j,k,QREINT) = eos_state % e * eos_state % rho    
             else
                ! add the transverse term to the p evolution eq here
                pnewr = qp(i,j,k,QPRES) - pxnew - pynew
                qpo(i,j,k,QPRES) = pnewr + hdt*srcQ(i,j,k,QPRES)
             endif

             qpo(i,j,k,QPRES) = max(qpo(i,j,k,QPRES),small_pres)

          else
             
             ! Update gammae with its transverse terms
             qpo(i,j,k,QGAME) = qp(i,j,k,QGAME) + gexnew + geynew

             ! and compute the p edge state from this and (rho e)
             qpo(i,j,k,QPRES) = qpo(i,j,k,QREINT)*(qpo(i,j,k,QGAME)-ONE)
             qpo(i,j,k,QPRES) = max(qpo(i,j,k,QPRES), small_pres)

          endif


          !-------------------------------------------------------------------
          ! qzmo state
          !-------------------------------------------------------------------          

          qmo(i,j,k+1,QRHO  ) = rrnewl        + hdt*srcQ(i,j,k,QRHO)
          qmo(i,j,k+1,QU    ) = runewl/rrnewl + hdt*srcQ(i,j,k,QU)
          qmo(i,j,k+1,QV    ) = rvnewl/rrnewl + hdt*srcQ(i,j,k,QV)
          qmo(i,j,k+1,QW    ) = rwnewl/rrnewl + hdt*srcQ(i,j,k,QW)

          ! note: we run the risk of (rho e) being negative here
          qmo(i,j,k+1,QREINT) = renewl - rhoekenl + hdt*srcQ(i,j,k,QREINT)

          if (transverse_reset_rhoe == 1) then
             ! If it is negative, reset the internal energy by using the discretized
             ! expression for updating (rho e).
             
             if (qmo(i,j,k+1,QREINT) .le. ZERO) then
                qmo(i,j,k+1,QREINT) = qm(i,j,k+1,QREINT) &
                     - hdtdx*(fxy(i+1,j,k,UEINT) - fxy(i,j,k,UEINT) + pxavm*duxm) &
                     - hdtdy*(fyx(i,j+1,k,UEINT) - fyx(i,j,k,UEINT) + pyavm*duym) &
                     + hdt*srcQ(i,j,k,QREINT)
                
                ! if we are still negative, then we need to reset
                if (qmo(i,j,k+1,QREINT) < ZERO) then
                   eos_state % rho = qmo(i,j,k+1,QRHO)
                   eos_state % T = small_temp
                   eos_state % xn(:) = qmo(i,j,k+1,QFS:QFS-1+nspec)
                   
                   call eos(eos_input_rt, eos_state)
                   
                   qmo(i,j,k+1,QREINT) = qmo(i,j,k+1,QRHO)*eos_state % e
                   qmo(i,j,k+1,QPRES) = eos_state % p
                endif
             endif
          endif

          if (ppm_predict_gammae == 0) then

             ! Optionally, use the EOS to calculate the pressure.

             if (transverse_use_eos .eq. 1) then
                eos_state % rho = qmo(i,j,k+1,QRHO)
                eos_state % e   = qmo(i,j,k+1,QREINT) / qmo(i,j,k+1,QRHO)
                eos_state % T   = small_temp
                eos_state % xn  = qmo(i,j,k+1,QFS:QFS+nspec-1)
                
                call eos(eos_input_re, eos_state)

                pnewl = eos_state % p
                qmo(i,j,k+1,QPRES ) = pnewl
                qmo(i,j,k+1,QREINT) = eos_state % e * eos_state % rho
             else
                ! add the transverse term to the p evolution eq here
                pnewl = qm(i,j,k+1,QPRES) - pxnewm - pynewm
                qmo(i,j,k+1,QPRES) = pnewl + hdt*srcQ(i,j,k,QPRES)
             endif

             qmo(i,j,k+1,QPRES) = max(qmo(i,j,k+1,QPRES),small_pres)
          
          else

             ! Update gammae with its transverse terms
             qmo(i,j,k+1,QGAME) = qm(i,j,k+1,QGAME) + gexnewm + geynewm
                 
             ! and compute the p edge state from this and (rho e)
             qmo(i,j,k+1,QPRES) = qmo(i,j,k+1,QREINT)*(qmo(i,j,k+1,QGAME)-ONE)
             qmo(i,j,k+1,QPRES) = max(qmo(i,j,k+1,QPRES), small_pres)

          endif

       enddo
       enddo
    enddo


    ! if ppm_trace_grav == 1, then we already added the piecewise parabolic traced
    ! gravity to the normal edge states
    if ((do_grav .eq. 1) .and. (ppm_trace_grav == 0 .or. ppm_type == 0)) then
       do    k = lo(3), hi(3)+1
          do j = lo(2), hi(2) 
          do i = lo(1), hi(1)
             qpo(i,j,k,QU    ) = qpo(i,j,k,QU    ) + hdt*grav(i,j,k,1)
             qpo(i,j,k,QV    ) = qpo(i,j,k,QV    ) + hdt*grav(i,j,k,2)
             qpo(i,j,k,QW    ) = qpo(i,j,k,QW    ) + hdt*grav(i,j,k,3)
          end do
          end do
       end do

       do    k = lo(3)-1, hi(3)
          do j = lo(2)  , hi(2) 
          do i = lo(1)  , hi(1)             
             qmo(i,j,k+1,QU    ) = qmo(i,j,k+1,QU    ) + hdt*grav(i,j,k,1)
             qmo(i,j,k+1,QV    ) = qmo(i,j,k+1,QV    ) + hdt*grav(i,j,k,2)
             qmo(i,j,k+1,QW    ) = qmo(i,j,k+1,QW    ) + hdt*grav(i,j,k,3)
          enddo
          enddo
       enddo
    endif

    ! if ppm_trace_rot == 1, then we already added the piecewise parabolic traced
    ! rotation to the normal edge states
    if ((do_rotation .eq. 1) .and. (ppm_trace_rot == 0 .or. ppm_type == 0)) then
       do    k = lo(3), hi(3)+1
          do j = lo(2), hi(2) 
          do i = lo(1), hi(1)
             qpo(i,j,k,QU    ) = qpo(i,j,k,QU    ) + hdt*rot(i,j,k,1)
             qpo(i,j,k,QV    ) = qpo(i,j,k,QV    ) + hdt*rot(i,j,k,2)
             qpo(i,j,k,QW    ) = qpo(i,j,k,QW    ) + hdt*rot(i,j,k,3)
          end do
          end do
       end do

       do    k = lo(3)-1, hi(3)
          do j = lo(2)  , hi(2) 
          do i = lo(1)  , hi(1)
             qmo(i,j,k+1,QU    ) = qmo(i,j,k+1,QU    ) + hdt*rot(i,j,k,1)
             qmo(i,j,k+1,QV    ) = qmo(i,j,k+1,QV    ) + hdt*rot(i,j,k,2)
             qmo(i,j,k+1,QW    ) = qmo(i,j,k+1,QW    ) + hdt*rot(i,j,k,3)
          enddo
          enddo
       enddo
    endif
    
  end subroutine transxy

  
  !===========================================================================
  ! transxz
  !===========================================================================
  subroutine transxz(qm,qmo,qp,qpo, qlo, qhi, &
                     fxz,fzx, flo, fhi, &
                     ugdnvx,pgdnvx,gegdnvx,ugdnvz,pgdnvz,gegdnvz, gdlo, gdhi, &
                     gamc,gclo,gchi, &
                     srcQ,slo,shi, &
                     grav,gvlo,gvhi, &
                     rot,rlo,rhi, &
                     hdt,hdtdx,hdtdz,lo,hi)
    
    use network, only : nspec, naux
    use meth_params_module, only : QVAR, NVAR, QRHO, QU, QV, QW, &
                                   QPRES, QREINT, QGAME, QESGS, QFA, QFS, &
                                   URHO, UMX, UMY, UMZ, UEDEN, UEINT, UESGS, UFA, UFS, &
                                   nadv, small_pres, small_temp, &
                                   npassive, upass_map, qpass_map, &
                                   ppm_predict_gammae, &
                                   transverse_use_eos, transverse_reset_density, transverse_reset_rhoe, &
                                   ppm_type, ppm_trace_grav, rot_period, ppm_trace_rot, &
                                   do_grav, do_rotation
    use eos_module

    implicit none

    integer, intent(in) :: qlo(3),qhi(3),flo(3),fhi(3),gdlo(3),gdhi(3),gclo(3),gchi(3), &
         slo(3),shi(3),gvlo(3),gvhi(3),rlo(3),rhi(3),lo(3),hi(3)
    double precision, intent(in) :: hdt, hdtdx, hdtdz
    double precision,intent(in )::     qm( qlo(1): qhi(1), qlo(2): qhi(2), qlo(3): qhi(3),QVAR)
    double precision,intent(out)::    qmo( qlo(1): qhi(1), qlo(2): qhi(2), qlo(3): qhi(3),QVAR)
    double precision,intent(in )::     qp( qlo(1): qhi(1), qlo(2): qhi(2), qlo(3): qhi(3),QVAR)
    double precision,intent(out)::    qpo( qlo(1): qhi(1), qlo(2): qhi(2), qlo(3): qhi(3),QVAR)
    double precision,intent(in )::    fxz( flo(1): fhi(1), flo(2): fhi(2), flo(3): fhi(3),NVAR)
    double precision,intent(in )::    fzx( flo(1): fhi(1), flo(2): fhi(2), flo(3): fhi(3),NVAR)
    double precision,intent(in ):: ugdnvx(gdlo(1):gdhi(1),gdlo(2):gdhi(2),gdlo(3):gdhi(3))
    double precision,intent(in ):: pgdnvx(gdlo(1):gdhi(1),gdlo(2):gdhi(2),gdlo(3):gdhi(3))
    double precision,intent(in )::gegdnvx(gdlo(1):gdhi(1),gdlo(2):gdhi(2),gdlo(3):gdhi(3))
    double precision,intent(in ):: ugdnvz(gdlo(1):gdhi(1),gdlo(2):gdhi(2),gdlo(3):gdhi(3))
    double precision,intent(in ):: pgdnvz(gdlo(1):gdhi(1),gdlo(2):gdhi(2),gdlo(3):gdhi(3))
    double precision,intent(in )::gegdnvz(gdlo(1):gdhi(1),gdlo(2):gdhi(2),gdlo(3):gdhi(3))
    double precision,intent(in )::   gamc(gclo(1):gchi(1),gclo(2):gchi(2),gclo(3):gchi(3))
    double precision,intent(in )::   srcQ( slo(1): shi(1), slo(2): shi(2), slo(3): shi(3),QVAR)
    double precision,intent(in )::   grav(gvlo(1):gvhi(1),gvlo(2):gvhi(2),gvlo(3):gvhi(3),3)
    double precision,intent(in )::    rot( rlo(1): rhi(1), rlo(2): rhi(2), rlo(3): rhi(3),3)

    integer i, j, k, n, nq, ipassive    
    
    double precision rrr, rur, rvr, rwr, rer, ekenr, rhoekenr
    double precision rrl, rul, rvl, rwl, rel, ekenl, rhoekenl
    double precision rrnewr, runewr, rvnewr, rwnewr, renewr
    double precision rrnewl, runewl, rvnewl, rwnewl, renewl
    double precision pnewr, pnewl
    double precision pgxp, pgxm, ugxp, ugxm, gegxp, gegxm, duxp, pxav, dux, pxnew, gexnew
    double precision pgzp, pgzm, ugzp, ugzm, gegzp, gegzm, duzp, pzav, duz, pznew, geznew
    double precision uxav, gexav, dgex, uzav, gezav, dgez
    double precision compr, compl, compnr, compnl

    type (eos_t) :: eos_state

    !-------------------------------------------------------------------------    
    ! update all of the passively-advected quantities with the
    ! transerse term and convert back to the primitive quantity
    !-------------------------------------------------------------------------

    do ipassive = 1,npassive
       n  = upass_map(ipassive)
       nq = qpass_map(ipassive)

       do    k = lo(3), hi(3)

          do j = lo(2), hi(2)+1 
          do i = lo(1), hi(1)

             rrr = qp(i,j,k,QRHO)
             compr = rrr*qp(i,j,k,nq)
             rrnewr = rrr - hdtdx*(fxz(i+1,j,k  ,URHO) - fxz(i,j,k,URHO)) &
                          - hdtdz*(fzx(i  ,j,k+1,URHO) - fzx(i,j,k,URHO))
             compnr = compr - hdtdx*(fxz(i+1,j,k  ,n) - fxz(i,j,k,n)) &
                            - hdtdz*(fzx(i  ,j,k+1,n) - fzx(i,j,k,n))
             qpo(i,j  ,k,nq) = compnr/rrnewr + hdt*srcQ(i,j,k,nq)

          enddo
          enddo

          do j = lo(2)-1, hi(2) 
          do i = lo(1)  , hi(1)

             rrl = qm(i,j+1,k,QRHO)             
             compl = rrl*qm(i,j+1,k,nq)
             rrnewl = rrl - hdtdx*(fxz(i+1,j,k  ,URHO) - fxz(i,j,k,URHO)) &
                          - hdtdz*(fzx(i  ,j,k+1,URHO) - fzx(i,j,k,URHO))
             compnl = compl - hdtdx*(fxz(i+1,j,k  ,n) - fxz(i,j,k,n)) &
                            - hdtdz*(fzx(i  ,j,k+1,n) - fzx(i,j,k,n))
             qmo(i,j+1,k,nq) = compnl/rrnewl + hdt*srcQ(i,j,k,nq)

          enddo
          enddo
       enddo
    enddo

    do    k = lo(3)  , hi(3)
       do j = lo(2)-1, hi(2)+1 
       do i = lo(1)  , hi(1)

          !-------------------------------------------------------------------
          ! add the transverse xz and zx differences to the y-states for the
          ! fluid variables
          !-------------------------------------------------------------------            
          pgxp = pgdnvx(i+1,j,k)
          pgxm = pgdnvx(i,j,k)
          ugxp = ugdnvx(i+1,j,k)
          ugxm = ugdnvx(i,j,k)
          gegxp = gegdnvx(i+1,j,k)
          gegxm = gegdnvx(i,j,k)
          
          pgzp = pgdnvz(i,j,k+1)
          pgzm = pgdnvz(i,j,k)
          ugzp = ugdnvz(i,j,k+1)
          ugzm = ugdnvz(i,j,k)
          gegzp = gegdnvz(i,j,k+1)
          gegzm = gegdnvz(i,j,k)

          ! Convert to conservation form
          rrr = qp(i,j,k,QRHO)
          rur = rrr*qp(i,j,k,QU)
          rvr = rrr*qp(i,j,k,QV)
          rwr = rrr*qp(i,j,k,QW)
          ekenr = HALF*rrr*(qp(i,j,k,QU)**2 + qp(i,j,k,QV)**2 + qp(i,j,k,QW)**2)
          rer = qp(i,j,k,QREINT) + ekenr
          
          rrl = qm(i,j+1,k,QRHO)
          rul = rrl*qm(i,j+1,k,QU)
          rvl = rrl*qm(i,j+1,k,QV)
          rwl = rrl*qm(i,j+1,k,QW)
          ekenl = HALF*rrl*(qm(i,j+1,k,QU)**2 + qm(i,j+1,k,QV)**2 + qm(i,j+1,k,QW)**2)
          rel = qm(i,j+1,k,QREINT) + ekenl
          
          ! Add transverse predictor
          rrnewr = rrr - hdtdx*(fxz(i+1,j,k,URHO) - fxz(i,j,k,URHO)) &
                       - hdtdz*(fzx(i,j,k+1,URHO) - fzx(i,j,k,URHO))
          runewr = rur - hdtdx*(fxz(i+1,j,k,UMX) - fxz(i,j,k,UMX)) &
                       - hdtdz*(fzx(i,j,k+1,UMX) - fzx(i,j,k,UMX))
          rvnewr = rvr - hdtdx*(fxz(i+1,j,k,UMY) - fxz(i,j,k,UMY)) &
                       - hdtdz*(fzx(i,j,k+1,UMY) - fzx(i,j,k,UMY))
          rwnewr = rwr - hdtdx*(fxz(i+1,j,k,UMZ) - fxz(i,j,k,UMZ)) &
                       - hdtdz*(fzx(i,j,k+1,UMZ) - fzx(i,j,k,UMZ))
          renewr = rer - hdtdx*(fxz(i+1,j,k,UEDEN) - fxz(i,j,k,UEDEN)) &
                       - hdtdz*(fzx(i,j,k+1,UEDEN) - fzx(i,j,k,UEDEN))

          rrnewl = rrl - hdtdx*(fxz(i+1,j,k,URHO) - fxz(i,j,k,URHO)) &
                       - hdtdz*(fzx(i,j,k+1,URHO) - fzx(i,j,k,URHO))
          runewl = rul - hdtdx*(fxz(i+1,j,k,UMX) - fxz(i,j,k,UMX)) &
                       - hdtdz*(fzx(i,j,k+1,UMX) - fzx(i,j,k,UMX))
          rvnewl = rvl - hdtdx*(fxz(i+1,j,k,UMY) - fxz(i,j,k,UMY)) &
                       - hdtdz*(fzx(i,j,k+1,UMY) - fzx(i,j,k,UMY))
          rwnewl = rwl - hdtdx*(fxz(i+1,j,k,UMZ) - fxz(i,j,k,UMZ)) &
                       - hdtdz*(fzx(i,j,k+1,UMZ) - fzx(i,j,k,UMZ))
          renewl = rel - hdtdx*(fxz(i+1,j,k,UEDEN) - fxz(i,j,k,UEDEN)) &
                       - hdtdz*(fzx(i,j,k+1,UEDEN) - fzx(i,j,k,UEDEN))

          ! Reset to original value if adding transverse terms made density negative
          if (transverse_reset_density == 1) then
             if (rrnewr .lt. ZERO) then
                rrnewr = rrr 
                runewr = rur
                rvnewr = rvr 
                rwnewr = rwr 
                renewr = rer 
             endif
             if (rrnewl .lt. ZERO) then
                rrnewl = rrl 
                runewl = rul 
                rvnewl = rvl 
                rwnewl = rwl 
                renewl = rel 
             endif
          endif

          rhoekenr = HALF*(runewr**2 + rvnewr**2 + rwnewr**2)/rrnewr
          rhoekenl = HALF*(runewl**2 + rvnewl**2 + rwnewl**2)/rrnewl

          duxp = pgxp*ugxp - pgxm*ugxm
          pxav = HALF*(pgxp+pgxm)
          uxav = HALF*(ugxp+ugxm)
          gexav = HALF*(gegxp+gegxm)
          dux = ugxp-ugxm
          dgex = gegxp-gegxm
          pxnew = hdtdx*(duxp + pxav*dux*(gamc(i,j,k)-ONE))
          gexnew = hdtdx*( (gexav-ONE)*(gexav-gamc(i,j,k))*dux - uxav*dgex )

          duzp = pgzp*ugzp - pgzm*ugzm
          pzav = HALF*(pgzp+pgzm)
          uzav = HALF*(ugzp+ugzm)
          gezav = HALF*(gegzp+gegzm)
          duz = ugzp-ugzm
          dgez = gegzp-gegzm
          pznew = hdtdz*(duzp + pzav*duz*(gamc(i,j,k)-ONE))
          geznew = hdtdz*( (gezav-ONE)*(gezav-gamc(i,j,k))*duz - uzav*dgez )


          !-------------------------------------------------------------------
          ! qypo state
          !-------------------------------------------------------------------

          ! Convert back to primitive form
          if (j.ge.lo(2)) then
             qpo(i,j,k,QRHO  ) = rrnewr        + hdt*srcQ(i,j,k,QRHO)
             qpo(i,j,k,QU    ) = runewr/rrnewr + hdt*srcQ(i,j,k,QU)
             qpo(i,j,k,QV    ) = rvnewr/rrnewr + hdt*srcQ(i,j,k,QV)
             qpo(i,j,k,QW    ) = rwnewr/rrnewr + hdt*srcQ(i,j,k,QW)

             ! note: we run the risk of (rho e) being negative here
             qpo(i,j,k,QREINT) = renewr - rhoekenr + hdt*srcQ(i,j,k,QREINT)

             if (transverse_reset_rhoe == 1) then
                ! If it is negative, reset the internal energy by using the discretized
                ! expression for updating (rho e).

                if (qpo(i,j,k,QREINT) .le. ZERO) then
                   qpo(i,j,k,QREINT) = qp(i,j,k,QREINT) &
                        - hdtdx*(fxz(i+1,j,k,UEINT) - fxz(i,j,k,UEINT) + pxav*dux) &
                        - hdtdz*(fzx(i  ,j,k+1,UEINT) - fzx(i,j,k,UEINT) + pzav*duz) &
                        + hdt*srcQ(i,j,k,QREINT)
                   
                   ! if we are still negative, then we need to reset
                   if (qpo(i,j,k,QREINT) < ZERO) then
                      eos_state % rho = qpo(i,j,k,QRHO)
                      eos_state % T = small_temp
                      eos_state % xn(:) = qpo(i,j,k,QFS:QFS-1+nspec)
                      
                      call eos(eos_input_rt, eos_state)
                      
                      qpo(i,j,k,QREINT) = qpo(i,j,k,QRHO)*eos_state % e
                      qpo(i,j,k,QPRES) = eos_state % p
                   endif
                endif
             endif

             if (ppm_predict_gammae == 0) then
                
                ! Optionally, use the EOS to calculate the pressure.

                if (transverse_use_eos .eq. 1) then

                   eos_state % rho = qpo(i,j,k,QRHO)
                   eos_state % e   = qpo(i,j,k,QREINT) / qpo(i,j,k,QRHO)
                   eos_state % T   = small_temp
                   eos_state % xn  = qpo(i,j,k,QFS:QFS+nspec-1)

                   call eos(eos_input_re, eos_state)

                   pnewr = eos_state % p
                   qpo(i,j,k,QPRES ) = pnewr
                   qpo(i,j,k,QREINT) = eos_state % e * eos_state % rho
                else
                   ! add the transverse term to the p evolution eq here
                   pnewr = qp(i,j,k,QPRES) - pxnew - pznew
                   qpo(i,j,k,QPRES) = pnewr + hdt*srcQ(i,j,k,QPRES)
                endif
                
                qpo(i,j,k,QPRES) = max(qpo(i,j,k,QPRES),small_pres)

             else
                
                ! Update gammae with its transverse terms
                qpo(i,j,k,QGAME) = qp(i,j,k,QGAME) + gexnew + geznew

                ! and compute the p edge state from this and (rho e)
                qpo(i,j,k,QPRES) = qpo(i,j,k,QREINT)*(qpo(i,j,k,QGAME)-ONE)
                qpo(i,j,k,QPRES) = max(qpo(i,j,k,QPRES), small_pres)

             endif

          end if


          !-------------------------------------------------------------------
          ! qymo state
          !-------------------------------------------------------------------          

          if (j.le.hi(2)) then
             qmo(i,j+1,k,QRHO  ) = rrnewl        + hdt*srcQ(i,j,k,QRHO)
             qmo(i,j+1,k,QU    ) = runewl/rrnewl + hdt*srcQ(i,j,k,QU)
             qmo(i,j+1,k,QV    ) = rvnewl/rrnewl + hdt*srcQ(i,j,k,QV)
             qmo(i,j+1,k,QW    ) = rwnewl/rrnewl + hdt*srcQ(i,j,k,QW)

             ! note: we run the risk of (rho e) being negative here
             qmo(i,j+1,k,QREINT) = renewl - rhoekenl + hdt*srcQ(i,j,k,QREINT)

             if (transverse_reset_rhoe == 1) then
                ! If it is negative, reset the internal energy by using the discretized
                ! expression for updating (rho e).
                
                if (qmo(i,j+1,k,QREINT) .le. ZERO) then
                   qmo(i,j+1,k,QREINT) = qm(i,j+1,k,QREINT) &
                        - hdtdx*(fxz(i+1,j,k,UEINT) - fxz(i,j,k,UEINT) + pxav*dux) &
                        - hdtdz*(fzx(i,j,k+1,UEINT) - fzx(i,j,k,UEINT) + pzav*duz) &
                        + hdt*srcQ(i,j,k,QREINT)
                   
                   ! if we are still negative, then we need to reset
                   if (qmo(i,j+1,k,QREINT) < ZERO) then
                      eos_state % rho = qmo(i,j+1,k,QRHO)
                      eos_state % T = small_temp
                      eos_state % xn(:) = qmo(i,j+1,k,QFS:QFS-1+nspec)
                      
                      call eos(eos_input_rt, eos_state)

                      qmo(i,j+1,k,QREINT) = qmo(i,j+1,k,QRHO)*eos_state % e
                      qmo(i,j+1,k,QPRES) = eos_state % p
                   endif
                endif
             endif
             
             
             if (ppm_predict_gammae == 0) then

                ! Optionally, use the EOS to calculate the pressure.

                if (transverse_use_eos .eq. 1) then
                   eos_state % rho = qmo(i,j+1,k,QRHO)
                   eos_state % e   = qmo(i,j+1,k,QREINT) / qmo(i,j+1,k,QRHO)
                   eos_state % T   = small_temp
                   eos_state % xn  = qmo(i,j+1,k,QFS:QFS+nspec-1)
                   
                   call eos(eos_input_re, eos_state)

                   pnewl = eos_state % p
                   qmo(i,j+1,k,QPRES ) = pnewl
                   qmo(i,j+1,k,QREINT) = eos_state % e * eos_state % rho
                else
                   pnewl = qm(i,j+1,k,QPRES) - pxnew - pznew
                   qmo(i,j+1,k,QPRES) = pnewl + hdt*srcQ(i,j,k,QPRES)
                endif

                qmo(i,j+1,k,QPRES) = max(qmo(i,j+1,k,QPRES),small_pres)
             
             else

                ! Update gammae with its transverse terms
                qmo(i,j+1,k,QGAME) = qm(i,j+1,k,QGAME) + gexnew + geznew

                ! and compute the p edge state from this and (rho e)
                qmo(i,j+1,k,QPRES) = qmo(i,j+1,k,QREINT)*(qmo(i,j+1,k,QGAME)-ONE)
                qmo(i,j+1,k,QPRES) = max(qmo(i,j+1,k,QPRES), small_pres)

             endif

          endif
          
       enddo
       enddo
    enddo

    ! if ppm_trace_grav == 1, then we already added the piecewise parabolic traced
    ! gravity to the normal edge states
    if ((do_grav .eq. 1) .and. (ppm_trace_grav == 0 .or. ppm_type == 0)) then
       do    k = lo(3), hi(3)
          do j = lo(2), hi(2)+1 
          do i = lo(1), hi(1)
             qpo(i,j,k,QU    ) = qpo(i,j,k,QU    ) + hdt*grav(i,j,k,1)
             qpo(i,j,k,QV    ) = qpo(i,j,k,QV    ) + hdt*grav(i,j,k,2)
             qpo(i,j,k,QW    ) = qpo(i,j,k,QW    ) + hdt*grav(i,j,k,3)
          end do
          end do
          do j = lo(2)-1, hi(2)
          do i = lo(1)  , hi(1)
             qmo(i,j+1,k,QU    ) = qmo(i,j+1,k,QU    ) + hdt*grav(i,j,k,1)
             qmo(i,j+1,k,QV    ) = qmo(i,j+1,k,QV    ) + hdt*grav(i,j,k,2)
             qmo(i,j+1,k,QW    ) = qmo(i,j+1,k,QW    ) + hdt*grav(i,j,k,3)
          enddo
          enddo
       enddo
    endif

    ! if ppm_trace_rot == 1, then we already added the piecewise parabolic traced
    ! rotation to the normal edge states
    if ((do_rotation .eq. 1) .and. (ppm_trace_rot == 0 .or. ppm_type == 0)) then
       do    k = lo(3), hi(3)
          do j = lo(2), hi(2)+1 
          do i = lo(1), hi(1)
             qpo(i,j,k,QU    ) = qpo(i,j,k,QU    ) + hdt*rot(i,j,k,1)
             qpo(i,j,k,QV    ) = qpo(i,j,k,QV    ) + hdt*rot(i,j,k,2)
             qpo(i,j,k,QW    ) = qpo(i,j,k,QW    ) + hdt*rot(i,j,k,3)
          end do
          end do
          do j = lo(2)-1, hi(2)
          do i = lo(1)  , hi(1)             
             qmo(i,j+1,k,QU    ) = qmo(i,j+1,k,QU    ) + hdt*rot(i,j,k,1)
             qmo(i,j+1,k,QV    ) = qmo(i,j+1,k,QV    ) + hdt*rot(i,j,k,2)
             qmo(i,j+1,k,QW    ) = qmo(i,j+1,k,QW    ) + hdt*rot(i,j,k,3)
          enddo
          enddo
       enddo
    endif
    
  end subroutine transxz


  !===========================================================================
  ! transyz
  !===========================================================================
  subroutine transyz(qm,qmo,qp,qpo, qlo, qhi, &
                     fyz,fzy, flo, fhi, &
                     ugdnvy,pgdnvy,gegdnvy,ugdnvz,pgdnvz,gegdnvz, gdlo, gdhi, &
                     gamc,gclo,gchi, &
                     srcQ,slo,shi, &
                     grav,gvlo,gvhi, &
                     rot,rlo,rhi, &
                     hdt,hdtdy,hdtdz,lo,hi)
    
    use network, only : nspec, naux
    use meth_params_module, only : QVAR, NVAR, QRHO, QU, QV, QW, &
                                   QPRES, QREINT, QGAME, QESGS, QFA, QFS, &
                                   URHO, UMX, UMY, UMZ, UEDEN, UEINT, UESGS, UFA, UFS, &
                                   nadv, small_pres, small_temp, &
                                   npassive, upass_map, qpass_map, &
                                   ppm_predict_gammae, &
                                   transverse_use_eos, transverse_reset_density, transverse_reset_rhoe, &
                                   ppm_type, ppm_trace_grav, rot_period, ppm_trace_rot, &
                                   do_grav, do_rotation
    use eos_module

    implicit none

    integer, intent(in) :: qlo(3),qhi(3),flo(3),fhi(3),gdlo(3),gdhi(3),gclo(3),gchi(3), &
         slo(3),shi(3),gvlo(3),gvhi(3),rlo(3),rhi(3),lo(3),hi(3)
    double precision, intent(in) :: hdt, hdtdy, hdtdz
    double precision,intent(in )::     qm( qlo(1): qhi(1), qlo(2): qhi(2), qlo(3): qhi(3),QVAR)
    double precision,intent(out)::    qmo( qlo(1): qhi(1), qlo(2): qhi(2), qlo(3): qhi(3),QVAR)
    double precision,intent(in )::     qp( qlo(1): qhi(1), qlo(2): qhi(2), qlo(3): qhi(3),QVAR)
    double precision,intent(out)::    qpo( qlo(1): qhi(1), qlo(2): qhi(2), qlo(3): qhi(3),QVAR)
    double precision,intent(in )::    fyz( flo(1): fhi(1), flo(2): fhi(2), flo(3): fhi(3),NVAR)
    double precision,intent(in )::    fzy( flo(1): fhi(1), flo(2): fhi(2), flo(3): fhi(3),NVAR)
    double precision,intent(in ):: ugdnvy(gdlo(1):gdhi(1),gdlo(2):gdhi(2),gdlo(3):gdhi(3))
    double precision,intent(in ):: pgdnvy(gdlo(1):gdhi(1),gdlo(2):gdhi(2),gdlo(3):gdhi(3))
    double precision,intent(in )::gegdnvy(gdlo(1):gdhi(1),gdlo(2):gdhi(2),gdlo(3):gdhi(3))
    double precision,intent(in ):: ugdnvz(gdlo(1):gdhi(1),gdlo(2):gdhi(2),gdlo(3):gdhi(3))
    double precision,intent(in ):: pgdnvz(gdlo(1):gdhi(1),gdlo(2):gdhi(2),gdlo(3):gdhi(3))
    double precision,intent(in )::gegdnvz(gdlo(1):gdhi(1),gdlo(2):gdhi(2),gdlo(3):gdhi(3))
    double precision,intent(in )::   gamc(gclo(1):gchi(1),gclo(2):gchi(2),gclo(3):gchi(3))
    double precision,intent(in )::   srcQ( slo(1): shi(1), slo(2): shi(2), slo(3): shi(3),QVAR)
    double precision,intent(in )::   grav(gvlo(1):gvhi(1),gvlo(2):gvhi(2),gvlo(3):gvhi(3),3)
    double precision,intent(in )::    rot( rlo(1): rhi(1), rlo(2): rhi(2), rlo(3): rhi(3),3)

    integer i, j, k, n, nq, ipassive
    
    double precision rrr, rur, rvr, rwr, rer, ekenr, rhoekenr
    double precision rrl, rul, rvl, rwl, rel, ekenl, rhoekenl
    double precision rrnewr, runewr, rvnewr, rwnewr, renewr
    double precision rrnewl, runewl, rvnewl, rwnewl, renewl
    double precision pnewr, pnewl
    double precision pgyp, pgym, ugyp, ugym, gegyp, gegym, duyp, pyav, duy, pynew, geynew
    double precision pgzp, pgzm, ugzp, ugzm, gegzp, gegzm, duzp, pzav, duz, pznew, geznew
    double precision uyav, geyav, dgey, uzav, gezav, dgez
    double precision compr, compl, compnr, compnl
    
    type (eos_t) :: eos_state

    !-------------------------------------------------------------------------
    ! update all of the passively-advected quantities with the
    ! transerse term and convert back to the primitive quantity
    !-------------------------------------------------------------------------

    do ipassive = 1,npassive
       n  = upass_map(ipassive)
       nq = qpass_map(ipassive)

       do    k = lo(3), hi(3)
       do    j = lo(2), hi(2) 
          do i = lo(1), hi(1)+1
             rrr = qp(i,j,k,QRHO)
             compr = rrr*qp(i,j,k,nq)
             rrnewr = rrr - hdtdy*(fyz(i,j+1,k  ,URHO) - fyz(i,j,k,URHO)) &
                          - hdtdz*(fzy(i,j  ,k+1,URHO) - fzy(i,j,k,URHO))
             compnr = compr - hdtdy*(fyz(i,j+1,k,n) - fyz(i,j,k,n)) &
                            - hdtdz*(fzy(i,j  ,k+1,n) - fzy(i,j,k,n))
             qpo(i  ,j,k,nq) = compnr/rrnewr + hdt*srcQ(i,j,k,nq)
          enddo
          do i = lo(1)-1, hi(1)
             rrl = qm(i+1,j,k,QRHO)
             compl = rrl*qm(i+1,j,k,nq)
             rrnewl = rrl - hdtdy*(fyz(i,j+1,k  ,URHO) - fyz(i,j,k,URHO)) &
                          - hdtdz*(fzy(i,j  ,k+1,URHO) - fzy(i,j,k,URHO))
             compnl = compl - hdtdy*(fyz(i,j+1,k,n) - fyz(i,j,k,n)) &
                            - hdtdz*(fzy(i,j  ,k+1,n) - fzy(i,j,k,n))
             qmo(i+1,j,k,nq) = compnl/rrnewl + hdt*srcQ(i,j,k,nq)             
          enddo
       enddo
       enddo
    enddo

    do    k = lo(3)  , hi(3)
       do j = lo(2)  , hi(2) 
       do i = lo(1)-1, hi(1)+1

          !-------------------------------------------------------------------
          ! add the transverse yz and zy differences to the x-states for the 
          ! fluid variables
          !-------------------------------------------------------------------          

          pgyp = pgdnvy(i,j+1,k)
          pgym = pgdnvy(i,j,k)
          ugyp = ugdnvy(i,j+1,k)
          ugym = ugdnvy(i,j,k)
          gegyp = gegdnvy(i,j+1,k)
          gegym = gegdnvy(i,j,k)
          
          pgzp = pgdnvz(i,j,k+1)
          pgzm = pgdnvz(i,j,k)
          ugzp = ugdnvz(i,j,k+1)
          ugzm = ugdnvz(i,j,k)
          gegzp = gegdnvz(i,j,k+1)
          gegzm = gegdnvz(i,j,k)
          
          ! Convert to conservation form
          rrr = qp(i,j,k,QRHO)
          rur = rrr*qp(i,j,k,QU)
          rvr = rrr*qp(i,j,k,QV)
          rwr = rrr*qp(i,j,k,QW)
          ekenr = HALF*rrr*(qp(i,j,k,QU)**2 + qp(i,j,k,QV)**2 + &
               qp(i,j,k,QW)**2)
          rer = qp(i,j,k,QREINT) + ekenr
          
          rrl = qm(i+1,j,k,QRHO)
          rul = rrl*qm(i+1,j,k,QU)
          rvl = rrl*qm(i+1,j,k,QV)
          rwl = rrl*qm(i+1,j,k,QW)
          ekenl = HALF*rrl*(qm(i+1,j,k,QU)**2 + qm(i+1,j,k,QV)**2 + &
               qm(i+1,j,k,QW)**2)
          rel = qm(i+1,j,k,QREINT) + ekenl
          
          ! Add transverse predictor
          rrnewr = rrr - hdtdy*(fyz(i,j+1,k,URHO) - fyz(i,j,k,URHO)) &
                       - hdtdz*(fzy(i,j,k+1,URHO) - fzy(i,j,k,URHO))
          runewr = rur - hdtdy*(fyz(i,j+1,k,UMX) - fyz(i,j,k,UMX)) &
                       - hdtdz*(fzy(i,j,k+1,UMX) - fzy(i,j,k,UMX))
          rvnewr = rvr - hdtdy*(fyz(i,j+1,k,UMY) - fyz(i,j,k,UMY)) &
                       - hdtdz*(fzy(i,j,k+1,UMY) - fzy(i,j,k,UMY))
          rwnewr = rwr - hdtdy*(fyz(i,j+1,k,UMZ) - fyz(i,j,k,UMZ)) &
                       - hdtdz*(fzy(i,j,k+1,UMZ) - fzy(i,j,k,UMZ))
          renewr = rer - hdtdy*(fyz(i,j+1,k,UEDEN) - fyz(i,j,k,UEDEN)) &
                       - hdtdz*(fzy(i,j,k+1,UEDEN) - fzy(i,j,k,UEDEN))

          rrnewl = rrl - hdtdy*(fyz(i,j+1,k,URHO) - fyz(i,j,k,URHO)) &
                       - hdtdz*(fzy(i,j,k+1,URHO) - fzy(i,j,k,URHO))
          runewl = rul - hdtdy*(fyz(i,j+1,k,UMX) - fyz(i,j,k,UMX)) &
                       - hdtdz*(fzy(i,j,k+1,UMX) - fzy(i,j,k,UMX))
          rvnewl = rvl - hdtdy*(fyz(i,j+1,k,UMY) - fyz(i,j,k,UMY)) &
                       - hdtdz*(fzy(i,j,k+1,UMY) - fzy(i,j,k,UMY))
          rwnewl = rwl - hdtdy*(fyz(i,j+1,k,UMZ) - fyz(i,j,k,UMZ)) &
                       - hdtdz*(fzy(i,j,k+1,UMZ) - fzy(i,j,k,UMZ))
          renewl = rel - hdtdy*(fyz(i,j+1,k,UEDEN) - fyz(i,j,k,UEDEN)) &
                       - hdtdz*(fzy(i,j,k+1,UEDEN) - fzy(i,j,k,UEDEN))

          ! Reset to original value if adding transverse terms made density negative
          if (transverse_reset_density == 1) then
             if (rrnewr .lt. ZERO) then
                rrnewr = rrr 
                runewr = rur 
                rvnewr = rvr 
                rwnewr = rwr 
                renewr = rer 
             endif
             if (rrnewl .lt. ZERO) then
                rrnewl = rrl 
                runewl = rul 
                rvnewl = rvl 
                rwnewl = rwl 
                renewl = rel 
             endif
          endif

          rhoekenr = HALF*(runewr**2 + rvnewr**2 + rwnewr**2)/rrnewr
          rhoekenl = HALF*(runewl**2 + rvnewl**2 + rwnewl**2)/rrnewl

          duyp = pgyp*ugyp - pgym*ugym
          pyav = HALF*(pgyp+pgym)
          uyav = HALF*(ugyp+ugym)
          geyav = HALF*(gegyp+gegym)
          duy = ugyp-ugym
          dgey = gegyp-gegym
          pynew = hdtdy*(duyp + pyav*duy*(gamc(i,j,k)-ONE))
          geynew = hdtdy*( (geyav-ONE)*(geyav-gamc(i,j,k))*duy - uyav*dgey )
          
          duzp = pgzp*ugzp - pgzm*ugzm
          pzav = HALF*(pgzp+pgzm)
          uzav = HALF*(ugzp+ugzm)
          gezav = HALF*(gegzp+gegzm)
          duz = ugzp-ugzm
          dgez = gegzp-gegzm
          pznew = hdtdz*(duzp + pzav*duz*(gamc(i,j,k)-ONE))
          geznew = hdtdz*( (gezav-ONE)*(gezav-gamc(i,j,k))*duz - uzav*dgez )

          !-------------------------------------------------------------------
          ! qxpo state
          !-------------------------------------------------------------------
          
          ! Convert back to primitive form
          if (i.ge.lo(1)) then
             qpo(i,j,k,QRHO  ) = rrnewr        + hdt*srcQ(i,j,k,QRHO)
             qpo(i,j,k,QU    ) = runewr/rrnewr + hdt*srcQ(i,j,k,QU)
             qpo(i,j,k,QV    ) = rvnewr/rrnewr + hdt*srcQ(i,j,k,QV)
             qpo(i,j,k,QW    ) = rwnewr/rrnewr + hdt*srcQ(i,j,k,QW)

             ! note: we run the risk of (rho e) being negative here
             qpo(i,j,k,QREINT) = renewr - rhoekenr + hdt*srcQ(i,j,k,QREINT)

             if (transverse_reset_rhoe == 1) then
                ! If it is negative, reset the internal energy by using the discretized
                ! expression for updating (rho e).
                
                if (qpo(i,j,k,QREINT) .le. ZERO) then
                   qpo(i,j,k,QREINT) = qp(i,j,k,QREINT) &
                        - hdtdy*(fyz(i,j+1,k,UEINT) - fyz(i,j,k,UEINT) + pyav*duy) &
                        - hdtdz*(fzy(i,j  ,k+1,UEINT) - fzy(i,j,k,UEINT) + pzav*duz) &
                        + hdt*srcQ(i,j,k,QREINT)
                   
                   ! if we are still negative, then we need to reset
                   if (qpo(i,j,k,QREINT) .le. ZERO) then
                      eos_state % rho = qpo(i,j,k,QRHO)
                      eos_state % T = small_temp
                      eos_state % xn(:) = qpo(i,j,k,QFS:QFS-1+nspec)
                      
                      call eos(eos_input_rt, eos_state)
                      
                      qpo(i,j,k,QREINT) = qpo(i,j,k,QRHO)*eos_state % e
                      qpo(i,j,k,QPRES) = eos_state % p
                   endif
                endif
             endif

             if (ppm_predict_gammae == 0) then

                ! Optionally, use the EOS to calculate the pressure.

                if (transverse_use_eos .eq. 1) then
                   eos_state % rho = qpo(i,j,k,QRHO)
                   eos_state % e   = qpo(i,j,k,QREINT) / qpo(i,j,k,QRHO)
                   eos_state % T   = small_temp
                   eos_state % xn  = qpo(i,j,k,QFS:QFS+nspec-1)

                   call eos(eos_input_re, eos_state)

                   pnewr = eos_state % p
                   qpo(i,j,k,QPRES ) = pnewr
                   qpo(i,j,k,QREINT) = eos_state % e * eos_state % rho
                else
                   ! add the transverse term to the p evolution eq here
                   pnewr = qp(i,j,k,QPRES) - pynew - pznew
                   qpo(i,j,k,QPRES) = pnewr + hdt*srcQ(i,j,k,QPRES)
                endif

                qpo(i,j,k,QPRES) = max(qpo(i,j,k,QPRES),small_pres)

             else

                ! Update gammae with its transverse terms
                qpo(i,j,k,QGAME) = qp(i,j,k,QGAME) + geynew + geznew

                ! and compute the p edge state from this and (rho e)
                qpo(i,j,k,QPRES) = qpo(i,j,k,QREINT)*(qpo(i,j,k,QGAME)-ONE)
                qpo(i,j,k,QPRES) = max(qpo(i,j,k,QPRES), small_pres)

             end if

          endif

          !-------------------------------------------------------------------
          ! qxmo state
          !-------------------------------------------------------------------

          if (i.le.hi(1)) then
             qmo(i+1,j,k,QRHO   ) = rrnewl        + hdt*srcQ(i,j,k,QRHO)
             qmo(i+1,j,k,QU     ) = runewl/rrnewl + hdt*srcQ(i,j,k,QU)
             qmo(i+1,j,k,QV     ) = rvnewl/rrnewl + hdt*srcQ(i,j,k,QV)
             qmo(i+1,j,k,QW     ) = rwnewl/rrnewl + hdt*srcQ(i,j,k,QW)

             ! note: we run the risk of (rho e) being negative here
             qmo(i+1,j,k,QREINT ) = renewl - rhoekenl + hdt*srcQ(i,j,k,QREINT)

             if (transverse_reset_rhoe == 1) then
                ! If it is negative, reset the internal energy by using the discretized
                ! expression for updating (rho e).
                
                if (qmo(i+1,j,k,QREINT) .le. ZERO) then
                   qmo(i+1,j,k,QREINT ) = qm(i+1,j,k,QREINT) &
                        - hdtdy*(fyz(i,j+1,k,UEINT) - fyz(i,j,k,UEINT) + pyav*duy) &
                        - hdtdz*(fzy(i,j  ,k+1,UEINT) - fzy(i,j,k,UEINT) + pzav*duz) &
                        + hdt*srcQ(i,j,k,QREINT)

                   ! if we are still negative, then we need to reset
                   if (qmo(i+1,j,k,QREINT) < ZERO) then
                      eos_state % rho = qmo(i+1,j,k,QRHO)
                      eos_state % T = small_temp
                      eos_state % xn(:) = qmo(i+1,j,k,QFS:QFS-1+nspec)
                      
                      call eos(eos_input_rt, eos_state)

                      qmo(i+1,j,k,QREINT) = qmo(i+1,j,k,QRHO)*eos_state % e
                      qmo(i+1,j,k,QPRES) = eos_state % p
                   endif
                endif
             endif

             if (ppm_predict_gammae == 0) then

                ! Optionally, use the EOS To calculate the pressure.

                if (transverse_use_eos .eq. 1) then
                   eos_state % rho = qmo(i+1,j,k,QRHO)
                   eos_state % e   = qmo(i+1,j,k,QREINT) / qmo(i+1,j,k,QRHO)
                   eos_state % T   = small_temp
                   eos_state % xn  = qmo(i+1,j,k,QFS:QFS+nspec-1)
                
                   call eos(eos_input_re, eos_state)
                   
                   pnewl = eos_state % p
                   qmo(i+1,j,k,QPRES ) = pnewl
                   qmo(i+1,j,k,QREINT) = eos_state % e * eos_state % rho
                else
                   ! add the transverse term to the p evolution eq here
                   pnewl = qm(i+1,j,k,QPRES) - pynew - pznew
                   qmo(i+1,j,k,QPRES  ) = pnewl + hdt*srcQ(i,j,k,QPRES)
                endif

                qmo(i+1,j,k,QPRES  ) = max(qmo(i+1,j,k,QPRES),small_pres)

             else
                
                ! Update gammae with its transverse terms
                qmo(i+1,j,k,QGAME) = qm(i+1,j,k,QGAME) + geynew + geznew

                ! and compute the p edge state from this and (rho e)
                qmo(i+1,j,k,QPRES) = qmo(i+1,j,k,QREINT)*(qmo(i+1,j,k,QGAME)-ONE)
                qmo(i+1,j,k,QPRES) = max(qmo(i+1,j,k,QPRES), small_pres)

             end if

          endif

       enddo
       enddo
    enddo
    
    ! if ppm_trace_grav == 1, then we already added the piecewise parabolic traced
    ! gravity to the normal edge states
    if ((do_grav .eq. 1) .and. (ppm_trace_grav == 0 .or. ppm_type == 0)) then
       do    k = lo(3), hi(3)
       do    j = lo(2), hi(2) 
          do i = lo(1), hi(1)+1
             qpo(i,j,k,QU    ) = qpo(i,j,k,QU    ) + hdt*grav(i,j,k,1)
             qpo(i,j,k,QV    ) = qpo(i,j,k,QV    ) + hdt*grav(i,j,k,2)
             qpo(i,j,k,QW    ) = qpo(i,j,k,QW    ) + hdt*grav(i,j,k,3)
          end do
          do i = lo(1)-1, hi(1)
             qmo(i+1,j,k,QU     ) = qmo(i+1,j,k,QU     ) + hdt*grav(i,j,k,1)
             qmo(i+1,j,k,QV     ) = qmo(i+1,j,k,QV     ) + hdt*grav(i,j,k,2)
             qmo(i+1,j,k,QW     ) = qmo(i+1,j,k,QW     ) + hdt*grav(i,j,k,3)
          enddo
       enddo
       enddo
    endif

    ! if ppm_trace_rot == 1, then we already added the piecewise parabolic traced
    ! rotation to the normal edge states
    if ((do_rotation .eq. 1) .and. (ppm_trace_rot == 0 .or. ppm_type == 0)) then
       do    k = lo(3), hi(3)
       do    j = lo(2), hi(2) 
          do i = lo(1), hi(1)+1
             qpo(i,j,k,QU    ) = qpo(i,j,k,QU    ) + hdt*rot(i,j,k,1)
             qpo(i,j,k,QV    ) = qpo(i,j,k,QV    ) + hdt*rot(i,j,k,2)
             qpo(i,j,k,QW    ) = qpo(i,j,k,QW    ) + hdt*rot(i,j,k,3)
          end do
          do i = lo(1)-1, hi(1)
             qmo(i+1,j,k,QU     ) = qmo(i+1,j,k,QU     ) + hdt*rot(i,j,k,1)
             qmo(i+1,j,k,QV     ) = qmo(i+1,j,k,QV     ) + hdt*rot(i,j,k,2)
             qmo(i+1,j,k,QW     ) = qmo(i+1,j,k,QW     ) + hdt*rot(i,j,k,3)
          enddo
       enddo
       enddo
    endif
    
  end subroutine transyz

end module transverse_module
