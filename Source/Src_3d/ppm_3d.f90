module ppm_module

  implicit none

  private

  public ppm

contains
  !
  ! characteristics based on u
  !
  subroutine ppm(s, slo, shi, &
                 u,cspd, ulo, uhi, &
                 flatn, flo, fhi, &
                 Ip,Im, Ilo, Ihi, &
                 lo,hi,dx,dy,dz,dt, &
                 force_type_in)

    use meth_params_module, only : ppm_type

    implicit none

    integer, intent(in) :: lo(3), hi(3), slo(3), shi(3), ulo(3), uhi(3), flo(3), fhi(3), Ilo(3), Ihi(3)
    double precision, intent(in) :: dx, dy, dz, dt
    double precision,intent(in )::    s(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
    double precision,intent(in )::    u(ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3),3)
    double precision,intent(in ):: cspd(ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3))
    double precision,intent(in )::flatn(flo(1):fhi(1),flo(2):fhi(2),flo(3):fhi(3))
    double precision,intent(out)::   Ip(Ilo(1):Ihi(1),Ilo(2):Ihi(2),Ilo(3):Ihi(3),3,3)
    double precision,intent(out)::   Im(Ilo(1):Ihi(1),Ilo(2):Ihi(2),Ilo(3):Ihi(3),3,3)
    integer, intent(in), optional :: force_type_in
    
    integer :: ppm_type_to_use

    ppm_type_to_use = ppm_type
    if (present(force_type_in)) ppm_type_to_use = force_type_in

    if (ppm_type_to_use == 1) then

       call ppm_type1(s, slo, shi, &
                      u,cspd, ulo, uhi, &
                      flatn, flo, fhi, &
                      Ip,Im, Ilo, Ihi, &
                      lo,hi,dx,dy,dz,dt)

    else if (ppm_type_to_use == 2) then

       call ppm_type2(s, slo, shi, &
                      u,cspd, ulo, uhi, &
                      flatn, flo, fhi, &
                      Ip,Im, Ilo, Ihi, &
                      lo,hi,dx,dy,dz,dt)

    end if

  end subroutine ppm

  ! :::
  ! ::: ----------------------------------------------------------------
  ! :::
  
  subroutine ppm_type1(s, slo, shi, &
                       u,cspd, ulo, uhi, &
                       flatn, flo, fhi, &
                       Ip,Im, Ilo, Ihi, &
                       lo,hi,dx,dy,dz,dt)

    use meth_params_module, only : ppm_type, ppm_flatten_before_integrals
    use bl_constants_module
  
    implicit none

    integer, intent(in) :: lo(3), hi(3), slo(3), shi(3), ulo(3), uhi(3), flo(3), fhi(3), Ilo(3), Ihi(3)
    double precision, intent(in) :: dx, dy, dz, dt
    double precision,intent(in )::    s(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
    double precision,intent(in )::    u(ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3),3)
    double precision,intent(in ):: cspd(ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3))
    double precision,intent(in )::flatn(flo(1):fhi(1),flo(2):fhi(2),flo(3):fhi(3))
    double precision,intent(out)::   Ip(Ilo(1):Ihi(1),Ilo(2):Ihi(2),Ilo(3):Ihi(3),3,3)
    double precision,intent(out)::   Im(Ilo(1):Ihi(1),Ilo(2):Ihi(2),Ilo(3):Ihi(3),3,3)

    ! local
    integer i,j,k,n

    double precision dsl, dsr, dsc
    double precision sigma, s6

    double precision :: sp, sm, dtdx, dtdy, dtdz
    
    double precision :: dsvlx(lo(1)-2:hi(1)+2)
    double precision :: sedgex(lo(1)-1:hi(1)+2)
    
    double precision, allocatable :: dsvly(:,:), sedgey(:,:)
    double precision, allocatable :: dsvlz(:,:,:), sedgez(:,:,:)
    
    if (ppm_type .ne. 1) &
         call bl_error("Should have ppm_type = 1 in ppm_type1")

    do n = 1, 3
       if (slo(n) .gt. lo(n)-3 .or. shi(n) .lt. ihi(n)+3) then
          print *, 'Direction: ', n
          print *, 'Bounds of array: ', slo(n), shi(n)
          print *, 'Bounds of  loop: ',  lo(n),  hi(n)
          call bl_error("Need more ghost cells on array in ppm_type1")
       end if
    end do

    dtdx = dt/dx
    dtdy = dt/dy
    dtdz = dt/dz    
    
    allocate(dsvly (lo(1)-1:hi(1)+1,lo(2)-2:hi(2)+2))
    allocate(sedgey(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+2))

    allocate(dsvlz (lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-2:hi(3)+2))
    allocate(sedgez(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+2))

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! x-direction
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    do    k = lo(3)-1, hi(3)+1
       do j = lo(2)-1, hi(2)+1

          ! compute s at x-edges

          ! compute van Leer slopes in x-direction

          do i = lo(1)-2, hi(1)+2
             dsc = HALF * (s(i+1,j,k) - s(i-1,j,k))
             dsl = TWO  * (s(i  ,j,k) - s(i-1,j,k))
             dsr = TWO  * (s(i+1,j,k) - s(i  ,j,k))
             if (dsl*dsr .gt. ZERO) then
                dsvlx(i) = sign(ONE,dsc)*min(abs(dsc),abs(dsl),abs(dsr))
             else
                dsvlx(i) = ZERO
             end if
          end do

          ! interpolate s to x-edges
          do i = lo(1)-1, hi(1)+2
             sedgex(i) = HALF*(s(i,j,k)+s(i-1,j,k)) - SIXTH*(dsvlx(i)-dsvlx(i-1))
             ! make sure sedge lies in between adjacent cell-centered values
             sedgex(i) = max(sedgex(i),min(s(i,j,k),s(i-1,j,k)))
             sedgex(i) = min(sedgex(i),max(s(i,j,k),s(i-1,j,k)))
          end do

          do i = lo(1)-1, hi(1)+1

             ! copy sedge into sp and sm
             sp = sedgex(i+1)
             sm = sedgex(i)

             if (ppm_flatten_before_integrals == 1) then
                ! flatten the parabola BEFORE doing the other                     
                ! monotonization -- this is the method that Flash does       
                sm = flatn(i,j,k)*sm + (ONE-flatn(i,j,k))*s(i,j,k)
                sp = flatn(i,j,k)*sp + (ONE-flatn(i,j,k))*s(i,j,k)
             endif

             ! modify using quadratic limiters -- note this version of the limiting comes
             ! from Colella and Sekora (2008), not the original PPM paper.
             if ((sp-s(i,j,k))*(s(i,j,k)-sm) .le. ZERO) then
                sp = s(i,j,k)
                sm = s(i,j,k)

             else if (abs(sp-s(i,j,k)) .ge. TWO*abs(sm-s(i,j,k))) then
                !else if (-(sp-sm)**2/SIX > &
                !     (sp - sm)*(s(i,j,k) - HALF*(sm + sp))) then
                sp = THREE*s(i,j,k) - TWO*sm

             else if (abs(sm-s(i,j,k)) .ge. TWO*abs(sp-s(i,j,k))) then
                !else if ((sp-sm)*(s(i,j,k) - HALF*(sm + sp)) > &
                !     (sp - sm)**2/SIX) then
                sm = THREE*s(i,j,k) - TWO*sp
             end if
             
             if (ppm_flatten_before_integrals == 2) then
                ! flatten the parabola AFTER doing the monotonization --
                ! this is the method that Miller & Colella do
                sm = flatn(i,j,k)*sm + (ONE-flatn(i,j,k))*s(i,j,k)
                sp = flatn(i,j,k)*sp + (ONE-flatn(i,j,k))*s(i,j,k)
             endif

             ! compute x-component of Ip and Im
             s6 = SIX*s(i,j,k) - THREE*(sm+sp)

             ! Ip/m is the integral under the parabola for the extent
             ! that a wave can travel over a timestep
             !
             ! Ip integrates to the right edge of a cell
             ! Im integrates to the left edge of a cell
             
             ! u-c wave
             sigma = abs(u(i,j,k,1)-cspd(i,j,k))*dtdx

             if (u(i,j,k,1)-cspd(i,j,k) <= ZERO) then
                Ip(i,j,k,1,1) = sp
             else
                Ip(i,j,k,1,1) = sp - &
                     HALF*sigma*(sp-sm-(ONE-TWO3RD*sigma)*s6)
             endif

             if (u(i,j,k,1)-cspd(i,j,k) >= ZERO) then
                Im(i,j,k,1,1) = sm 
             else
                Im(i,j,k,1,1) = sm + &
                     HALF*sigma*(sp-sm+(ONE-TWO3RD*sigma)*s6)
             endif
             
             ! u wave
             sigma = abs(u(i,j,k,1))*dtdx
             
             if (u(i,j,k,1) <= ZERO) then
                Ip(i,j,k,1,2) = sp 
             else
                Ip(i,j,k,1,2) = sp - &
                     HALF*sigma*(sp-sm-(ONE-TWO3RD*sigma)*s6)
             endif
             
             if (u(i,j,k,1) >= ZERO) then
                Im(i,j,k,1,2) = sm 
             else
                Im(i,j,k,1,2) = sm + &
                     HALF*sigma*(sp-sm+(ONE-TWO3RD*sigma)*s6)
             endif

             ! u+c wave
             sigma = abs(u(i,j,k,1)+cspd(i,j,k))*dtdx
             
             if (u(i,j,k,1)+cspd(i,j,k) <= ZERO) then
                Ip(i,j,k,1,3) = sp 
             else
                Ip(i,j,k,1,3) = sp - &
                     HALF*sigma*(sp-sm-(ONE-TWO3RD*sigma)*s6)
             endif
             
             if (u(i,j,k,1)+cspd(i,j,k) >= ZERO) then
                Im(i,j,k,1,3) = sm 
             else
                Im(i,j,k,1,3) = sm + &
                     HALF*sigma*(sp-sm+(ONE-TWO3RD*sigma)*s6)
             endif

          end do
       end do
    end do

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! y-direction
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    do k = lo(3)-1, hi(3)+1
    
       ! compute s at y-edges

       ! compute van Leer slopes in y-direction
       do    j = lo(2)-2, hi(2)+2
          do i = lo(1)-1, hi(1)+1
             dsc = HALF * (s(i,j+1,k) - s(i,j-1,k))
             dsl = TWO  * (s(i,j  ,k) - s(i,j-1,k))
             dsr = TWO  * (s(i,j+1,k) - s(i,j  ,k))
             if (dsl*dsr .gt. ZERO) then
                dsvly(i,j) = sign(ONE,dsc)*min(abs(dsc),abs(dsl),abs(dsr))
             else
                dsvly(i,j) = ZERO
             end if
          end do
       end do

       ! interpolate s to y-edges
       do    j = lo(2)-1, hi(2)+2
          do i = lo(1)-1, hi(1)+1
             sedgey(i,j) = HALF*(s(i,j,k)+s(i,j-1,k)) - SIXTH*(dsvly(i,j)-dsvly(i,j-1))
             ! make sure sedgey lies in between adjacent cell-centered values
             sedgey(i,j) = max(sedgey(i,j),min(s(i,j,k),s(i,j-1,k)))
             sedgey(i,j) = min(sedgey(i,j),max(s(i,j,k),s(i,j-1,k)))
          end do
       end do

       do    j = lo(2)-1, hi(2)+1
          do i = lo(1)-1, hi(1)+1

             ! copy sedgey into sp and sm
             sp = sedgey(i,j+1)
             sm = sedgey(i,j  )

             if (ppm_flatten_before_integrals == 1) then
                ! flatten the parabola BEFORE doing the other                     
                ! monotonization -- this is the method that Flash does       
                sm = flatn(i,j,k)*sm + (ONE-flatn(i,j,k))*s(i,j,k)
                sp = flatn(i,j,k)*sp + (ONE-flatn(i,j,k))*s(i,j,k)
             endif

             ! modify using quadratic limiters
             if ((sp-s(i,j,k))*(s(i,j,k)-sm) .le. ZERO) then
                sp = s(i,j,k)
                sm = s(i,j,k)
                
             else if (abs(sp-s(i,j,k)) .ge. TWO*abs(sm-s(i,j,k))) then
                !else if (-(sp-sm)**2/SIX > &
                !     (sp - sm)*(s(i,j,k) - HALF*(sm + sp))) then
                sp = THREE*s(i,j,k) - TWO*sm
                
             else if (abs(sm-s(i,j,k)) .ge. TWO*abs(sp-s(i,j,k))) then
                !else if ((sp-sm)*(s(i,j,k) - HALF*(sm + sp)) > &
                !     (sp - sm)**2/SIX) then
                sm = THREE*s(i,j,k) - TWO*sp
             end if
             
             if (ppm_flatten_before_integrals == 2) then
                ! flatten the parabola AFTER doing the monotonization --
                ! this is the method that Miller & Colella do
                sm = flatn(i,j,k)*sm + (ONE-flatn(i,j,k))*s(i,j,k)
                sp = flatn(i,j,k)*sp + (ONE-flatn(i,j,k))*s(i,j,k)
             endif
             
             ! compute y-component of Ip and Im
             s6 = SIX*s(i,j,k) - THREE*(sm+sp)
             
             ! v-c wave
             sigma = abs(u(i,j,k,2)-cspd(i,j,k))*dtdy
             
             if (u(i,j,k,2)-cspd(i,j,k) <= ZERO) then
                Ip(i,j,k,2,1) = sp
             else
                Ip(i,j,k,2,1) = sp - &
                     HALF*sigma*(sp-sm-(ONE-TWO3RD*sigma)*s6)
             endif
             
             if (u(i,j,k,2)-cspd(i,j,k) >= ZERO) then
                Im(i,j,k,2,1) = sm 
             else
                Im(i,j,k,2,1) = sm + &
                     HALF*sigma*(sp-sm+(ONE-TWO3RD*sigma)*s6)
             endif
             
             ! v wave
             sigma = abs(u(i,j,k,2))*dtdy
             
             if (u(i,j,k,2) <= ZERO) then
                Ip(i,j,k,2,2) = sp 
             else
                Ip(i,j,k,2,2) = sp - &
                     HALF*sigma*(sp-sm-(ONE-TWO3RD*sigma)*s6)
             endif

             if (u(i,j,k,2) >= ZERO) then
                Im(i,j,k,2,2) = sm 
             else
                Im(i,j,k,2,2) = sm + &
                     HALF*sigma*(sp-sm+(ONE-TWO3RD*sigma)*s6)
             endif
             
             ! v+c wave
             sigma = abs(u(i,j,k,2)+cspd(i,j,k))*dtdy
             
             if (u(i,j,k,2)+cspd(i,j,k) <= ZERO) then
                Ip(i,j,k,2,3) = sp 
             else
                Ip(i,j,k,2,3) = sp - &
                     HALF*sigma*(sp-sm-(ONE-TWO3RD*sigma)*s6)
             endif
             
             if (u(i,j,k,2)+cspd(i,j,k) >= ZERO) then
                Im(i,j,k,2,3) = sm 
             else
                Im(i,j,k,2,3) = sm + &
                     HALF*sigma*(sp-sm+(ONE-TWO3RD*sigma)*s6)
             endif
             
          end do
       end do
    end do

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! z-direction
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! compute s at z-edges

    ! compute van Leer slopes in z-direction

    do       k = lo(3)-2, hi(3)+2
       do    j = lo(2)-1, hi(2)+1
          do i = lo(1)-1, hi(1)+1

             dsc = HALF * (s(i,j,k+1) - s(i,j,k-1))
             dsl = TWO  * (s(i,j,k  ) - s(i,j,k-1))
             dsr = TWO  * (s(i,j,k+1) - s(i,j,k  ))
             if (dsl*dsr .gt. ZERO) then
                dsvlz(i,j,k) = sign(ONE,dsc)*min(abs(dsc),abs(dsl),abs(dsr))
             else
                dsvlz(i,j,k) = ZERO
             end if

          end do
       end do
    end do

    ! interpolate s to z-edges
    do       k = lo(3)-1, hi(3)+2
       do    j = lo(2)-1, hi(2)+1
          do i = lo(1)-1, hi(1)+1
             sedgez(i,j,k) = HALF*(s(i,j,k)+s(i,j,k-1)) - SIXTH*(dsvlz(i,j,k)-dsvlz(i,j,k-1))
             ! make sure sedgez lies in between adjacent cell-centered values
             sedgez(i,j,k) = max(sedgez(i,j,k),min(s(i,j,k),s(i,j,k-1)))
             sedgez(i,j,k) = min(sedgez(i,j,k),max(s(i,j,k),s(i,j,k-1)))
          end do
       end do
    end do
             

    do       k = lo(3)-1, hi(3)+2
       do    j = lo(2)-1, hi(2)+1
          do i = lo(1)-1, hi(1)+1
          
             ! copy sedgez into sp and sm
             sp = sedgez(i,j,k+1)
             sp = sedgez(i,j,k)

             if (ppm_flatten_before_integrals == 1) then
                ! flatten the parabola BEFORE doing the other                     
                ! monotonization -- this is the method that Flash does       
                sm = flatn(i,j,k)*sm + (ONE-flatn(i,j,k))*s(i,j,k)
                sp = flatn(i,j,k)*sp + (ONE-flatn(i,j,k))*s(i,j,k)
             endif

             ! modify using quadratic limiters
             if ((sp-s(i,j,k))*(s(i,j,k)-sm) .le. ZERO) then
                sp = s(i,j,k)
                sm = s(i,j,k)
                
             else if (abs(sp-s(i,j,k)) .ge. TWO*abs(sm-s(i,j,k))) then
                !else if (-(sp-sm)**2/SIX > &
                !     (sp - sm)*(s(i,j,k) - HALF*(sm + sp))) then
                sp = THREE*s(i,j,k) - TWO*sm
                
             else if (abs(sm-s(i,j,k)) .ge. TWO*abs(sp-s(i,j,k))) then
                !else if ((sp-sm)*(s(i,j,k) - HALF*(sm + sp)) > &
                !     (sp - sm)**2/SIX) then
                sm = THREE*s(i,j,k) - TWO*sp
             end if
             
             if (ppm_flatten_before_integrals == 2) then
                ! flatten the parabola AFTER doing the monotonization --
                ! this is the method that Miller & Colella do
                sm = flatn(i,j,k)*sm + (ONE-flatn(i,j,k))*s(i,j,k)
                sp = flatn(i,j,k)*sp + (ONE-flatn(i,j,k))*s(i,j,k)
             endif
             
             ! compute z-component of Ip and Im
             s6 = SIX*s(i,j,k) - THREE*(sm+sp)
             
             ! w-c wave
             sigma = abs(u(i,j,k,3)-cspd(i,j,k))*dtdz
             
             if (u(i,j,k,3)-cspd(i,j,k) <= ZERO) then
                Ip(i,j,k,3,1) = sp 
             else
                Ip(i,j,k,3,1) = sp - &
                     HALF*sigma*(sp-sm-(ONE-TWO3RD*sigma)*s6)
             endif
             
             if (u(i,j,k,3)-cspd(i,j,k) >= ZERO) then
                Im(i,j,k,3,1) = sm 
             else
                Im(i,j,k,3,1) = sm + &
                     HALF*sigma*(sp-sm+(ONE-TWO3RD*sigma)*s6)
             endif
             
             ! w wave
             sigma = abs(u(i,j,k,3))*dtdz
             
             if (u(i,j,k,3) <= ZERO) then
                Ip(i,j,k,3,2) = sp 
             else
                Ip(i,j,k,3,2) = sp - &
                     HALF*sigma*(sp-sm-(ONE-TWO3RD*sigma)*s6)
             endif
             
             if (u(i,j,k,3) >= ZERO) then
                Im(i,j,k,3,2) = sm 
             else
                Im(i,j,k,3,2) = sm + &
                     HALF*sigma*(sp-sm+(ONE-TWO3RD*sigma)*s6)
             endif
             
             ! w+c wave
             sigma = abs(u(i,j,k,3)+cspd(i,j,k))*dtdz
             
             if (u(i,j,k,3)+cspd(i,j,k) <= ZERO) then
                Ip(i,j,k,3,3) = sp 
             else
                Ip(i,j,k,3,3) = sp - &
                     HALF*sigma*(sp-sm-(ONE-TWO3RD*sigma)*s6)
             endif
             
             if (u(i,j,k,3)+cspd(i,j,k) >= ZERO) then
                Im(i,j,k,3,3) = sm 
             else
                Im(i,j,k,3,3) = sm + &
                     HALF*sigma*(sp-sm+(ONE-TWO3RD*sigma)*s6)
             endif

          end do
       end do
    end do

    deallocate(dsvly,sedgey,dsvlz,sedgez)

  end subroutine ppm_type1

  ! :::
  ! ::: ----------------------------------------------------------------
  ! :::

  subroutine ppm_type2(s, slo, shi, &
                       u,cspd, ulo, uhi, &
                       flatn, flo, fhi, &
                       Ip,Im, Ilo, Ihi, &
                       lo,hi,dx,dy,dz,dt)

    use meth_params_module, only : ppm_type, ppm_flatten_before_integrals
    use bl_constants_module
  
    implicit none

    integer, intent(in) :: lo(3), hi(3), slo(3), shi(3), ulo(3), uhi(3), flo(3), fhi(3), Ilo(3), Ihi(3)
    double precision, intent(in) :: dx, dy, dz, dt
    double precision,intent(in )::    s(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
    double precision,intent(in )::    u(ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3),3)
    double precision,intent(in ):: cspd(ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3))
    double precision,intent(in )::flatn(flo(1):fhi(1),flo(2):fhi(2),flo(3):fhi(3))
    double precision,intent(out)::   Ip(Ilo(1):Ihi(1),Ilo(2):Ihi(2),Ilo(3):Ihi(3),3,3)
    double precision,intent(out)::   Im(Ilo(1):Ihi(1),Ilo(2):Ihi(2),Ilo(3):Ihi(3),3,3)

    ! local
    integer i,j,k,n
    logical extremum, bigp, bigm

    double precision D2, D2C, D2L, D2R, D2LIM, alphap, alpham
    double precision sgn, sigma, s6
    double precision dafacem, dafacep, dabarm, dabarp, dafacemin, dabarmin
    double precision dachkm, dachkp
    double precision amax, delam, delap

    double precision :: sp, sm, dtdx, dtdy, dtdz
    
    double precision :: sedgex(lo(1)-2:hi(1)+3)    
    double precision, allocatable :: sedgey(:,:)
    double precision, allocatable :: sedgez(:,:,:)

    ! constant used in Colella 2008
    double precision, parameter :: C = 1.25d0

    if (ppm_type .ne. 2) &
         call bl_error("Should have ppm_type = 2 in ppm_type2")

        do n = 1, 3
       if (slo(n) .gt. lo(n)-3 .or. shi(n) .lt. ihi(n)+3) then
          print *, 'Direction: ', n
          print *, 'Bounds of array: ', slo(n), shi(n)
          print *, 'Bounds of  loop: ',  lo(n),  hi(n)
          call bl_error("Need more ghost cells on array in ppm_type1")
       end if
    end do

    dtdx = dt/dx
    dtdy = dt/dy
    dtdz = dt/dz    

    allocate(sedgey(lo(1)-1:hi(1)+1,lo(2)-2:hi(2)+3))
    allocate(sedgez(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-2:hi(3)+3))

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! x-direction
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    do    k = lo(3)-1, hi(3)+1
       do j = lo(2)-1, hi(2)+1
          
          ! compute s at x-edges
          do i = lo(1)-2, hi(1)+3          
             sedgex(i) = SEVEN12TH*(s(i-1,j,k)+s(i  ,j,k)) - TWELFTH*(s(i-2,j,k)+s(i+1,j,k))
             !
             ! limit sedge
             !
             if ((sedgex(i)-s(i-1,j,k))*(s(i,j,k)-sedgex(i)) .lt. ZERO) then
                D2  = THREE*(s(i-1,j,k)-TWO*sedgex(i)+s(i,j,k))
                D2L = s(i-2,j,k)-TWO*s(i-1,j,k)+s(i,j,k)
                D2R = s(i-1,j,k)-TWO*s(i,j,k)+s(i+1,j,k)
                sgn = sign(ONE,D2)
                D2LIM = sgn*max(min(C*sgn*D2L,C*sgn*D2R,sgn*D2),ZERO)
                sedgex(i) = HALF*(s(i-1,j,k)+s(i,j,k)) - SIXTH*D2LIM
             end if
          end do
          !
          ! Use Colella 2008 limiters.
          !
          ! This is a new version of the algorithm to eliminate sensitivity to roundoff.
          !
          do i = lo(1)-1, hi(1)+1

             alphap   = sedgex(i+1)-s(i,j,k)
             alpham   = sedgex(i  )-s(i,j,k)
             bigp     = abs(alphap).gt.TWO*abs(alpham)
             bigm     = abs(alpham).gt.TWO*abs(alphap)
             extremum = .false.
             
             if (alpham*alphap .ge. ZERO) then
                extremum = .true.
             else if (bigp .or. bigm) then
                !
                ! Possible extremum. We look at cell centered values and face
                ! centered values for a change in sign in the differences adjacent to
                ! the cell. We use the pair of differences whose minimum magnitude is the
                ! largest, and thus least susceptible to sensitivity to roundoff.
                !
                dafacem   = sedgex(i) - sedgex(i-1)
                dafacep   = sedgex(i+2) - sedgex(i+1)
                dabarm    = s(i,j,k) - s(i-1,j,k)
                dabarp    = s(i+1,j,k) - s(i,j,k)
                dafacemin = min(abs(dafacem),abs(dafacep))
                dabarmin  = min(abs(dabarm),abs(dabarp))
                if (dafacemin.ge.dabarmin) then
                   dachkm = dafacem
                   dachkp = dafacep
                else
                   dachkm = dabarm
                   dachkp = dabarp
                endif
                extremum = (dachkm*dachkp .le. ZERO)
             end if
             
             if (extremum) then
                D2     = SIX*(alpham + alphap)
                D2L    = s(i-2,j,k)-TWO*s(i-1,j,k)+s(i,j,k)
                D2R    = s(i,j,k)-TWO*s(i+1,j,k)+s(i+2,j,k)
                D2C    = s(i-1,j,k)-TWO*s(i,j,k)+s(i+1,j,k)
                sgn    = sign(ONE,D2)
                D2LIM  = max(min(sgn*D2,C*sgn*D2L,C*sgn*D2R,C*sgn*D2C),ZERO)
                alpham = alpham*D2LIM/max(abs(D2),1.d-10)
                alphap = alphap*D2LIM/max(abs(D2),1.d-10)
             else
                if (bigp) then
                   sgn   = sign(ONE,alpham)
                   amax  = -alphap**2 / (4*(alpham + alphap))
                   delam = s(i-1,j,k) - s(i,j,k)
                   if (sgn*amax .ge. sgn*delam) then
                      if (sgn*(delam - alpham).ge.1.d-10) then
                         alphap = (-TWO*delam - TWO*sgn*sqrt(delam**2 - delam*alpham))
                      else 
                         alphap = -TWO*alpham
                      endif
                   endif
                end if
                if (bigm) then
                   sgn   = sign(ONE,alphap)
                   amax  = -alpham**2 / (4*(alpham + alphap))
                   delap = s(i+1,j,k) - s(i,j,k)
                   if (sgn*amax .ge. sgn*delap) then
                      if (sgn*(delap - alphap).ge.1.d-10) then
                         alpham = (-TWO*delap - TWO*sgn*sqrt(delap**2 - delap*alphap))
                      else
                         alpham = -TWO*alphap
                      endif
                   endif
                end if
             end if

             sm = s(i,j,k) + alpham
             sp = s(i,j,k) + alphap

             if (ppm_flatten_before_integrals > 0) then
                ! flatten the parabola AFTER doing the monotonization 
                sm = flatn(i,j,k)*sm + (ONE-flatn(i,j,k))*s(i,j,k)
                sp = flatn(i,j,k)*sp + (ONE-flatn(i,j,k))*s(i,j,k)
             endif

             !
             ! Compute x-component of Ip and Im.
             !
             s6    = SIX*s(i,j,k) - THREE*(sm+sp)

             ! u-c wave
             sigma = abs(u(i,j,k,1)-cspd(i,j,k))*dt/dx

             if (u(i,j,k,1)-cspd(i,j,k) <= ZERO) then
                Ip(i,j,k,1,1) = sp
             else
                Ip(i,j,k,1,1) = sp - &
                     HALF*sigma*(sp-sm-(ONE-TWO3RD*sigma)*s6)
             endif
             
             if (u(i,j,k,1)-cspd(i,j,k) >= ZERO) then
                Im(i,j,k,1,1) = sm
             else
                Im(i,j,k,1,1) = sm + &
                     HALF*sigma*(sp-sm+(ONE-TWO3RD*sigma)*s6)
             endif
             
             ! u wave
             sigma = abs(u(i,j,k,1))*dt/dx
             
             if (u(i,j,k,1) <= ZERO) then
                Ip(i,j,k,1,2) = sp
             else
                Ip(i,j,k,1,2) = sp - &
                     HALF*sigma*(sp-sm-(ONE-TWO3RD*sigma)*s6)
             endif
             
             if (u(i,j,k,1) >= ZERO) then
                Im(i,j,k,1,2) = sm
             else
                Im(i,j,k,1,2) = sm + &
                     HALF*sigma*(sp-sm+(ONE-TWO3RD*sigma)*s6)
             endif
             
             ! u+c wave
             sigma = abs(u(i,j,k,1)+cspd(i,j,k))*dt/dx
             
             if (u(i,j,k,1)+cspd(i,j,k) <= ZERO) then
                Ip(i,j,k,1,3) = sp 
             else
                Ip(i,j,k,1,3) = sp - &
                     HALF*sigma*(sp-sm-(ONE-TWO3RD*sigma)*s6)
             endif
             
             if (u(i,j,k,1)+cspd(i,j,k) >= ZERO) then
                Im(i,j,k,1,3) = sm 
             else
                Im(i,j,k,1,3) = sm + &
                     HALF*sigma*(sp-sm+(ONE-TWO3RD*sigma)*s6)
             endif
             
          end do
       end do
    end do

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! y-direction
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    do    k = lo(3)-1, hi(3)+1

       ! compute s at y-edges

       ! interpolate s to y-edges
       do    j = lo(2)-2, hi(2)+3
          do i = lo(1)-1, hi(1)+1          
             sedgey(i,j) = SEVEN12TH*(s(i,j-1,k)+s(i,j,k)) - TWELFTH*(s(i,j-2,k)+s(i,j+1,k))
             !
             ! limit sedgey
             !
             if ((sedgey(i,j)-s(i,j-1,k))*(s(i,j,k)-sedgey(i,j)) .lt. ZERO) then
                D2  = THREE*(s(i,j-1,k)-TWO*sedgey(i,j)+s(i,j,k))
                D2L = s(i,j-2,k)-TWO*s(i,j-1,k)+s(i,j,k)
                D2R = s(i,j-1,k)-TWO*s(i,j,k)+s(i,j+1,k)
                sgn = sign(ONE,D2)
                D2LIM = sgn*max(min(C*sgn*D2L,C*sgn*D2R,sgn*D2),ZERO)
                sedgey(i,j) = HALF*(s(i,j-1,k)+s(i,j,k)) - SIXTH*D2LIM
             end if
          end do
       end do
       !
       ! Use Colella 2008 limiters.
       !
       ! This is a new version of the algorithm to eliminate sensitivity to roundoff.
       !
       do    j = lo(2)-1, hi(2)+1
          do i = lo(1)-1, hi(1)+1          

             alphap   = sedgey(i,j+1)-s(i,j,k)
             alpham   = sedgey(i,j  )-s(i,j,k)
             bigp     = abs(alphap).gt.TWO*abs(alpham)
             bigm     = abs(alpham).gt.TWO*abs(alphap)
             extremum = .false.

             if (alpham*alphap .ge. ZERO) then
                extremum = .true.
             else if (bigp .or. bigm) then
                !
                ! Possible extremum. We look at cell centered values and face
                ! centered values for a change in sign in the differences adjacent to
                ! the cell. We use the pair of differences whose minimum magnitude is the
                ! largest, and thus least susceptible to sensitivity to roundoff.
                !
                dafacem   = sedgey(i,j) - sedgey(i,j-1)
                dafacep   = sedgey(i,j+2) - sedgey(i,j+1)
                dabarm    = s(i,j,k) - s(i,j-1,k)
                dabarp    = s(i,j+1,k) - s(i,j,k)
                dafacemin = min(abs(dafacem),abs(dafacep))
                dabarmin  = min(abs(dabarm),abs(dabarp))
                if (dafacemin.ge.dabarmin) then
                   dachkm = dafacem
                   dachkp = dafacep
                else
                   dachkm = dabarm
                   dachkp = dabarp
                endif
                extremum = (dachkm*dachkp .le. ZERO)
             end if
             
             if (extremum) then
                D2     = SIX*(alpham + alphap)
                D2L    = s(i,j-2,k)-TWO*s(i,j-1,k)+s(i,j,k)
                D2R    = s(i,j,k)-TWO*s(i,j+1,k)+s(i,j+2,k)
                D2C    = s(i,j-1,k)-TWO*s(i,j,k)+s(i,j+1,k)
                sgn    = sign(ONE,D2)
                D2LIM  = max(min(sgn*D2,C*sgn*D2L,C*sgn*D2R,C*sgn*D2C),ZERO)
                alpham = alpham*D2LIM/max(abs(D2),1.d-10)
                alphap = alphap*D2LIM/max(abs(D2),1.d-10)
             else
                if (bigp) then
                   sgn   = sign(ONE,alpham)
                   amax  = -alphap**2 / (4*(alpham + alphap))
                   delam = s(i,j-1,k) - s(i,j,k)
                   if (sgn*amax .ge. sgn*delam) then
                      if (sgn*(delam - alpham).ge.1.d-10) then
                         alphap = (-TWO*delam - TWO*sgn*sqrt(delam**2 - delam*alpham))
                      else 
                         alphap = -TWO*alpham
                      endif
                   endif
                end if
                if (bigm) then
                   sgn   = sign(ONE,alphap)
                   amax  = -alpham**2 / (4*(alpham + alphap))
                   delap = s(i,j+1,k) - s(i,j,k)
                   if (sgn*amax .ge. sgn*delap) then
                      if (sgn*(delap - alphap).ge.1.d-10) then
                         alpham = (-TWO*delap - TWO*sgn*sqrt(delap**2 - delap*alphap))
                      else
                         alpham = -TWO*alphap
                      endif
                   endif
                end if
             end if

             sm = s(i,j,k) + alpham
             sp = s(i,j,k) + alphap

             if (ppm_flatten_before_integrals > 0) then
                ! flatten the parabola AFTER doing the monotonization 
                sm = flatn(i,j,k)*sm + (ONE-flatn(i,j,k))*s(i,j,k)
                sp = flatn(i,j,k)*sp + (ONE-flatn(i,j,k))*s(i,j,k)
             endif
             
             !
             ! Compute y-component of Ip and Im.
             !
             s6    = SIX*s(i,j,k) - THREE*(sm+sp)
             
             ! v-c wave
             sigma = abs(u(i,j,k,2)-cspd(i,j,k))*dt/dy
             
             if (u(i,j,k,2)-cspd(i,j,k) <= ZERO) then
                Ip(i,j,k,2,1) = sp 
             else
                Ip(i,j,k,2,1) = sp - &
                     HALF*sigma*(sp-sm-(ONE-TWO3RD*sigma)*s6)
             endif
             
             if (u(i,j,k,2)-cspd(i,j,k) >= ZERO) then
                Im(i,j,k,2,1) = sm 
             else
                Im(i,j,k,2,1) = sm + &
                     HALF*sigma*(sp-sm+(ONE-TWO3RD*sigma)*s6)
             endif
             
             ! v wave
             sigma = abs(u(i,j,k,2))*dt/dy
             
             if (u(i,j,k,2) <= ZERO) then
                Ip(i,j,k,2,2) = sp 
             else
                Ip(i,j,k,2,2) = sp - &
                     HALF*sigma*(sp-sm-(ONE-TWO3RD*sigma)*s6)
             endif
             
             if (u(i,j,k,2) >= ZERO) then
                Im(i,j,k,2,2) = sm 
             else
                Im(i,j,k,2,2) = sm + &
                     HALF*sigma*(sp-sm+(ONE-TWO3RD*sigma)*s6)
             endif
             
             ! v+c wave
             sigma = abs(u(i,j,k,2)+cspd(i,j,k))*dt/dy
             
             if (u(i,j,k,2)+cspd(i,j,k) <= ZERO) then
                Ip(i,j,k,2,3) = sp 
             else
                Ip(i,j,k,2,3) = sp - &
                     HALF*sigma*(sp-sm-(ONE-TWO3RD*sigma)*s6)
             endif
             
             if (u(i,j,k,2)+cspd(i,j,k) >= ZERO) then
                Im(i,j,k,2,3) = sm 
             else
                Im(i,j,k,2,3) = sm + &
                     HALF*sigma*(sp-sm+(ONE-TWO3RD*sigma)*s6)
             endif

          end do
       end do
    end do

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! z-direction
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! compute s at z-edges

    ! interpolate s to z-edges
    do       k = lo(3)-2, hi(3)+3
       do    j = lo(2)-1, hi(2)+1
          do i = lo(1)-1, hi(1)+1
             sedgez(i,j,k) = SEVEN12TH*(s(i,j,k-1)+s(i,j,k)) - TWELFTH*(s(i,j,k-2)+s(i,j,k+1))
             !
             ! limit sedgez
             !
             if ((sedgez(i,j,k)-s(i,j,k-1))*(s(i,j,k)-sedgez(i,j,k)) .lt. ZERO) then
                D2  = THREE*(s(i,j,k-1)-TWO*sedgez(i,j,k)+s(i,j,k))
                D2L = s(i,j,k-2)-TWO*s(i,j,k-1)+s(i,j,k)
                D2R = s(i,j,k-1)-TWO*s(i,j,k)+s(i,j,k+1)
                sgn = sign(ONE,D2)
                D2LIM = sgn*max(min(C*sgn*D2L,C*sgn*D2R,sgn*D2),ZERO)
                sedgez(i,j,k) = HALF*(s(i,j,k-1)+s(i,j,k)) - SIXTH*D2LIM
             end if
          end do
       end do
    end do
    !
    ! Use Colella 2008 limiters.
    !
    ! This is a new version of the algorithm to eliminate sensitivity to roundoff.
    !
    do       k = lo(3)-1, hi(3)+1
       do    j = lo(2)-1, hi(2)+1
          do i = lo(1)-1, hi(1)+1          

             alphap   = sedgez(i,j,k+1)-s(i,j,k)
             alpham   = sedgez(i,j,k  )-s(i,j,k)
             bigp     = abs(alphap).gt.TWO*abs(alpham)
             bigm     = abs(alpham).gt.TWO*abs(alphap)
             extremum = .false.
             
             if (alpham*alphap .ge. ZERO) then
                extremum = .true.
             else if (bigp .or. bigm) then
                !
                ! Possible extremum. We look at cell centered values and face
                ! centered values for a change in sign in the differences adjacent to
                ! the cell. We use the pair of differences whose minimum magnitude is the
                ! largest, and thus least susceptible to sensitivity to roundoff.
                !
                dafacem   = sedgez(i,j,k) - sedgez(i,j,k-1)
                dafacep   = sedgez(i,j,k+2) - sedgez(i,j,k+1)
                dabarm    = s(i,j,k) - s(i,j,k-1)
                dabarp    = s(i,j,k+1) - s(i,j,k)
                dafacemin = min(abs(dafacem),abs(dafacep))
                dabarmin  = min(abs(dabarm),abs(dabarp))
                if (dafacemin.ge.dabarmin) then
                   dachkm = dafacem
                   dachkp = dafacep
                else
                   dachkm = dabarm
                   dachkp = dabarp
                endif
                extremum = (dachkm*dachkp .le. ZERO)
             end if
             
             if (extremum) then
                D2     = SIX*(alpham + alphap)
                D2L    = s(i,j,k-2)-TWO*s(i,j,k-1)+s(i,j,k)
                D2R    = s(i,j,k)-TWO*s(i,j,k+1)+s(i,j,k+2)
                D2C    = s(i,j,k-1)-TWO*s(i,j,k)+s(i,j,k+1)
                sgn    = sign(ONE,D2)
                D2LIM  = max(min(sgn*D2,C*sgn*D2L,C*sgn*D2R,C*sgn*D2C),ZERO)
                alpham = alpham*D2LIM/max(abs(D2),1.d-10)
                alphap = alphap*D2LIM/max(abs(D2),1.d-10)
             else
                if (bigp) then
                   sgn   = sign(ONE,alpham)
                   amax  = -alphap**2 / (4*(alpham + alphap))
                   delam = s(i,j,k-1) - s(i,j,k)
                   if (sgn*amax .ge. sgn*delam) then
                      if (sgn*(delam - alpham).ge.1.d-10) then
                         alphap = (-TWO*delam - TWO*sgn*sqrt(delam**2 - delam*alpham))
                      else 
                         alphap = -TWO*alpham
                      endif
                   endif
                end if
                if (bigm) then
                   sgn   = sign(ONE,alphap)
                   amax  = -alpham**2 / (4*(alpham + alphap))
                   delap = s(i,j,k+1) - s(i,j,k)
                   if (sgn*amax .ge. sgn*delap) then
                      if (sgn*(delap - alphap).ge.1.d-10) then
                         alpham = (-TWO*delap - TWO*sgn*sqrt(delap**2 - delap*alphap))
                      else
                         alpham = -TWO*alphap
                      endif
                   endif
                end if
             end if
             
             sm = s(i,j,k) + alpham
             sp = s(i,j,k) + alphap
             
             if (ppm_flatten_before_integrals > 0) then
                ! flatten the parabola AFTER doing the monotonization (note k = k here)
                sm = flatn(i,j,k)*sm + (ONE-flatn(i,j,k))*s(i,j,k)
                sp = flatn(i,j,k)*sp + (ONE-flatn(i,j,k))*s(i,j,k)
             endif
             
             !
             ! Compute z-component of Ip and Im.
             !
             s6    = SIX*s(i,j,k) - THREE*(sm+sp)
             
             
             ! w-c wave
             sigma = abs(u(i,j,k,3)-cspd(i,j,k))*dt/dz
             
             if (u(i,j,k,3)-cspd(i,j,k) <= ZERO) then
                Ip(i,j,k,3,1) = sp 
             else
                Ip(i,j,k,3,1) = sp - &
                     HALF*sigma*(sp-sm-(ONE-TWO3RD*sigma)*s6)
             endif

             if (u(i,j,k,3)-cspd(i,j,k) >= ZERO) then
                Im(i,j,k,3,1) = sm 
             else
                Im(i,j,k,3,1) = sm + &
                     HALF*sigma*(sp-sm+(ONE-TWO3RD*sigma)*s6)
             endif
             
             ! w wave
             sigma = abs(u(i,j,k,3))*dt/dz
             
             if (u(i,j,k,3) <= ZERO) then
                Ip(i,j,k,3,2) = sp
             else
                Ip(i,j,k,3,2) = sp - &
                     HALF*sigma*(sp-sm-(ONE-TWO3RD*sigma)*s6)
             endif
             
             if (u(i,j,k,3) >= ZERO) then
                Im(i,j,k,3,2) = sm
             else
                Im(i,j,k,3,2) = sm + &
                     HALF*sigma*(sp-sm+(ONE-TWO3RD*sigma)*s6)
             endif
             
             ! w+c wave
             sigma = abs(u(i,j,k,3)+cspd(i,j,k))*dt/dz
             
             if (u(i,j,k,3)+cspd(i,j,k) <= ZERO) then
                Ip(i,j,k,3,3) = sp 
             else
                Ip(i,j,k,3,3) = sp - &
                     HALF*sigma*(sp-sm-(ONE-TWO3RD*sigma)*s6)
             endif
             
             if (u(i,j,k,3)+cspd(i,j,k) >= ZERO) then
                Im(i,j,k,3,3) = sm 
             else
                Im(i,j,k,3,3) = sm + &
                     HALF*sigma*(sp-sm+(ONE-TWO3RD*sigma)*s6)
             endif
             
          end do
       end do
    end do

    deallocate(sedgey,sedgez)
       
  end subroutine ppm_type2

end module ppm_module

