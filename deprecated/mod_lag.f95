
 
 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine INTERP_V(lag, XP, YP, ZP, UIN, VIN, WIN) ! OK
    ! linear interpolation of velocity field at particle positions
    use ALL_VARS, only : A1U, A2U, NBE, YC, XC, DZ, ZZ
    use LIMS, only : KBM1, N, KB

    class(LAG_OBJ), intent(inout) :: lag
    real(sp), intent(in), dimension(lag%ndrft) :: XP, YP, ZP ! ZP is sigma depth
    real(sp), intent(in), dimension(0:N, 1:KB) :: UIN, VIN, WIN

    integer, dimension(lag%ndrft) ::  INWATER
    integer :: ii, host, E1, E2, E3, K1, K2, K
    real(SP) :: DUDX, DUDY, DVDX, DVDY, DWDX, DWDY, UE01, UE02, VE01, VE02, WE01, WE02
    real(SP) :: ZF1, ZF2, X0C, Y0C
    logical :: ALL_FOUND

    INWATER(:) = 1
    lag%FOUND(:) = 0
    call lag%FHE_QUICK(XP, YP, ALL_FOUND) ! determine host element
    if (.not. ALL_FOUND) call lag%FHE_ROBUST(XP, YP, INWATER)

    particle_loop: do ii = 1, lag%ndrft
      if ( (lag%INDOMAIN(ii) .eq. 0) .or. (INWATER(ii)) .eq. 0) cycle ! skip particles outside domain
      host = lag%HOST(ii)
      E1  = NBE(host,1)
      E2  = NBE(host,2)
      E3  = NBE(host,3)
      X0C = XP(ii) - XC(host)
      Y0C = YP(ii) - YC(host)

      ! Determine sigma layers above and below particle
      if (ZP(ii) .gt. ZZ(1)) then ! Particle near surface
         K1  = 1
         K2  = 1
         ZF1 = 1.0_SP
         ZF2 = 0.0_SP
      else if (ZP(ii) .lt. ZZ(KBM1)) then ! Particle near bottom
         K1 = KBM1
         K2 = KBM1
         ZF1 = 1.0_SP
         ZF2 = 0.0_SP
      else
         K1 = int( (ZZ(1) - ZP(ii)) / DZ(1) ) + 1;
         K2 = K1 + 1
         ZF1 = (ZP(ii) - ZZ(K2)) / DZ(1)
         ZF2 = (ZZ(K1) - ZP(ii)) / DZ(1)
      end if


      ! Linear interpolation of velocity in sigma level above particle
      K = K1
      DUDX = A1U(host,1)*UIN(host,K)+A1U(host,2)*UIN(E1,K)+A1U(host,3)*UIN(E2,K)+A1U(host,4)*UIN(E3,K)
      DUDY = A2U(host,1)*UIN(host,K)+A2U(host,2)*UIN(E1,K)+A2U(host,3)*UIN(E2,K)+A2U(host,4)*UIN(E3,K)
      DVDX = A1U(host,1)*VIN(host,K)+A1U(host,2)*VIN(E1,K)+A1U(host,3)*VIN(E2,K)+A1U(host,4)*VIN(E3,K)
      DVDY = A2U(host,1)*VIN(host,K)+A2U(host,2)*VIN(E1,K)+A2U(host,3)*VIN(E2,K)+A2U(host,4)*VIN(E3,K)
      DWDX = A1U(host,1)*WIN(host,K)+A1U(host,2)*WIN(E1,K)+A1U(host,3)*WIN(E2,K)+A1U(host,4)*WIN(E3,K)
      DWDY = A2U(host,1)*WIN(host,K)+A2U(host,2)*WIN(E1,K)+A2U(host,3)*WIN(E2,K)+A2U(host,4)*WIN(E3,K)
      UE01 = UIN(host,K) + DUDX*X0C + DUDY*Y0C
      VE01 = VIN(host,K) + DVDX*X0C + DVDY*Y0C
      WE01 = WIN(host,K) + DWDX*X0C + DWDY*Y0C

      ! Linear interpolation of velocity in sigma level below particle
      K = K2
      DUDX = A1U(host,1)*UIN(host,K)+A1U(host,2)*UIN(E1,K)+A1U(host,3)*UIN(E2,K)+A1U(host,4)*UIN(E3,K)
      DUDY = A2U(host,1)*UIN(host,K)+A2U(host,2)*UIN(E1,K)+A2U(host,3)*UIN(E2,K)+A2U(host,4)*UIN(E3,K)
      DVDX = A1U(host,1)*VIN(host,K)+A1U(host,2)*VIN(E1,K)+A1U(host,3)*VIN(E2,K)+A1U(host,4)*VIN(E3,K)
      DVDY = A2U(host,1)*VIN(host,K)+A2U(host,2)*VIN(E1,K)+A2U(host,3)*VIN(E2,K)+A2U(host,4)*VIN(E3,K)
      DWDX = A1U(host,1)*WIN(host,K)+A1U(host,2)*WIN(E1,K)+A1U(host,3)*WIN(E2,K)+A1U(host,4)*WIN(E3,K)
      DWDY = A2U(host,1)*WIN(host,K)+A2U(host,2)*WIN(E1,K)+A2U(host,3)*WIN(E2,K)+A2U(host,4)*WIN(E3,K)
      UE02 = UIN(host,K) + DUDX*X0C + DUDY*Y0C
      VE02 = VIN(host,K) + DVDX*X0C + DVDY*Y0C
      WE02 = WIN(host,K) + DWDX*X0C + DWDY*Y0C

      ! Interpolate particle velocity between two sigma layers
      lag%UP(ii) = UE01*ZF1 + UE02*ZF2
      lag%VP(ii) = VE01*ZF1 + VE02*ZF2
      lag%WP(ii) = WE01*ZF1 + WE02*ZF2

    end do particle_loop
    return
  end subroutine INTERP_V


  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine INTERP_KH(self, ZKH, DZKH, KHOUT, DKHOUT)
  ! Obtain a spline interpolation on the vertical with the provided eddy diffusivity (ZKH) and its derivative (DZKH) at grid point points
  ! then linear interpolation on the horizontal
  !  RETURNS:   Both dkh/dz (DKHOUT) and kh (KHOUT)

    use ALL_VARS
    class(LAG_OBJ), intent(inout) :: self
    real(SP), intent(out), dimension(self%ndrft) :: DKHOUT, KHOUT
    real(SP), intent(in) ,dimension(0:M, KB) :: DZKH, ZKH

    real(SP) :: X0C, Y0C, COF1, COF2, COF3
    integer :: N1, N2, N3, ii
    real(SP) :: DKHR1, DKHR2, DKHR3 !, DKHR4
    real(SP) :: KHR1, KHR2, KHR3 !, KHR4
    !real(SP) :: DDKHR1, DDKHR2, DDKHR3, DDKHR4
    real(SP) :: DKHTMP, DZP
    integer, dimension(0:KB+1) :: NZRINDX
    real(SP) :: HK, HK2, AK, AK2, AK3, BK, BK2, BK3
    integer :: KLO, KHI, NZR, host

    ! Interpolate eddy diffusivity and its derivative
    KHOUT  = 0.0_SP
    DKHOUT = 0.0_SP
    NZRINDX(0) = 1
    do ii = 1, KB
      NZRINDX(ii) = ii
    end do
    NZRINDX(KB+1) = KB

    particle_loop: do ii = 1, self%ndrft
      if (self%INDOMAIN(ii) .eq. 0) cycle ! skip particles outisde domain
      host = self%host(ii) ! element containing particle
      N1 = NV(host, 1); N2 = NV(host, 2); N3 = NV(host, 3) ! get node indices of host element
      X0C = self%XP(ii) - XC(host); Y0C = self%YP(ii) - YC(host) ! distance from element center

      ! DERIVATIVE OF THE DIFFUSION
      ! find vertical location
      NZR = floor( float(KBM1)*abs(self%zp(ii)) ) ! guess value for hunt
      !call HUNT(Z, KB, self%ZP(ii), NZR) ! Z is sigma coordinate value from ALL_VARS, only necessary if sigma layers are different thicknesses

      KHI = NZRINDX(NZR)
      KLO = NZRINDX(NZR + 1)

      ! as k in z(k) increases, sigma decreases.
      HK=Z(KHI)-Z(KLO)
      AK=(Z(KHI)-self%ZP(ii))/HK
      AK2=AK**2
      BK=(self%ZP(ii)-Z(KLO))/HK
      BK2=BK**2

      DKHR1=(-ZKH(N1,KLO)+ZKH(N1,KHI))/HK+((-3*AK2+1)*DZKH(N1,KLO)+(3*BK2-1)*DZKH(N1,KHI))*HK/6
      DKHR2=(-ZKH(N2,KLO)+ZKH(N2,KHI))/HK+((-3*AK2+1)*DZKH(N2,KLO)+(3*BK2-1)*DZKH(N2,KHI))*HK/6
      DKHR3=(-ZKH(N3,KLO)+ZKH(N3,KHI))/HK+((-3*AK2+1)*DZKH(N3,KLO)+(3*BK2-1)*DZKH(N3,KHI))*HK/6

      COF1=AW0(self%HOST(ii),1)*DKHR1 +AW0(self%HOST(ii),2)*DKHR2+AW0(self%HOST(ii),3)*DKHR3
      COF2=AWX(self%HOST(ii),1)*DKHR1 +AWX(self%HOST(ii),2)*DKHR2+AWX(self%HOST(ii),3)*DKHR3
      COF3=AWY(self%HOST(ii),1)*DKHR1 +AWY(self%HOST(ii),2)*DKHR2+AWY(self%HOST(ii),3)*DKHR3

      DKHTMP = COF1 + COF2*X0C + COF3*Y0C
      DKHOUT(ii) = DKHTMP / (self%HP(ii) - self%EP(ii))     !--want the answer to be dkh/dz not dkh/dsigma, changed sign, keeney 1/7

      ! DIFFUSION ITSELF
      ! find z in grid again, as per visser, but in sigma
      DZP = self%ZP(ii) + 0.5*DKHOUT(ii)*DTRW/(self%HP(ii) - self%EP(ii)) ! changed sign, keeney 1/7
      ! adding 0.5*dkhtmp*dtrw, new z can be out of [0;-1]
      DZP = min(DZP, ZERO)
      DZP = max(DZP, -1.0_SP)

      ! find vertical location
      !call HUNT(Z, KB, DZP, NZR)
      NZR = floor( float(KBM1)*abs(DZP ) ) ! guess value for hunt
      KHI=NZRINDX(NZR)
      KLO=NZRINDX(NZR+1)

      HK = Z(KHI) - Z(KLO)
      HK2 = HK**2
      AK = (Z(KHI)-DZP)/HK
      AK3 = AK**3
      BK = (DZP-Z(KLO))/HK
      BK3 = BK**3

      KHR1 = AK*(ZKH(N1,KLO)) + BK*(ZKH(N1,KHI)) + ((AK3-AK)*(DZKH(N1,KLO))+(BK3-BK)*(DZKH(N1,KHI)))*HK2/6
      KHR2 = AK*(ZKH(N2,KLO)) + BK*(ZKH(N2,KHI)) + ((AK3-AK)*(DZKH(N2,KLO))+(BK3-BK)*(DZKH(N2,KHI)))*HK2/6
      KHR3 = AK*(ZKH(N3,KLO)) + BK*(ZKH(N3,KHI)) + ((AK3-AK)*(DZKH(N3,KLO))+(BK3-BK)*(DZKH(N3,KHI)))*HK2/6

      COF1=AW0(self%HOST(ii),1)*KHR1 + AW0(self%HOST(ii),2)*KHR2+AW0(self%HOST(ii),3)*KHR3
      COF2=AWX(self%HOST(ii),1)*KHR1 + AWX(self%HOST(ii),2)*KHR2+AWX(self%HOST(ii),3)*KHR3
      COF3=AWY(self%HOST(ii),1)*KHR1 + AWY(self%HOST(ii),2)*KHR2+AWY(self%HOST(ii),3)*KHR3

      KHOUT(ii) = COF1 + COF2*X0C + COF3*Y0C


