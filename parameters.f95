module parameters
  use MOD_PREC
  implicit none

  ! used in main
  logical, parameter :: includeAlgae = .true.
  logical, parameter :: includeFish = .true.
  logical, parameter :: use_ncd = .false.
  logical, parameter :: strict_integration = .false. ! set mass transfer
  logical, parameter :: continue_sim = .false.

  real(sp), parameter :: irradSurf = 650.0_SP ! W/M^2
  integer, parameter :: iocp = 101, iocs = 102, iotox = 103
  integer, parameter :: iofp = 201, iofs = 202, iophys = 301, iorun = 302, iovar=303

  ! physical and mathematical constants
  real(sp), parameter :: boltzmann = 1.3806488_SP * 10.0_SP**(-23.0_SP) ! m2 kg s-2 K-1
  real(sp), parameter :: microcystinRadius = 1.5_SP * 10.0_SP**(-9.0_SP) ! m
  real(sp), parameter :: avogadro = 6022.0_SP * 10.0_SP**(20.0_SP) ! per mol
  real(sp), parameter :: planckNumber = 663.0_SP * 10.0_SP**(-7.0_SP) ! Js
  real(sp), parameter :: lightSpeed = 2998.0_SP  * 10.0_SP**(5.0_SP) ! meters per second

  
  ! Runge-Kutta integration coefficients
  integer, parameter :: MSTAGE = 4 ! number of stages
  real(sp), parameter, dimension(4) :: A_RK = (/ 0.0_SP, 0.5_SP, 0.5_SP, 1.0_SP/) ! ERK coefficients (A)
  real(sp), parameter, dimension(4) :: B_RK = (/ 1.0_SP/6.0_SP, 1.0_SP/3.0_SP, 1.0_SP/3.0_SP, 1.0_SP/6.0_SP /) ! ERK coefficients (B)
  real(sp), parameter, dimension(4) :: C_RK = (/ 0.0_SP, 0.5_SP, 0.5_SP, 1.0_SP /)  ! ERK coefficients (C)



end module parameters
