module MOD_PREC ! Define floating point precision using kind
  implicit none
  !integer, parameter :: SP = SELECTED_REAL_KIND(6,30) ! single precision
  integer, parameter :: SP = SELECTED_REAL_KIND(12,300) ! double precision
end module MOD_PREC