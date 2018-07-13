! Recursive Fortran 95 quicksort routine
! sorts real numbers into ascending numerical order
! Author: Juli Rew, SCD Consulting (juliana@ucar.edu), 9/03
! Based on algorithm from Cormen et al., Introduction to Algorithms,
! 1997 printing

! Made F conformant by Walt Brainerd

module qsort_c_module

implicit none
public :: QsortC
private :: Partition

contains

recursive subroutine QsortC(A,B)
  real, intent(inout), dimension(:) :: A
  integer, intent(inout), dimension(:) :: B
  integer :: iq

  if(size(A) > 1) then
     call Partition(A, B, iq)
     call QsortC(A(:iq-1), B(:iq-1))
     call QsortC(A(iq:), B(iq:))
  endif
end subroutine QsortC

subroutine Partition(A, B, marker)
  real, intent(in out), dimension(:) :: A
  integer, intent(in out), dimension(:) :: B
  integer, intent(out) :: marker
  integer :: i, j
  integer :: tempA, tempB
  real :: x      ! pivot point
  x = A(1)
  i= 0
  j= size(A) + 1
  
  
  do
     j = j-1
     do
        if (A(j) <= x) exit
        j = j-1
     end do
     i = i+1
     do
        if (A(i) >= x) exit
        i = i+1
     end do
     if (i .lt. j) then
        ! exchange A(i) and A(j)
        tempA = A(i)
        A(i) = A(j)
        A(j) = tempA
        tempB = B(i)
        B(i) = B(j)
        B(j) = tempB
     elseif (i == j) then
        marker = i+1
        return
     else
        marker = i
        return
     endif
     
  end do


end subroutine Partition

end module qsort_c_module

program sortdriver
  use qsort_c_module
  implicit none
  integer, parameter :: r = 10
  real, dimension(1:r) :: myproxy
  real, dimension(1:r) :: myarray = &        ! (1:r)
     (/0, 50, 20, 25, 90, 10, 5, 1, 99, 75/)
  integer, dimension(1:r) :: myindices = &        ! (1:r)
     (/1, 2, 3, 4, 5, 6, 7, 8, 9, 10/)
     
  myproxy = myarray
  
  print *, "myproxy is ", myproxy
  print *, "myarray is ", myarray
  call QsortC(myproxy, myindices)
  print *, "myproxy is ", myproxy
  print *, "index order is ", myindices
end program sortdriver