subroutine hunt(sigma_nodes, KB, sigma_particle, jlo) ! Z, KB, self%ZP(ii), NZR
    ! from numerical recipies vol 2
    use parameters

    integer, intent(inout) :: jlo ! sigma layer below particle?
    integer, intent(in) :: KB ! number of sigma layers
    real(sp), dimension(KB), intent(in) :: sigma_nodes ! sigma coordinate value
    real(sp), intent(in) :: sigma_particle ! depth of particle

    integer :: inc, jhi, jm
    logical :: ascnd

    ascnd = (sigma_nodes(KB) > sigma_nodes(1)) ! bottom sigma layer greater than first layer
    if ((jlo <= 0) .or. (jlo > KB)) then
        jlo = 0
        jhi = KB + 1
        goto 3
    endif
    inc = 1

    if ((sigma_particle >= sigma_nodes(jlo)) .eqv. ascnd) then
        1    jhi = jlo + inc
        if(jhi > KB)then
            jhi=n+1
        else if ((sigma_particle >= sigma_nodes(jhi)) .eqv. ascnd) then
            jlo = jhi
            inc = inc + inc
            goto 1
        end if
    else
        jhi = jlo
        2    jlo = jhi - inc
        if (jlo < 1) then
            jlo=0
        else if ((sigma_particle < sigma_nodes(jlo)) .eqv. ascnd) then
            jhi = jlo
            inc = inc + inc
            go to 2
        end if
    end if
    3 if (jhi-jlo == 1) then
        if (sigma_particle == sigma_nodes(KB)) jlo = KB - 1
        if (sigma_particle == sigma_nodes(1)) jlo = 1
        return
    end if
    jm = (jhi + jlo)/2
    if (sigma_particle >= sigma_nodes(jm) .eqv. ascnd) then
        jlo = jm
    else
        jhi = jm
    end if
    go to 3

end subroutine hunt


subroutine spline(x, y, n2, yp1, ypn, y2)
    ! from numerical recipies vol 2, but modfied so that nmax=50
    use parameters, only : sp
    implicit none
    integer  :: n2
    real(sp), intent(in) :: x(n2), y(n2), yp1, ypn
    real(sp), intent(out), dimension(n2) :: y2

    integer  :: i, k
    integer, parameter :: nmax=50
    real(sp) :: p, qn, sig, un
    real(sp), dimension(nmax) :: u

    if (yp1 > 0.99e30) then ! force natural lower boundary
        y2(1) = 0.0_SP
        u(1) = 0.0_SP
    else ! or set specific values
        y2(1) = -0.5_SP
        u(1) = (3.0_SP/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
    end if
    do i = 2, n2-1 ! tridiagonal algorithm decomp
        sig = (x(i) - x(i-1)) / (x(i+1) - x(i-1))
        p = sig * y2(i-1) + 2.0_SP
        y2(i) = (sig-1.) / p
        u(i) = (6.0*((y(i+1) - y(i)) / (x(i+1) - x(i)) - (y(i) - y(i-1))/(x(i) - x(i-1)))/(x(i+1) - x(i-1)) - (- sig*u(i-1)))/p
    end do
    if (ypn > 0.99e30) then ! force natural upper boundary
        qn = 0.0_SP
        un = 0.0_SP
    else
        qn = 0.5_SP
        un = (3.0_SP/(x(n2)-x(n2-1)))*(ypn-(y(n2)-y(n2-1))/(x(n2)-x(n2-1)))
    end if
    y2(n2) = (un - qn * u(n2-1)) / (qn * y2(n2-1) + 1.0_SP)
    do k = n2-1, 1, -1
        y2(k) = y2(k) * y2(k+1) + u(k)
    end do

end subroutine spline
