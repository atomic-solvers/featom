module jacobi
use types, only: dp
use constants, only: pi
implicit none
contains

subroutine gauss_jacobi(n, a, b, x, w)
integer, intent(in) :: n
real(dp), intent(in) :: a, b
real(dp), intent(out) :: x(n), w(n)
real(dp), dimension(ceiling(n/2._dp)) :: x1, ders1
real(dp), dimension(n/2) :: x2, ders2
real(dp) :: ders(n), C
integer :: i
call recurrence(n, ceiling(n/2._dp), a, b, x1, ders1)
call recurrence(n, n/2, b, a, x2, ders2)
do i = 1, n/2
    x(i) = -x2(n/2-i+1)
    ders(i) = ders2(n/2-i+1)
end do
do i = 1, ceiling(n/2._dp)
    x(n/2+i) = x1(i)
    ders(n/2+i) = ders1(i)
end do
w = 1.0d0 / ((1.0d0 - x**2) * ders**2)
C = 2**(a+b+1) * exp( log_gamma(n+a+1) - log_gamma(n+a+b+1) + &
    log_gamma(n+b+1) - log_gamma(n+1._dp) );
w = w * C
end subroutine

subroutine recurrence(n, n2, a, b, x, PP)
integer, intent(in) :: n, n2
real(dp), intent(in) :: a, b
real(dp), intent(out) :: x(n2), PP(n2)
real(dp) :: dx(n2), P(n2)
integer :: r(n2), l, i
real(dp) :: C, T

do i = 1, n2
    r(i) = n2 - i + 1
end do

do i = 1, n2
    C = (2*r(i)+a-0.5d0)*pi/(2*n+a+b+1)
    T = C + 1/(2*n+a+b+1)**2 * ((0.25d0-a**2)/tan(0.5d0*C) - (0.25d0-b**2)*tan(0.5d0*C))
    x(i) = cos(T)
end do

dx = 1.0d0
l = 0

do while (maxval(abs(dx)) > sqrt(epsilon(1.0d0))/1000 .and. l < 10)
    l = l + 1
    call eval_jacobi_poly(x, n, a, b, P, PP)
    dx = -P / PP
    x = x + dx
end do

call eval_jacobi_poly(x, n, a, b, P, PP)
end subroutine

subroutine eval_jacobi_poly(x, n, a, b, P, Pp)
integer, intent(in) :: n
real(dp), intent(in) :: a, b
real(dp), intent(in) :: x(:)
real(dp), dimension(size(x)), intent(out) :: P, Pp
real(dp), dimension(size(x)) :: Pm1, Ppm1, Pa1, Ppa1
integer :: k, i
real(dp) :: A_val, B_val, C_val, D_val

P = 0.5d0 * (a - b + (a + b + 2) * x)
Pm1 = 1.0d0
Pp = 0.5d0 * (a + b + 2)
Ppm1 = 0.0d0

if (n == 0) then
    P = Pm1
    Pp = Ppm1
end if

do k = 1, n-1
    A_val = 2 * (k + 1) * (k + a + b + 1) * (2 * k + a + b)
    B_val = (2 * k + a + b + 1) * (a**2 - b**2)
    C_val = product([(2 * k + a + b + i, i = 0, 2)])
    D_val = 2 * (k + a) * (k + b) * (2 * k + a + b + 2)

    Pa1 = ((B_val + C_val * x) * P - D_val * Pm1) / A_val
    Ppa1 = ((B_val + C_val * x) * Pp + C_val * P - D_val * Ppm1) / A_val

    Pm1 = P
    P = Pa1
    Ppm1 = Pp
    Pp = Ppa1
end do
end subroutine

end module
