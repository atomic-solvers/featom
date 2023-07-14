module mixings

! This module contains SCF mixing algorithms.

use types, only: dp
use linalg, only: solve
implicit none
private
public mixing_linear, mixing_anderson, mixing_pulay

interface
    subroutine F_fn(x, y, energies)
    ! y = F(x), also return the calculated energies to converge
    import :: dp
    implicit none
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: y(:), energies(:)
    end subroutine

    real(dp) function integral_fn(x)
    ! Computes the integral of the vector 'x'
    import :: dp
    implicit none
    real(dp), intent(in) :: x(:)
    end function

    function matvec_fn(A, v) result(r)
    ! Computes the matmul(A, v)
    import :: dp
    implicit none
    real(dp), intent(in) :: A(:,:), v(:)
    real(dp) :: r(size(A,1))
    end function

    function matmat_fn(A, B) result(r)
    ! Computes the matmul(A, B)
    import :: dp
    implicit none
    real(dp), intent(in) :: A(:,:), B(:,:)
    real(dp) :: r(size(A,1),size(B,2))
    end function

end interface

contains

subroutine mixing_linear(F, integral, x0, nenergies, max_iter, alpha, &
        L2_eps, eig_eps, x_out)
! Finds "x" so that F(x) = x
procedure(F_fn) :: F
procedure(integral_fn) :: integral
real(dp), intent(in) :: x0(:)
integer, intent(in) :: nenergies, max_iter
real(dp), intent(in) :: alpha
real(dp), intent(in) :: L2_eps, eig_eps
real(dp), intent(out) :: x_out(:)

real(dp), dimension(size(x0)) :: x_i, y_i, R_i
real(dp) :: old_energies(nenergies), energies(nenergies)
real(dp) :: x_i_norm, R_i_norm
real(dp) :: err_old, err, L2_err
integer :: i
x_i = x0
err_old = 1e12_dp
old_energies = 1e12_dp
do i = 1, max_iter
    call F(x_i, y_i, energies)
    R_i = y_i-x_i

    ! L2 norm of the "input" potential:
    x_i_norm = sqrt(integral(x_i**2))
    ! L2 norm of the "output-input" potential:
    R_i_norm = sqrt(integral(R_i**2))
    if (x_i_norm < 1e-12_dp) x_i_norm = 1e-12_dp
    L2_err = R_i_norm / x_i_norm
    err = maxval(abs(energies - old_energies))
    ! Do at least 3 iterations
    if (i >= 3 .and. L2_err < L2_eps) then
        if (err < eig_eps .and. err_old < eig_eps) then
            x_out = x_i
            return
        end if
    end if
    old_energies = energies
    err_old = err

    x_i = x_i + alpha * R_i
end do
error stop "SCF didn't converge"
end subroutine

subroutine mixing_anderson(F, integral, x0, nenergies, max_iter, alpha, &
        L2_eps, eig_eps, x_out)
! Finds "x" so that F(x) = x, uses x0 as the initial estimate
procedure(F_fn) :: F
procedure(integral_fn) :: integral
real(dp), intent(in) :: x0(:)
integer, intent(in) :: nenergies, max_iter
real(dp), intent(in) :: alpha
real(dp), intent(in) :: L2_eps, eig_eps
real(dp), intent(out) :: x_out(:)

real(dp), dimension(size(x0)) :: x_i, y_i, x1_i, R_i, R1_i, delta_R, delta_x
real(dp) :: beta
real(dp) :: sn, sd
real(dp) :: old_energies(nenergies), energies(nenergies)
real(dp) :: x_i_norm, R_i_norm
real(dp) :: err_old, err, L2_err
integer :: i
x_i = x0
err_old = 1e12_dp
old_energies = 1e12_dp
do i = 1, max_iter
    call F(x_i, y_i, energies)
    R_i = y_i-x_i

    ! L2 norm of the "input" potential:
    x_i_norm = sqrt(integral(x_i**2))
    ! L2 norm of the "output-input" potential:
    R_i_norm = sqrt(integral(R_i**2))
    if (x_i_norm < 1e-12_dp) x_i_norm = 1e-12_dp
    L2_err = R_i_norm / x_i_norm
    err = maxval(abs(energies - old_energies))
    ! Do at least 3 iterations
    if (i >= 3 .and. L2_err < L2_eps) then
        if (err < eig_eps .and. err_old < eig_eps) then
            x_out = x_i
            return
        end if
    end if
    old_energies = energies
    err_old = err

    if (i > 1) then
        delta_x = x_i - x1_i
        delta_R = R_i - R1_i
    end if
    x1_i = x_i
    R1_i = R_i
    x_i = x_i + alpha * R_i
    if (i > 1) then
        sn = integral(R_i * delta_R)
        sd = integral(delta_R**2)
        beta = sn / sd
        x_i = x_i - beta * (delta_x + alpha * delta_R)
    end if
end do
error stop "SCF didn't converge"
end subroutine

subroutine mixing_pulay(g, integral, matvec, matmat, x0, nenergies, max_iter, &
        alpha, L2_eps, eig_eps, x_out, n, k)
! Finds "x" so that x = g(x)
! Implements Algorithm 1. from [1], we use the same notation, except we call
! the independent variable "x" instead of "rho":
! do i = 1, 2, ...
!   f_i = g(x_i) - x_i
!   if modulo(i+1, k) == 0
!     x_i = x_i + alpha*f_i - (R_i+alpha*F_i)*(F_i^T F_i)^-1 F_i^T f_i
!   else
!     x_i = x_i + alpha*f_i
! until |f_i| < tol
!
! [1] Banerjee, A. S., Suryanarayana, P., Pask, J. E. (2016). Periodic Pulay
! method for robust and efficient convergence acceleration of self-consistent
! field iterations. Chemical Physics Letters, 647, 31–35.
! http://doi.org/10.1016/j.cplett.2016.01.033
procedure(F_fn) :: g
procedure(integral_fn) :: integral
procedure(matvec_fn) :: matvec
procedure(matmat_fn) :: matmat
real(dp), intent(in) :: x0(:)
integer, intent(in) :: nenergies, max_iter, n, k
real(dp), intent(in) :: alpha
real(dp), intent(in) :: L2_eps, eig_eps
real(dp), intent(out) :: x_out(:)

real(dp), dimension(size(x0)) :: x_i, y_i, f_i
real(dp) :: old_energies(nenergies), energies(nenergies)
real(dp) :: x_i_norm, f_i_norm
real(dp) :: err_old, err, L2_err
integer :: i, j

! TODO: do not store the whole history in "f" and "x", but only the last "n"
! iterations:
real(dp) :: Ri(size(x0),n), Fi(size(x0),n), f(size(x0),max_iter), &
    x(size(x0),max_iter), FTF_inv(n,n), FTf(n)

x_i = x0
err_old = 1e12_dp
old_energies = 1e12_dp
do i = 1, max_iter
    call g(x_i, y_i, energies)
    f_i = y_i-x_i
    x(:,i) = x_i
    f(:,i) = f_i

    ! L2 norm of the "input" potential:
    x_i_norm = sqrt(integral(x_i**2))
    ! L2 norm of the "output-input" potential:
    f_i_norm = sqrt(integral(f_i**2))
    if (x_i_norm < 1e-12_dp) x_i_norm = 1e-12_dp
    L2_err = f_i_norm / x_i_norm
    err = maxval(abs(energies - old_energies))
    ! Do at least 3 iterations
    if (i >= 3 .and. L2_err < L2_eps) then
        print *, "SCF convergence error:", err
        if (err < eig_eps .and. err_old < eig_eps) then
            x_out = x_i
            return
        end if
    end if
    old_energies = energies
    err_old = err

    if (i > n .and. modulo(i, k) == 0) then
        do j = i-n+1, i
            Ri(:,j-i+n) = x(:,j)-x(:,j-1)
            Fi(:,j-i+n) = f(:,j)-f(:,j-1)
        end do
        FTF_inv = matmat(transpose(Fi), Fi)
        FTf = matvec(transpose(Fi), f_i)
        x_i = x_i + alpha * f_i - matmul((Ri+alpha*Fi), solve(FTF_inv, FTf))
    else
        x_i = x_i + alpha * f_i
    end if
end do
error stop "SCF didn't converge"
end subroutine

end module
