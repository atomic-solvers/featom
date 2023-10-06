! BEGIN_HEADER
! -----------------------------------------------------------------------------
! Gauss-Jacobi Quadrature Implementation
! Authors: Rohit Goswami <rgoswami[at]ieee.org>
! Source: GaussJacobiQuad Library
! License: MIT
! GitHub Repository: https://github.com/HaoZeke/GaussJacobiQuad
! Date: 2023-08-28
! Commit: 2a265ce
! -----------------------------------------------------------------------------
! This code is part of the GaussJacobiQuad library, providing an efficient
! implementation for Gauss-Jacobi quadrature nodes and weights computation.
! -----------------------------------------------------------------------------
! To cite this software:
! Rohit Goswami (2023). HaoZeke/GaussJacobiQuad: v0.1.0.
! Zenodo: https://doi.org/10.5281/ZENODO.8285112
! ---------------------------------------------------------------------
! END_HEADER

module gjp_types
implicit none
private
public sp, dp, hp, qp, gjp_sparse_matrix
integer, parameter :: dp = kind(0.d0), &
                      hp = selected_real_kind(15), &
                      qp = selected_real_kind(32), &
                      sp = kind(0.)
type gjp_sparse_matrix
    real(dp), allocatable :: diagonal(:)
    real(dp), allocatable :: off_diagonal(:)
end type gjp_sparse_matrix

end module gjp_types
module gjp_imtqlx
use gjp_types, only: dp
implicit none

contains
subroutine imtqlx(mat_size, diag, off_diag, sol_vec)
    implicit none
    integer, intent(in) :: mat_size
    real(dp), intent(inout) :: diag(mat_size)
    real(dp), intent(inout) :: off_diag(mat_size - 1)
    real(dp), intent(inout) :: sol_vec(mat_size)
    real(dp) :: precision, pivot_val, g_val, rot_val, scale_val, f_val, b_val, cos_val
    integer :: lower_bound, upper_bound, inner_i, i, iter_count
    integer, parameter :: max_iter = 30
    precision = epsilon(precision)
    off_diag(mat_size - 1) = 0.0_dp
    do lower_bound = 1, mat_size
        iter_count = 0
        do while (iter_count < max_iter)
            do upper_bound = lower_bound, mat_size
                if (upper_bound == mat_size) exit
                if (abs(off_diag(upper_bound)) <= precision * (abs(diag(upper_bound)) + abs(diag(upper_bound + 1)))) exit
            end do
            pivot_val = diag(lower_bound)
            if (upper_bound == lower_bound) exit
            if (iter_count > max_iter) then
                print*," "
                print*,"IMTQLX - Fatal error."
                print*,"Iteration limit exceeded."
                stop "Terminating due to iteration limit exceeded."
            end if

            iter_count = iter_count + 1
            g_val = (diag(lower_bound + 1) - pivot_val) / (2.0_dp * off_diag(lower_bound))
            rot_val = sqrt(g_val * g_val + 1.0_dp)
            g_val = diag(upper_bound) - pivot_val + off_diag(lower_bound) / (g_val + sign(rot_val, g_val))
            scale_val = 1.0_dp
            cos_val = 1.0_dp
            pivot_val = 0.0_dp
            do inner_i = 1, upper_bound - lower_bound
                i = upper_bound - inner_i
                f_val = scale_val * off_diag(i)
                b_val = cos_val * off_diag(i)
                if (abs(g_val) <= abs(f_val)) then
                    cos_val = g_val / f_val
                    rot_val = sqrt(cos_val * cos_val + 1.0_dp)
                    off_diag(i + 1) = f_val * rot_val
                    scale_val = 1.0_dp / rot_val
                    cos_val = cos_val * scale_val
                else
                    scale_val = f_val / g_val
                    rot_val = sqrt(scale_val * scale_val + 1.0_dp)
                    off_diag(i + 1) = g_val * rot_val
                    cos_val = 1.0_dp / rot_val
                    scale_val = scale_val * cos_val
                end if
                g_val = diag(i + 1) - pivot_val
                rot_val = (diag(i) - g_val) * scale_val + 2.0_dp * cos_val * b_val
                pivot_val = scale_val * rot_val
                diag(i + 1) = g_val + pivot_val
                g_val = cos_val * rot_val - b_val
                f_val = sol_vec(i + 1)
                sol_vec(i + 1) = scale_val * sol_vec(i) + cos_val * f_val
                sol_vec(i) = cos_val * sol_vec(i) - scale_val * f_val
            end do
            diag(lower_bound) = diag(lower_bound) - pivot_val
            off_diag(lower_bound) = g_val
            off_diag(upper_bound) = 0.0_dp
        end do
    end do
    call dsort2a(mat_size, diag, sol_vec)
end subroutine
subroutine dsort2a(n, x, y)
    integer, intent(in) :: n
    real(dp), intent(inout) :: x(n), y(n)
    integer :: i, j, min_idx
    real(dp) :: temp

    do i = 1, n - 1
        min_idx = i
        do j = i + 1, n
            if (x(j) < x(min_idx)) min_idx = j
        end do
        if (min_idx /= i) then
            temp = x(i)
            x(i) = x(min_idx)
            x(min_idx) = temp
            temp = y(i)
            y(i) = y(min_idx)
            y(min_idx) = temp
        end if
    end do
end subroutine dsort2a

end module gjp_imtqlx

module gjp_common
use gjp_types, only: dp, gjp_sparse_matrix
implicit none

contains
function jacobi_matrix(n, alpha, beta) result(jacmat)
    integer, intent(in) :: n
    real(dp), intent(in) :: alpha, beta
    type(gjp_sparse_matrix) :: jacmat
    integer :: idx
    real(dp) :: ab, abi, a2b2

    allocate (jacmat%diagonal(n))
    allocate (jacmat%off_diagonal(n - 1))

    ab = alpha + beta
    abi = 2.0_dp + ab
    jacmat%diagonal(1) = (beta - alpha) / abi
    jacmat%off_diagonal(1) = 4.0_dp * (1.0_dp + alpha) * (1.0_dp + beta) &
                             / ((abi + 1.0_dp) * abi * abi)
    a2b2 = beta * beta - alpha * alpha
    do idx = 2, n
        abi = 2.0_dp * idx + ab
        jacmat%diagonal(idx) = a2b2 / ((abi - 2.0_dp) * abi)
        abi = abi**2
        if (idx < n) then
            jacmat%off_diagonal(idx) = 4.0_dp * idx * (idx + alpha) * (idx + beta) &
                                       * (idx + ab) / ((abi - 1.0_dp) * abi)
        end if
    end do
    jacmat%off_diagonal(1:n - 1) = sqrt(jacmat%off_diagonal(1:n - 1))
end function jacobi_matrix
function jacobi_zeroeth_moment(alpha, beta) result(zmom)
    real(dp), intent(in) :: alpha, beta
    real(dp) :: zmom
    real(dp) :: ab, abi

    ab = alpha + beta
    abi = 2.0_dp + ab

    zmom = 2.0_dp**(alpha + beta + 1.0_dp) * exp(log_gamma(alpha + 1.0_dp) &
                                                 + log_gamma(beta + 1.0_dp) - log_gamma(abi))

    if (zmom <= 0.0_dp) then
        error stop "Zeroth moment is not positive but should be"
    end if
end function jacobi_zeroeth_moment

end module gjp_common

module gjp_algo665
use gjp_types, only: dp, gjp_sparse_matrix
use gjp_imtqlx, only: imtqlx
use gjp_common, only: jacobi_matrix, jacobi_zeroeth_moment
implicit none
contains
subroutine gauss_jacobi_algo665(npts, alpha, beta, x, wts)
    integer, intent(in) :: npts
    real(dp), intent(in) :: alpha, beta
    real(dp), intent(out) :: x(npts), wts(npts)
    real(dp) :: zeroeth_moment
    type(gjp_sparse_matrix) :: jacobi_mat
    real(dp) :: diagonal_elements(npts), &
                off_diagonal_elements(npts - 1)

    jacobi_mat = jacobi_matrix(npts, alpha, beta)
    zeroeth_moment = jacobi_zeroeth_moment(alpha, beta)
    diagonal_elements = jacobi_mat%diagonal(1:npts)
    off_diagonal_elements = jacobi_mat%off_diagonal(1:npts - 1)
    wts = 0.0_dp
    x = diagonal_elements
    wts(1) = sqrt(zeroeth_moment)
    call imtqlx(npts, x, off_diagonal_elements, wts)
    wts = wts**2

end subroutine gauss_jacobi_algo665

end module gjp_algo665
