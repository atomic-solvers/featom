! BEGIN_HEADER
! -----------------------------------------------------------------------------
! Gauss-Jacobi Quadrature Implementation
! Authors: Rohit Goswami <rgoswami[at]ieee.org>
! Source: GaussJacobiQuad Library
! License: MIT
! GitHub Repository: https://github.com/HaoZeke/GaussJacobiQuad
! Date: 2023-09-06
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

module gjp_lapack
implicit none

integer, parameter :: dp = kind(0.d0)
interface
    subroutine DSTEQR(COMPZ, N, D, E, Z, LDZ, WORK, INFO)
        import :: dp
        character :: COMPZ
        integer :: N, LDZ, INFO
        real(dp) :: D(*), E(*), Z(LDZ, *), WORK(*)
    end subroutine DSTEQR
end interface

contains

end module gjp_lapack

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

module gjp_gw
use gjp_types, only: dp, gjp_sparse_matrix
use gjp_lapack, only: DSTEQR
use gjp_common, only: jacobi_matrix, jacobi_zeroeth_moment
implicit none
contains
subroutine gauss_jacobi_gw(npts, alpha, beta, x, wts)
    integer, intent(in) :: npts
    real(dp), intent(in) :: alpha, beta
    real(dp), intent(out) :: x(npts), wts(npts)
    real(dp) :: zeroeth_moment
    type(gjp_sparse_matrix) :: jacobi_mat
    real(dp) :: diagonal_elements(npts), &
                off_diagonal_elements(npts - 1), &
                eigenvectors(npts, npts), &
                workspace(2 * npts - 2)
    integer :: computation_info, i

    jacobi_mat = jacobi_matrix(npts, alpha, beta)
    zeroeth_moment = jacobi_zeroeth_moment(alpha, beta)
    diagonal_elements = jacobi_mat%diagonal(1:npts)
    off_diagonal_elements = jacobi_mat%off_diagonal(1:npts - 1)
    eigenvectors = 0.0_dp
    do i = 1, npts
        eigenvectors(i, i) = 1.0_dp
    end do
    call DSTEQR('V', npts, diagonal_elements, off_diagonal_elements, &
                eigenvectors, npts, workspace, computation_info)

    if (computation_info /= 0) then
        write (*, *) 'Error in DSTEQR, info:', computation_info
        error stop
    end if
    x = diagonal_elements
    wts = eigenvectors(1, :)**2 * zeroeth_moment

end subroutine gauss_jacobi_gw

end module gjp_gw

