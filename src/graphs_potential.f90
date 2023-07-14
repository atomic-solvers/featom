module graphs_potential

use types, only: dp
use mesh, only: meshexp
use schroed_glob, only: solve_schroed
use dirac, only: solve_dirac
use feutils, only: define_connect, get_quad_pts, get_parent_quad_pts_wts, &
        get_parent_nodes, phih, dphih, phih_array, dphih_array
use fe, only: assemble_radial_S, assemble_radial_H, assemble_radial_dirac_SH
use linalg, only: eigh
use string_utils, only: str
use schroed_dirac_solver, only: total_energy
implicit none

contains

subroutine run_convergence_potential(study_type, dirac_int, p_or_Ne, &
     potential_type, alpha_int, directory)
! <study_type> can be,
!       0: error as p is varied
!       1: error as rmax is varied
!       2: error as Ne is varied
integer, intent(in) :: study_type
! <equation> can be,
!       0: Schroedinger
!       1: Dirac
integer, intent(in) :: dirac_int
! For <study_type>
!       0, 1: 3rd argument p_or_Ne = Ne (Number of elements)
!       2   : 3rd argument p_or_Ne = p  (Polynomial order)
integer, intent(in) :: p_or_Ne
! <potential_type> can be,
!       0: Coulomb
!       1: Harmonic
integer, intent(in) :: potential_type
! <alpha> can be 0, 1, -1 (-1 implies beta). used only for Dirac.
integer, intent(in) :: alpha_int

! The directory where to save the output files
character(len=*), intent(in) :: directory

integer :: p, Ne, Nq, Z, u, i, DOFs
real(dp) :: rmax, a, asympt, c
real(dp), allocatable :: alpha_j(:), alpha(:)
integer :: Lmax, Lmin, kappa
real(dp) :: optim_a(2:7)
integer :: Nes(11)
real(dp), allocatable :: lam(:), eigfn(:,:,:), xq(:,:)
character(len=:), allocatable :: filename
Z = 92
rmax = 50
a = 200
Ne = 4
Nq = 64
p = 25
c = 137.0359895_dp

if (study_type == 2) then
    p = p_or_Ne
else
    Ne = p_or_Ne
end if
filename = str(p_or_Ne)

if (alpha_int == 0) then
   filename = trim(filename) // "_0"
else if (alpha_int == 1) then
   filename = trim(filename) // "_1"
else
   filename = trim(filename) // "_beta"
end if

Lmax=6
Lmin=-7

allocate(alpha(Lmin:Lmax), alpha_j(Lmin:Lmax))
allocate(xq(Nq, Ne))
do kappa = Lmin, Lmax
    if (kappa == 0) cycle
    ! asymptotic at r = 0
    asympt = sqrt(kappa**2 - Z**2 / c**2)
    ! solve for P/r**alpha
    if (alpha_int == -1) then
        alpha(kappa) = asympt
        ! power of r for Gauss-Jacobi quadrature
        alpha_j(kappa) = 2*asympt - 2
    else
        alpha(kappa) = alpha_int
        ! don't use Gauss-Jacobi quadrature
        alpha_j(kappa) = -2
    end if
end do

optim_a = [58.985048682607555, 163.13530060338942, 340.82602039608668, &
            444.68894311026423, 591.72463734788732, 596.61404750045062]

if (Ne >= 2 .and. Ne <= 7) then
    a = optim_a(Ne)
end if

a = 100

if (study_type == 0) then
    filename = "conv_" // trim(filename) // ".txt"
else if (study_type == 1) then
    filename = "rmax_" // trim(filename) // ".txt"
else
    filename = "ne_" // trim(filename) // ".txt"
end if

if (dirac_int == 1) then
    filename = "dirac_" // trim(filename)
else
    filename = "schroed_" // trim(filename)
end if

if (potential_type == 0) then
    filename = "coulomb_" // trim(filename)
else
    filename = "harmonic_" // trim(filename)
end if

filename = directory // "/" // filename

open(newunit=u, file=filename, status="replace")
print "(a3,a6,a5,a8,a3,a3,a5)", "Z", "rmax", "Ne", "a", "p", "Nq", &
    "DOFs"

if (study_type == 0) then
    do i = 7, 63
        ! change p for p-conv study. p must be less than 31.
        if (dirac_int == 1 .and. i > 22) then
            exit
        end if
        p = i
        call total_energy(Z, rmax, Ne, a, p, Nq, DOFs, alpha_int, dirac_int, &
            c, potential_type, Lmin, alpha_j, alpha, lam, eigfn, xq)
        print "(i3, f6.1, i5, f8.1, i3, i3, i5, f22.12)", Z, rmax, Ne, a, p, Nq, &
            DOFs, sum(lam)
        write(u,*) Z, rmax, Ne, a, p, Nq, DOFs, sum(lam), lam
    end do
else if (study_type == 1) then
    rmax = 0.5_dp
    do while (rmax < 9)
        ! change rmax for rmax study.
        if (dirac_int == 1) then
            p = 22
        end if
        call total_energy(Z, rmax, Ne, a, p, Nq, DOFs, alpha_int, dirac_int, &
            c, potential_type, Lmin, alpha_j, alpha, lam, eigfn, xq)
        print "(i3, f6.1, i5, f8.1, i3, i3, i5, f22.12)", Z, rmax, Ne, a, p, Nq, &
            DOFs, sum(lam)
        write(u,*) Z, rmax, Ne, a, p, Nq, DOFs, sum(lam), lam
        rmax = rmax + 0.3_dp
    end do
else if (study_type == 2) then
    Nes = [3,4,5,6,7,8,9,10,20,25,30]
    a = 600
    do i = 1, 11
        Ne = Nes(i)
        deallocate(xq)
        allocate(xq(Nq, Ne))
        call total_energy(Z, rmax, Ne, a, p, Nq, DOFs, alpha_int, dirac_int, &
            c, potential_type, Lmin, alpha_j, alpha, lam, eigfn, xq)
        print "(i3, f6.1, i5, f8.1, i3, i3, i5, f22.12)", Z, rmax, Ne, a, p, Nq, &
            DOFs, sum(lam)
        write(u,*) Z, rmax, Ne, a, p, Nq, DOFs, sum(lam), lam
    end do
end if
close(u)
end subroutine

end module
