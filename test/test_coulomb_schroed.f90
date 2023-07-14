program test_coulomb_schroed

use types, only: dp
use mesh, only: meshexp
use schroed_glob, only: solve_schroed
use dirac, only: solve_dirac
use feutils, only: get_parent_quad_pts_wts
use graphs_potential, only: run_convergence_potential
use string_utils, only: str
use schroed_dirac_solver, only: total_energy
implicit none

real(dp), allocatable :: xe(:)              ! element coordinates
real(dp), allocatable :: xiq(:), wtq(:)     ! quadrature points and weights
integer :: p, Ne, Nq, Z, i, DOFs
real(dp) :: rmin, rmax, a, err, Etot
real(dp), allocatable :: energies(:), xq(:,:), eigfn(:,:,:)
real(dp), parameter :: Etot_ref = -10972.97142857139_dp
real(dp), parameter :: energies_ref(*) = [ &
                            -4232._dp, &
                            -1058._dp, &
                            -1058._dp, &
                            -470.2222222222222222_dp, &
                            -470.2222222222222222_dp, &
                            -470.2222222222222222_dp, &
                            -264.5_dp, &
                            -264.5_dp, &
                            -264.5_dp, &
                            -264.5_dp, &
                            -169.28_dp, &
                            -169.28_dp, &
                            -169.28_dp, &
                            -169.28_dp, &
                            -169.28_dp, &
                            -117.55555555555555_dp, &
                            -117.55555555555555_dp, &
                            -117.55555555555555_dp, &
                            -117.55555555555555_dp, &
                            -117.55555555555555_dp, &
                            -117.55555555555555_dp, &
                            -86.36734693877_dp, &
                            -86.36734693877_dp, &
                            -86.36734693877_dp, &
                            -86.36734693877_dp, &
                            -86.36734693877_dp, &
                            -86.36734693877_dp, &
                            -86.36734693877_dp &
                            ]


! <study_type> can be,
!       0: error as p is varied
!       1: error as rmax is varied
!       2: error as Ne is varied
integer :: study_type
! <equation> can be,
!       0: Schroedinger
!       1: Dirac
integer :: dirac_int
! For <study_type>
!       0, 1: 3rd argument p_or_Ne = Ne (Number of elements)
!       2   : 3rd argument p_or_Ne = p  (Polynomial order)
integer :: p_or_Ne
! <potential_type> can be,
!       0: Coulomb
!       1: Harmonic
integer :: potential_type
! <alpha> can be 0, 1, -1 (-1 implies beta). used only for Dirac.
integer :: alpha_int

! The directory where to save the output files
character(len=:), allocatable :: directory

integer :: u, ind
real(dp) :: asympt, c
real(dp), allocatable :: alpha_j(:), alpha(:)
integer :: Lmax, Lmin, kappa
real(dp) :: optim_a(2:7)
integer :: Nes(11)
character(len=128) :: arg

Z = 92
rmin = 0
rmax = 50
a = 200
Ne = 4
Nq = 53
p = 26

i = 7
!call run_convergence_potential(0, 0, i, &
!    0, -1, ".")


study_type = 0
dirac_int = 0
p_or_Ne = i
potential_type = 0
alpha_int = -1
directory = "."

c = 137.0359895_dp

if (study_type == 2) then
    p = p_or_Ne
else
    Ne = p_or_Ne
end if

Lmax=6
Lmin=-7

allocate(alpha(Lmin:Lmax), alpha_j(Lmin:Lmax))
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

print "(a3,a6,a5,a8,a3,a3,a5)", "Z", "rmax", "Ne", "a", "p", "Nq", &
    "DOFs"

i = 31
! change p for p-conv study. p must be less than 31.
!if (dirac_int == 1 .and. i > 23) then
!    exit
!end if
p = i
allocate(xq(Nq, Ne))
call total_energy(Z, rmax, Ne, a, p, Nq, DOFs, alpha_int, dirac_int, &
    c, potential_type, Lmin, alpha_j, alpha, energies, eigfn, xq)
Etot = sum(energies)
print "(i3, f6.1, i5, f8.1, i3, i3, i5, f22.12)", Z, rmax, Ne, a, p, Nq, &
    DOFs, Etot
print *
print *, "Comparison of calculated and reference energies"
print *
print *, "Total energy:"
print "(a20,a20,a10)", "E", "E_ref", "error"
err = abs(Etot - Etot_ref)
print "(f20.12, f20.12, es10.2)", Etot, Etot_ref, err
if ( .not. (err < 1e-9_dp)) then
   error stop 'assert failed'
end if
print *, "Eigenvalues:"
print "(a4,a20,a20,a10)", "n", "E", "E_ref", "error"
do i = 1, size(energies)
    err = abs(energies(i) - energies_ref(i))
    print "(i4, f20.12, f20.12, es10.2)", i, energies(i), energies_ref(i), err
    if ( .not. (err < 5e-10_dp)) then
       error stop 'assert failed'
    end if
end do

print *, "Eigenfunctions saved in data_coulomb_schroed.txt"
open(newunit=u, file="data_coulomb_schroed.txt", status="replace")
write(u, *) xq
do i = 1, size(energies)
    write(u, *) eigfn(:,:,i)
end do
close(u)

end program
