program test_coulomb_dirac
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
real(dp) :: Etot_ref
real(dp), allocatable  :: energies_ref(:)
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
integer :: n, l, relat, relat_max

Z = 92
rmin = 0
rmax = 50
a = 200
Ne = 4
Nq = 53
p = 26

i = 7
!call run_convergence_potential(0, 1, i, &
!    0, -1, ".")


study_type = 0
dirac_int = 1
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

Lmax=5
Lmin=-6

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

i = 23
! change p for p-conv study. p must be less than 31.
!if (dirac_int == 1 .and. i > 23) then
!    exit
!end if
p = i
allocate(xq(Nq, Ne))
call total_energy(Z, rmax, Ne, a, p, Nq, DOFs, alpha_int, dirac_int, &
    c, potential_type, Lmin, alpha_j, alpha, energies, eigfn, xq)
Etot = sum(energies)

allocate(energies_ref(size(energies)))
i = 1
outer: do n = 1, 7
    do l = 0, n-1
        if (l == 0) then
            relat_max = 2
        else
            relat_max = 3
        end if
        do relat = 2, relat_max
            energies_ref(i) = E_nl(c, n, l, Z, relat)
            if (i == size(energies_ref)) exit outer
            i = i + 1
        end do
    end do
end do outer
Etot_ref = sum(energies_ref)

print "(i3, f6.1, i5, f8.1, i3, i3, i5, f22.12)", Z, rmax, Ne, a, p, Nq, &
    DOFs, Etot
print *
print *, "Comparison of calculated and reference energies"
print *
print *, "Total energy:"
print "(a20,a20,a10)", "E", "E_ref", "error"
err = abs(Etot - Etot_ref)
print "(f20.12, f20.12, es10.2)", Etot, Etot_ref, err
if ( .not. (err < 5e-9_dp)) then
   error stop 'assert failed'
end if

print *, "Eigenvalues:"
print "(a4,a20,a20,a10)", "n", "E", "E_ref", "error"
do i = 1, size(energies)
    err = abs(energies(i) - energies_ref(i))
    print "(i4, f20.12, f20.12, es10.2)", i, energies(i), energies_ref(i), err
    if ( .not. (err < 5e-9_dp)) then
       error stop 'assert failed'
    end if
end do

print *, "Eigenfunctions saved in data_coulomb_dirac.txt"
open(newunit=u, file="data_coulomb_dirac.txt", status="replace")
write(u, *) xq
do i = 1, size(energies)
    write(u, *) eigfn(:,:,i)
end do
close(u)

contains

    real(dp) function E_nl(c, n, l, Z, relat)
    ! Calculates exact energy for the radial Schroedinger/Dirac equations
    real(dp), intent(in) :: c ! speed of light in atomic units
    integer, intent(in) :: n, l, Z, relat
    ! quantum numbers (n, l), atomic number (z)
    ! relat == 0 ... Schroedinger equation
    ! relat == 2 ... Dirac equation, spin up
    ! relat == 3 ... Dirac equation, spin down

    integer :: kappa
    real(dp) :: beta
    if (.not. (l >= 0)) error stop "'l' must be positive or zero"
    if (.not. (n > l)) error stop "'n' must be greater than 'l'"
    if (l == 0 .and. relat == 3) error stop "Spin must be up for l==0."
    if (relat == 0) then
        E_nl = - Z**2 / (2.0_dp * n**2)
    else
        if (relat == 2) then
            kappa = -l - 1
        else
            kappa = l
        end if
        beta = sqrt(kappa**2 - (Z/c)**2)
        E_nl = c**2/sqrt(1 + (Z/c)**2/(n - abs(kappa) + beta)**2) - c**2
    end if
    end function

end program
