program test_harmonic_dirac
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
real(dp), allocatable :: energies_ref(:)
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
real(dp) :: E_exact(7, -7:7) ! E_exact(n, kappa)
E_exact = 0
E_exact(7, -7) = 7.4996755269_dp
E_exact(6:7, -6) = [ &
    6.4997620499_dp, &
    8.4994625812_dp &
]
E_exact(5:7, -5) = [ &
    5.4998352632_dp, &
    7.4995757160_dp, &
    9.4992363376_dp &
]
E_exact(4:7, -4) = [ &
    4.4998951661_dp, &
    6.4996755428_dp, &
    8.4993760822_dp, &
    10.4989967965_dp &
]
E_exact(3:7, -3) = [ &
    3.4999417582_dp, &
    5.4997620612_dp, &
    7.4995025208_dp, &
    9.4991631492_dp, &
    11.4987439585_dp &
]
E_exact(2:7, -2) = [ &
    2.4999750389_dp, &
    4.4998352706_dp, &
    6.4996156527_dp, &
    8.4993161976_dp, &
    10.4989369173_dp, &
    12.4984778242_dp &
]
E_exact(1:7, -1) = [ &
    1.4999950078_dp, &
    3.4998951705_dp, &
    5.4997154776_dp, &
    7.4994559413_dp, &
    9.4991165738_dp, &
    11.4986973873_dp, &
    13.4981983940_dp &
]
E_exact(2:7, 1) = [ &
    2.4999351051_dp, &
    4.4997953424_dp, &
    6.4995757249_dp, &
    8.4992762722_dp, &
    10.4988969952_dp, &
    12.4984379048_dp &
]
E_exact(3:7, 2) = [ &
    3.4998752033_dp, &
    5.4996955116_dp, &
    7.4994359765_dp, &
    9.4990966102_dp, &
    11.4986774249_dp &
]
E_exact(4:7, 3) = [ &
    4.4998019930_dp, &
    6.4995823772_dp, &
    8.4992829240_dp, &
    10.4989036457_dp &
]
E_exact(5:7, 4) = [ &
    5.4997154739_dp, &
    7.4994559363_dp, &
    9.4991165675_dp &
]
E_exact(6:7, 5) = [ &
    6.4996156467_dp, &
    8.4993161897_dp &
]
E_exact(7, 6) = 7.4995025118_dp


Z = 92
rmin = 0
rmax = 50
a = 200
Ne = 4
Nq = 53
p = 26

i = 7
!call run_convergence_potential(0, 1, i, &
!    1, 0, ".")

study_type = 0
dirac_int = 1
p_or_Ne = i
potential_type = 1
alpha_int = 0
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
            if (relat == 2) then
                kappa = -l - 1
            else
                kappa = l
            end if
            energies_ref(i) = E_exact(n, kappa)
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
if ( .not. (err < 1e-8_dp)) then
   error stop 'assert failed'
end if
print *, "Eigenvalues:"
print "(a4,a15,a15,a10)", "n", "E", "E_ref", "error"
do i = 1, size(energies)
    err = abs(energies(i) - energies_ref(i))
    print "(i4, f15.8, f15.8, es10.2)", i, energies(i), energies_ref(i), err
    if ( .not. (err < 1e-8_dp)) then
       error stop 'assert failed'
    end if
end do

print *, "Eigenfunctions saved in data_harmonic_dirac.txt"
open(newunit=u, file="data_harmonic_dirac.txt", status="replace")
write(u, *) xq
do i = 1, size(energies)
    write(u, *) eigfn(:,:,i)
end do
close(u)

end program
