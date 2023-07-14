module graphs

use types, only: dp
use mesh, only: meshexp
use schroed_glob, only: solve_schroed
use dirac, only: solve_dirac
use feutils, only: get_parent_quad_pts_wts
use string_utils, only: str
implicit none

contains

subroutine run_convergence(study_type, dirac_int, p_or_Ne, directory)
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

! The directory where to save the output files
character(len=*), intent(in) :: directory

integer :: p, Ne, Nq, Z, u, i, DOFs
real(dp) :: rmax, a, Etot
real(dp) :: optim_a(2:7)
integer, allocatable :: Nes(:), rmax_values(:)
real(dp), allocatable :: energies(:)
character(len=:), allocatable :: filename
Z = 92
rmax = 50
a = 200
Ne = 4
Nq = 64
p = 25

if (study_type == 2) then
    p = p_or_Ne
else
    Ne = p_or_Ne
end if
filename = str(p_or_Ne)


optim_a = [58.985048682607555, 163.13530060338942, 340.82602039608668, &
            444.68894311026423, 591.72463734788732, 596.61404750045062]

if (Ne >= 2 .and. Ne <= 7) then
    a = optim_a(Ne)
end if
if (dirac_int == 1 .and. study_type == 0) then
    a = 600
    rmax = 30
end if
if (dirac_int == 1 .and. study_type == 1) then
    a = 600
    p = 25
end if
if (dirac_int == 1 .and. study_type == 2) then
    a = 600
    rmax = 30
end if

if (study_type == 0) then
    filename = "conv_" // trim(filename) // ".txt"
else if (study_type == 1) then
    filename = "rmax_" // trim(filename) // ".txt"
else
    filename = "ne_" // trim(filename) // ".txt"
end if

if (dirac_int == 1) then
    filename = "dft_dirac_" // trim(filename)
else
    filename = "dft_schroed_" // trim(filename)
end if

filename = directory // "/" // filename

open(newunit=u, file=filename, status="replace")
print "(a3,a6,a5,a8,a3,a3,a5,a22)", "Z", "rmax", "Ne", "a", "p", "Nq", &
    "DOFs", "Etot"


if (study_type == 0) then
    do i = 7, 31
        ! change p for p-conv study. p must be less than 31.
        p = i
        call total_energy(Z, rmax, Ne, a, p, Nq, DOFs, Etot)
        print "(i3, f6.1, i5, f8.1, i3, i3, i5, f22.12)", Z, rmax, Ne, a, p, Nq, &
            DOFs, Etot
        write(u,*) Z, rmax, Ne, a, p, Nq, DOFs, Etot
    end do
else if (study_type == 1) then
    allocate(rmax_values(5))
    rmax_values = [5, 10, 20, 30, 40]
    do i = 1, size(rmax_values)
        ! change rmax for rmax study.
        rmax = rmax_values(i)
        call total_energy(Z, rmax, Ne, a, p, Nq, DOFs, Etot)
        print "(i3, f6.1, i5, f8.1, i3, i3, i5, 40f22.12)", Z, rmax, Ne, a, p, Nq, &
            DOFs, Etot, energies
        write(u,*) Z, rmax, Ne, a, p, Nq, DOFs, Etot, energies
    end do
else if (study_type == 2) then
    allocate(Nes(3))
    Nes = [6,10,20]
    do i = 1, size(Nes)
        Ne = Nes(i)
        call total_energy(Z, rmax, Ne, a, p, Nq, DOFs, Etot)
        print "(i3, f6.1, i5, f8.1, i3, i3, i5, f22.12)", Z, rmax, Ne, a, p, Nq, &
            DOFs, Etot
        write(u,*) Z, rmax, Ne, a, p, Nq, DOFs, Etot, energies
    end do
end if
close(u)

contains

    subroutine total_energy(Z, rmax, Ne, a, p, Nq, DOFs, Etot)
    integer, intent(in) :: Z, Ne, Nq, p
    real(dp), intent(in) :: rmax, a
    integer, intent(out) :: DOFs
    real(dp), intent(out) :: Etot
    real(dp), allocatable :: xe(:), xiq(:), wtq(:), V(:,:)
    real(dp) :: rmin, r0
    rmin = 0
    allocate(xe(Ne+1), xiq(Nq), wtq(Nq), V(Nq, Ne))
    if (dirac_int == 1 .and. study_type == 0) then
        r0 = 0.005_dp
        xe(1) = rmin
        xe(2:) = meshexp(r0, rmax, a, Ne-1)
    else if (dirac_int == 1 .and. study_type == 1) then
        r0 = 0.005_dp
        xe(1) = rmin
        xe(2:) = meshexp(r0, rmax, a, Ne-1)
    else if (dirac_int == 1 .and. study_type == 2) then
        r0 = 0.005_dp
        xe(1) = rmin
        xe(2:) = meshexp(r0, rmax, a, Ne-1)
    else
        xe = meshexp(rmin, rmax, a, Ne)
    end if
    call get_parent_quad_pts_wts(1, Nq, xiq, wtq)
    if (dirac_int == 1) then
        call solve_dirac(Z, p, xiq, wtq, xe, 1e-9_dp, energies, Etot, V, DOFs)
    else
        call solve_schroed(Z, p, xiq, wtq, xe, 1e-9_dp, energies, Etot, V, DOFs)
    end if
    end subroutine

end subroutine

end module
