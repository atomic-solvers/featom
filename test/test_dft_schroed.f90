program test_dft_schroed

use types, only: dp
use mesh, only: meshexp
use schroed_glob, only: solve_schroed
use dirac, only: solve_dirac
use feutils, only: get_parent_quad_pts_wts
implicit none

real(dp), allocatable :: xe(:)              ! element coordinates
real(dp), allocatable :: xiq(:), wtq(:)     ! quadrature points and weights
integer :: p, Ne, Nq, Z, i, DOFs
real(dp) :: rmin, rmax, a, err, Etot
real(dp), allocatable :: energies(:), V(:,:)
real(dp), parameter :: Etot_ref = -25658.4178888534_dp
real(dp), parameter :: energies_ref(*) = [ &
                -3689.3551398369_dp, &
                -639.7787280866_dp, &
                -619.1085501807_dp, &
                -161.1180732100_dp, &
                -150.9789801633_dp, &
                -131.9773582831_dp, &
                -40.5280842452_dp, &
                -35.8533208325_dp, &
                -27.1232122996_dp, &
                -15.0274600691_dp, &
                -8.8240894015_dp, &
                -7.0180922045_dp, &
                -3.8661751349_dp, &
                -0.3665433531_dp, &
                -1.3259763180_dp, &
                -0.8225379709_dp, &
                -0.1431901813_dp, &
                -0.1309478622_dp ]


Z = 92
rmin = 0
rmax = 50
a = 200
Ne = 4
Nq = 53
p = 26

allocate(xe(Ne+1), xiq(Nq), wtq(Nq), V(Nq, Ne))
xe = meshexp(rmin, rmax, a, Ne)
call get_parent_quad_pts_wts(1, Nq, xiq, wtq)

call solve_schroed(Z, p, xiq, wtq, xe, 1e-8_dp, energies, Etot, V, DOFs)
if ( .not. (size(energies) == size(energies_ref))) then
   error stop 'assert failed'
end if
print *, "Comparison of calculated and reference energies"
print *
print *, "Total energy:"
print "(a16,a16,a10)", "E", "E_ref", "error"
err = abs(Etot - Etot_ref)
print "(f16.8, f16.8, es10.2)", Etot, Etot_ref, err
if ( .not. (err < 1e-8_dp)) then
   error stop 'assert failed'
end if
print *
print *, "Eigenvalues:"
print "(a4,a16,a16,a10)", "n", "E", "E_ref", "error"
do i = 1, size(energies)
    err = abs(energies(i) - energies_ref(i))
    print "(i4, f16.8, f16.8, es10.2)", i, energies(i), energies_ref(i), err
    if ( .not. (err < 1e-8_dp)) then
       error stop 'assert failed'
    end if
end do
end program
