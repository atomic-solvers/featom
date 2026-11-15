program test_dft_schroed_cation

use types, only: dp
use mesh, only: meshexp
use schroed_glob, only: solve_schroed
use feutils, only: get_parent_quad_pts_wts
implicit none

real(dp), allocatable :: xe(:)              ! element coordinates
real(dp), allocatable :: xiq(:), wtq(:)     ! quadrature points and weights
integer :: p, Ne, Nq, Z, i, DOFs
real(dp) :: rmin, rmax, a, err, Etot
real(dp), allocatable :: energies(:), V(:,:)
integer :: no(1), lo(1)
real(dp) :: fo(1)
! Li+ 1s2, NIST Atomic Reference Data (DB 141), LDA.
real(dp), parameter :: Etot_ref = -7.142818_dp
real(dp), parameter :: energies_ref(*) = [ &
                -2.190276_dp ]


Z = 3
no = [ 1 ]
lo = [ 0 ]
fo = [ 2 ]
rmin = 0
rmax = 50
a = 200
Ne = 4
Nq = 53
p = 26

allocate(xe(Ne+1), xiq(Nq), wtq(Nq), V(Nq, Ne))
xe = meshexp(rmin, rmax, a, Ne)
call get_parent_quad_pts_wts(1, Nq, xiq, wtq)

call solve_schroed(Z, p, xiq, wtq, xe, 1e-8_dp, energies, Etot, V, DOFs, &
        no, lo, fo)
if ( .not. (size(energies) == size(energies_ref))) then
   error stop 'assert failed'
end if
print *, "Comparison of calculated and reference energies"
print *
print *, "Total energy:"
print "(a16,a16,a10)", "E", "E_ref", "error"
err = abs(Etot - Etot_ref)
print "(f16.8, f16.8, es10.2)", Etot, Etot_ref, err
if ( .not. (err < 1e-6_dp)) then
   error stop 'assert failed'
end if
print *
print *, "Eigenvalues:"
print "(a4,a16,a16,a10)", "n", "E", "E_ref", "error"
do i = 1, size(energies)
    err = abs(energies(i) - energies_ref(i))
    print "(i4, f16.8, f16.8, es10.2)", i, energies(i), energies_ref(i), err
    if ( .not. (err < 1e-6_dp)) then
       error stop 'assert failed'
    end if
end do
end program
