program test_dft_dirac

use types, only: dp
use mesh, only: meshexp
use dirac, only: solve_dirac
use feutils, only: get_parent_quad_pts_wts
implicit none

real(dp), allocatable :: xe(:)        ! element coordinates
real(dp), allocatable :: xiq(:), wtq(:)     ! quadrature points and weights
integer :: p, Ne, Nq, Z, i, DOFs
real(dp) :: rmin, rmax, a, c, err, Etot, r0
real(dp), allocatable :: energies(:), V(:,:)
real(dp), parameter :: Etot_ref = -28001.1323254868_dp
real(dp), parameter :: energies_ref(*) = [ &
                -4223.4190204552_dp, &
                -789.4897823303_dp, &
                -761.3744759730_dp, &
                -622.8480945649_dp, &
                -199.4298056450_dp, &
                -186.6637131249_dp, &
                -154.7010266741_dp, &
                -134.5411802896_dp, &
                -128.0166573820_dp, &
                -50.7889480646_dp, &
                -45.0371712884_dp, &
                -36.6886104859_dp, &
                -27.5293062430_dp, &
                -25.9854289064_dp, &
                -13.8895142333_dp, &
                -13.4854696912_dp, &
                -11.2955870987_dp, &
                -9.0579642498_dp, &
                -7.0692956350_dp, &
                -3.7974162278_dp, &
                -3.5012171832_dp, &
                -0.1467883850_dp, &
                -0.1160471651_dp, &
                -1.7480399541_dp, &
                -1.1011189998_dp, &
                -0.7757841787_dp, &
                -0.1030408153_dp, &
                -0.0848020246_dp, &
                -0.1609472826_dp ]


Z = 92
c = 137.0359895_dp
rmin = 0
r0 = 0.005_dp
rmax = 30
a = 600
Ne = 6
Nq = 64
p = 25

allocate(xe(Ne+1), xiq(Nq), wtq(Nq), V(Nq, Ne))
xe(1) = rmin
xe(2:) = meshexp(r0, rmax, a, Ne-1)
call get_parent_quad_pts_wts(1, Nq, xiq, wtq)

call solve_dirac(Z, p, xiq, wtq, xe, 1e-9_dp, energies, Etot, V, DOFs)

if ( .not. (size(energies) == size(energies_ref))) then
   error stop 'assert failed'
end if
print *, "Comparison of calculated and reference energies"
print *
print *, "Total energy:"
print "(a16,a16,a10)", "E", "E_ref", "error"
err = abs(Etot - Etot_ref)
print "(f16.8, f16.8, es10.2)", Etot, Etot_ref, err
if ( .not. (err < 2e-8_dp)) then
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
