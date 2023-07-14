module schroed_glob

use types, only: dp
use mesh, only: meshexp
use feutils, only: define_connect, get_quad_pts, get_parent_quad_pts_wts, &
        get_parent_nodes, phih, dphih, c2fullc2, fe2quad_core, get_nodes, &
        integrate, proj_fn, phih_array, dphih_array
use linalg, only: eigh
use fe, only: assemble_radial_H, assemble_radial_S, assemble_radial_H_setup, &
        assemble_radial_H_complete
use constants, only: pi
use hartree_screening, only: assemble_poisson_A, hartree_potential3
use xc, only: xc_vwn3
use mixings, only: mixing_linear, mixing_pulay
use states, only: get_atomic_states_nonrel_focc, get_atomic_states_rel_focc, &
        nlsf2focc, get_atomic_states_rel, nlf2focc, get_atomic_states_nonrel
use energies, only: thomas_fermi_potential
use iso_c_binding, only: c_double, c_int
use solvers, only: solve_eig_irange, solve_sym_setup
use lapack, only: dsytrf, dsytrs
implicit none
private
public solve_schroed!, csolve_schroed, total_energy

contains

subroutine solve_schroed(Z, p, xiq, wtq, xe, eps, energies, Etot, V, DOFs)

integer, intent(in) :: Z, p
real(dp), intent(in) :: xe(:)        ! element coordinates
real(dp), intent(in) :: xiq(:)       ! quadrature points
real(dp), intent(in) :: wtq(:)       ! quadrature weights
real(dp), intent(in) :: eps
real(dp), allocatable, intent(out) :: energies(:)
real(dp), intent(out) :: Etot
real(dp), intent(out) :: V(:,:)       ! SCF potential
integer, intent(out) :: DOFs
integer :: n, Nq
real(dp), allocatable :: H(:,:), S(:), D(:,:), lam(:)

real(dp), allocatable :: xin(:)       ! parent basis nodes
integer, allocatable :: ib(:, :)       ! basis connectivity: ib(i,j) = index of
   ! basis function associated with local basis function i of element j. 0 = no
   ! associated basis fn.
integer, allocatable :: in(:, :)
real(dp), allocatable :: xq(:, :), fullc(:), uq(:,:), rho(:,:), Vee(:,:), &
    Vxc(:,:), exc(:,:), Vin(:,:), Vout(:,:), xn(:), un(:), xq1(:,:), Am_p(:,:), &
    bv_p(:), Hl(:,:,:)
integer, allocatable :: ipiv(:)
real(dp), allocatable :: phihq(:,:)   ! parent basis at quadrature points
real(dp), allocatable :: dphihq(:,:)  ! parent basis derivative at quadrature points

real(dp), allocatable :: xin_p(:)       ! parent basis nodes for poisson grid
real(dp), allocatable :: xn_p(:)
real(dp), allocatable :: phihq_p(:,:)   ! parent basis at quadrature points
real(dp), allocatable :: dphihq_p(:,:)  ! parent basis derivative at quadrature points
integer, allocatable :: in_p(:, :), ib_p(:, :)

integer :: Ne, Nb, Nn, Nb_p, Nn_p, pp
real(dp) :: rmin, rmax
integer :: l, i, j, Lmax, al, eimin, eimax
real(dp) :: c
real(dp), allocatable :: focc(:,:), focc_idx_r(:,:), tmp(:)
real(dp) :: scf_alpha, scf_L2_eps, scf_eig_eps
integer, parameter :: mixing_scheme_linear = 1, mixing_scheme_pulay = 3
integer :: mixing_scheme
integer, allocatable :: no(:), lo(:), focc_idx(:,:), eirange(:,:)
real(dp), allocatable :: fo_idx(:), fo(:)
real(dp) :: T_s, E_ee, E_en, EE_xc
real(dp) :: xiq_lob(p+1), wtq_lob(p+1)


integer :: nband, scf_max_iter, iter

Nq = size(xiq)
c = 137.0359895_dp
rmin = 0
rmax = 50
Ne = size(xe)-1
mixing_scheme = mixing_scheme_pulay

call get_parent_quad_pts_wts(2, p+1, xiq_lob, wtq_lob)


call get_atomic_states_nonrel_focc(Z, focc)
call get_atomic_states_nonrel(Z, no, lo, fo)
allocate(fo_idx(size(fo)))
do i = 1, size(fo_idx)
    fo_idx(i) = i
end do
call nlf2focc(no, lo, fo_idx, focc_idx_r)
allocate(focc_idx(size(focc_idx_r,1),lbound(focc,2):ubound(focc,2)))
focc_idx = int(focc_idx_r)

Lmax = ubound(focc,2)

Nn = Ne*p+1

allocate(xin(p+1))
call get_parent_nodes(2, p, xin)
allocate(in(p+1, Ne), ib(p+1, Ne))
call define_connect(1, 1, Ne, p, in, ib)
Nb = maxval(ib)
if ( .not. (Nn == maxval(in)) ) then
   error stop 'Wrong size for Nn'
end if
DOFs = Nb
allocate(xq(Nq, Ne), xq1(Nq, 1), xn(Nn), un(Nn))
call get_quad_pts(xe, xiq, xq)
call get_nodes(xe, xin, xn)

pp = 2*p
Nn_p = Ne*pp+1
allocate(xin_p(pp+1))
call get_parent_nodes(2, pp, xin_p)
allocate(in_p(pp+1, Ne), ib_p(pp+1, Ne))
call define_connect(2, 1, Ne, pp, in_p, ib_p)
Nb_p = maxval(ib_p)
if ( .not. (Nn_p == maxval(in_p)) ) then
   error stop 'Wrong size for Nn_p'
end if
allocate(xn_p(Nn_p))
call get_nodes(xe, xin_p, xn_p)


allocate(phihq(Nq, p+1))
allocate(dphihq(Nq, p+1))
allocate(phihq_p(Nq, pp+1))
allocate(dphihq_p(Nq, pp+1))
! tabulate parent basis at quadrature points
do al = 1, p+1
    call phih_array(xin, al, xiq, phihq(:, al))
    call dphih_array(xin, al, xiq, dphihq(:, al))
end do
do al = 1, pp+1
    call phih_array(xin_p, al, xiq, phihq_p(:, al))
    call dphih_array(xin_p, al, xiq, dphihq_p(:, al))
end do

n = Nb
allocate(H(n, n), S(n))
allocate(D(n, n), lam(n), fullc(Nn), uq(Nq,Ne), rho(Nq,Ne), Vee(Nq,Ne), &
    Vxc(Nq,Ne), exc(Nq,Ne), Vin(Nq,Ne), Vout(Nq,Ne))

nband = count(focc > 0)
scf_max_iter = 100
scf_alpha = 0.7_dp
scf_L2_eps = 1e-4_dp
scf_eig_eps = eps

allocate(tmp(Nq*Ne))
allocate(energies(nband))
Vin = reshape(thomas_fermi_potential(reshape(xq, [Nq*Ne]), Z), [Nq, Ne]) + &
    Z / xq
iter = 0

call assemble_radial_S(xin, xe, ib, wtq_lob, S)
do i = 1, size(S)
    S(i) = 1/sqrt(S(i))
end do

allocate(Am_p(Nb_p, Nb_p), bv_p(Nb_p), ipiv(Nb_p))
Am_p = 0
call assemble_poisson_A(xin_p, xe, ib_p, xiq, wtq, dphihq_p, Am_p)
call solve_sym_setup(Am_p, ipiv)

allocate(Hl(n,n,0:Lmax))
call assemble_radial_H_setup(0, Lmax, xin, xe, ib, xiq, wtq, phihq, dphihq, Hl)

allocate(eirange(2,0:Lmax))

do l = 0, Lmax
    eimax = 0
    do i = 1, size(focc,1)
        if (focc(i,l) < tiny(1._dp)) cycle
        eimax = i
    end do
    eirange(2, l) = eimax

    eimin = 1
    do i = 1, size(focc,1)
        if (focc(i,l) >= tiny(1._dp)) exit
        eimin = i
    end do
    eirange(1, l) = eimin
end do

select case (mixing_scheme)
    case (mixing_scheme_linear)
        call mixing_linear &
            (Ffunc, integral, reshape(Vin, [Nq*Ne]), &
            nband, scf_max_iter, scf_alpha, scf_L2_eps, scf_eig_eps, tmp)
    case (mixing_scheme_pulay)
        call mixing_pulay &
            (Ffunc, integral, matvec, matmat, reshape(Vin, [Nq*Ne]), &
            nband, scf_max_iter, scf_alpha, scf_L2_eps, scf_eig_eps, tmp, 5, 3)
    case default
        error stop "Type of mixing not implemented."
end select

contains

    subroutine Ffunc(x, y, eng)
    ! Converge Vee+Vxc only (the other components are constant)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: y(:), eng(:)
    integer :: idx
    iter = iter + 1
    Vin = reshape(x, shape(Vin))
    rho = 0
    idx = 0
    V = Vin - Z/xq
    do l = 0, Lmax

        call assemble_radial_H_complete(V, xin, xe, ib, xiq, wtq, phihq, Hl(:,:,l), H)

        do concurrent (i = 1:size(S), j = 1:size(S), i>=j)
            H(i, j) = H(i, j)*S(i)*S(j)
        end do

        eimin = eirange(1, l)
        eimax = eirange(2, l)

        call solve_eig_irange(H, eimin, eimax, lam, D)

        do i = 1, size(S)
            D(i, eimin:eimax) = D(i, eimin:eimax)*S(i)
        end do

        do i = eimin,eimax
            if (focc(i,l) < tiny(1._dp)) cycle

            call c2fullc2(in, ib, D(:Nb,i), fullc)
            call fe2quad_core(xe, xin, in, fullc, phihq, uq)
            rho = rho - focc(i,l)*(uq/xq)**2 / (4*pi)

            idx = idx + 1
            eng(focc_idx(i,l)) = lam(i)
        end do
    end do
    if ( .not. (size(eng) == idx) ) then
       error stop 'Energy size mismatch'
    end if
    energies = eng

    !print *, "Energies:"
    !do i = 1, size(energies)
    !    print *, i, no(i), lo(i), energies(i)
    !end do

    Vee = hartree_potential3(real(Z, dp), pp, xe, xin_p, in_p, ib_p, xiq, wtq, phihq_p, Am_p, ipiv, -4*pi*rho)
    call xc_vwn3(size(rho), -rho, .FALSE., c, exc, Vxc)
    call total_energy(xe, wtq, fo, energies, Vin-Z/xq, Vee, -Z/xq, exc, &
        xq, -rho, T_s, E_ee, E_en, EE_xc, Etot)
    Vout = Vee + Vxc ! This term is added later: -Z/xq
    y = reshape(Vout, shape(y))
    end subroutine

    real(dp) function integral(x)
    ! Computes the integral of the vector 'x'
    real(dp), intent(in) :: x(:)
    integral = integrate(xe, wtq, reshape(x, shape(uq)))
    end function

    function matvec(A, b) result(r)
    real(dp), intent(in) :: A(:,:), b(:)
    real(dp) :: r(size(A,1))
    r = matmul(A, b)
    end function

    function matmat(A, B) result(r)
    real(dp), intent(in) :: A(:,:), B(:,:)
    real(dp) :: r(size(A,1), size(B,2))
    r = matmul(A, B)
    end function

end subroutine


subroutine total_energy(xe, wtq, fo, ks_energies, V_in, V_h, V_coulomb, e_xc, &
    R, n, T_s, E_ee, E_en, EE_xc, Etot)
! This is a variational, quadratically convergent form of total energy
real(dp), intent(in) :: xe(:), wtq(:)
real(dp), intent(in) :: R(:,:) ! Function 'r' on quadrature grid
real(dp), intent(in) :: fo(:), ks_energies(:) ! occupations, energies
real(dp), intent(in) :: V_in(:,:) ! Total input effective potential
real(dp), intent(in) :: V_h(:,:) ! Hartree energy, solution of Poiss. eq.
real(dp), intent(in) :: V_coulomb(:,:) ! Coulomb inter. -Z/r  (negative)
real(dp), intent(in) :: e_xc(:,:) ! XC density
real(dp), intent(in) :: n(:,:) ! number density (positive)
real(dp), intent(out) :: Etot ! Total energy
real(dp), intent(out) :: T_s, E_ee, E_en, EE_xc ! Parts of the total energy

real(dp) :: rho(size(n,1), size(n,2))
real(dp) :: E_c, E_band
rho = -n

E_band = sum(fo * ks_energies)
T_s = E_band + 4*pi * integrate(xe, wtq, V_in * rho * R**2)

E_ee = -2*pi * integrate(xe, wtq, V_h * rho * R**2)
E_en =  4*pi * integrate(xe, wtq, (-V_coulomb) * rho * R**2)
E_c = E_ee + E_en

EE_xc = -4*pi * integrate(xe, wtq, e_xc * rho * R**2)

Etot = T_s + E_c + EE_xc
end subroutine


subroutine csolve_schroed(Z, p, xiq, wtq, xe, eps, Nq, Ne, &
        nenergies, energies, Etot, V) bind(c)
integer(c_int), intent(in) :: Z, p, Nq, Ne, nenergies
real(c_double), intent(in) :: xe(Ne+1)        ! element coordinates
real(c_double), intent(in) :: xiq(Nq)       ! quadrature points
real(c_double), intent(in) :: wtq(Nq)       ! quadrature weights
real(c_double), intent(in) :: eps
real(c_double), intent(out) :: energies(nenergies)
real(c_double), intent(out) :: Etot
real(c_double), intent(out) :: V(Nq, Ne)    ! SCF potential
integer :: DOFs
real(dp), allocatable :: energies2(:)

call solve_schroed(Z, p, xiq, wtq, xe, eps, energies2, Etot, V, DOFs)
energies = energies2
end subroutine

end module
