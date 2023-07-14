module dirac

use types, only: dp
use mesh, only: meshexp
use feutils, only: define_connect, get_quad_pts, get_parent_quad_pts_wts, &
        get_parent_nodes, phih, dphih, c2fullc2, fe2quad_core, get_nodes, &
        integrate, proj_fn, phih_array, integrate2, fe2quad
use linalg, only: eigh
use jacobi, only: gauss_jacobi
use fe, only: assemble_radial_SH, assemble_radial_dirac_SH
use constants, only: pi, c => c_1986
use hartree_screening, only: hartree_potential_gj
use xc, only: xc_vwn3
use mixings, only: mixing_linear, mixing_pulay
use states, only: get_atomic_states_nonrel_focc, get_atomic_states_rel_focc, &
    nlsf2focc, get_atomic_states_rel, nlf2focc, get_atomic_states_nonrel
use energies, only: thomas_fermi_potential
use iso_c_binding, only: c_double, c_int
implicit none
private
public solve_dirac, csolve_dirac

contains

subroutine solve_dirac(Z, p, xiq, wtq, xe, eps, energies, Etot, V, DOFs)

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
real(dp), allocatable :: H(:,:), S(:,:), D(:,:), lam(:), lam_tmp(:)
real(dp), allocatable :: xiq_gj(:,:)      ! quadrature points first element
real(dp), allocatable :: wtq_gj(:,:)      ! quadrature weights first element
real(dp), allocatable :: xin(:)       ! parent basis nodes
integer, allocatable :: ib(:, :)       ! basis connectivity: ib(i,j) = index of
   ! basis function associated with local basis function i of element j. 0 = no
   ! associated basis fn.
integer, allocatable :: in(:, :)
real(dp), allocatable :: xq(:, :), fullc(:), uq(:,:), rho(:,:), Vee(:,:), &
    Vxc(:,:), exc(:,:), Vin(:,:), Vout(:,:), xn(:), un(:), xq1(:,:), &
    phihq(:,:), xq2(:,:), rho0(:,:), rho1(:,:)
integer :: Ne, Nb, Nn
real(dp) :: rmin, rmax, asympt
real(dp), allocatable :: alpha_j(:), alpha(:)
integer :: kappa, i, Lmin, Lmax, al
real(dp), allocatable :: focc(:,:), focc_idx_r(:,:), tmp(:)
real(dp) :: scf_alpha, scf_L2_eps, scf_eig_eps
integer, parameter :: mixing_scheme_linear = 1, mixing_scheme_pulay = 3
integer :: mixing_scheme
integer, allocatable :: no(:), lo(:), so(:), focc_idx(:,:)
real(dp), allocatable :: fo_idx(:), fo(:)
real(dp) :: T_s, E_ee, E_en, EE_xc


integer :: nband, scf_max_iter, iter

Nq = size(xiq)
rmin = 0
rmax = 50
Ne = size(xe)-1
mixing_scheme = mixing_scheme_pulay

call get_atomic_states_rel_focc(Z, focc)
call get_atomic_states_rel(Z, no, lo, so, fo)
allocate(fo_idx(size(fo)))
do i = 1, size(fo_idx)
    fo_idx(i) = i
end do
call nlsf2focc(no, lo, so, fo_idx, focc_idx_r)
allocate(focc_idx(size(focc_idx_r,1),lbound(focc,2):ubound(focc,2)))
focc_idx = int(focc_idx_r)

Lmax = ubound(focc,2)
Lmin = lbound(focc,2)

allocate(alpha(Lmin:Lmax), alpha_j(Lmin:Lmax))

do kappa = Lmin, Lmax
    if (kappa == 0) cycle
    ! asymptotic at r = 0
    asympt = sqrt(kappa**2 - Z**2 / c**2)
    ! alpha can be [0, asympt]
    alpha(kappa) = asympt
    ! power of r for Gauss-Jacobi quadrature
    alpha_j(kappa) = 2*alpha(kappa) - 2
end do

Nn = Ne*p+1

allocate(xin(p+1))
call get_parent_nodes(2, p, xin)
allocate(in(p+1, Ne), ib(p+1, Ne))
call define_connect(2, 1, Ne, p, in, ib)
Nb = maxval(ib)
if ( .not. (Nn == maxval(in)) ) then
   error stop 'Size mismatch'
end if
allocate(xq(Nq, Ne), xq1(Nq, 1), xn(Nn), un(Nn))
call get_quad_pts(xe, xiq, xq)
call get_nodes(xe, xin, xn)

allocate(xiq_gj(Nq, Lmin:Lmax), wtq_gj(Nq, Lmin:Lmax))
if (any(alpha_j > - 1)) then
    !print *, "Using Gauss-Jacobi quadrature"
    do kappa = Lmin, Lmax
        if (kappa == 0) cycle
        if (alpha_j(kappa) > -1) then
            call gauss_jacobi(Nq, 0.0_dp, alpha_j(kappa), xiq_gj(:, kappa), wtq_gj(:, kappa))
        end if
    end do
else
    !print *, "Using Gauss-Legendre quadrature"
end if

allocate(phihq(Nq, p+1))
! tabulate parent basis at quadrature points
do al = 1, p+1
    call phih_array(xin, al, xiq, phihq(:, al))
end do

n = Nb * 2
DOFs = n
allocate(H(n, n), S(n, n))
allocate(D(n, n), lam(n), lam_tmp(n), fullc(Nn), uq(Nq,Ne), rho(Nq,Ne), Vee(Nq,Ne), &
    Vxc(Nq,Ne), exc(Nq,Ne), Vin(Nq,Ne), Vout(Nq,Ne), xq2(Nq, Ne), &
    rho0(Nq,Ne), rho1(Nq,Ne))

nband = count(focc > 0)
scf_max_iter = 100
scf_alpha = 0.4_dp
scf_L2_eps = 1e-4_dp
scf_eig_eps = eps

allocate(tmp(Nq*Ne))
allocate(energies(nband))
xq2=xq
call get_quad_pts(xe(:2), xiq_gj(:, -1), xq2)
Vin = reshape(thomas_fermi_potential(reshape(xq2, [Nq*Ne]), Z), [Nq, Ne]) + &
    Z / xq2
iter = 0

select case (mixing_scheme)
    case (mixing_scheme_linear)
        call mixing_linear &
            (Ffunc, integral, reshape(Vin, [Nq*Ne]), &
            nband, scf_max_iter, scf_alpha, scf_L2_eps, scf_eig_eps, tmp)
    case (mixing_scheme_pulay)
        call mixing_pulay &
            (Ffunc, integral, matvec, matmat, reshape(Vin, [Nq*Ne]), &
            nband, scf_max_iter, scf_alpha, scf_L2_eps, scf_eig_eps, tmp, 5, 1)
    case default
        error stop "Type of mixing not implemented."
end select

contains

    subroutine Ffunc(x, y, eng)
    ! Converge Vee+Vxc only (the other components are constant)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: y(:), eng(:)
    integer :: idx
    logical :: accurate_eigensolver
    accurate_eigensolver = .true.
    iter = iter + 1
    print *, "SCF iteration:", iter
    Vin = reshape(x, shape(Vin))
    rho = 0
    idx = 0
    Vee = 0
    rho0 = 0
    rho1 = 0
    do kappa = Lmin, Lmax
        if (kappa == 0) cycle
        !print *, "Calculating kappa =", kappa
        if (alpha_j(kappa) > -1) then
            call get_quad_pts(xe(:2), xiq_gj(:, kappa), xq1)
            call proj_fn(Nq-1, xe(:2), xiq_gj(:,-1), wtq_gj(:,-1), xiq_gj(:, kappa), Vin, V(:,:1))
            V(:,1) = V(:,1) - Z/xq1(:,1)
            V(:,2:) = Vin(:,2:) - Z/xq(:,2:)
        else
            V = Vin - Z/xq
        endif

        call assemble_radial_dirac_SH(V, kappa, xin, xe, ib, xiq, wtq, &
            xiq_gj(:, kappa), wtq_gj(:, kappa), alpha(kappa), alpha_j(kappa), c, S, H)

        ! One can enforce symmetry using the following lines, it helps but
        ! we still need two seperate eigensolves for 1e-8 accuracy:
        !H = (H + transpose(H))/2
        !S = (S + transpose(S))/2

        if (accurate_eigensolver) then
            call eigh(H, S, lam)
            call eigh(H, S, lam_tmp, D)
        else
            call eigh(H, S, lam, D)
        end if

        do i = 1, size(focc,1)
            if (focc(i,kappa) < tiny(1._dp)) cycle

            call c2fullc2(in, ib, D(:Nb,i), fullc)
            call fe2quad(xe, xin, xiq, in, fullc, uq)
            rho0 = rho0 - focc(i,kappa)*uq**2 * xq**alpha_j(kappa)
            call fe2quad(xe, xin, xiq_gj(:,-1), in, fullc, uq)
            rho1(:,1) = rho1(:,1) - focc(i,kappa)*uq(:,1)**2 * xq2(:,1)**alpha_j(kappa)

            call c2fullc2(in, ib, D(Nb+1:,i), fullc)
            call fe2quad(xe, xin, xiq, in, fullc, uq)
            rho0 = rho0 - focc(i,kappa)*uq**2 * xq**alpha_j(kappa)
            call fe2quad(xe, xin, xiq_gj(:,-1), in, fullc, uq)
            rho1(:,1) = rho1(:,1) - focc(i,kappa)*uq(:,1)**2 * xq2(:,1)**alpha_j(kappa)

            idx = idx + 1
            eng(focc_idx(i,kappa)) = sqrt(lam(i)) - c**2
        end do

    end do

    if ( .not. (size(eng) == idx) ) then
       error stop 'Size mismatch in energy array'
    end if
    energies = eng

    rho0(:,1) = rho1(:,1)
    rho = rho0  / (4*pi)
    Vee = Vee + hartree_potential_gj(Nq-1, xe, xiq, wtq, &
        xiq_gj(:,-1), wtq_gj(:,-1), alpha_j(-1), -rho0)

    !print *, "Energies:"
    !do i = 1, size(eng)
    !    print *, i, no(i), lo(i), so(i), energies(i)
    !end do

    call xc_vwn3(size(rho), -rho, .TRUE. , c, exc, Vxc)
    Vout = Vee + Vxc ! This term is added later: -Z/xq
    y = reshape(Vout, shape(y))

    call total_energy(xe, xiq, wtq, xiq_gj(:,-1), wtq_gj(:,-1), alpha_j(-1), fo, energies, Vin-Z/xq2, Vee, -Z/xq2, exc, &
        xq2, -rho, T_s, E_ee, E_en, EE_xc, Etot)
    !print *, Etot
    end subroutine

    real(dp) function integral(x)
    ! Computes the integral of the vector 'x'
    real(dp), intent(in) :: x(:)
    integral = integrate2(xe, xiq, xiq_gj(:,-1), wtq, wtq_gj(:,-1), alpha_j(-1), reshape(x, shape(uq)))
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

    subroutine total_energy(xe, xiq, wtq, xiq1, wtq1, alpha, fo, ks_energies, V_in, V_h, V_coulomb, e_xc, &
    R, n, T_s, E_ee, E_en, EE_xc, Etot)
    ! This is a variational, quadratically convergent form of total energy
    real(dp), intent(in) :: xe(:), xiq(:), xiq1(:), wtq(:), wtq1(:), alpha
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
    T_s = E_band + 4*pi * integrate2(xe, xiq, xiq1, wtq, wtq1, alpha, V_in * rho * R**2)

    E_ee = -2*pi * integrate2(xe, xiq, xiq1, wtq, wtq1, alpha, V_h * rho * R**2)
    E_en =  4*pi * integrate2(xe, xiq, xiq1, wtq, wtq1, alpha, (-V_coulomb) * rho * R**2)
    E_c = E_ee + E_en

    EE_xc = -4*pi * integrate2(xe, xiq, xiq1, wtq, wtq1, alpha, e_xc * rho * R**2)

    Etot = T_s + E_c + EE_xc
    end subroutine

end subroutine

subroutine csolve_dirac(Z, p, xiq, wtq, xe, eps, Nq, Ne, nenergies, energies, &
        Etot, V) bind(c)
integer(c_int), intent(in) :: Z, p, Nq, Ne, nenergies
real(c_double), intent(in) :: xe(Ne)        ! element coordinates
real(c_double), intent(in) :: xiq(Nq)       ! quadrature points
real(c_double), intent(in) :: wtq(Nq)       ! quadrature weights
real(c_double), intent(in) :: eps
real(c_double), intent(out) :: energies(nenergies)
real(c_double), intent(out) :: Etot
real(c_double), intent(out) :: V(size(xiq),size(xe)-1)    ! SCF potential
integer :: DOFs
real(dp), allocatable :: energies2(:)

call solve_dirac(Z, p, xiq, wtq, xe, eps, energies2, Etot, V, DOFs)
energies = energies2
end subroutine

end module
