module schroed_dirac_solver
use types, only: dp
use constants, only: pi
use mesh, only: meshexp
use schroed_glob, only: solve_schroed
use dirac, only: solve_dirac, solve_dirac_eigenproblem
use feutils, only: define_connect, get_quad_pts, get_parent_quad_pts_wts, &
        get_parent_nodes, phih, dphih, phih_array, dphih_array, c2fullc2, &
        fe2quad, fe2quad_core
use fe, only: assemble_radial_S, assemble_radial_H, assemble_radial_dirac_SH
use gjp_gw, only: gauss_jacobi_gw
use linalg, only: eigh
use string_utils, only: str
implicit none

contains

    subroutine total_energy(Z, rmax, Ne, a, p, Nq, DOFs, alpha_int, dirac_int, &
            c, potential_type, Lmin, alpha_j, alpha, lam, eigfn, xq)
    integer, intent(in) :: Z, Ne, Nq, p, alpha_int, dirac_int, potential_type, &
                Lmin
    real(dp), intent(in) :: rmax, a, c, alpha_j(Lmin:), alpha(Lmin:)
    integer, intent(out) :: DOFs
    real(dp), allocatable, intent(out) :: lam(:), eigfn(:,:,:)
    real(dp), intent(out) :: xq(Nq, Ne)
    real(dp), allocatable :: xe(:), xiq(:), wtq(:), V(:,:), Vin(:,:), xin(:), &
        H(:,:), S(:), S2(:,:), DSQ(:), phihq(:,:), dphihq(:,:), xiq_lob(:), wtq_lob(:), &
        D(:,:), lam2(:), xiq1(:), wtq1(:), xq1(:,:), fullc(:), uq(:,:), rho(:,:), &
        focc(:,:), lam_tmp(:), rho1(:,:), xiq_gj(:,:), wtq_gj(:,:)
    integer, allocatable :: in(:,:), ib(:,:), focc_idx(:,:)
    real(dp) :: rmin
    integer :: al, i, j, l, k, ind, kappa, Nb, Nn, idx, Lmin2, Lmax
    real(dp) :: E_dirac_shift
    rmin = 0
    allocate(xe(Ne+1), xiq(Nq), xiq1(Nq), wtq(Nq), wtq1(Nq), V(Nq, Ne), xiq_lob(p+1), wtq_lob(p+1))
    allocate(phihq(Nq, p+1), dphihq(Nq, p+1), in(p+1, Ne), ib(p+1, Ne), xin(p+1), xq1(Nq,1))

    xe = meshexp(rmin, rmax, a, Ne)
    call get_parent_quad_pts_wts(1, Nq, xiq, wtq)
    call get_parent_quad_pts_wts(2, p+1, xiq_lob, wtq_lob)
    call get_parent_nodes(2, p, xin)
    call get_quad_pts(xe, xiq, xq)

    ! tabulate parent basis at quadrature points
    do al = 1, p+1
        call phih_array(xin, al, xiq, phihq(:, al))
        call dphih_array(xin, al, xiq, dphihq(:, al))
    end do

    call define_connect(1, 1, Ne, p, in, ib)
    if (dirac_int == 1 .and. alpha_int == -1) then
        call define_connect(2, 1, Ne, p, in, ib)
    end if
    Nb = maxval(ib)
    DOFS = Nb
    if (dirac_int == 1) then
        DOFS = 2 * Nb
        allocate(S2(DOFS, DOFS))
    end if
    Nn = Ne*p+1
    if ( .not. (Nn == maxval(in)) ) then
        error stop 'Wrong size for Nn'
    end if

    allocate(H(DOFS, DOFS), S(DOFS), DSQ(DOFS), uq(Nq,Ne))
    allocate(D(DOFS, DOFS), lam2(DOFS), rho(Nq,Ne), fullc(Nn))
    if (dirac_int == 1) then
        allocate(lam(47))
        allocate(eigfn(Nq, Ne, 49))
    else
        allocate(lam(28))
        allocate(eigfn(Nq, Ne, 28))
    end if

    if (potential_type == 0) then
        V = -Z/xq
        E_dirac_shift = 500
    else
        V = xq**2/2
        E_dirac_shift = 1000
    end if
    if (dirac_int == 1) then
        V = V - E_dirac_shift
    end if

    if (dirac_int == 0) then
        do l = 0, 6
            call assemble_radial_H(V, l, xin, xe, ib, xiq, wtq, phihq, dphihq, H)
            call assemble_radial_S(xin, xe, ib, wtq_lob, S)
            do i = 1, Nb
                DSQ(i) = 1/sqrt(S(i))
            end do
            do concurrent (i = 1:Nb, j = 1:Nb)
                H(i, j) = H(i, j)*DSQ(i)*DSQ(j)
            end do
            call eigh(H, lam2, D)
            do i = 1, size(S)
                D(i,:) = D(i,:)*DSQ(i)
            end do
            do k = 0, 6-l
                ind = (k+l)*(k+l+1)/2+l+1
                if (k+1 > size(lam2)) then
                    lam(ind) = 0
                    eigfn(:,:,ind) = 0
                else
                    lam(ind) = lam2(k+1)
                    call c2fullc2(in, ib, D(:Nb,k+1), fullc)
                    call fe2quad_core(xe, xin, in, fullc, phihq, uq)
                    eigfn(:,:,ind) = uq/xq
                end if
            end do
        end do
    else
        Lmin2 = -6
        Lmax = 5
        allocate(Vin(size(xq,1), size(xq,2)))
        Vin = 0
        allocate(xiq_gj(size(xiq1),Lmin:Lmax))
        allocate(wtq_gj(size(wtq1),Lmin:Lmax))
        allocate(rho1(Nq,Ne))

        ! Initialize focc and focc_idx
        allocate(focc(max(Lmax,abs(Lmin2))+1,Lmin2:Lmax))
        allocate(focc_idx(max(Lmax,abs(Lmin2))+1,Lmin2:Lmax))
        focc = 0
        focc_idx = 0
        do kappa = Lmin2, Lmax
            if (kappa == 0) cycle
            if (alpha_j(kappa) > -1) then
                call gauss_jacobi_gw(Nq, 0.0_dp, alpha_j(kappa), xiq1, wtq1)
                xiq_gj(:,kappa) = xiq1
                wtq_gj(:,kappa) = wtq1
            end if
            if (kappa < 0) then
                l = -kappa-1
            else
                l = kappa
            end if
            do k = 1, 7-l
                ind = (kappa+1)*(kappa+1)+(k-1)*(2*kappa+1)+(k-1)*(k-2)
                if (kappa < 0) then
                    ind = ind + 2*k-2+(4*k-2)*(-1-kappa)
                end if
                if (kappa == -1) then
                    ind = ind + 1
                end if
                focc_idx(k,kappa) = ind
                focc(k,kappa) = 1
            end do
        end do
        call solve_dirac_eigenproblem(Nb, Nq, Lmin2, Lmax, alpha, alpha_j, xe, xiq_gj, &
            xq, xq1, wtq_gj, V, Z, Vin, D, S2, H, lam2, rho, rho1, .false., fullc, &
            ib, in, idx, lam_tmp, uq, wtq, xin, xiq, focc, focc_idx, lam, xq, &
            E_dirac_shift)
    end if
    end subroutine total_energy

end module
