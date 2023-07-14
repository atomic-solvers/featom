module schroed_dirac_solver
use types, only: dp
use constants, only: pi
use mesh, only: meshexp
use schroed_glob, only: solve_schroed
use dirac, only: solve_dirac
use feutils, only: define_connect, get_quad_pts, get_parent_quad_pts_wts, &
        get_parent_nodes, phih, dphih, phih_array, dphih_array, c2fullc2, &
        fe2quad, fe2quad_core
use fe, only: assemble_radial_S, assemble_radial_H, assemble_radial_dirac_SH
use jacobi, only: gauss_jacobi
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
    real(dp), allocatable :: xe(:), xiq(:), wtq(:), V(:,:), xin(:), &
        H(:,:), S(:), S2(:,:), DSQ(:), phihq(:,:), dphihq(:,:), xiq_lob(:), wtq_lob(:), &
        D(:,:), lam2(:), xiq1(:), wtq1(:), xq1(:,:), fullc(:), uq(:,:), rho(:,:)
    integer, allocatable :: in(:,:), ib(:,:)
    real(dp) :: rmin
    integer :: al, i, j, l, k, ind, kappa, Nb, Nn
    real(dp) :: E_dirac_shift
    rmin = 0
    allocate(xe(Ne+1), xiq(Nq), xiq1(Nq), wtq(Nq), wtq1(Nq), V(Nq, Ne), xiq_lob(p+1), wtq_lob(p+1))
    allocate(phihq(Nq, p+1), dphihq(Nq, p+1), in(p+1, Ne), ib(p+1, Ne), xin(p+1), xq1(Nq, Ne))

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
        do kappa = -6, 5
            if (kappa == 0) cycle
            if (alpha_j(kappa) > -1) then
                call gauss_jacobi(Nq, 0.0_dp, alpha_j(kappa), xiq1, wtq1)
                call get_quad_pts(xe(:2), xiq1, xq1)
                V(:,:1) = -Z/xq1(:,:1)-E_dirac_shift
            end if
            call assemble_radial_dirac_SH(V, kappa, xin, xe, ib, xiq, wtq, xiq1, wtq1, &
                                            alpha(kappa), alpha_j(kappa), c, S2, H)
            call eigh(H, S2, lam2, D)
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
                lam(ind) = sqrt(lam2(k)) - c**2
                lam(ind) = lam(ind) + E_dirac_shift

                call c2fullc2(in, ib, D(:Nb,k), fullc)
                call fe2quad(xe, xin, xiq, in, fullc, uq)
                eigfn(:,:,ind) = uq * xq**(alpha_j(kappa)/2)
                if (alpha_j(kappa) > -1) then
                    call fe2quad(xe, xin, xiq1, in, fullc, uq)
                    eigfn(:,1,ind) = uq(:,1) * xq1(:,1)**(alpha_j(kappa)/2)
                end if
            end do
        end do
    end if
    end subroutine total_energy

end module
