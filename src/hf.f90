module hf

use types, only: dp
use mesh,  only: meshexp
use feutils, only: define_connect, get_quad_pts, get_parent_quad_pts_wts, &
                   get_parent_nodes, phih, dphih, c2fullc2, fe2quad_core, get_nodes, &
                   integrate, proj_fn, phih_array, dphih_array
use linalg, only: eigh
use fe,     only: assemble_radial_H, assemble_radial_S, assemble_radial_H_setup, &
                   assemble_radial_H_complete
use constants, only: pi
use hartree_screening, only: assemble_poisson_A, hartree_potential3
use mixings, only: mixing_anderson
use states,  only: get_atomic_states_nonrel_focc, nlf2focc, get_atomic_states_nonrel
use energies, only: thomas_fermi_potential
use iso_c_binding, only: c_double, c_int
use solvers, only: solve_eig_irange, solve_sym_setup
implicit none
private
public solve_hf_schroed

contains

subroutine solve_hf_schroed(Z, p, xiq, wtq, xe, eps, energies, Etot, V, DOFs)
    integer, intent(in) :: Z, p
    real(dp), intent(in) :: xe(:), xiq(:), wtq(:)
    real(dp), intent(in) :: eps
    real(dp), allocatable, intent(out) :: energies(:)
    real(dp), intent(out) :: Etot
    real(dp), allocatable, intent(out) :: V(:,:)
    integer, intent(out) :: DOFs

    integer :: Nq, Ne, Nb, Nn, pp, Nb_p, Nn_p
    real(dp), allocatable :: xin(:), xn(:), xin_p(:), xn_p(:)
    integer, allocatable :: in(:, :), ib(:, :), in_p(:, :), ib_p(:, :)
    real(dp), allocatable :: xq(:, :), phihq(:,:), dphihq(:,:)
    real(dp), allocatable :: phihq_p(:,:), dphihq_p(:,:)

    integer :: n, Lmax, l, i, j, al, eimax, iter
    integer, allocatable :: no(:), lo(:), focc_idx(:,:)
    real(dp), allocatable :: focc(:,:), fo(:)

    real(dp), allocatable :: H(:,:), S(:), D(:,:), lam(:), fullc(:)
    real(dp), allocatable :: Hl(:,:,:)      ! per-l base H (no V)
    real(dp), allocatable :: Am_p(:,:), bv_p(:)
    integer, allocatable :: ipiv(:)

    real(dp), allocatable :: uq(:,:), rho(:,:), n_r(:,:), Vee(:,:), Vx(:,:), Vin(:,:), Vout(:,:), Veff(:,:)
    real(dp), allocatable :: orbitals(:,:,:), orbital_densities(:,:,:)

    real(dp) :: scf_alpha, scf_L2_eps, scf_eig_eps
    integer :: nband, scf_max_iter
    real(dp), allocatable :: tmp(:)
    real(dp) :: Etot_Slater, Etot_old
    real(dp) :: E_x_LDA, E_H, E_band, En_coupling, E_x_exact

    real(dp) :: xiq_lob(p+1), wtq_lob(p+1)

    Nq = size(xiq)
    Ne = size(xe) - 1

    call get_parent_quad_pts_wts(2, p+1, xiq_lob, wtq_lob)

    call get_atomic_states_nonrel_focc(Z, focc)
    call get_atomic_states_nonrel(Z, no, lo, fo)
    Lmax = ubound(focc, 2)

    Nn = Ne*p + 1
    allocate(xin(p+1));  call get_parent_nodes(2, p, xin)
    allocate(in(p+1,Ne), ib(p+1,Ne)); call define_connect(1, 1, Ne, p, in, ib)
    Nb = maxval(ib)
    if (Nn /= maxval(in)) error stop 'Wrong size for Nn'
    DOFs = Nb

    allocate(xq(Nq,Ne)); call get_quad_pts(xe, xiq, xq)
    allocate(xn(Nn));    call get_nodes(xe, xin, xn)

    pp = 2*p
    Nn_p = Ne*pp + 1
    allocate(xin_p(pp+1));  call get_parent_nodes(2, pp, xin_p)
    allocate(in_p(pp+1,Ne), ib_p(pp+1,Ne)); call define_connect(2, 1, Ne, pp, in_p, ib_p)
    Nb_p = maxval(ib_p)
    if (Nn_p /= maxval(in_p)) error stop 'Wrong size for Nn_p'
    allocate(xn_p(Nn_p));   call get_nodes(xe, xin_p, xn_p)

    allocate(phihq(Nq,p+1), dphihq(Nq,p+1))
    allocate(phihq_p(Nq,pp+1), dphihq_p(Nq,pp+1))
    do al=1, p+1
        call phih_array(xin,  al, xiq, phihq(:,al))
        call dphih_array(xin, al, xiq, dphihq(:,al))
    end do
    do al=1, pp+1
        call phih_array(xin_p,  al, xiq, phihq_p(:,al))
        call dphih_array(xin_p, al, xiq, dphihq_p(:,al))
    end do

    n = Nb
    allocate(H(n,n), S(n))
    allocate(D(n,n), lam(n), fullc(Nn))
    allocate(uq(Nq,Ne), rho(Nq,Ne), n_r(Nq,Ne))
    allocate(Vee(Nq,Ne), Vx(Nq,Ne), Vin(Nq,Ne), Vout(Nq,Ne), Veff(Nq,Ne))

    nband        = count(focc > 0.0_dp)
    scf_max_iter = 400
    scf_alpha    = 0.20_dp
    scf_L2_eps   = 1.0e-9_dp
    scf_eig_eps  = eps*0.001
    allocate(tmp(Nq*Ne), energies(nband))

    Vin  = reshape(thomas_fermi_potential(reshape(xq, [Nq*Ne]), Z), [Nq,Ne]) + real(Z,dp) / xq
    iter = 0
    Etot = 0.0_dp

    print *, "=============================================================="
    print *, "Restricted Hartree-Fock: Slater-SCF + exact exchange (all l via multipoles)"
    print *, "=============================================================="
    print *, "Iter |   E_Slater (Ha)  |  dE      |  ||ΔV||_L2"

    call assemble_radial_S(xin, xe, ib, wtq_lob, S)
    do i=1, size(S)
        S(i) = 1.0_dp/sqrt(S(i))
    end do

    allocate(Am_p(Nb_p,Nb_p), bv_p(Nb_p), ipiv(Nb_p))
    Am_p = 0.0_dp
    call assemble_poisson_A(xin_p, xe, ib_p, xiq, wtq, dphihq_p, Am_p)
    call solve_sym_setup(Am_p, ipiv)

    allocate(Hl(n,n,0:Lmax))
    call assemble_radial_H_setup(0, Lmax, xin, xe, ib, xiq, wtq, phihq, dphihq, Hl)

    allocate(focc_idx(size(focc,1), 0:Lmax)); focc_idx = 0
    j = 1
    do l=0, Lmax
        do i=1, size(focc,1)
            if (focc(i,l) > 0.0_dp) then
                focc_idx(i,l) = j
                j = j + 1
            end if
        end do
    end do

    allocate(orbitals(Nq,Ne,nband), orbital_densities(Nq,Ne,nband))

    call mixing_anderson(Ffunc, integral, reshape(Vin, [Nq*Ne]), &
                         nband, scf_max_iter, scf_alpha, scf_L2_eps, scf_eig_eps, tmp)

    call compute_exact_hf_exchange_all_l(xe, xiq, wtq, xq, orbitals, focc_idx, focc, E_x_exact)

    E_band      = sum(fo * energies)
    n_r         = -rho
    En_coupling = 4.0_dp*pi * integrate(xe, wtq, n_r * (Vee + Vx) * xq**2)
    E_H         = 0.5_dp * 4.0_dp*pi * integrate(xe, wtq, n_r * Vee * xq**2)
    Etot        = E_band - En_coupling + E_H + E_x_exact
    V           = -real(Z,dp)/xq + Vee + Vx

    print *, ""
    print '(a, f16.10, a)', "Post-SCF exact exchange (all l): ", E_x_exact, " Ha"
    print '(a, f16.10, a)', "Final RHF total energy:         ", Etot,      " Ha"
    print *, "=============================================================="

contains

    subroutine Ffunc(x, y, eng)
        real(dp), intent(in)  :: x(:)
        real(dp), intent(out) :: y(:), eng(:)
        integer :: i_local, j_local, idx

        Etot_old = Etot
        iter     = iter + 1
        Vin      = reshape(x, shape(Vin))
        rho      = 0.0_dp
        idx      = 0

        Veff = Vin - real(Z,dp)/xq

        do l=0, Lmax
            call assemble_radial_H_complete(Veff, xin, xe, ib, xiq, wtq, phihq, Hl(:,:,l), H)

            do concurrent(i_local=1:size(S), j_local=1:size(S), i_local >= j_local)
                H(i_local,j_local) = H(i_local,j_local) * S(i_local) * S(j_local)
            end do

            eimax = 0
            do i_local=1, size(focc,1)
                if (focc(i_local,l) > 0.0_dp) eimax = i_local
            end do
            if (eimax == 0) cycle

            call solve_eig_irange(H, 1, eimax, lam, D)

            do i_local=1, size(S)
                D(i_local,1:eimax) = D(i_local,1:eimax) * S(i_local)
            end do

            do i_local=1, eimax
                if (focc(i_local,l) <= tiny(1.0_dp)) cycle

                call c2fullc2(in, ib, D(:Nb,i_local), fullc)
                call fe2quad_core(xe, xin, in, fullc, phihq, uq)

                idx = focc_idx(i_local,l)
                eng(idx)                  = lam(i_local)
                orbitals(:,:,idx)         = uq
                orbital_densities(:,:,idx)= (uq/xq)**2 / (4.0_dp*pi)

                rho = rho - focc(i_local,l) * orbital_densities(:,:,idx)
            end do
        end do
        energies = eng

        Vee = hartree_potential3(sum(fo), pp, xe, xin_p, in_p, ib_p, xiq, wtq, &
                                 phihq_p, Am_p, ipiv, -4.0_dp*pi*rho)

        call slater_exchange_potential(-rho, Vx)
        Vout = Vee + Vx

        E_x_LDA     = 0.75_dp * 4.0_dp*pi * integrate(xe, wtq, (-rho) * Vx  * xq**2)
        E_band      = sum(fo * energies)
        E_H         = 0.5_dp * 4.0_dp*pi * integrate(xe, wtq, (-rho) * Vee * xq**2)
        En_coupling = 4.0_dp*pi * integrate(xe, wtq, (-rho) * (Vee + Vx) * xq**2)
        Etot_Slater = E_band - En_coupling + E_H + E_x_LDA

        if (iter > 1) then
            print '(i4, " | ", f16.8, " | ", es9.2, " | ", es9.2)', &
                   iter, Etot_Slater, Etot_Slater-Etot_old, &
                   sqrt(integrate(xe, wtq, (Vout - (Vin - real(Z,dp)/xq))**2))
        else
            print '(i4, " | ", f16.8, " | ", a9, " | ", a9)', iter, Etot_Slater, '---', '---'
        end if

        y   = reshape(Vout, shape(y))
        Etot= Etot_Slater
    end subroutine Ffunc

    real(dp) function integral(x)
        real(dp), intent(in) :: x(:)
        integral = integrate(xe, wtq, reshape(x, [size(xq,1), size(xq,2)]))
    end function integral

end subroutine solve_hf_schroed


subroutine slater_exchange_potential(n, Vx)
    real(dp), intent(in)  :: n(:,:)
    real(dp), intent(out) :: Vx(:,:)
    real(dp), parameter :: Cx = - (3.0_dp/pi)**(1.0_dp/3.0_dp)

    Vx = 0.0_dp
    where (n > 1.0e-18_dp)
        Vx = Cx * n**(1.0_dp/3.0_dp)
    end where
end subroutine slater_exchange_potential


subroutine compute_exact_hf_exchange_all_l(xe, xiq, wtq, xq, orbitals, focc_idx, focc, E_x)
    ! Exact HF exchange for all occupied (n,l) pairs via multipole expansion.
    ! Angular factor: (2l_i+1)(2l_j+1) * [ ( l_i  k  l_j ; 0 0 0 ) ]^2
    ! Radial integral R_k uses P(r)=u(r) with ∫ u^2 dr = 1 (no r^2 factor).

    real(dp), intent(in) :: xe(:), xiq(:), wtq(:), xq(:,:)
    real(dp), intent(in) :: orbitals(:,:,:)
    integer,  intent(in) :: focc_idx(:,0:)
    real(dp), intent(in) :: focc(:,0:)
    real(dp), intent(out):: E_x

    integer :: Nq, Ne, nbocc, a, b, n_i, l_i, n_j, l_j, k, kmin, kmax
    integer, allocatable :: bands(:), nlist(:), llist(:)
    real(dp) :: threej, angfac, Rk
    real(dp), allocatable :: A_t(:,:), r(:), wJ(:), Avals(:), rpk(:), rpn(:)

    Nq = size(xq,1); Ne = size(xq,2)

    ! Build occupied (n,l) list
    nbocc = 0
    do l_i = 0, ubound(focc,2)
        do n_i = 1, size(focc,1)
            if (focc(n_i,l_i) > 0.0_dp) nbocc = nbocc + 1
        end do
    end do
    if (nbocc == 0) then
        E_x = 0.0_dp; return
    end if
    allocate(bands(nbocc), nlist(nbocc), llist(nbocc))
    nbocc = 0
    do l_i = 0, ubound(focc,2)
        do n_i = 1, size(focc,1)
            if (focc(n_i,l_i) > 0.0_dp) then
                nbocc = nbocc + 1
                bands(nbocc) = focc_idx(n_i,l_i)
                nlist(nbocc) = n_i
                llist(nbocc) = l_i
            end if
        end do
    end do

    allocate(A_t(Nq,Ne))
    allocate(r(Nq*Ne), wJ(Nq*Ne), Avals(Nq*Ne), rpk(Nq*Ne), rpn(Nq*Ne))

    E_x = 0.0_dp

    do a = 1, nbocc
        l_i = llist(a)
        n_i = nlist(a)
        do b = 1, nbocc
            l_j = llist(b)
            n_j = nlist(b)

            ! Pair product P_i P_j = u_i u_j
            A_t = orbitals(:,:,bands(a)) * orbitals(:,:,bands(b))

            kmin = abs(l_i - l_j)
            kmax = l_i + l_j

            do k = kmin, kmax
                if (mod(l_i + l_j + k, 2) /= 0) cycle  ! parity

                threej = wigner3j_000(l_i, k, l_j)
                if (abs(threej) < 1.0e-30_dp) cycle

                angfac = real((2*l_i+1)*(2*l_j+1), dp) * (threej*threej)

                call radial_slater_integral_k(xe, xq, wtq, A_t, k, r, wJ, Avals, rpk, rpn, Rk)

                ! No extra 1/2 or occupancy weights here (closed-shell spatial set)
                E_x = E_x - angfac * Rk
            end do
        end do
    end do

    deallocate(bands, nlist, llist, A_t, r, wJ, Avals, rpk, rpn)
end subroutine compute_exact_hf_exchange_all_l


subroutine radial_slater_integral_k(xe, xq, wtq, A, k, r, wJ, Avals, rpk, rpn, Rk)
    ! Compute R_k = ∬ A(r) A(r') r_<^k / r_>^{k+1} dr dr'
    ! Numerically robust:
    !  - global argsort of all (r, wJ, A) points across elements
    !  - Kahan compensated summation for prefix/suffix and final accumulation
    !  - strict < and > split, plus single diagonal add-back: Σ wJ^2 A^2 / r
    !
    ! Inputs:
    !   xe(:)      element nodes
    !   xq(:,:)    quadrature radii (Nq, Ne)
    !   wtq(:)     parent quadrature weights (length Nq)
    !   A(:,:)     pair radial product P_i P_j = u_i * u_j on (Nq, Ne)
    !   k          multipole order (integer >= 0)
    !
    ! Work arrays (intent(inout) as per original signature, reused internally):
    !   r(:), wJ(:), Avals(:), rpk(:), rpn(:) must have length Nq*Ne in caller
    !
    ! Output:
    !   Rk         Slater-type radial integral
    !
    integer,  intent(in)    :: k
    real(dp), intent(in)    :: xe(:), xq(:,:), wtq(:), A(:,:)
    real(dp), intent(inout) :: r(:), wJ(:), Avals(:), rpk(:), rpn(:)
    real(dp), intent(out)   :: Rk

    integer :: Nq, Ne, M, e, q, m_idx, ofs
    integer, allocatable :: ord_global(:)
    real(dp), allocatable :: rS(:), wJS(:), AS(:), rpkS(:), rpnS(:)
    real(dp), allocatable :: prefix(:), suffix(:)
    real(dp) :: he, rr, wq, tinyr
    real(dp) :: s, c, y, t, self_pk, self_pn
    real(dp) :: Rk_acc, Rk_c, diag_acc, diag_c

    Nq   = size(xq,1)
    Ne   = size(xq,2)
    M    = Nq*Ne
    tinyr = 1.0e-24_dp

    ! 1) Flatten all element quadrature samples into (r, wJ, Avals)
    ofs = 0
    do e = 1, Ne
        he = 0.5_dp * (xe(e+1) - xe(e))
        do q = 1, Nq
            ofs        = ofs + 1
            rr         = max(xq(q,e), tinyr)
            wq         = wtq(q) * he
            r(ofs)     = rr
            wJ(ofs)    = wq
            Avals(ofs) = A(q,e)
        end do
    end do

    ! 2) Global argsort by ascending r
    allocate(ord_global(M))
    call argsort_ascending_global(r, ord_global)

    allocate(rS(M), wJS(M), AS(M), rpkS(M), rpnS(M), prefix(M), suffix(M))

    do m_idx = 1, M
        rS(m_idx)   = r(   ord_global(m_idx))
        wJS(m_idx)  = wJ(  ord_global(m_idx))
        AS(m_idx)   = Avals(ord_global(m_idx))
    end do

    ! 3) Precompute r^k and r^(-k-1)
    do m_idx = 1, M
        rpkS(m_idx) = rS(m_idx)**k
        rpnS(m_idx) = rS(m_idx)**(-k-1)
    end do

    ! 4) Kahan-compensated prefix: prefix(m) = Σ_{t<=m} AS(t)*rpkS(t)*wJS(t)
    s = 0.0_dp; c = 0.0_dp
    do m_idx = 1, M
        y = AS(m_idx)*rpkS(m_idx)*wJS(m_idx) - c
        t = s + y
        c = (t - s) - y
        s = t
        prefix(m_idx) = s
    end do

    ! 5) Kahan-compensated suffix: suffix(m) = Σ_{t>=m} AS(t)*rpnS(t)*wJS(t)
    s = 0.0_dp; c = 0.0_dp
    do m_idx = M, 1, -1
        y = AS(m_idx)*rpnS(m_idx)*wJS(m_idx) - c
        t = s + y
        c = (t - s) - y
        s = t
        suffix(m_idx) = s
    end do

    ! 6) Assemble Rk with strict < and > split (exclude diagonal once), Kahan-accumulated
    Rk_acc = 0.0_dp
    Rk_c   = 0.0_dp
    do m_idx = 1, M
        self_pk = AS(m_idx)*rpkS(m_idx)*wJS(m_idx)
        self_pn = AS(m_idx)*rpnS(m_idx)*wJS(m_idx)
        ! sum over t < m: use prefix(m_idx) - self_pk
        ! sum over t > m: use suffix(m_idx) - self_pn
        y = wJS(m_idx)*AS(m_idx) * ( rpnS(m_idx) * (prefix(m_idx) - self_pk) + rpkS(m_idx) * (suffix(m_idx) - self_pn) )
        y = y - Rk_c
        t = Rk_acc + y
        Rk_c   = (t - Rk_acc) - y
        Rk_acc = t
    end do

    ! 7) Add diagonal exactly once: Σ wJ^2 A^2 / r
    diag_acc = 0.0_dp
    diag_c   = 0.0_dp
    do m_idx = 1, M
        y = (wJS(m_idx)*wJS(m_idx)) * (AS(m_idx)*AS(m_idx)) / rS(m_idx) - diag_c
        t = diag_acc + y
        diag_c   = (t - diag_acc) - y
        diag_acc = t
    end do

    Rk = Rk_acc + diag_acc

    deallocate(ord_global, rS, wJS, AS, rpkS, rpnS, prefix, suffix)
end subroutine radial_slater_integral_k


subroutine argsort_ascending_global(v, idx)
    ! Global argsort by ascending v(:)
    real(dp), intent(in)  :: v(:)
    integer, intent(out)  :: idx(:)
    integer :: n, i, j, imin, tmp

    n = size(v)
    do i = 1, n
        idx(i) = i
    end do
    ! Simple selection sort is fine for typical M (<= a few thousand)
    do i = 1, n-1
        imin = i
        do j = i+1, n
            if (v(idx(j)) < v(idx(imin))) imin = j
        end do
        if (imin /= i) then
            tmp     = idx(i)
            idx(i)  = idx(imin)
            idx(imin) = tmp
        end if
    end do
end subroutine argsort_ascending_global


subroutine argsort_ascending(v, idx)
    real(dp), intent(in) :: v(:)
    integer, intent(out) :: idx(:)
    integer :: i, j, n, imin, tmp
    n = size(v)
    do i = 1, n
        idx(i) = i
    end do
    do i = 1, n-1
        imin = i
        do j = i+1, n
            if (v(idx(j)) < v(idx(imin))) imin = j
        end do
        if (imin /= i) then
            tmp = idx(i); idx(i) = idx(imin); idx(imin) = tmp
        end if
    end do
end subroutine argsort_ascending


pure real(dp) function wigner3j_000(l1, l2, l3) result(val)
    integer, intent(in) :: l1, l2, l3
    integer :: L, sgn
    real(dp) :: lnval

    if ( (abs(l1-l2) > l3) .or. (l1 + l2 < l3) ) then
        val = 0.0_dp; return
    end if
    L = l1 + l2 + l3
    if (mod(L,2) /= 0) then
        val = 0.0_dp; return
    end if

    sgn = 1
    if (mod(L/2,2) /= 0) sgn = -1

    lnval = 0.5_dp*( lnfac(L - 2*l1) + lnfac(L - 2*l2) + lnfac(L - 2*l3) - lnfac(L + 1) ) &
            + lnfac(L/2) - ( lnfac(L/2 - l1) + lnfac(L/2 - l2) + lnfac(L/2 - l3) )

    val = sgn * exp(lnval)
end function wigner3j_000


pure real(dp) function lnfac(n) result(l)
    integer, intent(in) :: n
    if (n < 0) then
        l = -huge(1.0_dp)
    else
        l = log_gamma(real(n+1, dp))
    end if
end function lnfac

end module
