module fe
use types, only: dp
use constants, only: pi
use feutils, only: get_quad_pts, phih, dphih, phih_array, dphih_array
implicit none
private
public assemble_radial_SH, &
    assemble_radial_dirac_SH, assemble_radial_S, assemble_radial_H, &
    assemble_radial_H_setup, assemble_radial_H_complete

contains

subroutine assemble_radial_SH(V, l, xin, xe, ib, xiq, wtq, alpha, S, H)
real(dp), intent(in) :: V(:,:) ! V(:,e) potential at quadrature grid for
                               ! element `e`
integer, intent(in) :: l
real(dp), intent(in) :: xin(:)       ! parent basis nodes
real(dp), intent(in) :: xe(:)        ! element coordinates
integer, intent(in) :: ib(:, :)       ! basis connectivity: ib(i,j) = index of
   ! basis function associated with local basis function i of element j. 0 = no
   ! associated basis fn.
real(dp), intent(in) :: xiq(:)       ! quadrature points
real(dp), intent(in) :: wtq(:)       ! quadrature weights
real(dp), intent(in) :: alpha
real(dp), intent(out) :: S(:,:), H(:,:)
integer Ne, Nb, p, e, i, j, al, be
real(dp) xa, xb, jac
real(dp), dimension(size(xiq), size(xin)) :: phihq, dphihq
real(dp), dimension(size(xiq)) :: hq, x, Bi, Bj, Bip, Bjp, m

Ne = size(xe)-1
Nb = maxval(ib)
p = size(xin)-1
! tabulate parent basis and derivatives at quadrature points
do al = 1, p+1
    call phih_array(xin, al, xiq, phihq(:, al))
    call dphih_array(xin, al, xiq, dphihq(:, al))
end do
S = 0
H = 0
do e = 1, Ne
    xa = xe(e)
    xb = xe(e+1)
    jac = (xb-xa)/2  ! affine mapping
    x = (xiq+1)/2*(xb-xa)+xa
    m = x**(2*alpha) * jac*wtq
    do be = 1, p+1
        j = ib(be, e)
        if (j==0) cycle              ! omit boundary basis fns for Dirichlet BCs
        do al = 1, p+1
            i = ib(al, e)
            if (i==0) cycle          ! omit boundary basis fns for Dirichlet BCs
            if (j>i) cycle           ! compute only lower triangles

            Bi = phihq(:,al)
            Bj = phihq(:,be)
            Bip = dphihq(:,al)/jac
            Bjp = dphihq(:,be)/jac

            hq = Bi * Bj
            S(i,j) = S(i,j) + sum(hq*m)
            hq = 0.5_dp * Bip*Bjp &
                + Bi * (l*(l+1) - alpha*(alpha-1))/(2*x**2) * Bj &
                + Bi * V(:,e) * Bj
            H(i,j) = H(i,j) + sum(hq*m)
        end do
    end do
end do
! fill in upper triangles
do j = 1, Nb
    do i = 1, j-1
        S(i,j) = S(j,i)
        H(i,j) = H(j,i)
    end do
end do
end subroutine

subroutine assemble_radial_H_setup(Lmin, Lmax, xin, xe, ib, xiq, wtq, phihq, dphihq, H)
integer, intent(in) :: Lmax, Lmin
real(dp), intent(in) :: xin(:)       ! parent basis nodes
real(dp), intent(in) :: xe(:)        ! element coordinates
integer, intent(in) :: ib(:, :)       ! basis connectivity: ib(i,j) = index of
   ! basis function associated with local basis function i of element j. 0 = no
   ! associated basis fn.
real(dp), intent(in) :: xiq(:)       ! quadrature points
real(dp), intent(in) :: wtq(:)       ! quadrature weights
real(dp), intent(in) :: phihq(:,:)   ! parent basis at quadrature points
real(dp), intent(in) :: dphihq(:,:)  ! parent basis derivative at quadrature points
real(dp), intent(out) :: H(maxval(ib),maxval(ib),Lmin:Lmax)
integer Ne, Nb, p, e, i, j, al, be, l
real(dp) xa, xb, jac, h0
real(dp), dimension(size(xiq)) :: x, Bi, Bj, Bip, Bjp, m, m2

Ne = size(xe)-1
Nb = maxval(ib)
p = size(xin)-1
H = 0
do e = 1, Ne
    xa = xe(e)
    xb = xe(e+1)
    jac = (xb-xa)/2  ! affine mapping
    x = (xiq+1)/2*(xb-xa)+xa
    m  = jac*wtq * 0.5_dp / x**2
    m2 = jac*wtq * 0.5_dp / jac**2
    do be = 1, p+1
        j = ib(be, e)
        if (j==0) cycle              ! omit boundary basis fns for Dirichlet BCs
        do al = 1, p+1
            i = ib(al, e)
            if (i==0) cycle          ! omit boundary basis fns for Dirichlet BCs
            if (j>i) cycle           ! compute only lower triangles

            Bi = phihq(:,al)
            Bj = phihq(:,be)
            Bip = dphihq(:,al)
            Bjp = dphihq(:,be)

            h0 = sum(m2*Bip*Bjp)
            do l = Lmin, Lmax
                H(i, j, l) = H(i, j, l) + l*(l+1)*sum(Bi*Bj*m) + h0
            end do
        end do
    end do
end do

end subroutine

subroutine assemble_radial_H_complete(V, xin, xe, ib, xiq, wtq, phihq, H0, H)
real(dp), intent(in) :: V(:,:) ! V(:,e) potential at quadrature grid for
                               ! element `e`
real(dp), intent(in) :: xin(:)       ! parent basis nodes
real(dp), intent(in) :: xe(:)        ! element coordinates
integer, intent(in) :: ib(:, :)       ! basis connectivity: ib(i,j) = index of
   ! basis function associated with local basis function i of element j. 0 = no
   ! associated basis fn.
real(dp), intent(in) :: xiq(:)       ! quadrature points
real(dp), intent(in) :: wtq(:)       ! quadrature weights
real(dp), intent(in) :: phihq(:,:)   ! parent basis at quadrature points
real(dp), intent(in) :: H0(:,:)
real(dp), intent(out) :: H(maxval(ib),maxval(ib))
integer Ne, Nb, p, e, i, j, al, be
real(dp) xa, xb, jac
real(dp), dimension(size(xiq)) :: Bi, Bj, m

Ne = size(xe)-1
Nb = maxval(ib)
p = size(xin)-1
H = 0
do e = 1, Ne
    xa = xe(e)
    xb = xe(e+1)
    jac = (xb-xa)/2  ! affine mapping
    m = jac*wtq * V(:,e)
    do be = 1, p+1
        j = ib(be, e)
        if (j==0) cycle              ! omit boundary basis fns for Dirichlet BCs
        do al = 1, p+1
            i = ib(al, e)
            if (i==0) cycle          ! omit boundary basis fns for Dirichlet BCs
            if (j>i) cycle           ! compute only lower triangles
            Bi = phihq(:,al)
            Bj = phihq(:,be)
            H(i,j) = H(i,j) + sum(Bi * Bj * m)
        end do
    end do
end do

do concurrent (i = 1:Nb, j = 1:Nb, i>=j)
    H(i, j) = H(i, j) + H0(i, j)
end do

end subroutine


subroutine assemble_radial_H(V, l, xin, xe, ib, xiq, wtq, phihq, dphihq, H)
real(dp), intent(in) :: V(:,:) ! V(:,e) potential at quadrature grid for
                               ! element `e`
integer, intent(in) :: l
real(dp), intent(in) :: xin(:)       ! parent basis nodes
real(dp), intent(in) :: xe(:)        ! element coordinates
integer, intent(in) :: ib(:, :)       ! basis connectivity: ib(i,j) = index of
   ! basis function associated with local basis function i of element j. 0 = no
   ! associated basis fn.
real(dp), intent(in) :: xiq(:)       ! quadrature points
real(dp), intent(in) :: wtq(:)       ! quadrature weights
real(dp), intent(in) :: phihq(:,:)   ! parent basis at quadrature points
real(dp), intent(in) :: dphihq(:,:)  ! parent basis derivative at quadrature points
real(dp), intent(out) :: H(:,:)
integer Ne, Nb, p, e, i, j, al, be
real(dp) xa, xb, jac
real(dp), dimension(size(xiq)) :: hq, x, Bi, Bj, Bip, Bjp, m

Ne = size(xe)-1
Nb = maxval(ib)
p = size(xin)-1
H = 0
do e = 1, Ne
    xa = xe(e)
    xb = xe(e+1)
    jac = (xb-xa)/2  ! affine mapping
    x = (xiq+1)/2*(xb-xa)+xa
    m = jac*wtq
    do be = 1, p+1
        j = ib(be, e)
        if (j==0) cycle              ! omit boundary basis fns for Dirichlet BCs
        do al = 1, p+1
            i = ib(al, e)
            if (i==0) cycle          ! omit boundary basis fns for Dirichlet BCs
            if (j>i) cycle           ! compute only lower triangles

            Bi = phihq(:,al)
            Bj = phihq(:,be)
            Bip = dphihq(:,al)/jac
            Bjp = dphihq(:,be)/jac
            hq = 0.5_dp * Bip*Bjp &
                + Bi * (l*(l+1))/(2*x**2) * Bj &
                + Bi * V(:,e) * Bj
            H(i,j) = H(i,j) + sum(hq*m)
        end do
    end do
end do
end subroutine

subroutine assemble_radial_S(xin, xe, ib, wtq, S)
real(dp), intent(in) :: xin(:)       ! parent basis nodes
real(dp), intent(in) :: xe(:)        ! element coordinates
integer, intent(in) :: ib(:, :)       ! basis connectivity: ib(i,j) = index of
   ! basis function associated with local basis function i of element j. 0 = no
   ! associated basis fn.
real(dp), intent(in) :: wtq(:)       ! quadrature weights
real(dp), intent(out) :: S(:)
integer Ne, p, e, j, be
real(dp) xa, xb, jac

Ne = size(xe)-1
p = size(xin)-1
S = 0
do e = 1, Ne
    xa = xe(e)
    xb = xe(e+1)
    jac = (xb-xa)/2  ! affine mapping
    do be = 1, p+1
        j = ib(be, e)
        if (j==0) cycle              ! omit boundary basis fns for Dirichlet BCs
        S(j) = S(j) + jac*wtq(be)
    end do
end do
end subroutine

subroutine assemble_radial_dirac_SH(V, kappa, xin, xe, ib, xiq, wtq, xiq1, wtq1, &
        alpha, alpha_j, c, S, H)
real(dp), intent(in) :: V(:,:) ! V(:,e) potential at quadrature grid for
                               ! element `e
integer, intent(in) :: kappa
real(dp), intent(in) :: xin(:)       ! parent basis nodes
real(dp), intent(in) :: xe(:)        ! element coordinates
integer, intent(in) :: ib(:, :)       ! basis connectivity: ib(i,j) = index of
   ! basis function associated with local basis function i of element j. 0 = no
   ! associated basis fn.
real(dp), intent(in) :: xiq(:)       ! quadrature points
real(dp), intent(in) :: wtq(:)       ! quadrature weights
real(dp), intent(in) :: xiq1(:)       ! quadrature points first element
real(dp), intent(in) :: wtq1(:)       ! quadrature weights first element
real(dp), intent(in) :: alpha ! alpha in r^(2*alpha) in Dirac equations
real(dp), intent(in) :: alpha_j ! The alpha to use in Gauss-Jacobi quadrature
                                ! On the first element. Set alpha_j = -2 to use
                                ! Gauss-Legendre quadrature on the first element
real(dp), intent(in) :: c ! Speed of light in atomic units
real(dp), intent(out) :: S(:,:), H(:,:)
integer Ne, Nb, p, e, i, j, al, be
real(dp) xa, xb, jac
real(dp), dimension(size(xiq), size(xin)) :: phihq, dphihq, phihq1, dphihq1
real(dp), dimension(size(xiq)) :: hq, x, Bi, Bj, Bip, Bjp, Vq, m

Ne = size(xe)-1
Nb = maxval(ib)
p = size(xin)-1
! tabulate parent basis and derivatives at quadrature points
do al = 1, p+1
    call phih_array(xin, al, xiq, phihq(:, al))
    call phih_array(xin, al, xiq1, phihq1(:, al))
    call dphih_array(xin, al, xiq, dphihq(:, al))
    call dphih_array(xin, al, xiq1, dphihq1(:, al))
end do

S = 0
H = 0
do e = 1, Ne
    xa = xe(e)
    xb = xe(e+1)
    jac = (xb-xa)/2  ! affine mapping
    ! m contains the r^(2*alpha) term, either directly, or as part
    ! of the Gauss-Jacobi quadrature on the first element
    if (e == 1 .and. alpha_j > -1) then
        x = (xiq1+1)/2*(xb-xa)+xa
        m = x**(2*alpha-alpha_j) * jac**(alpha_j) * jac*wtq1
    else
        x = (xiq+1)/2*(xb-xa)+xa
        m = x**(2*alpha-alpha_j) *   x**(alpha_j) * jac*wtq
        ! Note: The last line is equivalent to:
        !m = x**(2*alpha) * jac*wtq
    end if
    Vq = V(:, e)
    do be = 1, p+1
        j = ib(be, e)
        if (j==0) cycle              ! omit boundary basis fns for Dirichlet BCs
        do al = 1, p+1
            i = ib(al, e)
            if (i==0) cycle          ! omit boundary basis fns for Dirichlet BCs
            if (e == 1 .and. alpha_j > -1) then
                Bi = phihq1(:,al)
                Bj = phihq1(:,be)
                Bip = dphihq1(:,al)/jac
                Bjp = dphihq1(:,be)/jac
            else
                Bi = phihq(:,al)
                Bj = phihq(:,be)
                Bip = dphihq(:,al)/jac
                Bjp = dphihq(:,be)/jac
            end if
            ! The diagonal terms are given by the equations
            ! A11
            hq = c**2*Bip*Bjp + ((Vq + c**2)**2 + &
                    c**2*(kappa*(kappa+1)-alpha*(alpha-1))/x**2)*Bi*Bj
            H(i,j) = H(i,j) + sum(m*hq)
            ! A22
            hq = c**2*Bip*Bjp + ((Vq - c**2)**2 + &
                    c**2*(kappa*(kappa-1)-alpha*(alpha-1))/x**2)*Bi*Bj
            H(i+Nb,j+Nb) = H(i+Nb,j+Nb) + sum(m*hq)


            !hq = -c*2*Vq*Bi*Bjp + c*(2*(kappa-alpha)/x*Vq)*Bi*Bj &
            !       - c*Vqp*Bi*Bj
            ! We expand the last term using:

            ! V' * Bi * Bj * r^(2*alpha)
            ! = - V * (Bi * Bj * r^(2*alpha))'
            ! = - V * (Bi' * Bj * r^(2*alpha) + Bi * Bj' * r^(2*alpha)
            !       + Bi * Bj * 2*alpha*r^(2*alpha-1))
            ! = - V * (Bi'*Bj + Bi*Bj' + Bi*Bj * 2*alpha/r) * r^(2*alpha)

            ! to get:
            !        + c*Vq*(Bip*Bj + Bi*Bjp + Bi*Bj*2*alpha/x)
            ! and add things together:
            !hq = -c*Vq*Bi*Bjp + c*(2*(kappa-alpha)/x*Vq)*Bi*Bj &
            !       + c*Vq*Bip*Bj + c*Vq*Bi*Bj*2*alpha/x
            ! The alpha terms cancel:
            !hq = -c*Vq*Bi*Bjp + c*(2*kappa/x*Vq)*Bi*Bj &
            !       + c*Vq*Bip*Bj
            ! We get c*Vq out:
            !hq = c*Vq(-Bi*Bjp + (2*kappa/x)*Bi*Bj + Bip*Bj)

            ! If you exchange i <-> j, the A12 term becomes A21, i.e.:
            ! A^{12}_{ij} = A^{21}_{ji}
            ! A12
            hq = c*Vq*(-Bi*Bjp + 2*kappa/x*Bi*Bj + Bip*Bj)
            H(i,j+Nb) = H(i,j+Nb) + sum(m*hq)
            ! A21
            hq = c*Vq*( Bi*Bjp + 2*kappa/x*Bi*Bj - Bip*Bj)
            H(i+Nb,j) = H(i+Nb,j) + sum(m*hq)

            ! B11
            hq = Bi*Bj
            S(i,j) = S(i,j) + sum(m*hq)
            ! B22
            S(i+Nb,j+Nb) = S(i,j)
        end do
    end do
end do

!print *, "Checking symmetry"
do j = 1, 2*Nb
    do i = 1, j-1
        if (abs(H(i,j) - H(j,i)) > 1e-8_dp) then
            if (abs(H(i,j)-H(j,i)) / max(abs(H(i,j)), abs(H(j,i))) &
                    > 1e-8_dp) then
                print *, i, j, H(i,j)-H(j,i), H(i,j), H(j,i)
                error stop "H not symmetric"
            end if
        end if
        if (abs(S(i,j)-S(j,i)) > 1e-12_dp) error stop "S not symmetric"
   end do
end do
end subroutine


end module
