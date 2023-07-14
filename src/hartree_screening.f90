module hartree_screening
! Hartree screening function
use types, only: dp
use feutils, only: define_connect, c2fullc2, fe2quad, get_parent_nodes, &
        phih, dphih, phih_array, dphih_array, get_quad_pts, &
        fe2quad_core, fe2quad_gj
use solvers, only: solve_sym, solve_sym2
implicit none
private
public assemble_poisson_A, hartree_potential3, hartree_potential_gj

contains

subroutine assemble_poisson_A(xin, xe, ib, xiq, wtq, dphihq, Am)
! forms system equation matrices corresponding to the problem
! -u''(r) = f(x)*r
! subject to boundary conditions consistent with basis specified by ib
! The weak formulation is:
! \int u'(r)*v'(r) \d r = \int f(x) * r * v(x) \d r
real(dp), intent(in) :: xin(:)       ! parent basis nodes
real(dp), intent(in) :: xe(:)        ! element coordinates
integer, intent(in) :: ib(:,:)       ! basis connectivity: ib(i,j) = index of
    ! basis function associated with local basis function i of element j.
    ! 0 = no associated basis fn.
real(dp), intent(in) :: xiq(:)       ! quadrature points
real(dp), intent(in) :: wtq(:)       ! quadrature weights
real(dp), intent(in) :: dphihq(:,:)  ! parent basis derivative at quadrature points
real(dp), intent(out) :: Am(:,:)     ! system matrix: Am c = bv
integer :: Ne, Nb                ! number of elements, basis functions
integer :: p                     ! order of FE/SE basis
integer :: e                     ! element index
integer :: i,j                   ! basis fn indices
integer :: al,be                 ! "alpha", "beta": local basis fn indices
real(dp) :: xa,xb                ! element boundary node coordinates
real(dp) :: jac                  ! Jacobian of transformation from parent
    ! coords xi in [-1,1] to coords x in [xa,xb]: x = (xb-xa)/2 xi + (xb+xa)/2
real(dp), dimension(size(xiq)) :: intq, xq, m

! initializations
Ne = size(ib, 2)
Nb = maxval(ib)
p = size(xin)-1
if (size(xin) /= size(ib,1)) &
    error stop "Error: inconsistent parent node and connectivity dimensions."
if (size(Am,1) /= Nb) &
    error stop "Error: size of Am inconsistent with Nb."
! accumulate Am matrix and bv vector
! compute lower triangle
Am = 0
do e = 1, Ne
    xa = xe(e)
    xb = xe(e+1)
    jac = (xb - xa)/2  ! affine mapping
    xq = (xiq+1)/2*(xb-xa)+xa
    m = wtq*xq**2/jac
    ! compute matrix/vector elements (integrals transformed to [-1,1])
    do al = 1, p+1
        i = ib(al, e)
        if (i == 0) cycle              ! omit boundary basis fns for Dirichlet BCs
        do be = 1, p+1
            j = ib(be, e)
            if (j == 0) cycle           ! omit boundary basis fns for Dirichlet BCs
            if (j > i) cycle            ! compute only lower triangles
            intq = dphihq(:, al) * dphihq(:, be)
            Am(i, j) = Am(i, j) + sum(m*intq)
        end do
    end do
end do
! fill in upper triangle
do j = 1, Nb
    do i = 1, j-1
        Am(i, j) = Am(j, i)
    end do
end do
end subroutine

subroutine assemble_poisson_b(f_vals, xin, xe, ib, xiq, wtq, phihq, bv)
! forms system equation matrices corresponding to the problem
! -u''(r) = f(x)*r
! subject to boundary conditions consistent with basis specified by ib
! The weak formulation is:
! \int u'(r)*v'(r) \d r = \int f(x) * r * v(x) \d r
real(dp), intent(in) :: f_vals(:,:)   ! f(x) at quadrature points:
   ! f_vals(i,j) = value at ith point in jth element
real(dp), intent(in) :: xin(:)       ! parent basis nodes
real(dp), intent(in) :: xe(:)        ! element coordinates
integer, intent(in) :: ib(:,:)       ! basis connectivity: ib(i,j) = index of
    ! basis function associated with local basis function i of element j.
    ! 0 = no associated basis fn.
real(dp), intent(in) :: xiq(:)       ! quadrature points
real(dp), intent(in) :: wtq(:)       ! quadrature weights
real(dp), intent(in) :: phihq(:,:)   ! parent basis at quadrature points
real(dp), intent(out) :: bv(:)       ! source vector: Am c = bv
integer :: Ne, Nb                ! number of elements, basis functions
integer :: p                     ! order of FE/SE basis
integer :: e                     ! element index
integer :: i                     ! basis fn indices
integer :: al                    ! "alpha", local basis fn indices
real(dp) :: xa,xb                ! element boundary node coordinates
real(dp) :: jac                  ! Jacobian of transformation from parent
    ! coords xi in [-1,1] to coords x in [xa,xb]: x = (xb-xa)/2 xi + (xb+xa)/2
real(dp), dimension(size(xiq)) :: xq, bq, m

! initializations
Ne = size(ib, 2)
Nb = maxval(ib)
p = size(xin)-1
if (size(xin) /= size(ib,1)) &
    error stop "Error: inconsistent parent node and connectivity dimensions."
if (size(bv,1) /= Nb) &
    error stop "Error: size of bv inconsistent with Nb."
! accumulate Am matrix and bv vector
! compute lower triangle
bv = 0
do e = 1, Ne
    xa = xe(e)
    xb = xe(e+1)
    jac = (xb - xa)/2  ! affine mapping
    xq = (xiq+1)/2*(xb-xa)+xa
    m = wtq*xq**2*jac
    ! compute matrix/vector elements (integrals transformed to [-1,1])
    do al = 1, p+1
        i = ib(al, e)
        if (i == 0) cycle              ! omit boundary basis fns for Dirichlet BCs
        bq = f_vals(:,e)*phihq(:,al)
        bv(i) = bv(i) + sum(m * bq)
    end do
end do
end subroutine

function hartree_potential3(Z, p, xe, xin, in, ib, xiq, wtq, phihq, Am, ipiv, f) result(uq)
! Equivalently, the V(r) is the solution of the radial Poisson equation:
!   V''(r) + (2/r)*V'(r) = -f(r)
! Note that in DFT, f(r) = 4*pi*n(r)
! We rewrite the equation as:
!   -u''(r) = f(r)*r
! where u(r) = V(r)*r.
! The weak form is:
!   \int u'(r)*v'(r) \d r = \int f(x) * r * v(x) \d r
! And the boundary conditions are u(0) = 0 and u'(rmax) = 0. Finally, the
! potential is calculated as V(r) = u(r) / r.
!
! The input f is the function f(r) defined on Gauss-Legendre (GL) quadrature
! grid.
!
! The V(r) is returned on the GL grid.
!
real(dp), intent(in) :: Z
integer, intent(in) :: p ! polynomial order to use for Y^k
real(dp), intent(in) :: xe(:), xin(:), xiq(:), wtq(:), f(:,:), phihq(:,:), Am(:,:)
integer, intent(in) :: in(:,:), ib(:,:), ipiv(:)
real(dp), dimension(size(xiq), size(xe)-1) :: uq, xq
integer :: Ne
real(dp) :: fullc((size(xe-1)-1)*p+1)
real(dp), allocatable :: bv(:), u(:)
integer :: Nb
Ne = size(xe)-1
! Make sure the fullc has the right size:
if ( .not. (size(fullc) == maxval(in)) ) then
   error stop 'Wrong size for fullc'
end if
call get_quad_pts(xe, xiq, xq)

Nb = maxval(ib)
allocate(bv(Nb), u(Nb))
call assemble_poisson_b(f, xin, xe, ib, xiq, wtq, phihq, bv)
! solve
u = solve_sym2(Am, bv, ipiv)
! Transform solution from FE coefficient vector to full coefficient vector
call c2fullc2(in, ib, u, fullc)
! transform solution to quadrature grid
call fe2quad_core(xe, xin, in, fullc, phihq, uq)
uq = uq + Z/xe(Ne+1)
end function


subroutine assemble_poisson_gj(f_vals, xin, xe, ib, xiq, wtq, xiq1, wtq1, &
        zeta, Am, bv)
! forms system equation matrices corresponding to the problem
! -u''(r) = f(x)*r
! subject to boundary conditions consistent with basis specified by ib
! The weak formulation is:
! \int u'(r)*v'(r) \d r = \int f(x) * r * v(x) \d r
! The f(r) is rewritten as f(r) = r^zeta * (f(r)/r^zeta) and the GJ quadrature
! integrates f(r)/r^zeta. The zeta and xiq1/wtq1 must be consistent and
! containt he GJ quadrature.
real(dp), intent(in) :: f_vals(:,:)   ! f(x) at quadrature points:
   ! f_vals(i,j) = value at ith point in jth element
real(dp), intent(in) :: xin(:)       ! parent basis nodes
real(dp), intent(in) :: xe(:)        ! element coordinates
integer, intent(in) :: ib(:,:)       ! basis connectivity: ib(i,j) = index of
    ! basis function associated with local basis function i of element j.
    ! 0 = no associated basis fn.
real(dp), intent(in) :: xiq(:), xiq1(:)       ! quadrature points
real(dp), intent(in) :: wtq(:), wtq1(:)       ! quadrature weights
real(dp), intent(in) :: zeta
real(dp), intent(out) :: Am(:,:)     ! system matrix: Am c = bv
real(dp), intent(out) :: bv(:)       ! source vector: Am c = bv
integer :: Ne, Nb                ! number of elements, basis functions
integer :: p                     ! order of FE/SE basis
integer :: e                     ! element index
integer :: i,j                   ! basis fn indices
integer :: al,be                 ! "alpha", "beta": local basis fn indices
integer :: iq                    ! quadrature point index
real(dp) :: xa,xb                ! element boundary node coordinates
real(dp) :: jac                  ! Jacobian of transformation from parent
    ! coords xi in [-1,1] to coords x in [xa,xb]: x = (xb-xa)/2 xi + (xb+xa)/2
real(dp), dimension(size(xiq),size(xin)) :: phihq, dphihq, phihq1   ! parent
    ! basis fns and derivs at quadrature points:
    ! phihq(i,j) = value of jth function at ith quadrature point
real(dp), dimension(size(xiq)) :: intq, bq, xq

! initializations
Ne = size(ib, 2)
Nb = maxval(ib)
p = size(xin)-1
if (size(xin) /= size(ib,1)) &
    error stop "Error: inconsistent parent node and connectivity dimensions."
if (size(Am,1) /= Nb .or. size(bv,1) /= Nb) &
    error stop "Error: size of Am and/or bv inconsistent with Nb."
! tabulate parent basis and derivatives at quadrature points
do al = 1, p+1
    do iq = 1, size(xiq)
        phihq(iq, al) = phih(xin, al, xiq(iq))
        phihq1(iq, al) = phih(xin, al, xiq1(iq))
        dphihq(iq, al) = dphih(xin, al, xiq(iq))
    end do
end do
! accumulate Am matrix and bv vector
! compute lower triangle
Am = 0; bv = 0
do e = 1, Ne
    xa = xe(e)
    xb = xe(e+1)
    jac = (xb - xa)/2  ! affine mapping
    if (e == 1 .and. zeta > -1) then
        xq = (xiq1+1)/2*(xb-xa)+xa
        bq = xq * jac**(zeta) * jac*wtq1 / xq**(zeta)
    else
        xq = (xiq+1)/2*(xb-xa)+xa
        bq = xq * jac*wtq
    end if
    ! compute matrix/vector elements (integrals transformed to [-1,1])
    do al = 1, p+1
        i = ib(al, e)
        if (i == 0) cycle              ! omit boundary basis fns for Dirichlet BCs
        if (e == 1 .and. zeta > -1) then
            bv(i) = bv(i) + sum(bq*f_vals(:, e)*phihq1(:,al))
        else
            bv(i) = bv(i) + sum(bq*f_vals(:, e)*phihq(:,al))
        end if
        do be = 1, p+1
            j = ib(be, e)
            if (j == 0) cycle           ! omit boundary basis fns for Dirichlet BCs
            if (j > i) cycle            ! compute only lower triangles
            intq = dphihq(:, al) * dphihq(:, be) / jac**2
            Am(i, j) = Am(i, j) + sum(wtq*intq*jac)
        end do
    end do
end do
! fill in upper triangle
do j = 1, Nb
    do i = 1, j-1
        Am(i, j) = Am(j, i)
    end do
end do
end subroutine

function hartree_potential_gj(p, xe, xiq, wtq, xiq1, wtq1, zeta, f0) result(uq)
! Equivalently, the V(r) is the solution of the radial Poisson equation:
!   V''(r) + (2/r)*V'(r) = -f(r)
! Note that in DFT, f(r) = 4*pi*n(r)
! We rewrite the equation as:
!   -u''(r) = f(r)*r
! where u(r) = V(r)*r.
! The weak form is:
!   \int u'(r)*v'(r) \d r = \int f(x) * r * v(x) \d r
! And the boundary conditions are u(0) = 0 and u'(rmax) = 0. Finally, the
! potential is calculated as V(r) = u(r) / r.
!
! The input f0 is the function f(r) defined on Gauss-Legendre (GL) quadrature
! grid for all the elements except for the first element where f(r) is defined
! on a Gauss-Jacobi (GJ) quadrature grid with exponent zeta at the origin.
!
! The solution V(r) is returned on the same grid.
!
integer, intent(in) :: p ! polynomial order to use for Y^k
real(dp), intent(in) :: xe(:), xiq(:), wtq(:), xiq1(:), wtq1(:), f0(:,:)
real(dp), dimension(size(xiq), size(xe)-1) :: uq, xq
real(dp), intent(in) :: zeta
real(dp) :: Yq(size(xiq), size(xe)-1)
!real(dp) :: xn(p+1, size(xe)-1)
real(dp), allocatable :: xin(:)
integer, allocatable :: ib(:, :), in(:, :)
integer :: Ne
real(dp) :: fullc((size(xe-1)-1)*p+1), f(size(f0,1), size(f0,2))
real(dp), allocatable :: Am(:, :), bv(:), u(:)
integer :: Nb
Ne = size(xe)-1
allocate(xin(p+1))
call get_parent_nodes(2, p, xin)
allocate(in(p+1, Ne), ib(p+1, Ne))
call define_connect(1, 2, Ne, p, in, ib)
! Make sure the fullc has the right size:
if ( .not. (size(fullc) == maxval(in)) ) then
   error stop 'Wrong size for fullc'
end if

call get_quad_pts(xe, xiq, xq)

!print *, "Poisson: zeta =", zeta
!print *, "Poisson: using GJ"
f = f0

Nb = maxval(ib)
allocate(Am(Nb, Nb), bv(Nb), u(Nb))
call assemble_poisson_gj(f, xin, xe, ib, xiq, wtq, xiq1, wtq1, zeta, Am, bv)

! solve
u = solve_sym(Am, bv)
! Transform solution from FE coefficient vector to full coefficient vector
call c2fullc2(in, ib, u, fullc)
! transform solution to quadrature grid
call fe2quad(xe, xin, xiq, in, fullc, Yq)

call get_quad_pts(xe(:2), xiq1, xq(:,:1))
call fe2quad_gj(xe, xin, xiq, xiq1, in, fullc, Yq)
uq = Yq/xq
end function

end module
