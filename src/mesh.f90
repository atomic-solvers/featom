!> @brief This module is used to generate meshes and derivatives
!> @details The supported meshes are currently:
!>   - Linearly spaced
!>   - Exponential meshes
!> The time step used for derivatives is a uniformly spaced series of integers from 1 to N+1.
module mesh

use types, only: dp

implicit none

private
public meshexp, meshexp_der, get_meshexp_pars, meshexp_der2, &
    linspace, meshgrid

contains

!> @brief This function generates an exponential mesh of N elements on @f$ [r_{\mathrm{min}}, r_{\mathrm{max}}] @f$
!> @param[in] rmin Inclusive endpoint
!> @param[in] rmax Inclusive endpoint
!> @param[in] a Ratio of the rightmost and leftmost elements in the mesh,
!>  controls the mesh spacing. Must satisfy @f$ a > 0 @f$. For @f$ a == 1 @f$ a uniform mesh is
!>  returned. For @f$ a > 1 @f$ this is the largest/smallest ratio.
!> @param[in] N The number of elements in the mesh.
!> @returns mesh The generated mesh as an array of N+1 one elements
!>
!> @details Every exponential mesh is fully determined by the set of parameters
!> `(rmin, rmax, a, N)`. Use the @ref get_meshexp_pars() subroutine to obtain them
!> from the given mesh.
!>
!> @b Example:
!>
!>     real(dp) :: r(11)
!>     r = meshexp(0._dp, 50._dp, 1e9_dp, 10)
function meshexp(rmin, rmax, a, N) result(mesh)
! The domain [rmin, rmax], the mesh will contain both endpoints.
real(dp), intent(in) :: rmin, rmax
! The ratio of the rightmost to leftmost element lengths in the mesh
real(dp), intent(in) :: a
! The number of elements in the mesh:
integer, intent(in) :: N
! The generated mesh:
real(dp) :: mesh(N+1)

integer :: i
real(dp) :: alpha, beta
if (a < 0) then
    error stop "meshexp: a > 0 required"
else if (abs(a - 1) < tiny(1._dp)) then
    alpha = (rmax - rmin) / N
    do i = 1, N+1
        mesh(i) = alpha * (i-1.0_dp) + rmin
    end do
else
    if (N > 1) then
        beta = log(a) / (N-1)
        alpha = (rmax - rmin) / (exp(beta*N) - 1)
        do i = 1, N+1
            mesh(i) = alpha * (exp(beta*(i-1)) - 1) + rmin
        end do
    else if (N == 1) then
        mesh(1) = rmin
        mesh(2) = rmax
    else
        error stop "meshexp: N >= 1 required"
    end if
end if
end function

!> @brief Generates the first derivative dR/dt where R(t) is the mesh returned by @ref meshexp()
!> @param[in] rmin Inclusive endpoint
!> @param[in] rmax Inclusive endpoint
!> @param[in] a Ratio of the rightmost and leftmost elements in the mesh,
!>  controls the mesh spacing. Must satisfy @f$ a > 0 @f$. For @f$ a == 1 @f$ a uniform mesh is
!>  returned. For @f$ a > 1 @f$ this is the largest/smallest ratio.
!> @param[in] N The number of elements in the mesh.
!> @returns Rp(N+1) The first derivative w.r.t time for the mesh as an array of N+1 one elements
!>
!> @details The input parameters are the same as for @ref meshexp().
!> The variable "t" is defined by:
!> t = 1, 2, ..., N+1
!> So it describes a uniform mesh, with a step size 1, and the corresponding
!> physical points are given by the R(t) array.
function meshexp_der(rmin, rmax, a, N) result(Rp)
real(dp), intent(in) :: rmin
real(dp), intent(in) :: rmax
real(dp), intent(in) :: a
integer, intent(in) :: N
real(dp) :: Rp(N+1)

integer :: i
real(dp) :: alpha, beta
if (a < 0) then
    error stop "meshexp_der: a > 0 required"
else if (abs(a - 1) < tiny(1._dp)) then
    error stop "meshexp_der: a == 1 not implemented"
else
    if (N > 1) then
        beta = log(a)/(N-1)
        alpha = (rmax - rmin) / (exp(beta*N) - 1)
        do i = 1, N+1
            Rp(i) = alpha * beta * exp(beta*(i-1))
        end do
    else
        error stop "meshexp_der: N > 1 required"
    end if
end if
end function

!> @brief Generates the second derivative d^R/dt^2 where R(t) is the mesh returned by @ref meshexp()
!> @param[in] rmin Inclusive endpoint
!> @param[in] rmax Inclusive endpoint
!> @param[in] a Ratio of the rightmost and leftmost elements in the mesh,
!>  controls the mesh spacing. Must satisfy @f$ a > 0 @f$. For @f$ a == 1 @f$ a uniform mesh is
!>  returned. For @f$ a > 1 @f$ this is the largest/smallest ratio.
!> @param[in] N The number of elements in the mesh.
!> @returns Rpp(N+1) The second derivative w.r.t time for the mesh as an array of N+1 one elements
!>
!> @details The input parameters are the same as for @ref meshexp().
!> The variable "t" is defined by:
!> t = 1, 2, ..., N+1
!> So it describes a uniform mesh, with a step size 1, and the corresponding
!> physical points are given by the R(t) array.
function meshexp_der2(rmin, rmax, a, N) result(Rpp)
real(dp), intent(in) :: rmin
real(dp), intent(in) :: rmax
real(dp), intent(in) :: a
integer, intent(in) :: N
real(dp) :: Rpp(N+1)

integer :: i
real(dp) :: alpha, beta
if (a < 0) then
    error stop "meshexp_der2: a > 0 required"
else if (abs(a - 1) < tiny(1._dp)) then
    error stop "meshexp_der2: a == 1 not implemented"
else
    if (N > 1) then
        beta = log(a)/(N-1)
        alpha = (rmax - rmin) / (exp(beta*N) - 1)
        do i = 1, N+1
            Rpp(i) = alpha * beta**2 * exp(beta*(i-1))
        end do
    else
        error stop "meshexp_der2: N > 1 required"
    end if
end if
end function

!> @brief Given any exponential mesh R, it determines the parameters
!> @param[in] R The input mesh array
!> @param[out] rmin Inclusive endpoint
!> @param[out] rmax Inclusive endpoint
!> @param[out] a Ratio of the rightmost and leftmost elements in the mesh,
!>  controls the mesh spacing. Must satisfy @f$ a > 0 @f$. For @f$ a == 1 @f$ a uniform mesh is
!>  returned. For @f$ a > 1 @f$ this is the largest/smallest ratio.
!> @param[out] N The number of elements in the mesh.
!>
!> @details This only looks at the number of elements, the leftmost and the rightmost
!>  elements (so the middle elements are not checked/taken into account).
subroutine get_meshexp_pars(R, rmin, rmax, a, N)
real(dp), intent(in) :: R(:)
real(dp), intent(out) :: rmin, rmax, a
integer, intent(out) :: N
rmin = R(1)
rmax = R(size(R))
a = (R(size(R)) - R(size(R)-1)) / (R(2) - R(1))
N = size(R) - 1
end subroutine

!> @brief A helper function to generate a linearly spaced mesh
!> @param[in] a Inclusive endpoint
!> @param[in] b Inclusive endpoint
!> @param[in] n The number of elements in the mesh.
!> @returns s(n) The generated mesh as an array of N+1 one elements
!> @see meshexp()
!>
!> @details This calls @ref meshexp()
function linspace(a, b, n) result(s)
real(dp), intent(in) :: a, b
integer, intent(in) :: n
real(dp) :: s(n)
s = meshexp(a, b, 1.0_dp, n-1)
end function

!> @brief A helper subroutine to generate a two dimensional mesh
!> @param[in] x Mesh along dimension 1
!> @param[in] y Mesh along dimension 2
!> @param[out] x2 Mesh elements replicated y times along dimension 1
!> @param[out] y2 Mesh elements replicated x times along dimension 2
subroutine meshgrid(x, y, x2, y2)
real(dp), intent(in) :: x(:), y(:)
real(dp), intent(out) :: x2(:, :), y2(:, :)
x2 = spread(x, 1, size(y))
y2 = spread(y, 2, size(x))
end subroutine

end module
