module solvers
use types, only: dp
use lapack, only: dlamch, dsysv, dsyevx, dsytrf, dsytrs
implicit none
private
public solve_sym, solve_eig_irange, solve_sym2, solve_sym_setup

contains

subroutine solve_eig_irange(Am, l, h, lam, c)
! solve standard eigenvalue problem and return l-th to h-th eigenvalues
! and eigenvectors
real(dp), intent(in) :: Am(:,:)   ! LHS matrix: Am c = lam c
integer, intent(in) :: l, h
real(dp), intent(out) :: lam(:)   ! eigenvalues: Am c = lam c
real(dp), intent(out) :: c(:,:)   ! eigenvectors: Am c = lam c; c(i,j) = ith component of jth vec.
integer n
! lapack variables
integer lwork, info, m
integer, allocatable:: iwork(:), ifail(:)
real(dp) abstol
real(dp), allocatable:: Amt(:,:), work(:)

! solve
n=size(Am,1)
lwork=8*n
if ( .not. (h >= l) ) then
   error stop 'h must be greater than l'
end if
if ( .not. (size(c,1) == n) ) then
   error stop 'Wrong size for the first eigenvector'
end if
if ( .not. (size(c,2) >= h-l+1) ) then
   error stop 'Wrong size for the second eigenvector'
end if
if ( .not. (size(lam) == n) ) then
   error stop 'Wrong size for the eigenvalues'
end if
allocate(Amt(n,n),work(lwork),iwork(5*n),ifail(n))
Amt=Am; ! Amt temporaries overwritten by dsyevx
abstol=2*dlamch('S')
call dsyevx('V','I','L',n,Amt,n,0.0_dp,0.0_dp,l,h,abstol,m, &
    lam,c,n,work,lwork,iwork,ifail,info)
if (info/=0) then
    print *, "dsyevx returned info =", info
    if (info > 0) then
        print *,  "algorithm failed to compute an eigenvalue while working ", &
        "on the submatrix lying in rows and columns ", info/(n+1), &
        " through ", mod(info, n+1)
    else
        print *, -info, "-th argument had an illegal value"
    end if
    error stop 'DSYEVX ERROR'
end if
if ( .not. (m == h-l+1) ) then
   error stop 'Wrong size for the lapack helper'
end if
end subroutine

function solve_sym(Am, bv) result(c)
! solves symmetric dense system of equations
real(dp), intent(in) :: Am(:,:)   ! system matrix: Am c = bv
real(dp), intent(in) :: bv(:)     ! source vector: Am c = bv
real(dp) :: c(size(bv))     ! solution vector: Am c = bv
!real(dp) :: r(size(bv))
integer :: n
! lapack variables
integer :: lwork, info
integer, allocatable :: ipiv(:)
real(dp), allocatable :: Amt(:,:),bm(:,:),work(:)

n = size(c)
lwork = n
allocate(Amt(n,n), bm(n,1), ipiv(n), work(lwork))
Amt=Am; bm(:,1)=bv   ! temporaries for dsysv
call dsysv('L', n, 1, Amt, n, ipiv, bm, n, work, lwork, info)
if (info < 0) then
    print *, "The", -info, "-th argument had illegal value"
    error stop 'DSYSV ERROR.'
end if
if (info > 0) then
    print *, "D(", info, ",", info, ") is exactly zero."
    print *, "The factorization has been completed, but the block diagonal"
    print *, "matrix D is exactly singular, so the solution could not be"
    print *, "computed."
    error stop 'DSYSV ERROR.'
end if
c = bm(:, 1)
! error
!r=matmul(Am, c)-bv
!write(*,'(1x,a,es18.11)') "Solution vector residual ||Am c - bv||/||bv||: ", &
!   sqrt(dot_product(r,r)/dot_product(bv,bv))
end function

subroutine solve_sym_setup(Am, ipiv)
! factorizes symmetric matrix Am, such that solve_sym2 can use the factorization
real(dp), intent(inout) :: Am(:,:)   ! system matrix: Am c = bv
integer, intent(out) :: ipiv(size(Am, 1))   ! internal details
integer :: n
! lapack variables
integer :: lwork, info
real(dp), allocatable :: work(:)

n = size(Am, 1)
lwork = n
allocate(work(lwork))
call dsytrf('L', n, Am, n, ipiv, work, lwork, info)
if (info < 0) then
    print *, "The", -info, "-th argument had illegal value"
    error stop 'DSYTRF ERROR.'
end if
if (info > 0) then
    print *, "D(", info, ",", info, ") is exactly zero."
    print *, "The factorization has been completed, but the block diagonal"
    print *, "matrix D is exactly singular, and division by zero will occur"
    print *, "if it is used to solve a system of equations."
    error stop 'DSYTRF ERROR.'
end if
end subroutine

function solve_sym2(Am, bv, ipiv) result(c)
! Uses factorization of Am from solve_sym_setup to solve a symmetric system
real(dp), intent(in) :: Am(:,:)   ! system matrix: Am c = bv
real(dp), intent(in) :: bv(:)     ! source vector: Am c = bv
integer, intent(in) :: ipiv(size(Am, 1))   ! internal details
real(dp) :: c(size(bv))     ! solution vector: Am c = bv
integer :: n
! lapack variables
integer :: info

n = size(Am, 1)
c = bv
call dsytrs('L', n, 1, Am, n, ipiv, c, n, info)
if (info < 0) then
    print *, "The", -info, "-th argument had illegal value"
    error stop 'DSYTRS ERROR.'
end if
end function

end module
