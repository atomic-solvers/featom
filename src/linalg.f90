module linalg
  use types, only: dp
  use lapack, only: dsyevd, dsygvd, dgesv
  implicit none
  private
  public eigh, solve

  interface eigh
     module procedure deigh_generalized
     module procedure deigh_generalized_values
     module procedure deigh_simple
  end interface eigh

contains

  subroutine deigh_generalized(Am, Bm, lam, c)
    ! solves generalized eigen value problem for all eigenvalues and eigenvectors
    ! Am must by symmetric, Bm symmetric positive definite.
    ! Only the lower triangular part of Am and Bm is used.
    real(dp), intent(in) :: Am(:,:)   ! LHS matrix: Am c = lam Bm c
    real(dp), intent(in) :: Bm(:,:)   ! RHS matrix: Am c = lam Bm c
    real(dp), intent(out) :: lam(:)   ! eigenvalues: Am c = lam Bm c
    real(dp), intent(out) :: c(:,:)   ! eigenvectors: Am c = lam Bm c; c(i,j) = ith component of jth vec.
    integer :: n
    ! lapack variables
    integer :: lwork, liwork, info
    integer, allocatable :: iwork(:)
    real(dp), allocatable :: Bmt(:,:), work(:)

    ! solve
    n = size(Am,1)
    call assert_shape(Am, [n, n], "eigh", "Am")
    call assert_shape(Bm, [n, n], "eigh", "B")
    call assert_shape(c, [n, n], "eigh", "c")
    lwork = 1 + 6*n + 2*n**2
    liwork = 3 + 5*n
    allocate(Bmt(n,n), work(lwork), iwork(liwork))
    c = Am; Bmt = Bm  ! Bmt temporaries overwritten by dsygvd
    call dsygvd(1,'V','L',n,c,n,Bmt,n,lam,work,lwork,iwork,liwork,info)
    if (info /= 0) then
       print *, "dsygvd returned info =", info
       if (info < 0) then
          print *, "the", -info, "-th argument had an illegal value"
       else if (info <= n) then
          print *, "the algorithm failed to compute an eigenvalue while working"
          print *, "on the submatrix lying in rows and columns", 1.0_dp*info/(n+1)
          print *, "through", mod(info, n+1)
       else
          print *, "The leading minor of order ", info-n, &
               "of B is not positive definite. The factorization of B could ", &
               "not be completed and no eigenvalues or eigenvectors were computed."
       end if
       error stop 'eigh: dsygvd error'
    end if
  end subroutine deigh_generalized

  subroutine deigh_generalized_values(Am, Bm, lam)
    ! solves generalized eigen value problem for all eigenvalues
    ! Am must by symmetric, Bm symmetric positive definite.
    ! Only the upper triangular part of Am and Bm is used.
    real(dp), intent(in) :: Am(:,:)   ! LHS matrix: Am c = lam Bm c
    real(dp), intent(in) :: Bm(:,:)   ! RHS matrix: Am c = lam Bm c
    real(dp), intent(out) :: lam(:)   ! eigenvalues: Am c = lam Bm c
    integer :: n
    ! lapack variables
    integer :: lwork, liwork, info
    integer, allocatable :: iwork(:)
    real(dp), allocatable :: work(:)
    real(dp) :: c(size(Am, 1), size(Am, 2)), Bmt(size(Bm, 1), size(Bm, 2))

    ! solve
    n = size(Am,1)
    call assert_shape(Am, [n, n], "eigh", "Am")
    call assert_shape(Bm, [n, n], "eigh", "B")
    lwork = 1 + 2*n
    liwork = 1
    allocate(work(lwork), iwork(liwork))
    c = Am; Bmt = Bm  ! Bmt temporaries overwritten by dsygvd
    call dsygvd(1,'N','U',n,c,n,Bmt,n,lam,work,lwork,iwork,liwork,info)
    if (info /= 0) then
       print *, "dsygvd returned info =", info
       if (info < 0) then
          print *, "the", -info, "-th argument had an illegal value"
       else if (info <= n) then
          print *, " the algorithm failed to converge; "
          print *, info, " off-diagonal elements of an intermediate tridiagonal form "
          print *, "did not converge to zero"
       else
          print *, "The leading minor of order ", info-n, &
               "of B is not positive definite. The factorization of B could ", &
               "not be completed and no eigenvalues or eigenvectors were computed."
       end if
       error stop 'eigh: dsygvd error'
    end if
  end subroutine deigh_generalized_values

  subroutine deigh_simple(Am, lam, c)
    ! solves eigen value problem for all eigenvalues and eigenvectors
    ! Am must by symmetric
    ! Only the lower triangular part of Am is used.
    real(dp), intent(in) :: Am(:,:)   ! LHS matrix: Am c = lam c
    real(dp), intent(out) :: lam(:)   ! eigenvalues: Am c = lam c
    real(dp), intent(out) :: c(:,:)   ! eigenvectors: Am c = lam c; c(i,j) = ith component of jth vec.
    integer :: n
    ! lapack variables
    integer :: lwork, liwork, info
    integer, allocatable :: iwork(:)
    real(dp), allocatable :: work(:)

    ! solve
    n = size(Am,1)
    call assert_shape(Am, [n, n], "eigh", "Am")
    call assert_shape(c, [n, n], "eigh", "c")
    lwork = 1 + 6*n + 2*n**2
    liwork = 3 + 5*n
    allocate(work(lwork), iwork(liwork))
    c = Am
    call dsyevd('V','L',n,c,n,lam,work,lwork,iwork,liwork,info)
    if (info /= 0) then
       print *, "dsyevd returned info =", info
       if (info < 0) then
          print *, "the", -info, "-th argument had an illegal value"
       else
          print *, "the algorithm failed to compute an eigenvalue while working"
          print *, "on the submatrix lying in rows and columns", 1.0_dp*info/(n+1)
          print *, "through", mod(info, n+1)
       end if
       error stop 'eigh: dsyevd error'
    end if
  end subroutine deigh_simple

  function solve(A, b) result(x)
    ! solves a system of equations A x = b with one right hand side
    real(dp), intent(in) :: A(:,:)  ! coefficient matrix A
    real(dp), intent(in) :: b(:)  ! right-hand-side A x = b
    real(dp), allocatable :: x(:)
    ! LAPACK variables:
    real(dp), allocatable :: At(:,:), bt(:,:)
    integer :: n, info, lda
    integer, allocatable :: ipiv(:)

    n = size(A(1,:))
    lda = size(A(:, 1))  ! TODO: remove lda (which is = n!)
    call assert_shape(A, [n, n], "solve", "A")
    allocate(At(lda,n), bt(n,1), ipiv(n), x(n))
    At = A
    bt(:,1) = b(:)
    call dgesv(n, 1, At, lda, ipiv, bt, n, info)
    if(info /= 0) then
       print *, "dgesv returned info =", info
       if (info < 0) then
          print *, "the", -info, "-th argument had an illegal value"
       else
          print *, "U(", info, ",", info, ") is exactly zero; The factorization"
          print *, "has been completed, but the factor U is exactly"
          print *, "singular, so the solution could not be computed."
       end if
       error stop 'inv: dgesv error'
    endif
    x = bt(:,1)
  end function solve

  subroutine assert_shape(A, shap, routine, matname)
    ! make sure a given real matrix has a given shape
    real(dp), intent(in) :: A(:,:)
    integer, intent(in) :: shap(:)
    character(len=*) :: routine, matname

    if(any(shape(A) /= shap)) then
       print *, "In routine " // routine // " matrix " // matname // " has illegal shape ", shape(A)
       print *, "Shape should be ", shap
       error stop "Aborting due to illegal matrix operation"
    end if
  end subroutine assert_shape

end module linalg
