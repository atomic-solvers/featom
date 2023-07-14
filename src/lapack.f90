!> @brief Double precision LAPACK subroutines
!> @details This is the precision that LAPACK "d" routines were compiled with (typically
!> double precision, unless a special compiler option was used while compiling
!> LAPACK). This "dp" is only used in lapack.f90
!> The "d" routines data type is defined as "double precision", so
!> we make "dp" the same kind as 0.d0 ("double precision"), so
!> as long as LAPACK and this file were compiled with the same compiler options,
!> it will be consistent. (If for example all double precision is promoted to
!> quadruple precision, it will be promoted both in LAPACK and here.)
!>
!> The remaining documentation is taken from [LAPACK](http://www.netlib.org/lapack/explore-html/index.html) and can be augmented by the [quick reference guide](https://www.maths.tcd.ie/~domijank/lapack.pdf)
module lapack
implicit none

integer, parameter:: dp=kind(0.d0)

interface

    SUBROUTINE DGESV( N, NRHS, A, LDA, IPIV, B, LDB, INFO )
    import :: dp
    INTEGER            INFO, LDA, LDB, N, NRHS
    INTEGER            IPIV( * )
    REAL(dp)           A( LDA, * ), B( LDB, * )
    END SUBROUTINE

    SUBROUTINE DSYSV( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, WORK, LWORK, INFO )
    import :: dp
    CHARACTER          UPLO
    INTEGER            INFO, LDA, LDB, LWORK, N, NRHS
    INTEGER            IPIV( * )
    REAL(dp)           A( LDA, * ), B( LDB, * ), WORK( * )
    END SUBROUTINE

    SUBROUTINE DSYTRS( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, INFO )
    import :: dp
    CHARACTER          UPLO
    INTEGER            INFO, LDA, LDB, N, NRHS
    INTEGER            IPIV( * )
    REAL(dp)           A( LDA, * ), B( LDB, * )
    END SUBROUTINE

    SUBROUTINE DSYTRF( UPLO, N, A, LDA, IPIV, WORK, LWORK, INFO )
    import :: dp
    CHARACTER          UPLO
    INTEGER            INFO, LDA, LWORK, N
    INTEGER            IPIV( * )
    REAL(dp)           A( LDA, * ), WORK( * )
    END SUBROUTINE

    SUBROUTINE DSYEVD( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, IWORK, &
                       LIWORK, INFO )
    import :: dp
    CHARACTER          JOBZ, UPLO
    INTEGER            INFO, LDA, LIWORK, LWORK, N
    INTEGER            IWORK( * )
    REAL(dp)           A( LDA, * ), W( * ), WORK( * )
    END SUBROUTINE

    SUBROUTINE DSYEVX( JOBZ, RANGE, UPLO, N, A, LDA, &
                       VL, VU, IL, IU, ABSTOL, M, W, Z, LDZ, WORK, &
                       LWORK, IWORK, IFAIL, INFO )
    import :: dp
    CHARACTER          JOBZ, RANGE, UPLO
    INTEGER            IL, INFO, IU, LDA, LDZ, LWORK, M, N
    REAL(dp)           ABSTOL, VL, VU
    INTEGER            IFAIL( * ), IWORK( * )
    REAL(dp)           A( LDA, * ), W( * ), WORK( * ), &
                       Z( LDZ, * )
    END SUBROUTINE

    SUBROUTINE DSYGVD( ITYPE, JOBZ, UPLO, N, A, LDA, B, LDB, W, WORK, &
                       LWORK, IWORK, LIWORK, INFO )
    import :: dp
    CHARACTER          JOBZ, UPLO
    INTEGER            INFO, ITYPE, LDA, LDB, LIWORK, LWORK, N
    INTEGER            IWORK( * )
    REAL(dp)           A( LDA, * ), B( LDB, * ), W( * ), WORK( * )
    END SUBROUTINE

    REAL(dp) FUNCTION DLAMCH( CMACH )
    import :: dp
    CHARACTER          CMACH
    END FUNCTION

end interface

contains

end module
