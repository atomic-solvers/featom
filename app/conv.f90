program conv

! Compute convergence studies for Schroedinger or Dirac.
! One can do convergence with respect to any parameter.
! run ./conv <study_type> <equation> <Ne/p>
!
! <study_type> can be,
!       0: error as p is varied
!       1: error as rmax is varied
!       2: error as Ne is varied
!
! <equation> can be,
!       0: Schroedinger
!       1: Dirac
!
! For <study_type>
!       0, 1: 3rd parameter = Ne (Number of elements)
!       2   : 3rd parameter = p  (Polynomial order)

use graphs, only: run_convergence
implicit none

integer :: dirac_int, study_type, p_or_Ne
character(len=128) :: arg

if (command_argument_count() /= 3) then
    print *, "./conv <study_type> <equation> <Ne/p>"
    error stop "Must supply 3 arguments"
end if

call get_command_argument(1, arg)
read(arg, '(i4)') study_type
call get_command_argument(2, arg)
read(arg, '(i4)') dirac_int
call get_command_argument(3, arg)
read(arg, '(i4)') p_or_Ne

call run_convergence(study_type, dirac_int, p_or_Ne, ".")

end program
