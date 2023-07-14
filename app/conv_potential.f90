program conv_potential

! Compute convergence studies for Schroedinger or Dirac.
! One can do convergence with respect to any parameter.
! run ./conv_potential <study_type> <equation> <Ne/p> <potential_type> <alpha>
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
!
! <potential_type> can be,
!       0: Coulomb
!       1: Harmonic
!
! <alpha> can be 0, 1, -1 (-1 implies beta). used only for Dirac.

use graphs_potential, only: run_convergence_potential
implicit none

integer :: dirac_int, study_type, p_or_Ne, potential_type, alpha_int
character(len=128) :: arg

if (command_argument_count() /= 5) then
    print *, "./conv_potential <study_type> <equation> <Ne/p> <potential_type> <alpha>"
    error stop "Must supply 5 arguments"
end if

call get_command_argument(1, arg)
read(arg, '(i4)') study_type
call get_command_argument(2, arg)
read(arg, '(i4)') dirac_int
call get_command_argument(3, arg)
read(arg, '(i4)') p_or_Ne
call get_command_argument(4, arg)
read(arg, '(i4)') potential_type
call get_command_argument(5, arg)
read(arg, '(i4)') alpha_int

call run_convergence_potential(study_type, dirac_int, p_or_Ne, &
     potential_type, alpha_int, ".")

end program
