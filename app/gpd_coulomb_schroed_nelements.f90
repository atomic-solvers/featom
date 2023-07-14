program gpd_coulomb_schroed_nelements
use graphs_potential, only: run_convergence_potential
implicit none
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
integer, parameter, dimension(*) :: p_values = [2, 3, 4, 8, 16, 32]
integer :: idx


do idx = 1, size(p_values)
   print*, "Running Schroedinger Coulomb number of elements ", &
        2, 0, p_values(idx), 0, -1
   call run_convergence_potential(2, 0, p_values(idx), &
        0, -1, "app")
end do

end program gpd_coulomb_schroed_nelements
