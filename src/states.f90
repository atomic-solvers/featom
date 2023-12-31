module states

! This module lists nonrelativistic and relativistic atomic configurations.
! The nonrelativistic configurations are the same as at NIST [1] and are simply
! hardcoded in the subroutine for each atom. The relativistic configuration is
! then calculated from the nonrelativistic by splitting the occupancy according
! to the degeneracy (see the soubroutine for more info).
!
! [1] "Local-density-functional calculations of the energy of atoms,"
! S. Kotochigova, Z.H. Levine, E.L. Shirley, M.D. Stiles, and C.W. Clark,
! Phys. Rev.  A 55, 191-199 (1997).

use types, only: dp
use iso_c_binding, only: c_int
implicit none
private
public get_atomic_states_nonrel, get_atomic_states_rel, nlf2focc, &
    get_atomic_states_nonrel_focc, focc2nlf, nlsf2focc, focc2nlsf, &
    get_atomic_states_rel_focc, get_num_atomic_states_rel, &
    get_num_atomic_states_nonrel

contains

subroutine get_atomic_states_nonrel(Z, no, lo, fo)
! Returns all electron states of the form (n, l, f) for a given Z
integer, intent(in) :: Z ! atomic number
integer, intent(out), allocatable :: no(:), lo(:) ! quantum numbers "n" and "l"
real(dp), intent(out), allocatable :: fo(:) ! occupancy of the (n, l) state
! Note: sum(fo) == Z

integer :: n

select case (Z)

    case (1)
        n = 1
        allocate(no(n), lo(n), fo(n))
        no = [ 1 ]
        lo = [ 0 ]
        fo = [ 1 ]

    case (2)
        n = 1
        allocate(no(n), lo(n), fo(n))
        no = [ 1 ]
        lo = [ 0 ]
        fo = [ 2 ]

    case (3)
        n = 2
        allocate(no(n), lo(n), fo(n))
        no = [ 1, 2 ]
        lo = [ 0, 0 ]
        fo = [ 2, 1 ]

    case (4)
        n = 2
        allocate(no(n), lo(n), fo(n))
        no = [ 1, 2 ]
        lo = [ 0, 0 ]
        fo = [ 2, 2 ]

    case (5)
        n = 3
        allocate(no(n), lo(n), fo(n))
        no = [ 1, 2, 2 ]
        lo = [ 0, 0, 1 ]
        fo = [ 2, 2, 1 ]

    case (6)
        n = 3
        allocate(no(n), lo(n), fo(n))
        no = [ 1, 2, 2 ]
        lo = [ 0, 0, 1 ]
        fo = [ 2, 2, 2 ]

    case (7)
        n = 3
        allocate(no(n), lo(n), fo(n))
        no = [ 1, 2, 2 ]
        lo = [ 0, 0, 1 ]
        fo = [ 2, 2, 3 ]

    case (8)
        n = 3
        allocate(no(n), lo(n), fo(n))
        no = [ 1, 2, 2 ]
        lo = [ 0, 0, 1 ]
        fo = [ 2, 2, 4 ]

    case (9)
        n = 3
        allocate(no(n), lo(n), fo(n))
        no = [ 1, 2, 2 ]
        lo = [ 0, 0, 1 ]
        fo = [ 2, 2, 5 ]

    case (10)
        n = 3
        allocate(no(n), lo(n), fo(n))
        no = [ 1, 2, 2 ]
        lo = [ 0, 0, 1 ]
        fo = [ 2, 2, 6 ]

    case (11)
        n = 4
        allocate(no(n), lo(n), fo(n))
        no = [ 1, 2, 2, 3 ]
        lo = [ 0, 0, 1, 0 ]
        fo = [ 2, 2, 6, 1 ]

    case (12)
        n = 4
        allocate(no(n), lo(n), fo(n))
        no = [ 1, 2, 2, 3 ]
        lo = [ 0, 0, 1, 0 ]
        fo = [ 2, 2, 6, 2 ]

    case (13)
        n = 5
        allocate(no(n), lo(n), fo(n))
        no = [ 1, 2, 2, 3, 3 ]
        lo = [ 0, 0, 1, 0, 1 ]
        fo = [ 2, 2, 6, 2, 1 ]

    case (14)
        n = 5
        allocate(no(n), lo(n), fo(n))
        no = [ 1, 2, 2, 3, 3 ]
        lo = [ 0, 0, 1, 0, 1 ]
        fo = [ 2, 2, 6, 2, 2 ]

    case (15)
        n = 5
        allocate(no(n), lo(n), fo(n))
        no = [ 1, 2, 2, 3, 3 ]
        lo = [ 0, 0, 1, 0, 1 ]
        fo = [ 2, 2, 6, 2, 3 ]

    case (16)
        n = 5
        allocate(no(n), lo(n), fo(n))
        no = [ 1, 2, 2, 3, 3 ]
        lo = [ 0, 0, 1, 0, 1 ]
        fo = [ 2, 2, 6, 2, 4 ]

    case (17)
        n = 5
        allocate(no(n), lo(n), fo(n))
        no = [ 1, 2, 2, 3, 3 ]
        lo = [ 0, 0, 1, 0, 1 ]
        fo = [ 2, 2, 6, 2, 5 ]

    case (18)
        n = 5
        allocate(no(n), lo(n), fo(n))
        no = [ 1, 2, 2, 3, 3 ]
        lo = [ 0, 0, 1, 0, 1 ]
        fo = [ 2, 2, 6, 2, 6 ]

    case (19)
        n = 6
        allocate(no(n), lo(n), fo(n))
        no = [ 1, 2, 2, 3, 3, 4 ]
        lo = [ 0, 0, 1, 0, 1, 0 ]
        fo = [ 2, 2, 6, 2, 6, 1 ]

    case (20)
        n = 6
        allocate(no(n), lo(n), fo(n))
        no = [ 1, 2, 2, 3, 3, 4 ]
        lo = [ 0, 0, 1, 0, 1, 0 ]
        fo = [ 2, 2, 6, 2, 6, 2 ]

    case (21)
        n = 7
        allocate(no(n), lo(n), fo(n))
        no = [ 1, 2, 2, 3, 3, 3, 4 ]
        lo = [ 0, 0, 1, 0, 1, 2, 0 ]
        fo = [ 2, 2, 6, 2, 6, 1, 2 ]

    case (22)
        n = 7
        allocate(no(n), lo(n), fo(n))
        no = [ 1, 2, 2, 3, 3, 3, 4 ]
        lo = [ 0, 0, 1, 0, 1, 2, 0 ]
        fo = [ 2, 2, 6, 2, 6, 2, 2 ]

    case (23)
        n = 7
        allocate(no(n), lo(n), fo(n))
        no = [ 1, 2, 2, 3, 3, 3, 4 ]
        lo = [ 0, 0, 1, 0, 1, 2, 0 ]
        fo = [ 2, 2, 6, 2, 6, 3, 2 ]

    case (24)
        n = 7
        allocate(no(n), lo(n), fo(n))
        no = [ 1, 2, 2, 3, 3, 3, 4 ]
        lo = [ 0, 0, 1, 0, 1, 2, 0 ]
        fo = [ 2, 2, 6, 2, 6, 5, 1 ]

    case (25)
        n = 7
        allocate(no(n), lo(n), fo(n))
        no = [ 1, 2, 2, 3, 3, 3, 4 ]
        lo = [ 0, 0, 1, 0, 1, 2, 0 ]
        fo = [ 2, 2, 6, 2, 6, 5, 2 ]

    case (26)
        n = 7
        allocate(no(n), lo(n), fo(n))
        no = [ 1, 2, 2, 3, 3, 3, 4 ]
        lo = [ 0, 0, 1, 0, 1, 2, 0 ]
        fo = [ 2, 2, 6, 2, 6, 6, 2 ]

    case (27)
        n = 7
        allocate(no(n), lo(n), fo(n))
        no = [ 1, 2, 2, 3, 3, 3, 4 ]
        lo = [ 0, 0, 1, 0, 1, 2, 0 ]
        fo = [ 2, 2, 6, 2, 6, 7, 2 ]

    case (28)
        n = 7
        allocate(no(n), lo(n), fo(n))
        no = [ 1, 2, 2, 3, 3, 3, 4 ]
        lo = [ 0, 0, 1, 0, 1, 2, 0 ]
        fo = [ 2, 2, 6, 2, 6, 8, 2 ]

    case (29)
        n = 7
        allocate(no(n), lo(n), fo(n))
        no = [ 1, 2, 2, 3, 3, 3, 4 ]
        lo = [ 0, 0, 1, 0, 1, 2, 0 ]
        fo = [ 2, 2, 6, 2, 6, 10, 1 ]

    case (30)
        n = 7
        allocate(no(n), lo(n), fo(n))
        no = [ 1, 2, 2, 3, 3, 3, 4 ]
        lo = [ 0, 0, 1, 0, 1, 2, 0 ]
        fo = [ 2, 2, 6, 2, 6, 10, 2 ]

    case (31)
        n = 8
        allocate(no(n), lo(n), fo(n))
        no = [ 1, 2, 2, 3, 3, 3, 4, 4 ]
        lo = [ 0, 0, 1, 0, 1, 2, 0, 1 ]
        fo = [ 2, 2, 6, 2, 6, 10, 2, 1 ]

    case (32)
        n = 8
        allocate(no(n), lo(n), fo(n))
        no = [ 1, 2, 2, 3, 3, 3, 4, 4 ]
        lo = [ 0, 0, 1, 0, 1, 2, 0, 1 ]
        fo = [ 2, 2, 6, 2, 6, 10, 2, 2 ]

    case (33)
        n = 8
        allocate(no(n), lo(n), fo(n))
        no = [ 1, 2, 2, 3, 3, 3, 4, 4 ]
        lo = [ 0, 0, 1, 0, 1, 2, 0, 1 ]
        fo = [ 2, 2, 6, 2, 6, 10, 2, 3 ]

    case (34)
        n = 8
        allocate(no(n), lo(n), fo(n))
        no = [ 1, 2, 2, 3, 3, 3, 4, 4 ]
        lo = [ 0, 0, 1, 0, 1, 2, 0, 1 ]
        fo = [ 2, 2, 6, 2, 6, 10, 2, 4 ]

    case (35)
        n = 8
        allocate(no(n), lo(n), fo(n))
        no = [ 1, 2, 2, 3, 3, 3, 4, 4 ]
        lo = [ 0, 0, 1, 0, 1, 2, 0, 1 ]
        fo = [ 2, 2, 6, 2, 6, 10, 2, 5 ]

    case (36)
        n = 8
        allocate(no(n), lo(n), fo(n))
        no = [ 1, 2, 2, 3, 3, 3, 4, 4 ]
        lo = [ 0, 0, 1, 0, 1, 2, 0, 1 ]
        fo = [ 2, 2, 6, 2, 6, 10, 2, 6 ]

    case (37)
        n = 9
        allocate(no(n), lo(n), fo(n))
        no = [ 1, 2, 2, 3, 3, 3, 4, 4, 5 ]
        lo = [ 0, 0, 1, 0, 1, 2, 0, 1, 0 ]
        fo = [ 2, 2, 6, 2, 6, 10, 2, 6, 1 ]

    case (38)
        n = 9
        allocate(no(n), lo(n), fo(n))
        no = [ 1, 2, 2, 3, 3, 3, 4, 4, 5 ]
        lo = [ 0, 0, 1, 0, 1, 2, 0, 1, 0 ]
        fo = [ 2, 2, 6, 2, 6, 10, 2, 6, 2 ]

    case (39)
        n = 10
        allocate(no(n), lo(n), fo(n))
        no = [ 1, 2, 2, 3, 3, 3, 4, 4, 4, 5 ]
        lo = [ 0, 0, 1, 0, 1, 2, 0, 1, 2, 0 ]
        fo = [ 2, 2, 6, 2, 6, 10, 2, 6, 1, 2 ]

    case (40)
        n = 10
        allocate(no(n), lo(n), fo(n))
        no = [ 1, 2, 2, 3, 3, 3, 4, 4, 4, 5 ]
        lo = [ 0, 0, 1, 0, 1, 2, 0, 1, 2, 0 ]
        fo = [ 2, 2, 6, 2, 6, 10, 2, 6, 2, 2 ]

    case (41)
        n = 10
        allocate(no(n), lo(n), fo(n))
        no = [ 1, 2, 2, 3, 3, 3, 4, 4, 4, 5 ]
        lo = [ 0, 0, 1, 0, 1, 2, 0, 1, 2, 0 ]
        fo = [ 2, 2, 6, 2, 6, 10, 2, 6, 4, 1 ]

    case (42)
        n = 10
        allocate(no(n), lo(n), fo(n))
        no = [ 1, 2, 2, 3, 3, 3, 4, 4, 4, 5 ]
        lo = [ 0, 0, 1, 0, 1, 2, 0, 1, 2, 0 ]
        fo = [ 2, 2, 6, 2, 6, 10, 2, 6, 5, 1 ]

    case (43)
        n = 10
        allocate(no(n), lo(n), fo(n))
        no = [ 1, 2, 2, 3, 3, 3, 4, 4, 4, 5 ]
        lo = [ 0, 0, 1, 0, 1, 2, 0, 1, 2, 0 ]
        fo = [ 2, 2, 6, 2, 6, 10, 2, 6, 5, 2 ]

    case (44)
        n = 10
        allocate(no(n), lo(n), fo(n))
        no = [ 1, 2, 2, 3, 3, 3, 4, 4, 4, 5 ]
        lo = [ 0, 0, 1, 0, 1, 2, 0, 1, 2, 0 ]
        fo = [ 2, 2, 6, 2, 6, 10, 2, 6, 7, 1 ]

    case (45)
        n = 10
        allocate(no(n), lo(n), fo(n))
        no = [ 1, 2, 2, 3, 3, 3, 4, 4, 4, 5 ]
        lo = [ 0, 0, 1, 0, 1, 2, 0, 1, 2, 0 ]
        fo = [ 2, 2, 6, 2, 6, 10, 2, 6, 8, 1 ]

    case (46)
        n = 9
        allocate(no(n), lo(n), fo(n))
        no = [ 1, 2, 2, 3, 3, 3, 4, 4, 4 ]
        lo = [ 0, 0, 1, 0, 1, 2, 0, 1, 2 ]
        fo = [ 2, 2, 6, 2, 6, 10, 2, 6, 10 ]

    case (47)
        n = 10
        allocate(no(n), lo(n), fo(n))
        no = [ 1, 2, 2, 3, 3, 3, 4, 4, 4, 5 ]
        lo = [ 0, 0, 1, 0, 1, 2, 0, 1, 2, 0 ]
        fo = [ 2, 2, 6, 2, 6, 10, 2, 6, 10, 1 ]

    case (48)
        n = 10
        allocate(no(n), lo(n), fo(n))
        no = [ 1, 2, 2, 3, 3, 3, 4, 4, 4, 5 ]
        lo = [ 0, 0, 1, 0, 1, 2, 0, 1, 2, 0 ]
        fo = [ 2, 2, 6, 2, 6, 10, 2, 6, 10, 2 ]

    case (49)
        n = 11
        allocate(no(n), lo(n), fo(n))
        no = [ 1, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5 ]
        lo = [ 0, 0, 1, 0, 1, 2, 0, 1, 2, 0, 1 ]
        fo = [ 2, 2, 6, 2, 6, 10, 2, 6, 10, 2, 1 ]

    case (50)
        n = 11
        allocate(no(n), lo(n), fo(n))
        no = [ 1, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5 ]
        lo = [ 0, 0, 1, 0, 1, 2, 0, 1, 2, 0, 1 ]
        fo = [ 2, 2, 6, 2, 6, 10, 2, 6, 10, 2, 2 ]

    case (51)
        n = 11
        allocate(no(n), lo(n), fo(n))
        no = [ 1, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5 ]
        lo = [ 0, 0, 1, 0, 1, 2, 0, 1, 2, 0, 1 ]
        fo = [ 2, 2, 6, 2, 6, 10, 2, 6, 10, 2, 3 ]

    case (52)
        n = 11
        allocate(no(n), lo(n), fo(n))
        no = [ 1, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5 ]
        lo = [ 0, 0, 1, 0, 1, 2, 0, 1, 2, 0, 1 ]
        fo = [ 2, 2, 6, 2, 6, 10, 2, 6, 10, 2, 4 ]

    case (53)
        n = 11
        allocate(no(n), lo(n), fo(n))
        no = [ 1, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5 ]
        lo = [ 0, 0, 1, 0, 1, 2, 0, 1, 2, 0, 1 ]
        fo = [ 2, 2, 6, 2, 6, 10, 2, 6, 10, 2, 5 ]

    case (54)
        n = 11
        allocate(no(n), lo(n), fo(n))
        no = [ 1, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5 ]
        lo = [ 0, 0, 1, 0, 1, 2, 0, 1, 2, 0, 1 ]
        fo = [ 2, 2, 6, 2, 6, 10, 2, 6, 10, 2, 6 ]

    case (55)
        n = 12
        allocate(no(n), lo(n), fo(n))
        no = [ 1, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 6 ]
        lo = [ 0, 0, 1, 0, 1, 2, 0, 1, 2, 0, 1, 0 ]
        fo = [ 2, 2, 6, 2, 6, 10, 2, 6, 10, 2, 6, 1 ]

    case (56)
        n = 12
        allocate(no(n), lo(n), fo(n))
        no = [ 1, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 6 ]
        lo = [ 0, 0, 1, 0, 1, 2, 0, 1, 2, 0, 1, 0 ]
        fo = [ 2, 2, 6, 2, 6, 10, 2, 6, 10, 2, 6, 2 ]

    case (57)
        n = 13
        allocate(no(n), lo(n), fo(n))
        no = [ 1, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5, 6 ]
        lo = [ 0, 0, 1, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0 ]
        fo = [ 2, 2, 6, 2, 6, 10, 2, 6, 10, 2, 6, 1, 2 ]

    case (58)
        n = 14
        allocate(no(n), lo(n), fo(n))
        no = [ 1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 6 ]
        lo = [ 0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 0 ]
        fo = [ 2, 2, 6, 2, 6, 10, 2, 6, 10, 1, 2, 6, 1, 2 ]

    case (59)
        n = 13
        allocate(no(n), lo(n), fo(n))
        no = [ 1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 6 ]
        lo = [ 0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 0 ]
        fo = [ 2, 2, 6, 2, 6, 10, 2, 6, 10, 3, 2, 6, 2 ]

    case (60)
        n = 13
        allocate(no(n), lo(n), fo(n))
        no = [ 1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 6 ]
        lo = [ 0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 0 ]
        fo = [ 2, 2, 6, 2, 6, 10, 2, 6, 10, 4, 2, 6, 2 ]

    case (61)
        n = 13
        allocate(no(n), lo(n), fo(n))
        no = [ 1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 6 ]
        lo = [ 0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 0 ]
        fo = [ 2, 2, 6, 2, 6, 10, 2, 6, 10, 5, 2, 6, 2 ]

    case (62)
        n = 13
        allocate(no(n), lo(n), fo(n))
        no = [ 1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 6 ]
        lo = [ 0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 0 ]
        fo = [ 2, 2, 6, 2, 6, 10, 2, 6, 10, 6, 2, 6, 2 ]

    case (63)
        n = 13
        allocate(no(n), lo(n), fo(n))
        no = [ 1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 6 ]
        lo = [ 0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 0 ]
        fo = [ 2, 2, 6, 2, 6, 10, 2, 6, 10, 7, 2, 6, 2 ]

    case (64)
        n = 14
        allocate(no(n), lo(n), fo(n))
        no = [ 1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 6 ]
        lo = [ 0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 0 ]
        fo = [ 2, 2, 6, 2, 6, 10, 2, 6, 10, 7, 2, 6, 1, 2 ]

    case (65)
        n = 13
        allocate(no(n), lo(n), fo(n))
        no = [ 1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 6 ]
        lo = [ 0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 0 ]
        fo = [ 2, 2, 6, 2, 6, 10, 2, 6, 10, 9, 2, 6, 2 ]

    case (66)
        n = 13
        allocate(no(n), lo(n), fo(n))
        no = [ 1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 6 ]
        lo = [ 0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 0 ]
        fo = [ 2, 2, 6, 2, 6, 10, 2, 6, 10, 10, 2, 6, 2 ]

    case (67)
        n = 13
        allocate(no(n), lo(n), fo(n))
        no = [ 1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 6 ]
        lo = [ 0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 0 ]
        fo = [ 2, 2, 6, 2, 6, 10, 2, 6, 10, 11, 2, 6, 2 ]

    case (68)
        n = 13
        allocate(no(n), lo(n), fo(n))
        no = [ 1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 6 ]
        lo = [ 0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 0 ]
        fo = [ 2, 2, 6, 2, 6, 10, 2, 6, 10, 12, 2, 6, 2 ]

    case (69)
        n = 13
        allocate(no(n), lo(n), fo(n))
        no = [ 1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 6 ]
        lo = [ 0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 0 ]
        fo = [ 2, 2, 6, 2, 6, 10, 2, 6, 10, 13, 2, 6, 2 ]

    case (70)
        n = 13
        allocate(no(n), lo(n), fo(n))
        no = [ 1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 6 ]
        lo = [ 0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 0 ]
        fo = [ 2, 2, 6, 2, 6, 10, 2, 6, 10, 14, 2, 6, 2 ]

    case (71)
        n = 14
        allocate(no(n), lo(n), fo(n))
        no = [ 1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 6 ]
        lo = [ 0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 0 ]
        fo = [ 2, 2, 6, 2, 6, 10, 2, 6, 10, 14, 2, 6, 1, 2 ]

    case (72)
        n = 14
        allocate(no(n), lo(n), fo(n))
        no = [ 1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 6 ]
        lo = [ 0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 0 ]
        fo = [ 2, 2, 6, 2, 6, 10, 2, 6, 10, 14, 2, 6, 2, 2 ]

    case (73)
        n = 14
        allocate(no(n), lo(n), fo(n))
        no = [ 1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 6 ]
        lo = [ 0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 0 ]
        fo = [ 2, 2, 6, 2, 6, 10, 2, 6, 10, 14, 2, 6, 3, 2 ]

    case (74)
        n = 14
        allocate(no(n), lo(n), fo(n))
        no = [ 1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 6 ]
        lo = [ 0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 0 ]
        fo = [ 2, 2, 6, 2, 6, 10, 2, 6, 10, 14, 2, 6, 4, 2 ]

    case (75)
        n = 14
        allocate(no(n), lo(n), fo(n))
        no = [ 1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 6 ]
        lo = [ 0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 0 ]
        fo = [ 2, 2, 6, 2, 6, 10, 2, 6, 10, 14, 2, 6, 5, 2 ]

    case (76)
        n = 14
        allocate(no(n), lo(n), fo(n))
        no = [ 1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 6 ]
        lo = [ 0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 0 ]
        fo = [ 2, 2, 6, 2, 6, 10, 2, 6, 10, 14, 2, 6, 6, 2 ]

    case (77)
        n = 14
        allocate(no(n), lo(n), fo(n))
        no = [ 1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 6 ]
        lo = [ 0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 0 ]
        fo = [ 2, 2, 6, 2, 6, 10, 2, 6, 10, 14, 2, 6, 7, 2 ]

    case (78)
        n = 14
        allocate(no(n), lo(n), fo(n))
        no = [ 1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 6 ]
        lo = [ 0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 0 ]
        fo = [ 2, 2, 6, 2, 6, 10, 2, 6, 10, 14, 2, 6, 9, 1 ]

    case (79)
        n = 14
        allocate(no(n), lo(n), fo(n))
        no = [ 1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 6 ]
        lo = [ 0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 0 ]
        fo = [ 2, 2, 6, 2, 6, 10, 2, 6, 10, 14, 2, 6, 10, 1 ]

    case (80)
        n = 14
        allocate(no(n), lo(n), fo(n))
        no = [ 1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 6 ]
        lo = [ 0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 0 ]
        fo = [ 2, 2, 6, 2, 6, 10, 2, 6, 10, 14, 2, 6, 10, 2 ]

    case (81)
        n = 15
        allocate(no(n), lo(n), fo(n))
        no = [ 1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 6, 6 ]
        lo = [ 0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 0, 1 ]
        fo = [ 2, 2, 6, 2, 6, 10, 2, 6, 10, 14, 2, 6, 10, 2, 1 ]

    case (82)
        n = 15
        allocate(no(n), lo(n), fo(n))
        no = [ 1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 6, 6 ]
        lo = [ 0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 0, 1 ]
        fo = [ 2, 2, 6, 2, 6, 10, 2, 6, 10, 14, 2, 6, 10, 2, 2 ]

    case (83)
        n = 15
        allocate(no(n), lo(n), fo(n))
        no = [ 1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 6, 6 ]
        lo = [ 0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 0, 1 ]
        fo = [ 2, 2, 6, 2, 6, 10, 2, 6, 10, 14, 2, 6, 10, 2, 3 ]

    case (84)
        n = 15
        allocate(no(n), lo(n), fo(n))
        no = [ 1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 6, 6 ]
        lo = [ 0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 0, 1 ]
        fo = [ 2, 2, 6, 2, 6, 10, 2, 6, 10, 14, 2, 6, 10, 2, 4 ]

    case (85)
        n = 15
        allocate(no(n), lo(n), fo(n))
        no = [ 1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 6, 6 ]
        lo = [ 0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 0, 1 ]
        fo = [ 2, 2, 6, 2, 6, 10, 2, 6, 10, 14, 2, 6, 10, 2, 5 ]

    case (86)
        n = 15
        allocate(no(n), lo(n), fo(n))
        no = [ 1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 6, 6 ]
        lo = [ 0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 0, 1 ]
        fo = [ 2, 2, 6, 2, 6, 10, 2, 6, 10, 14, 2, 6, 10, 2, 6 ]

    case (87)
        n = 16
        allocate(no(n), lo(n), fo(n))
        no = [ 1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 6, 6, 7 ]
        lo = [ 0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 0, 1, 0 ]
        fo = [ 2, 2, 6, 2, 6, 10, 2, 6, 10, 14, 2, 6, 10, 2, 6, 1 ]

    case (88)
        n = 16
        allocate(no(n), lo(n), fo(n))
        no = [ 1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 6, 6, 7 ]
        lo = [ 0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 0, 1, 0 ]
        fo = [ 2, 2, 6, 2, 6, 10, 2, 6, 10, 14, 2, 6, 10, 2, 6, 2 ]

    case (89)
        n = 17
        allocate(no(n), lo(n), fo(n))
        no = [ 1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 6, 6, 6, 7 ]
        lo = [ 0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 0, 1, 2, 0 ]
        fo = [ 2, 2, 6, 2, 6, 10, 2, 6, 10, 14, 2, 6, 10, 2, 6, 1, 2 ]

    case (90)
        n = 17
        allocate(no(n), lo(n), fo(n))
        no = [ 1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 6, 6, 6, 7 ]
        lo = [ 0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 0, 1, 2, 0 ]
        fo = [ 2, 2, 6, 2, 6, 10, 2, 6, 10, 14, 2, 6, 10, 2, 6, 2, 2 ]

    case (91)
        n = 18
        allocate(no(n), lo(n), fo(n))
        no = [ 1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 6, 6, 6, 7 ]
        lo = [ 0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 0 ]
        fo = [ 2, 2, 6, 2, 6, 10, 2, 6, 10, 14, 2, 6, 10, 2, 2, 6, 1, 2 ]

    case (92)
        n = 18
        allocate(no(n), lo(n), fo(n))
        no = [ 1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 6, 6, 6, 7 ]
        lo = [ 0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 0 ]
        fo = [ 2, 2, 6, 2, 6, 10, 2, 6, 10, 14, 2, 6, 10, 3, 2, 6, 1, 2 ]

    case default
        print *, "Z =", Z
        error stop "Z not supported."
end select
end subroutine

subroutine get_atomic_states_rel(Z, no, lo, so, fo)
! Returns all relativistic electron states of the form
! (n, l, s, f) for a given Z
! We just take the non-relativistic occupation fo(:) for each l-level and
! in the relativistic case, each l-level for l > 0 is split into j=l+1/2 (so=1,
! degeneracy 2l+2) and j=l-1/2 (so=0, degeneracy 2l) sublevels. The occupancy
! is split according to the degeneracy, so fo*(2l+2)/(2l+2+2l) goes into the
! so=1 sublevel and fo*2l/(2l+2+2l) goes into the so=0 sublevel.
integer, intent(in) :: Z ! atomic number
integer, intent(out), allocatable :: no(:), lo(:) ! quantum numbers "n" and "l"
integer, intent(out), allocatable :: so(:) ! spin (1..up or 0..down)
real(dp), intent(out), allocatable :: fo(:) ! occupancy of the (n, l, s) state
! Note: sum(fo) == Z

integer, allocatable :: no2(:), lo2(:)
real(dp), allocatable :: fo2(:)
! This needs to hold all atomic states (29 for Uranium), so 100 should be a
! very safe value for all purposes (the subroutine will ask to increase this
! value if it is too low):
integer, parameter :: MAX_STATES=100
integer :: no_tmp(MAX_STATES), lo_tmp(MAX_STATES), so_tmp(MAX_STATES)
real(dp) :: fo_tmp(MAX_STATES)
integer :: i, j
call get_atomic_states_nonrel(Z, no2, lo2, fo2)
j = 1
do i = 1, size(no2)
    ! In the l > 0 case below, we are filling "j" and "j+1" elements of the
    ! arrays, so we need to make sure we fit into the preallocated arrays:
    if (j + 1 > MAX_STATES) then
        error stop "get_atomic_states_rel: increase MAX_STATES"
    end if
    if (lo2(i) == 0) then
        no_tmp(j) = no2(i)
        lo_tmp(j) = lo2(i)
        fo_tmp(j) = fo2(i)
        so_tmp(j) = 1
    else
        no_tmp(j) = no2(i)
        lo_tmp(j) = lo2(i)
        fo_tmp(j) = fo2(i) * (2*lo2(i)) / (2*(2*lo2(i)+1))
        so_tmp(j) = 0
        j = j + 1
        no_tmp(j) = no2(i)
        lo_tmp(j) = lo2(i)
        fo_tmp(j) = fo2(i) * (2*lo2(i)+2) / (2*(2*lo2(i)+1))
        so_tmp(j) = 1
    end if
    j = j + 1
end do
j = j - 1
allocate(no(j), lo(j), fo(j), so(j))
no(:) = no_tmp(:j)
lo(:) = lo_tmp(:j)
fo(:) = fo_tmp(:j)
so(:) = so_tmp(:j)
end subroutine

subroutine nlf2focc(no, lo, fo, focc)
! Converts the (no, lo, fo) triplet into a single focc(:,:) array
integer, intent(in) :: no(:), lo(:)
real(dp), intent(in) :: fo(:)
real(dp), intent(out), allocatable :: focc(:,:)
integer :: Lmax, Nmax, i
Lmax = maxval(lo)
Nmax = maxval(no)
allocate(focc(Nmax, 0:Lmax))
focc = 0
do i = 1, size(no)
    focc(no(i)-lo(i),lo(i)) = fo(i)
end do
end subroutine

subroutine focc2nlf(focc, no, lo, fo)
! Converts the focc(:,:) array into the (no, lo, fo) triplet
real(dp), intent(in) :: focc(:,0:)
integer, intent(out), allocatable :: no(:), lo(:)
real(dp), intent(out), allocatable :: fo(:)
integer :: Lmax, Nmax, n, l, idx
Lmax = ubound(focc,2)
Nmax = size(focc,1)
n = count(focc > 0)
allocate(no(n), lo(n), fo(n))
idx = 0
do n = 1, Nmax
    do l = 0, min(n-1, Lmax)
        if (focc(n-l,l) > 0) then
            idx = idx + 1
            no(idx) = n
            lo(idx) = l
            fo(idx) = focc(n-l,l)
        end if
    end do
end do
end subroutine

subroutine nlsf2focc(no, lo, so, fo, focc)
! Converts the (no, lo, fo) triplet into a single focc(:,:) array
integer, intent(in) :: no(:), lo(:), so(:)
real(dp), intent(in) :: fo(:)
real(dp), intent(out), allocatable :: focc(:,:)
integer :: kappa(size(no))
integer :: Nmax, i
Nmax = maxval(no)
where (so > 0)
    kappa = -lo-1
else where
    kappa = lo
end where
allocate(focc(Nmax, minval(kappa):maxval(kappa)))
focc = 0
do i = 1, size(no)
    focc(no(i)-lo(i),kappa(i)) = fo(i)
end do
end subroutine

subroutine focc2nlsf(kappa_min, focc, no, lo, so, fo)
! Converts the focc(:,:) array into the (no, lo, fo) triplet
integer, intent(in) :: kappa_min
real(dp), intent(in) :: focc(:,kappa_min:)
integer, intent(out), allocatable :: no(:), lo(:), so(:)
real(dp), intent(out), allocatable :: fo(:)
integer :: Nmax, n, l, idx, kappa
Nmax = size(focc,1)
n = count(focc > 0)
allocate(no(n), lo(n), so(n), fo(n))
idx = 0
do n = 1, Nmax
    do l = 0, n-1
        kappa = l
        if (kappa > 0 .and. kappa <= ubound(focc,2)) then
            if (focc(n-l,kappa) > 0) then
                idx = idx + 1
                no(idx) = n
                lo(idx) = l
                so(idx) = 0
                fo(idx) = focc(n-l,kappa)
            end if
        end if

        kappa = -l-1
        if (kappa >= lbound(focc,2)) then
            if (focc(n-l,kappa) > 0) then
                idx = idx + 1
                no(idx) = n
                lo(idx) = l
                so(idx) = 1
                fo(idx) = focc(n-l,kappa)
            end if
        end if
    end do
end do
end subroutine

subroutine get_atomic_states_nonrel_focc(Z, focc)
integer, intent(in) :: Z ! atomic number
real(dp), intent(out), allocatable :: focc(:,:)
integer, allocatable :: no(:), lo(:)
real(dp), allocatable :: fo(:)
call get_atomic_states_nonrel(Z, no, lo, fo)
call nlf2focc(no, lo, fo, focc)
end subroutine

subroutine get_atomic_states_rel_focc(Z, focc)
integer, intent(in) :: Z ! atomic number
real(dp), intent(out), allocatable :: focc(:,:)
integer, allocatable :: no(:), lo(:), so(:)
real(dp), allocatable :: fo(:)
call get_atomic_states_rel(Z, no, lo, so, fo)
call nlsf2focc(no, lo, so, fo, focc)
end subroutine

integer function get_atom_orb_nonrel(Z) result(n_orb)
! Returns the number of orbitals for the given atom
integer, intent(in) :: Z
integer, allocatable :: no(:), lo(:)
real(dp), allocatable :: fo(:)
call get_atomic_states_nonrel(Z, no, lo, fo)
n_orb = size(no)
end function

integer function get_atom_orb_rel(Z) result(n_orb)
! Returns the number of orbitals for the given atom
integer, intent(in) :: Z
integer, allocatable :: no(:), lo(:), so(:)
real(dp), allocatable :: fo(:)
call get_atomic_states_rel(Z, no, lo, so, fo)
n_orb = size(no)
end function


integer(c_int) function get_num_atomic_states_rel(Z) bind(c)
integer(c_int), intent(in) :: Z
get_num_atomic_states_rel = get_atom_orb_rel(Z)
end function

integer(c_int) function get_num_atomic_states_nonrel(Z) bind(c)
integer(c_int), intent(in) :: Z
get_num_atomic_states_nonrel = get_atom_orb_nonrel(Z)
end function

end module
