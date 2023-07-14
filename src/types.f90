module types
implicit none
private
public sp, dp, hp, qp, ivector, dvector, zvector

integer, parameter :: dp=kind(0.d0), &             ! double precision
                      hp=selected_real_kind(15), & ! high precision
                      qp=selected_real_kind(32), & ! quadruple precision
                      sp = kind(0.)                ! single precision

type ivector                       ! allocatable integer vector
   integer, pointer :: vec(:) => null()
end type

type dvector                       ! allocatable real double precision vector
   real(dp), pointer :: vec(:) => null()
end type

type zvector                       ! allocatable complex double precision vector
   complex(dp), pointer :: vec(:) => null()
end type

end module
