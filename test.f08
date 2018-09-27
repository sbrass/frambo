program main
  use phs
  use phs_rambo

  implicit none

  class(phs_t), allocatable :: phs_obj
  real(default), dimension(3 * 4 - 4) :: r
  allocate (phs_rambo_t :: phs_obj)
  phs_obj = phs_rambo_t (4, [500._default, 0._default, 0._default, 0._default], [0._default, 0._default, 0._default, 0._default])

  call random_number (r)

  call phs_obj%generate (r)
  call phs_obj%write ()
end program main
