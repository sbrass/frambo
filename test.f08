program main
  use ISO_FORTRAN_ENV, only: ERROR_UNIT
  use phs
  use phs_rambo

  implicit none

  class(phs_t), allocatable :: phs_obj
  real(default), dimension(3 * 3 - 4) :: r, rp
  real(default), dimension(4) :: p
  integer :: i, j, k
  allocate (phs_rambo_t :: phs_obj)
  phs_obj = phs_rambo_t (3, [500._default, 0._default, 0._default, 0._default], [0._default, 0._default, 0._default])

  ! call random_number (r)

  ! call phs_obj%generate (r)
  ! call phs_obj%write ()

  ! select type (phs_obj)
  ! type is (phs_rambo_t)
  !    call phs_obj%invert(rp)
  ! end select
  ! print *, r
  ! print *, rp
  ! print *, r - rp

  i = 1; k = 1
  do while (i <= 100000)
     call random_number (r)
     call phs_obj%generate (r)
     k = k + 1
     if (.not. phs_obj%is_valid ())  then
        ! call phs_obj%write (ERROR_UNIT)
        cycle
     end if
     do j = 1, 3
        p = phs_obj%get_p (j)
        print *, p
     end do
     i = i + 1
  end do
  write (ERROR_UNIT, *) i, k
end program main
