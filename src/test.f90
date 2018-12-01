program main
  use ISO_FORTRAN_ENV, only: ERROR_UNIT
  use phs
  use phs_rambo

  implicit none

  class(phs_t), allocatable :: phs_obj
  integer, parameter :: n = 4
  real(default), dimension(3 * n - 4) :: r, rp
  real(default), dimension(n) :: p
  integer :: i, j, k
  allocate (phs_rambo_t :: phs_obj)
  ! phs_obj = phs_rambo_t (3, [500._default, 0._default, 0._default, 0._default], [0._default, 0._default, 0._default])
  phs_obj = phs_rambo_t (n, [500._default, 0._default, 0._default, 0._default], [0._default, 0._default, 0._default, 0._default])

  call phs_obj%write ()

  call random_number (r)

  call phs_obj%generate (r)
  call phs_obj%write ()

  select type (phs_obj)
  type is (phs_rambo_t)
     call phs_obj%invert(rp)
  end select
  call phs_obj%write ()
  print *, r
  print *, rp
  print *, r - rp

  ! i = 1; k = 1
  ! do while (i <= 10000)
  !    call random_number (r)
  !    call phs_obj%generate (r)
  !    k = k + 1
  !    if (.not. phs_obj%is_valid ())  then
  !       call phs_obj%write (ERROR_UNIT)
  !       cycle
  !    end if
  !    do j = 1, 3
  !       p = phs_obj%get_p (j)
  !       print *, p
  !    end do
  !    i = i + 1
  ! end do
  ! write (ERROR_UNIT, *) i, k
end program main
