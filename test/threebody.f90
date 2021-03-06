!> \file threebody.f90
!! \brief Test three-body decay at fixed energy for massive particles.
!! \details We test different setups:
!! 1. Equal mass for all particles,
!! 2. One heavy, three massless particles,
!! 3. Completely different masses (through different scales).
program main
  use phs_rambo

  implicit none

  type(phs_rambo_t) :: phs_obj

  real(default), dimension(4) :: Q
  real(default), dimension(3) :: m
  real(default), dimension(6) :: r

  integer, allocatable :: seed(:)
  integer :: n

  call random_seed(size = n)
  allocate(seed(n))
  seed = 1961991
  call random_seed(put = seed)

  Q = [500, 0, 0, 0]
  m = [100, 100, 100]
  call random_number(r)

  phs_obj = phs_rambo_t (3, Q, m)

  call phs_obj%generate (r)

  call phs_obj%write ()

  ! 2. case
  m = [100, 0, 0]
  call random_number(r)

  phs_obj = phs_rambo_t (3, Q, m)

  call phs_obj%generate (r)

  call phs_obj%write ()

  ! 3. case
  m = [1, 10, 100]
  call random_number(r)

  phs_obj = phs_rambo_t (3, Q, m)

  call phs_obj%generate (r)

  call phs_obj%write ()
end program main
