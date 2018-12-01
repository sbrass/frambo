!> \file twobody.f90
!! \brief Test two-body decay at fixed energy for massless particles.
!! \details We test the two-body decay which only free parameters are the two
!! angles between the particles. (A massive case is redundant.)
program main
  use phs_rambo

  implicit none

  type(phs_rambo_t) :: phs_obj

  real(default), dimension(4) :: Q
  real(default), dimension(2) :: m
  real(default), dimension(2) :: r

  integer, allocatable :: seed(:)
  integer :: n

  call random_seed(size = n)
  allocate(seed(n))
  seed = 1961991
  call random_seed(put = seed)

  Q = [500, 0, 0, 0]
  m = [0, 0]

  phs_obj = phs_rambo_t (2, Q, m)

  call random_number(r)
  call phs_obj%generate (r)

  call phs_obj%write ()
end program main
