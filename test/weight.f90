!> \file weight.f90
!! \brief Test weight distribution of energy for three massless particles.
program main
  use phs_rambo

  implicit none

  type(phs_rambo_t) :: phs_obj

  real(default), dimension(4) :: Q
  real(default), dimension(3) :: m
  real(default), dimension(6) :: r

  real(default), dimension(4) :: p
  real(default), dimension(3) :: energy
  real(default) :: weight

  integer, allocatable :: seed(:)
  integer :: n, i

  call random_seed(size = n)
  allocate(seed(n))
  seed = 1961991
  call random_seed(put = seed)

  Q = [500, 0, 0, 0]
  m = [0, 0, 0]

  phs_obj = phs_rambo_t (3, Q, m)

  do i = 1, 10000
     call random_number (r)
     call phs_obj%generate (r)
     weight = phs_obj%weight ()
     p = phs_obj%p(1)
     energy(1) = p(1)
     p = phs_obj%p(2)
     energy(2) = p(1)
     p = phs_obj%p(3)
     energy(3) = p(1)
     write (*, *) weight, energy
  end do
end program main
