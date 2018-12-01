!> \file particle.f90
!! \brief Definition of a particle and four-vectors.

!> Module: particle
!> \author Simon Braß
!> \date 9.11.2018
module particle
  use ISO_FORTRAN_ENV, only: OUTPUT_UNIT
  use constants, only: default
  use vector

  implicit none

  private

  !> \brief Slim particle extending the vector4_t type.
  !! \details Missing.
  !! \note Simple wrapper around four-vector with mass.
  type, extends(vector4_t) :: particle_t
     private
     real(default) :: m = 0
   contains
     procedure :: write => particle_write
     procedure :: mass => particle_mass
     procedure :: p => particle_get_momentum !> Abbrev. for get_p.
     generic :: set_p => set_momentum, set_four_momentum
     procedure, private :: set_momentum => particle_set_momentum
     procedure, private :: set_four_momentum => particle_set_four_momentum
  end type particle_t

  interface particle_t
     module procedure :: particle_init, particle_init_at_rest
  end interface particle_t

  public :: particle_t
contains
  !> \brief Construct particle_t with mass and momentum.
  !! \details We set the mass and the momentum, but we do not check the mass shell condition.
  !! @param mass Mass of particle.
  !! @param momentum Momentum of particle.
  type(particle_t) function particle_init (mass, p) result (prt)
    real(default), intent(in) :: mass
    real(default), dimension(4), intent(in) :: p
    prt%m = mass
    call prt%vector4_t%set_time (p(1))
    call prt%vector4_t%set_space (p(2:4))
  end function particle_init

  !> Construct particle_t with mass at rest.
  !! @param mass On-shell mass of particle.
  type(particle_t) function particle_init_at_rest (mass) result (prt)
    real(default), intent(in) :: mass
    prt%m = mass
    prt%vector4_t = zero_vector
    call prt%vector4_t%set_time (mass)
  end function particle_init_at_rest

  !> \brief Write particle_t to unit.
  !! @param unit Output unit, defaults to OUTPUT_UNIT.
  subroutine particle_write (v, unit)
    class(particle_t), intent(in) :: v
    integer, intent(in), optional :: unit
    integer :: u
    u = OUTPUT_UNIT; if (present (unit)) u = unit
    write (u, "(A,1X,ES12.5)") "M =", v%m
    call v%vector4_t%write (unit)
  end subroutine particle_write

  !> \brief Get particle's mass.
  real(default) function particle_mass (prt) result (m)
    class(particle_t), intent(in) :: prt
    m = prt%m
  end function particle_mass

  !> \brief Get particle's momentum as 4-dim. array.
  !! \return p vector4_t object.
  function particle_get_momentum (prt) result (p)
    class(particle_t), intent(in) :: prt
    type(vector4_t) :: p
    p = prt%vector4_t
  end function particle_get_momentum

  !> \brief Set particle's momentum.
  !! \param x0 Time component of momentum vector.
  !! \param x  Space component of momentum vector.
  subroutine particle_set_momentum (prt, x0, x)
    class(particle_t), intent(inout) :: prt
    real(default), intent(in) :: x0
    real(default), dimension(3) :: x
    call prt%set_time (x0)
    call prt%set_space (x)
  end subroutine particle_set_momentum

  !> \brief Set particle's momentum.
  !! \param p Four-vector (vector4_t).
  subroutine particle_set_four_momentum (prt, p)
    class(particle_t), intent(inout) :: prt
    type(vector4_t), intent(in) :: p
    prt%vector4_t = p
  end subroutine particle_set_four_momentum
end module particle
