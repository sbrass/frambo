!> \file phs.f90
!! \brief Basic phase-space type with number of particles, total momentum, particle masses.

!> Module: phs
!> \author Simon Braß
!> \date 07.11.2018
module phs
  use ISO_FORTRAN_ENV, only: OUTPUT_UNIT
  use constants
  use vector
  use particle

  implicit none

  !> \brief Basic phase-space type.
  !! \details We store the number of particles, the total momentum of the system, the masses of the particles.
  !! We split the computation of the overall phase-space weight into
  !! - energy-dependent part, flux and jacobian
  !! - energy-independent part, volume
  !!
  !! The overall weight is the product of flux, volume and jacobian.
  !!
  !! The flux \f$F\f$ is defined as
  !! \f[
  !!    F = 4 \sqrt{\left( p_a \cdot p_b \right)^2 - m_a^2 m_b^2}.
  !! \f]
  !! We absorb an additional prefactor of \f$(2\pi)^4\f$ coming from the definition of lorentz-invariant matrix-element into the flux.
  !! The energy-independent part of the phase-space volume is given as for RAMBO as
  !! \f[
  !!   \Phi_n = \frac{1}{(2\pi)^{3n}} \left( \frac{\pi}{2} \right)^{n - 1} \frac{1}{(n - 1)!(n - 2)!}.
  !! \f]
  type, abstract :: phs_t
     integer :: n_particles = 0
     type(vector4_t) :: total_p = zero_vector
     real(default) :: flux = 0
     real(default) :: volume = 0
     real(default) :: jacobian = 0
     type(particle_t), dimension(:), allocatable :: prt
   contains
     procedure :: init => phs_init
     procedure(phs_write), deferred :: write
     procedure :: base_write => phs_base_write
     procedure :: p => phs_p
     procedure(phs_generate), deferred :: generate
     procedure :: weight => phs_get_weight
  end type phs_t

  abstract interface
     !> Generate phase-space using random numbers.
     !! @param r_in \f$3n - 4\f$ random numbers.
     subroutine phs_generate (phs, r_in)
       import
       class(phs_t), intent(inout) :: phs
       real(default), dimension(3 * phs%n_particles - 4), intent(in) :: r_in
     end subroutine phs_generate

     !> Write phase-space object to unit.
     !! @param unit Fortran file unit, default shall be OUTPUT_UNIT.
     subroutine phs_write (phs, unit)
       import
       class(phs_t), intent(in) :: phs
       integer, intent(in), optional :: unit
     end subroutine phs_write
  end interface

contains
  !> Construct phs_t requiring number of particles, total momenta and particles masses.
  !! @param n_particles Number of phase-space particles.
  !! @param total_momentum Total momentum \f$Q\f$.
  !! @param mass Masses of phase-space particles.
  subroutine phs_init (phs, n_particles, total_momentum, mass)
    class(phs_t), intent(inout) :: phs
    integer, intent(in) :: n_particles
    real(default), dimension(4) :: total_momentum
    real(default), dimension(n_particles) :: mass
    integer :: i
    phs%n_particles = n_particles
    phs%total_p = vector4_t (total_momentum)
    allocate (phs%prt(n_particles))
    do i = 1, n_particles
       phs%prt(i) = particle_t (mass(i))
    end do
  end subroutine phs_init

  !> Write basic output.
  !! @param unit
  subroutine phs_base_write (phs, unit)
    class(phs_t), intent(in) :: phs
    integer, intent(in), optional :: unit
    integer :: u, i
    type(vector4_t) :: total_p
    u = OUTPUT_UNIT; if (present (unit)) u = unit
    write (u, "(80(A))") "================================================================================"
    write (u, "(A)") "Phasespace basic type:"
    write (u, "(A,I0)") "Number of particles: ", phs%n_particles
    write (u, "(A,ES12.5)") "Q**2: ", squared (phs%total_p)
    write (u, "(A,ES12.5)") "Volume: ", phs%volume
    write (u, "(A,ES12.5)") "Jacobian: ", phs%jacobian
    call phs%total_p%write (unit)
    total_p = phs%total_p
    do i = 1, phs%n_particles
       associate (prt => phs%prt(i))
         call phs%prt(i)%write (unit)
         total_p = total_p - prt%p ()
       end associate
    end do
    write (u, "(A,1X,74(A))") "TOTAL", "=========================================================================="
    call total_p%write (unit)
  end subroutine phs_base_write

  !> Get i-th particle momentum.
  !! @param i Index of i-th particle.
  !! @return p Momentum of i-th particle.
  function phs_p (phs, i) result(p)
    class(phs_t), intent(in) :: phs
    integer, intent(in) :: i
    real(default), dimension(4) :: p
    associate (prt => phs%prt(i))
     p = prt%get_x ()
    end associate
  end function phs_p

  !> Get phase-space weight.
  !! @return w Phase-space weight.
  real(default) function phs_get_weight (phs) result (w)
    class(phs_t), intent(in) :: phs
    ! w = phs%jacobian * phs%volume * phs%flux
    w = phs%jacobian * phs%volume
  end function phs_get_weight
end module phs
