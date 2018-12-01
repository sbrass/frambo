!> \file rambo.f90
!! \brief Flat phase-space with the RAMBO algorithm.
!! \details The RAMBO algorithm equally populates the phase-space volume  with each point having the constant weight
!! \f[
!!    \Phi_n(0) = \frac{1}{8^{n - 1}(2\pi)^{3n}} \frac{Q^{n - 2}}{(n - 1)!(n - 2)!}
!! \f]
!! for massless particles at total momentum \f$Q\f$.
!! The relation to the massive case is given by the reweighting factor:
!! \f[
!!    \Phi_n(m) = \prod_{i = 2}^{n} \frac{ρ(M_{i - 1}, M_i, m_{i - 1})}{ρ(K_{i - 1}, K_i, 0)} \prod_{i = 2}^{n - 1} \frac{M_i}{K_i} \Phi_n(0).
!! \f]

!> Module: phs_rambo
!> \author Simon Braß
!> \date 07.11.2018
module phs_rambo
  use ISO_FORTRAN_ENV, only: OUTPUT_UNIT
  use constants
  use numerics
  use vector
  use particle
  use phs

  implicit none

  !> \brief Flat phase-space type RAMBO
  type, extends(phs_t) :: phs_rambo_t
     private
     real(default), dimension(:), allocatable :: K, M
   contains
     procedure :: write => phs_rambo_write
     procedure :: generate => phs_rambo_generate
     procedure :: invert => phs_rambo_invert
     procedure, private :: generate_intermediates => phs_rambo_generate_intermediates
     procedure, private :: invert_intermediates => phs_rambo_invert_intermediates
  end type phs_rambo_t

  !> \brief Overload type constructor with custom phs_rambo_init.
  interface phs_rambo_t
     module procedure :: phs_rambo_init
  end interface phs_rambo_t
contains
  !> \brief Initialize phs_rambo and calculate the energy-independent phase-space volume.
  type(phs_rambo_t) function phs_rambo_init (n_particles, total_momentum, masses) result (phs)
    integer, intent(in) :: n_particles
    real(default), dimension(4) :: total_momentum
    real(default), dimension(:), intent(in) :: masses
    call phs%init (n_particles, total_momentum, masses)
    allocate (phs%K(phs%n_particles), source=0._default)
    allocate (phs%M(phs%n_particles), source=0._default)
    phs%volume = 1. / (2. * PI)**(3. * n_particles) &
         * (PI / 2.)**(n_particles - 1) &
         / (faculty (n_particles - 1) * faculty (n_particles - 2))
  contains
  end function phs_rambo_init

  !> \brief Generate flat phase-space using \f$3n - 4\f$ random numbers.
  !! \param r_in \f$3n - 4\f$ random numbers in \f$(0, 1)\f$.
  subroutine phs_rambo_generate (phs, r_in)
    class(phs_rambo_t), intent(inout) :: phs
    real(default), dimension(3 * phs%n_particles - 4), intent(in) :: r_in
    real(default), dimension(2) :: r
    type(vector4_t), dimension(2) :: p
    type(vector4_t) :: Q, QNext
    integer :: i
    associate (n => phs%n_particles, M => phs%M)
      call phs%generate_intermediates (r_in(1:n - 2))
      Q = phs%total_p
      do i = 2, n
         associate (prt => phs%prt(i - 1))
           r(1) = r_in (n - 5 + 2 * i)
           r(2) = r_in (n - 4 + 2 * i)
           call decay_two_body_system ( &
                r, M(i - 1), M(i), prt%mass (), p)
           QNext = p(1)
           call prt%set_p (p(2))
           call prt%boost (Q)
           call QNext%boost (Q)
           Q = QNext
         end associate
      end do
      call phs%prt(n)%set_p (Q)
    end associate
  contains
    !> \brief Decay mass \f$M\f$ into two-body systems with masses \f$m_1,
    !! m_2\f$ and angles given by \f$r\f$.
    !! \details The magic happens mostly here...
    !! \params r Two random numbers in \f$(0, 1)\f$. First: \f$\theta\f$.
    !! Second: \f$\varphi\f$.
    subroutine decay_two_body_system (r, m, m1, m2, p)
      real(default), dimension(2), intent(in) :: r
      real(default), intent(in) :: m, m1, m2
      type(vector4_t), dimension(2), intent(out) :: p
      real(default) :: p_space_abs
      real(default), dimension(3) :: p_space
      real(default) :: phi, cos_theta
      cos_theta = 2. * r(1) - 1.
      phi = 2. * PI * r(2)
      if (phi > PI) phi = phi - 2. * PI
      p_space_abs = 4. * m * rho (m, m1, m2)
      p_space = [cos(phi) * sqrt(1. - cos_theta**2), &
                 sin(phi) * sqrt(1. - cos_theta**2), &
                 cos_theta]
      p_space = p_space_abs * p_space
      p(1) = vector4_t ([sqrt (p_space_abs**2 + m1**2), p_space])
      p(2) = vector4_t ([sqrt (p_space_abs**2 + m2**2), -p_space])
    end subroutine decay_two_body_system
  end subroutine phs_rambo_generate

  !> \brief Generate intermediate mass systems.
  !! \details The intermediate mass systems are defined recursively by...
  subroutine phs_rambo_generate_intermediates (phs, r)
    class(phs_rambo_t), intent(inout) :: phs
    real(default), dimension(phs%n_particles - 2), intent(in) :: r
    integer :: i, j
    associate (n => phs%n_particles, K => phs%K, M => phs%M)
      M(1) = mass (phs%total_p)
      M(n) = phs%prt(n)%mass ()
      call calculate_k (r)
      do i = 1, n - 1
         M(i) = K(i)
         do j = i, n
            M(i) = M(i) + phs%prt(j)%mass ()
         end do
      end do
      phs%jacobian = K(1)**(2 * n - 4) &
           * 8. * rho(M(n - 1), phs%prt(n)%mass (), phs%prt(n - 1)%mass ())
      do i = 2, n - 1
         phs%jacobian = phs%jacobian * &
              rho(M(i - 1), M(i), phs%prt(i - 1)%mass ()) / &
              rho(K(i - 1), K(i), 0._default) * &
              M(i) / K(i)
      end do
    end associate
  contains
    subroutine calculate_k (r)
      real(default), dimension(phs%n_particles - 2), intent(in) :: r
      real(default), dimension(:), allocatable :: u
      integer :: i
      associate (n => phs%n_particles, K => phs%K, M => phs%M)
        K = 0; K(1) = M(1)
        do i = 1, n
           K(1) = K(1) - phs%prt(i)%mass()
        end do
        allocate (u(2:n - 1), source=0._default)
        call solve_for_u (r, u)
        do i = 2, n - 1
           K(i) = sqrt (u(i) * K(i - 1)**2)
        end do
      end associate
    end subroutine calculate_k

    !> Solve for integration variable \f$u_i\f$.
    !!
    !! Search for one roots in \f$[0, 1]\f$ for the function
    !! \f[
    !!    r_i = (n + 1)u_i^n - n u_i^{n + 1}.
    !! \f]
    !! Bisection method to find a root in a given interval.
    !! The solution is hard-coded for the interval \f$[0,1]\f$.
    subroutine solve_for_u (r, u)
      real(default), dimension(phs%n_particles - 2), intent(in) :: r
      real(default), dimension(2:phs%n_particles - 1), intent(out) :: u
      integer :: i, j
      real(default) :: f, f_mid, xl, xr, xmid
      associate (n => phs%n_particles)
        do i = 2, n - 1
           xl = 0
           xr = 1
           if (r(i - 1) == 1 .or. r(i - 1) == 0) then
              u(i) = r(i - 1)
           else
              do j = 1, 100
                 xmid = (xl + xr) / 2.
                 f = f_rambo (xl, n - i) - r(i - 1)
                 f_mid = f_rambo (xmid, n - i) - r(i - 1)
                 if (f * f_mid > 0) then
                    xl = xmid
                 else
                    xr = xmid
                 end if
                 if (abs(xl - xr) < 1e-15) exit
              end do
              u(i) = xmid
           end if
        end do
      end associate
    end subroutine solve_for_u

    real(default) function f_rambo(u, n)
      real(default), intent(in) :: u
      integer, intent(in) :: n
      f_rambo = (n + 1) * u**n - n * u**(n + 1)
    end function f_rambo
  end subroutine phs_rambo_generate_intermediates

  !> \brief Invert momenta back to random numbers and compute phase-space weight (cross-check).
  subroutine phs_rambo_invert (phs, r)
    class(phs_rambo_t), intent(inout) :: phs
    real(default), dimension(3 * phs%n_particles - 4), intent(out) :: r
    type(vector4_t), dimension(:), allocatable :: Q
    type(vector4_t) :: p
    real(default) :: cos_theta, phi
    integer :: i, j
    associate (n => phs%n_particles, M => phs%M)
      allocate (Q(n), source=zero_vector)
      M(1) = mass (phs%total_p)
      Q(1) = vector4_t (M(1))
      Q(n) = phs%prt(n)%p ()
      do i = 2, n - 1
         do j = i, n
            Q(i) = Q(i) + phs%prt(j)%p ()
         end do
         M(i) = mass (Q(i))
      end do
      call phs%invert_intermediates (r(1:n - 2))
      do i = 2, n
         p = phs%prt(i - 1)%p ()
         call p%boost (Q(i - 1), invert = .true.)
         cos_theta = p%cos_theta ()
         phi = p%phi (); if (phi < 0.) phi = phi + 2. * PI
         r(n - 5 + 2 * i) = (cos_theta + 1.) / 2.
         r(n - 4 + 2 * i) = phi / (2. * PI)
      end do
    end associate
  end subroutine phs_rambo_invert

  subroutine phs_rambo_invert_intermediates (phs, r)
    class(phs_rambo_t), intent(inout) :: phs
    real(default), dimension(phs%n_particles - 2), intent(out) :: r
    integer :: i, j
    associate (n => phs%n_particles, K => phs%K, M => phs%M)
      K = M
      do i = 1, n - 1
         do j = i, n
            K(i) = K(i) - phs%prt(j)%mass ()
         end do
      end do
      call solve_for_r (r)
      phs%jacobian = K(1)**(2 * n - 4) &
           * 8. * rho(M(n - 1), phs%prt(n)%mass (), phs%prt(n - 1)%mass ())
      do i = 2, n - 1
         phs%jacobian = phs%jacobian * &
              rho(M(i - 1), M(i), phs%prt(i - 1)%mass ()) / &
              rho(K(i - 1), K(i), 0._default) * &
              M(i) / K(i)
      end do
    end associate
  contains
    subroutine solve_for_r (r)
      real(default), dimension(phs%n_particles - 2), intent(out) :: r
      integer :: i
      real(default) :: u
      associate (n => phs%n_particles)
        do i = 2, n - 1
           u = (phs%K(i) / phs%K(i - 1))**2
           r(i - 1) = (n + 1 - i) * u**(n - i) &
                - (n - i) * u**(n + 1 - i)
        end do
      end associate
    end subroutine solve_for_r
  end subroutine phs_rambo_invert_intermediates

  !> Write phs_rambo_t to unit.
  subroutine phs_rambo_write (phs, unit)
    class(phs_rambo_t), intent(in) :: phs
    integer, intent(in), optional :: unit
    integer :: u, i
    u = OUTPUT_UNIT; if (present (unit)) u = unit
    write (u, "(80(A))") "================================================================================"
    write (u, "(A)") "Phasespace RAMBO type:"
    write (u, "(11X,A,11X,A,1X,A)") "K", "M"
    do i = 1, phs%n_particles
       write (u, "(2(ES12.5))") phs%K(i), phs%M(i)
    end do
    call phs%base_write (unit)
  end subroutine phs_rambo_write
end module phs_rambo
