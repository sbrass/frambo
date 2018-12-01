!> \file vector.f90
!! \brief Definition of a particle and four-vectors.

!> Module: vector
!> \author Simon Braß
!> \date 9.11.2018
module vector
  use ISO_FORTRAN_ENV, only: OUTPUT_UNIT
  use constants, only: default

  implicit none

  private

  !> \brief Four-vector type.
  type :: vector4_t
     private
     real(default), dimension(4) :: x = &
          [0._default, 0._default, 0._default, 0._default]
   contains
     procedure :: write => vector4_write
     procedure :: get_x => vector4_get_x
     procedure :: time => vector4_time
     procedure :: space => vector4_space
     procedure :: cos_theta => vector4_cos_theta
     procedure :: phi => vector4_phi
     generic :: operator(+) => add
     generic :: operator(-) => subtract
     procedure, private :: add => vector4_add
     procedure, private :: subtract => vector4_subtract
     procedure :: set_time => vector4_set_time
     procedure :: set_space => vector4_set_space
     procedure :: boost => vector4_boost
  end type vector4_t
  !! \details We store the four vector as 4-dim. array.
  !! We supply additional functionality.

  type(vector4_t), parameter :: zero_vector = vector4_t ([0._default, 0._default, 0._default, 0._default])

  interface mass
     module procedure :: mass_vector4, mass_array4
  end interface mass

  interface squared
     module procedure :: squared_vector4, squared_array4
  end interface squared

  interface vector4_t
     module procedure :: vector4_init, vector4_init_at_rest
  end interface vector4_t

  public :: vector4_t, zero_vector, mass, squared, rho
contains
  real(default) function rho (M1, M2, m)
    real(default), intent(in) :: M1, M2, m
    real(default) :: MP, MM
    ! rho = sqrt ((M1**2 - (M2 + m)**2) * (M1**2 - (M2 - m)**2))
    MP = (M1 - (M2 + m)) * (M1 + (M2 + m))
    MM = (M1 - (M2 - m)) * (M1 + (M2 - m))
    rho = sqrt (MP * MM)
    rho = rho / (8._default * M1**2)
  end function rho

  !> \brief Compute lorentz-invariant \f$s^2\f$ of four-vector.
  !! @param v Four-vector.
  real(default) function squared_vector4 (v) result (v_sq)
    class(vector4_t), intent(in) :: v
    v_sq = v%time ()**2 - dot_product (v%space (), v%space ())
  end function squared_vector4

  !> \brief Compute lorentz-invariant \f$s^2\f$ of 4-dim. array.
  !! @param v 4-dim. array (representing a four-vector).
  real(default) function squared_array4 (v) result (v_sq)
    real(default), dimension(4), intent(in) :: v
    v_sq = v(1) * v(1) - dot_product (v(2:), v(2:))
  end function squared_array4

  !> \brief Compute mass of four-vector.
  !! @param
  real(default) function mass_vector4 (v) result (m)
    class(vector4_t), intent(in) :: v
    m = mass (v%x)
  end function mass_vector4

  !> Compute mass of 4-dim. array.
  !! @param v 4-dim. array.
  real(default) function mass_array4 (v) result (m)
    real(default), dimension(4), intent(in) :: v
    real(default) :: v_sq
    v_sq = squared (v)
    if (v_sq >= 0) then
       m = sqrt(v_sq)
    else
       m = sqrt (-v_sq)
    end if
  end function mass_array4

  !> Construct vector4_t with momentum.
  type(vector4_t) function vector4_init (x) result (v)
    real(default), dimension(4), intent(in) :: x
    v%x = x
  end function vector4_init

  !> Construct vector4_t at rest.
  !! \param x0 Time-component of four-vector.
  !! \return v Four-vector at rest.
  type(vector4_t) function vector4_init_at_rest (x0) result (v)
    real(default), intent(in) :: x0
    v%x = 0; v%x(1) = x0
  end function vector4_init_at_rest

  !> Write vector4_t to unit.
  !! @param unit Output unit, defaults to OUTPUT_UNIT.
  subroutine vector4_write (v, unit)
    class(vector4_t), intent(in) :: v
    integer, intent(in), optional :: unit
    integer :: u
    u = OUTPUT_UNIT; if (present (unit)) u = unit
    write (u, "(A)") "P (P^2)"
    write (u, "(4(1X,ES12.5),1X,A,ES12.5,A)") v%x, "(", squared (v) ,")"
  end subroutine vector4_write

  !> \brief Get 4-dim. array repressentation of four-vector.
  !! \return x 4-dim. array with x(1) time component and x(2:4) space component.
  function vector4_get_x (v) result (x)
    class(vector4_t), intent(in) :: v
    real(default), dimension(4) :: x
    ! Retrieve time and space component implementation independently.
    x(1) = v%time ()
    x(2:4) = v%space ()
  end function vector4_get_x

  !> Get time-component of four-vector.
  real(default) function vector4_time (v) result (time)
    class(vector4_t), intent(in) :: v
    time = v%x(1)
  end function vector4_time

  !> Get space-component of four-vector.
  function vector4_space (v) result (vector)
    class(vector4_t), intent(in) :: v
    real(default), dimension(3) :: vector
    vector = v%x(2:4)
  end function vector4_space

  !> Compute \f$\cos \theta = \frac{z}{|\vec{x}|}\f$.
  function vector4_cos_theta (v) result (ctheta)
    class(vector4_t), intent(in) :: v
    real(default) :: ctheta
    ctheta = v%x(4) / sqrt (dot_product (v%space (), v%space ()))
  end function vector4_cos_theta

  !> Compute \f$\varphi = \, \mathrm{atan} \, (\frac{z}{y})\f$.
  function vector4_phi (v) result (phi)
    class(vector4_t), intent(in) :: v
    real(default) :: phi
    phi = atan2 (v%x(3), v%x(2))
  end function vector4_phi

  !> Add v1 and v2.
  type(vector4_t) function vector4_add (v1, v2) result (v)
    class(vector4_t), intent(in) :: v1
    type(vector4_t), intent(in) :: v2
    v%x = v1%x + v2%x
  end function vector4_add

  !> Subtract v2 from v1.
  type(vector4_t) function vector4_subtract (v1, v2) result (v)
    class(vector4_t), intent(in) :: v1
    type(vector4_t), intent(in) :: v2
    v%x = v1%x - v2%x
  end function vector4_subtract

  !> Set time-component of four-vector.
  !! @param time New time component.
  subroutine vector4_set_time (v, time)
    class(vector4_t), intent(inout) :: v
    real(default), intent(in) :: time
    v%x(1) = time
  end subroutine vector4_set_time

  !> Set space-component of four-vector.
  !! @param space New 3-space component.
  subroutine vector4_set_space (v, space)
    class(vector4_t), intent(inout) :: v
    real(default), dimension(3) :: space
    v%x(2:4) = space
  end subroutine vector4_set_space

  !> Boost four-vector 4 in direction of boost (or opposite if invert).
  !! \details Add Formula here.
  !! @param boost Four-vector defining direction and gamma of boost.
  !! @param invert Invert direction of boost vector.
  subroutine vector4_boost (v, boost, invert)
    class(vector4_t), intent(inout) :: v
    type(vector4_t), intent(in) :: boost
    logical, intent(in), optional :: invert
    real(default) :: g, g_sq, b_sq, bv
    real(default), dimension(3) :: b
    if (all (boost%space () == 0)) return
    b = boost%space () / boost%time ()
    if (present(invert)) then
       if (invert) b = - b
    end if
    b_sq = dot_product (b, b)
    g = 1. / sqrt (1. - b_sq)
    g_sq = (g - 1.) / b_sq
    bv = dot_product (b, v%space ())
    call v%set_space (v%space () + g_sq * bv * b + g * b * v%time ())
    call v%set_time (g * v%time () + g * bv)
  end subroutine vector4_boost
end module vector
