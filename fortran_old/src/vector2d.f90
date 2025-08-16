module vector2d
  !use iso_fortran_env, only: real32,real64
  use nrtype
  implicit none
  !integer, parameter :: wp = real64 ! working precision
  
  ! 2D-vector type(class) declaration
  type, public :: vec2d
    private
    real(wp) :: vx=0._wp, vy=0._wp
    real(wp) :: mag=0._wp, phi=0._wp
  contains
    private
    procedure, public :: setcrt => set_vector_cartesian
    procedure, public :: setpol => set_vector_polar  
    procedure, public :: m => get_mag
    procedure, public :: m2 => get_magsq
    procedure, public :: f => get_phi
    procedure, public :: x => get_vx
    procedure, public :: y => get_vy
    generic, public :: assignment(=) => copy
    procedure, pass(this) :: copy
    generic, public :: operator(+) => add
    procedure, pass(this) :: add
    generic, public :: operator(-) => minus
    procedure, pass(this) :: minus
    generic, public :: operator(*) => scalar_L,scalar_R,dot
    procedure, pass(this) :: scalar_L,scalar_R,dot
    final :: destructor
  end type vec2d

  interface vec2d
    module procedure :: constructor
  end interface vec2d
  
  interface abs
    module procedure :: magnitude
  end interface abs

contains

  ! default constructors and destructor
  pure function constructor() result(this)
  type(vec2d) :: this
  call this%setpol(0._wp,0._wp)
  end function constructor

  pure subroutine destructor(this)
  type(vec2d), intent(inout) :: this
  end subroutine destructor

  ! setters
  pure subroutine set_vector_cartesian(this,x,y)
  class(vec2d), intent(inout) :: this
  real(wp), intent(in) :: x,y
  this%vx=x
  this%vy=y
  this%mag=hypot(x,y)
  this%phi=merge(0._wp,atan2(y,x),x.eq.0._wp.and.y.eq.0._wp)
  end subroutine set_vector_cartesian

  pure subroutine set_vector_polar(this,m,f)
  class(vec2d), intent(inout) :: this
  real(wp), intent(in) :: m,f
  this%mag=m
  this%phi=f
  this%vx=m*cos(f)
  this%vy=m*sin(f)
  end subroutine set_vector_polar

  ! getters
  pure function get_mag(this) result(res)
  class(vec2d), intent(in) :: this
  real(wp) :: res
  res=this%mag
  end function get_mag

  pure function get_magsq(this) result(res)
  class(vec2d), intent(in) :: this
  real(wp) :: res
  res=this%mag*this%mag
  end function get_magsq

  pure function get_phi(this) result(res)
  class(vec2d), intent(in) :: this
  real(wp) :: res
  res=this%phi
  end function get_phi

  pure function get_vx(this) result(res)
  class(vec2d), intent(in) :: this
  real(wp) :: res
  res=this%vx
  end function get_vx

  pure function get_vy(this) result(res)
  class(vec2d), intent(in) :: this
  real(wp) :: res
  res=this%vy
  end function get_vy

  ! intrinsic operator/function overloading
  pure subroutine copy(this,from)
  class(vec2d), intent(inout) :: this
  class(vec2d), intent(in) :: from
  call this%setcrt(from%x(),from%y())
  end subroutine copy

  pure function add(this,that) result(res)
  class(vec2d), intent(in) :: this,that
  type(vec2d) :: res
  call res%setcrt(this%x()+that%x(),this%y()+that%y())
  end function add

  pure function minus(this,that) result(res)
  class(vec2d), intent(in) :: this,that
  type(vec2d) :: res
  call res%setcrt(this%x()-that%x(),this%y()-that%y())
  end function minus

  pure function scalar_L(const,this) result(that)
  real(wp), intent(in) :: const
  class(vec2d), intent(in) :: this
  type(vec2d) :: that
  call that%setpol(const*this%m(), this%f())
  end function scalar_L

  pure function scalar_R(this,const) result(that)
  real(wp), intent(in) :: const
  class(vec2d), intent(in) :: this
  type(vec2d) :: that
  call that%setpol(const*this%m(), this%f())
  end function scalar_R

  pure function dot(this,that) result(res)
  class(vec2d), intent(in) :: this,that
  real(wp) :: res
  res = this%x()*that%x() + this%y()*that%y()
  end function dot

  pure function magnitude(this) result(res)
  class(vec2d), intent(in) :: this
  real(wp) :: res
  res = this%m()
  end function magnitude

  ! other functions

end module vector2d
