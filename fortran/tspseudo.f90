module pseudo_transient
  use iso_fortran_env, only: wp => real64
  implicit none
  private

  public :: PTCSolver, wp

  type :: PTCSolver
    integer :: neq !! number of ODEs
    procedure(lsoda_rhs_fcn), pointer :: f => NULL() !! right-hand-side of ODEs
    integer :: jt !! Jacobian type indicator
    procedure(lsoda_jac_fcn), pointer :: jac => NULL() !! jacobian of ODEs

    ...
  contains
    procedure :: initialize => PTCSolver_initialize
    procedure :: step => PTCSolver_step
    procedure :: solve => PTCSolver_solve
    procedure :: check_convergence => PTCSolver_check_convergence
  end type

  abstract interface

    function rhs_fcn(self, u, udot) result(ierr)
      import :: wp, PTCSolver
      implicit none
      class(PTCSolver), intent(inout) :: self
      real(dp), target, intent(in) :: u(:) !! State vector
      real(dp), target, intent(out) :: udot(:) !! Derivative vector
      integer :: ierr !! Set to 0 if successful.
                      !! Set to > 0 if there was a recoverable error. Code should re-try step with smaller timestep.
                      !! Set to < 0 to terminate the integration.
    end function

    function jac_fcn(self, u, jac) result(ierr)
      import :: wp, PTCSolver
      implicit none
      class(PTCSolver), intent(inout) :: self
      real(dp), target, intent(in) :: u(:) !! State vector
      real(dp), target, intent(out) :: jac(:,:)
      !! The Jacobian of the system. Can be dense or banded.
      integer :: ierr !! Set to 0 if successful.
                      !! Set to > 0 if there was a recoverable error. Code should re-try step with smaller timestep.
                      !! Set to < 0 to terminate the integration.
    end function

  end interface

  interface
    ! Interface to LAPACK routines here.
  end interface

contains

  subroutine PTCSolver_initialize(self, f, jac, ...)
    class(PTCSolver), intent(inout) :: self
    procedure(rhs_fcn), intent(in) :: f
    procedure(jac_fcn), intent(in) :: jac

  end subroutine

  subroutine PTCSolver_step(self, ...)
    class(PTCSolver), intent(inout) :: self
    ! Takes a single internal step
  end subroutine

  subroutine PTCSolver_solve(self, ...)
    class(PTCSolver), intent(inout) :: self
    ! Steps until convergence is reached
  end subroutine

  subroutine PTCSolver_check_convergence(self, ...)
    class(PTCSolver), intent(inout) :: self
    ! Checks for convergence
  end subroutine

end module