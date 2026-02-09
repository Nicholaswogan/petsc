module pseudo_transient
  use iso_fortran_env, only: wp => real64
  implicit none
  private

  integer, parameter, public :: PTC_JAC_DENSE = 1
  integer, parameter, public :: PTC_JAC_BAND  = 2

  integer, parameter, public :: PTC_REASON_NONE                  = 0
  integer, parameter, public :: PTC_CONVERGED_PSEUDO_FATOL       = 1
  integer, parameter, public :: PTC_CONVERGED_PSEUDO_FRTOL       = 2
  integer, parameter, public :: PTC_DIVERGED_STEP_REJECTED       = -1
  integer, parameter, public :: PTC_DIVERGED_CALLBACK_FATAL      = -2
  integer, parameter, public :: PTC_DIVERGED_NOT_INITIALIZED     = -3
  integer, parameter, public :: PTC_DIVERGED_INVALID_INPUT       = -4
  integer, parameter, public :: PTC_DIVERGED_MAX_STEPS           = -5

  public :: PTCSolver, wp

  abstract interface
    subroutine rhs_fcn(u, udot, ierr)
      import :: wp
      implicit none
      real(wp), intent(in) :: u(:)
      real(wp), intent(out) :: udot(:)
      integer, intent(out) :: ierr
    end subroutine rhs_fcn

    subroutine jac_fcn(u, jac, ierr)
      import :: wp
      implicit none
      real(wp), intent(in) :: u(:)
      real(wp), intent(out) :: jac(:, :)
      integer, intent(out) :: ierr
    end subroutine jac_fcn

    subroutine verify_step_fcn(x_new, dt, accept, ierr)
      import :: wp
      implicit none
      real(wp), intent(in) :: x_new(:)
      real(wp), intent(inout) :: dt
      logical, intent(out) :: accept
      integer, intent(out) :: ierr
    end subroutine verify_step_fcn

    subroutine timestep_fcn(fnorm, fnorm_initial, fnorm_previous, dt, dt_initial, dt_increment, increment_from_initial_dt, dt_max, new_dt, ierr)
      import :: wp
      implicit none
      real(wp), intent(in) :: fnorm, fnorm_initial, fnorm_previous
      real(wp), intent(in) :: dt, dt_initial, dt_increment, dt_max
      logical, intent(in) :: increment_from_initial_dt
      real(wp), intent(out) :: new_dt
      integer, intent(out) :: ierr
    end subroutine timestep_fcn
  end interface

  interface
    subroutine dgesv(n, nrhs, a, lda, ipiv, b, ldb, info)
      import :: wp
      implicit none
      integer, intent(in) :: n, nrhs, lda, ldb
      integer, intent(out) :: ipiv(*)
      real(wp), intent(inout) :: a(lda, *)
      real(wp), intent(inout) :: b(ldb, *)
      integer, intent(out) :: info
    end subroutine dgesv

    subroutine dgbsv(n, kl, ku, nrhs, ab, ldab, ipiv, b, ldb, info)
      import :: wp
      implicit none
      integer, intent(in) :: n, kl, ku, nrhs, ldab, ldb
      integer, intent(out) :: ipiv(*)
      real(wp), intent(inout) :: ab(ldab, *)
      real(wp), intent(inout) :: b(ldb, *)
      integer, intent(out) :: info
    end subroutine dgbsv
  end interface

  type :: PTCSolver
    integer :: neq = 0
    integer :: jacobian_type = 0
    integer :: kl = 0
    integer :: ku = 0
    integer :: ldab = 0

    procedure(rhs_fcn), pointer, nopass :: f => null()
    procedure(jac_fcn), pointer, nopass :: jac => null()
    procedure(verify_step_fcn), pointer, nopass :: verify => null()
    procedure(timestep_fcn), pointer, nopass :: compute_dt => null()

    real(wp) :: dt = 0.0_wp
    real(wp) :: dt_initial = 0.0_wp
    real(wp) :: dt_increment = 1.1_wp
    real(wp) :: dt_max = 0.0_wp
    logical :: increment_dt_from_initial_dt = .false.

    real(wp) :: fatol = 1.0e-50_wp
    real(wp) :: frtol = 1.0e-12_wp

    real(wp) :: fnorm = -1.0_wp
    real(wp) :: fnorm_initial = -1.0_wp
    real(wp) :: fnorm_previous = -1.0_wp

    integer :: steps = 0
    integer :: rejects_total = 0
    integer :: max_reject = 10
    integer :: max_steps = 10000
    integer :: reason = PTC_REASON_NONE

    logical :: initialized = .false.

    real(wp), allocatable :: x(:)
    real(wp), allocatable :: x_old(:)
    real(wp), allocatable :: fvec(:)
    real(wp), allocatable :: step_vec(:)
    real(wp), allocatable :: rhs_mat(:, :)

    real(wp), allocatable :: jac_mat(:, :)
    real(wp), allocatable :: a_dense(:, :)

    ! Banded Jacobian compact storage (LAPACK standard):
    ! jac_mat(ku+1+i-j, j) = J(i,j), for max(1,j-ku) <= i <= min(n,j+kl)
    real(wp), allocatable :: a_band(:, :)

    integer, allocatable :: ipiv(:)
  contains
    procedure :: initialize => PTCSolver_initialize
    procedure :: step => PTCSolver_step
    procedure :: solve => PTCSolver_solve
    procedure :: check_convergence => PTCSolver_check_convergence

    procedure :: set_verify_timestep => PTCSolver_set_verify_timestep
    procedure :: set_compute_timestep => PTCSolver_set_compute_timestep
  end type PTCSolver

contains

  subroutine PTCSolver_initialize(self, x0, f, jacobian_type, dt0, jac, kl, ku, fatol, frtol, dt_increment, dt_max, increment_dt_from_initial_dt, max_reject, max_steps)
    class(PTCSolver), intent(inout) :: self
    real(wp), intent(in) :: x0(:)
    procedure(rhs_fcn) :: f
    integer, intent(in) :: jacobian_type
    real(wp), intent(in) :: dt0
    procedure(jac_fcn) :: jac
    integer, intent(in), optional :: kl, ku, max_reject, max_steps
    real(wp), intent(in), optional :: fatol, frtol, dt_increment, dt_max
    logical, intent(in), optional :: increment_dt_from_initial_dt

    call reset_storage(self)

    self%reason = PTC_REASON_NONE
    if (size(x0) <= 0 .or. dt0 <= 0.0_wp) then
      self%reason = PTC_DIVERGED_INVALID_INPUT
      return
    end if

    self%neq = size(x0)
    self%dt = dt0
    self%dt_initial = dt0
    self%steps = 0
    self%rejects_total = 0

    self%f => f

    if (present(dt_increment)) self%dt_increment = dt_increment
    if (present(dt_max)) self%dt_max = dt_max
    if (present(fatol)) self%fatol = fatol
    if (present(frtol)) self%frtol = frtol
    if (present(increment_dt_from_initial_dt)) self%increment_dt_from_initial_dt = increment_dt_from_initial_dt
    if (present(max_reject)) self%max_reject = max_reject
    if (present(max_steps)) self%max_steps = max_steps

    self%fnorm = -1.0_wp
    self%fnorm_initial = -1.0_wp
    self%fnorm_previous = -1.0_wp

    allocate(self%x(self%neq), self%x_old(self%neq), self%fvec(self%neq), self%step_vec(self%neq), self%rhs_mat(self%neq, 1), self%ipiv(self%neq))
    self%x = x0

    self%jacobian_type = jacobian_type
    self%jac => jac

    select case (self%jacobian_type)
    case (PTC_JAC_DENSE)
      allocate(self%jac_mat(self%neq, self%neq), self%a_dense(self%neq, self%neq))

    case (PTC_JAC_BAND)
      if (.not. present(kl) .or. .not. present(ku)) then
        self%reason = PTC_DIVERGED_INVALID_INPUT
        return
      end if
      if (kl < 0 .or. ku < 0) then
        self%reason = PTC_DIVERGED_INVALID_INPUT
        return
      end if
      self%kl = kl
      self%ku = ku
      self%ldab = 2 * self%kl + self%ku + 1
      allocate(self%jac_mat(self%kl + self%ku + 1, self%neq), self%a_band(self%ldab, self%neq))

    case default
      self%reason = PTC_DIVERGED_INVALID_INPUT
      return
    end select

    self%initialized = .true.
  end subroutine PTCSolver_initialize

  subroutine PTCSolver_set_verify_timestep(self, verify)
    class(PTCSolver), intent(inout) :: self
    procedure(verify_step_fcn) :: verify

    self%verify => verify
  end subroutine PTCSolver_set_verify_timestep

  subroutine PTCSolver_set_compute_timestep(self, compute_dt)
    class(PTCSolver), intent(inout) :: self
    procedure(timestep_fcn) :: compute_dt

    self%compute_dt => compute_dt
  end subroutine PTCSolver_set_compute_timestep

  subroutine PTCSolver_step(self)
    class(PTCSolver), intent(inout) :: self

    integer :: ierr, rejections
    real(wp) :: next_dt, reject_dt
    logical :: accept

    if (.not. self%initialized) then
      self%reason = PTC_DIVERGED_NOT_INITIALIZED
      return
    end if

    if (self%reason /= PTC_REASON_NONE) return

    if (self%steps == 0) self%dt_initial = self%dt

    rejections = 0

    do
      self%x_old = self%x
      call PTCSolver_take_newton_update(self, ierr)
      if (ierr < 0) then
        self%reason = PTC_DIVERGED_CALLBACK_FATAL
        self%x = self%x_old
        return
      end if

      if (ierr > 0) then
        reject_dt = max(0.5_wp * self%dt, tiny(1.0_wp))
        call PTCSolver_reject_step(self, rejections, reject_dt)
        if (self%reason /= PTC_REASON_NONE) return
        cycle
      end if

      accept = .true.
      reject_dt = self%dt
      if (associated(self%verify)) then
        call self%verify(self%x, reject_dt, accept, ierr)
        if (ierr < 0) then
          self%reason = PTC_DIVERGED_CALLBACK_FATAL
          self%x = self%x_old
          return
        end if
        if (ierr > 0) accept = .false.
      end if

      if (.not. accept) then
        call PTCSolver_reject_step(self, rejections, max(reject_dt, tiny(1.0_wp)))
        if (self%reason /= PTC_REASON_NONE) return
        cycle
      end if

      call PTCSolver_compute_residual(self, self%x, self%fvec, self%fnorm, ierr)
      if (ierr < 0) then
        self%reason = PTC_DIVERGED_CALLBACK_FATAL
        self%x = self%x_old
        return
      end if
      if (ierr > 0) then
        reject_dt = max(0.5_wp * self%dt, tiny(1.0_wp))
        call PTCSolver_reject_step(self, rejections, reject_dt)
        if (self%reason /= PTC_REASON_NONE) return
        cycle
      end if

      if (self%fnorm_initial < 0.0_wp) then
        self%fnorm_initial = self%fnorm
        self%fnorm_previous = self%fnorm
      end if

      call PTCSolver_compute_next_dt(self, next_dt, ierr)
      if (ierr < 0) then
        self%reason = PTC_DIVERGED_CALLBACK_FATAL
        self%x = self%x_old
        return
      end if
      if (ierr > 0) then
        reject_dt = max(0.5_wp * self%dt, tiny(1.0_wp))
        call PTCSolver_reject_step(self, rejections, reject_dt)
        if (self%reason /= PTC_REASON_NONE) return
        cycle
      end if

      self%dt = next_dt
      self%fnorm_previous = self%fnorm
      self%steps = self%steps + 1

      call PTCSolver_check_convergence(self)
      return
    end do
  end subroutine PTCSolver_step

  subroutine PTCSolver_solve(self, max_steps)
    class(PTCSolver), intent(inout) :: self
    integer, intent(in), optional :: max_steps

    integer :: step_limit

    if (.not. self%initialized) then
      self%reason = PTC_DIVERGED_NOT_INITIALIZED
      return
    end if

    step_limit = self%max_steps
    if (present(max_steps)) step_limit = max_steps

    do while (self%reason == PTC_REASON_NONE)
      if (self%steps >= step_limit) then
        self%reason = PTC_DIVERGED_MAX_STEPS
        exit
      end if
      call self%step()
    end do
  end subroutine PTCSolver_solve

  subroutine PTCSolver_check_convergence(self)
    class(PTCSolver), intent(inout) :: self

    if (self%fnorm < self%fatol) then
      self%reason = PTC_CONVERGED_PSEUDO_FATOL
      return
    end if

    if (self%fnorm_initial > 0.0_wp) then
      if ((self%fnorm / self%fnorm_initial) < self%frtol) then
        self%reason = PTC_CONVERGED_PSEUDO_FRTOL
        return
      end if
    end if

    self%reason = PTC_REASON_NONE
  end subroutine PTCSolver_check_convergence

  subroutine PTCSolver_compute_residual(self, x, fvec, fnorm, ierr)
    class(PTCSolver), intent(inout) :: self
    real(wp), intent(in) :: x(:)
    real(wp), intent(out) :: fvec(:)
    real(wp), intent(out) :: fnorm
    integer, intent(out) :: ierr

    call self%f(x, fvec, ierr)
    if (ierr /= 0) then
      fnorm = -1.0_wp
      return
    end if

    fnorm = norm2(fvec)
  end subroutine PTCSolver_compute_residual

  subroutine PTCSolver_take_newton_update(self, ierr)
    class(PTCSolver), intent(inout) :: self
    integer, intent(out) :: ierr

    integer :: info, i, j, row_compact, row_solve
    real(wp) :: inv_dt

    ierr = 0
    if (self%dt <= 0.0_wp) then
      ierr = 1
      return
    end if

    call PTCSolver_compute_residual(self, self%x, self%fvec, self%fnorm, ierr)
    if (ierr /= 0) return

    inv_dt = 1.0_wp / self%dt

    select case (self%jacobian_type)
    case (PTC_JAC_DENSE)
      call self%jac(self%x, self%jac_mat, ierr)
      if (ierr /= 0) return

      self%a_dense = -self%jac_mat
      do i = 1, self%neq
        self%a_dense(i, i) = self%a_dense(i, i) + inv_dt
      end do

      self%rhs_mat(:, 1) = self%fvec
      call dgesv(self%neq, 1, self%a_dense, self%neq, self%ipiv, self%rhs_mat, self%neq, info)
      if (info /= 0) then
        ierr = 1
        return
      end if

      self%step_vec = self%rhs_mat(:, 1)
      self%x = self%x + self%step_vec

    case (PTC_JAC_BAND)
      call self%jac(self%x, self%jac_mat, ierr)
      if (ierr /= 0) return

      self%a_band = 0.0_wp
      do j = 1, self%neq
        do i = max(1, j - self%ku), min(self%neq, j + self%kl)
          row_compact = self%ku + 1 + i - j
          row_solve = self%kl + row_compact
          self%a_band(row_solve, j) = -self%jac_mat(row_compact, j)
        end do
        self%a_band(self%kl + self%ku + 1, j) = self%a_band(self%kl + self%ku + 1, j) + inv_dt
      end do

      self%rhs_mat(:, 1) = self%fvec
      call dgbsv(self%neq, self%kl, self%ku, 1, self%a_band, self%ldab, self%ipiv, self%rhs_mat, self%neq, info)
      if (info /= 0) then
        ierr = 1
        return
      end if

      self%step_vec = self%rhs_mat(:, 1)
      self%x = self%x + self%step_vec

    case default
      ierr = -1
    end select
  end subroutine PTCSolver_take_newton_update

  subroutine PTCSolver_compute_next_dt(self, next_dt, ierr)
    class(PTCSolver), intent(inout) :: self
    real(wp), intent(out) :: next_dt
    integer, intent(out) :: ierr

    ierr = 0

    if (associated(self%compute_dt)) then
      call self%compute_dt(self%fnorm, self%fnorm_initial, self%fnorm_previous, self%dt, self%dt_initial, self%dt_increment, self%increment_dt_from_initial_dt, self%dt_max, next_dt, ierr)
      if (ierr /= 0) return
    else
      if (self%fnorm == 0.0_wp) then
        next_dt = 1.0e12_wp * self%dt_increment * self%dt
      else if (self%increment_dt_from_initial_dt) then
        next_dt = self%dt_increment * self%dt_initial * self%fnorm_initial / self%fnorm
      else
        next_dt = self%dt_increment * self%dt * self%fnorm_previous / self%fnorm
      end if
      if (self%dt_max > 0.0_wp) next_dt = min(next_dt, self%dt_max)
    end if

    if (next_dt <= 0.0_wp) then
      ierr = 1
      return
    end if
  end subroutine PTCSolver_compute_next_dt

  subroutine PTCSolver_reject_step(self, rejections, new_dt)
    class(PTCSolver), intent(inout) :: self
    integer, intent(inout) :: rejections
    real(wp), intent(in) :: new_dt

    self%x = self%x_old
    self%dt = max(new_dt, tiny(1.0_wp))
    self%rejects_total = self%rejects_total + 1

    rejections = rejections + 1
    if (rejections > self%max_reject) then
      self%reason = PTC_DIVERGED_STEP_REJECTED
    end if
  end subroutine PTCSolver_reject_step

  subroutine reset_storage(self)
    class(PTCSolver), intent(inout) :: self

    if (allocated(self%x)) deallocate(self%x)
    if (allocated(self%x_old)) deallocate(self%x_old)
    if (allocated(self%fvec)) deallocate(self%fvec)
    if (allocated(self%step_vec)) deallocate(self%step_vec)
    if (allocated(self%rhs_mat)) deallocate(self%rhs_mat)
    if (allocated(self%jac_mat)) deallocate(self%jac_mat)
    if (allocated(self%a_dense)) deallocate(self%a_dense)
    if (allocated(self%a_band)) deallocate(self%a_band)
    if (allocated(self%ipiv)) deallocate(self%ipiv)

    self%initialized = .false.
    self%f => null()
    self%jac => null()
    self%verify => null()
    self%compute_dt => null()
  end subroutine reset_storage

end module pseudo_transient
