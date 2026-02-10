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
      real(wp), intent(in) :: u(:)  !! Input state vector.
      real(wp), intent(out) :: udot(:)  !! Residual/right-hand-side vector `f(u)`.
      integer, intent(out) :: ierr  !! Callback status (`0` success, nonzero failure).
    end subroutine rhs_fcn

    subroutine jac_fcn(u, jac, ierr)
      import :: wp
      implicit none
      real(wp), intent(in) :: u(:)  !! Input state vector.
      real(wp), intent(out) :: jac(:, :)  !! Jacobian in configured dense or compact-banded layout.
      integer, intent(out) :: ierr  !! Callback status (`0` success, nonzero failure).
    end subroutine jac_fcn

    subroutine verify_step_fcn(x_new, dt, accept, ierr)
      import :: wp
      implicit none
      real(wp), intent(in) :: x_new(:)  !! Candidate updated state.
      real(wp), intent(inout) :: dt  !! Proposed timestep (may be modified by callback).
      logical, intent(out) :: accept  !! Set true to accept candidate step.
      integer, intent(out) :: ierr  !! Callback status (`0` success, nonzero failure).
    end subroutine verify_step_fcn

    subroutine timestep_fcn(fnorm, fnorm_initial, fnorm_previous, dt, dt_initial, dt_increment, increment_from_initial_dt, dt_max, new_dt, ierr)
      import :: wp
      implicit none
      real(wp), intent(in) :: fnorm  !! Current residual norm.
      real(wp), intent(in) :: fnorm_initial  !! Residual norm from first accepted step.
      real(wp), intent(in) :: fnorm_previous  !! Residual norm from previous accepted step.
      real(wp), intent(in) :: dt  !! Current pseudo-time step.
      real(wp), intent(in) :: dt_initial  !! Initial pseudo-time step.
      real(wp), intent(in) :: dt_increment  !! Timestep growth factor.
      logical, intent(in) :: increment_from_initial_dt  !! Select initial-reference vs previous-step adaptation formula.
      real(wp), intent(in) :: dt_max  !! Maximum allowed timestep (non-positive means no cap).
      real(wp), intent(out) :: new_dt  !! Computed timestep for the next accepted step.
      integer, intent(out) :: ierr  !! Callback status (`0` success, nonzero failure).
    end subroutine timestep_fcn
  end interface

  interface
    subroutine dgesv(n, nrhs, a, lda, ipiv, b, ldb, info)
      import :: wp
      implicit none
      integer, intent(in) :: n  !! System size.
      integer, intent(in) :: nrhs  !! Number of right-hand sides.
      integer, intent(in) :: lda  !! Leading dimension of `a`.
      integer, intent(in) :: ldb  !! Leading dimension of `b`.
      integer, intent(out) :: ipiv(*)  !! LAPACK pivot indices.
      real(wp), intent(inout) :: a(lda, *)  !! Coefficient matrix (overwritten by LU factors).
      real(wp), intent(inout) :: b(ldb, *)  !! Right-hand side(s), overwritten by solution(s).
      integer, intent(out) :: info  !! LAPACK status code.
    end subroutine dgesv

    subroutine dgbsv(n, kl, ku, nrhs, ab, ldab, ipiv, b, ldb, info)
      import :: wp
      implicit none
      integer, intent(in) :: n  !! System size.
      integer, intent(in) :: kl  !! Number of sub-diagonals.
      integer, intent(in) :: ku  !! Number of super-diagonals.
      integer, intent(in) :: nrhs  !! Number of right-hand sides.
      integer, intent(in) :: ldab  !! Leading dimension of `ab`.
      integer, intent(in) :: ldb  !! Leading dimension of `b`.
      integer, intent(out) :: ipiv(*)  !! LAPACK pivot indices.
      real(wp), intent(inout) :: ab(ldab, *)  !! Banded coefficient matrix (overwritten by LU factors).
      real(wp), intent(inout) :: b(ldb, *)  !! Right-hand side(s), overwritten by solution(s).
      integer, intent(out) :: info  !! LAPACK status code.
    end subroutine dgbsv
  end interface

  type :: PTCSolver
    integer :: neq = 0  !! Number of unknowns in the state vector.
    integer :: jacobian_type = 0  !! Jacobian mode (`PTC_JAC_DENSE` or `PTC_JAC_BAND`).
    integer :: kl = 0  !! Number of sub-diagonals for banded Jacobians.
    integer :: ku = 0  !! Number of super-diagonals for banded Jacobians.
    integer :: ldab = 0  !! Leading dimension for LAPACK banded matrix storage.

    procedure(rhs_fcn), pointer, nopass :: f => null()  !! User residual callback computing `f(x)`.
    procedure(jac_fcn), pointer, nopass :: jac => null()  !! User Jacobian callback (dense or compact banded layout).
    procedure(verify_step_fcn), pointer, nopass :: verify => null()  !! Optional step verification callback.
    procedure(timestep_fcn), pointer, nopass :: compute_dt => null()  !! Optional timestep update callback.

    real(wp) :: dt = 0.0_wp  !! Current pseudo-time step size.
    real(wp) :: dt_initial = 0.0_wp  !! Initial pseudo-time step size.
    real(wp) :: dt_increment = 1.1_wp  !! Growth factor used by default timestep adaptation.
    real(wp) :: dt_max = 0.0_wp  !! Maximum pseudo-time step (`<=0` disables cap).
    logical :: increment_dt_from_initial_dt = .false.  !! If true, adapt from initial `(dt, fnorm)` pair.

    real(wp) :: fatol = 1.0e-50_wp  !! Absolute convergence tolerance on residual metric.
    real(wp) :: frtol = 1.0e-12_wp  !! Relative convergence tolerance on residual metric.
    logical :: use_weighted_norm = .false.  !! If true, use WRMS norm as the residual metric.
    real(wp) :: weighted_rtol = 0.0_wp  !! Relative tolerance used in WRMS weights.

    real(wp) :: fnorm = -1.0_wp  !! Current residual metric (2-norm or WRMS norm).
    real(wp) :: fnorm_initial = -1.0_wp  !! Residual metric at first accepted step.
    real(wp) :: fnorm_previous = -1.0_wp  !! Residual metric from previous accepted step.
    real(wp) :: fnorm_l2 = -1.0_wp  !! Current unweighted residual 2-norm `||f(x)||_2`.
    real(wp) :: fnorm_wrms = -1.0_wp  !! Current weighted residual WRMS norm.

    integer :: steps = 0  !! Number of accepted pseudo-steps.
    integer :: rejects_total = 0  !! Total number of rejected step attempts.
    integer :: max_reject = 10  !! Maximum rejects allowed per `step()` before failure.
    integer :: max_steps = 10000  !! Maximum accepted steps allowed in `solve()`.
    integer :: reason = PTC_REASON_NONE  !! Solver state/reason code.

    logical :: initialized = .false.  !! True once arrays and callbacks are configured.

    real(wp), allocatable :: x(:)  !! Current solution iterate.
    real(wp), allocatable :: x_old(:)  !! Backup iterate for rollback on rejected steps.
    real(wp), allocatable :: fvec(:)  !! Residual workspace.
    real(wp), allocatable :: step_vec(:)  !! Linear correction vector workspace.
    real(wp), allocatable :: rhs_mat(:, :)  !! Right-hand-side workspace passed to LAPACK solvers.

    real(wp), allocatable :: jac_mat(:, :)  !! Jacobian workspace (dense or compact banded).
    real(wp), allocatable :: a_dense(:, :)  !! Dense system matrix workspace for `(I/dt - J)`.

    ! Banded Jacobian compact storage (LAPACK standard):
    ! jac_mat(ku+1+i-j, j) = J(i,j), for max(1,j-ku) <= i <= min(n,j+kl)
    real(wp), allocatable :: a_band(:, :)  !! Banded system matrix workspace for LAPACK `dgbsv`.
    real(wp), allocatable :: weighted_atol(:)  !! Per-component absolute tolerances used in WRMS weighting.

    integer, allocatable :: ipiv(:)  !! Pivot indices returned by LAPACK factorizations.
  contains
    procedure :: initialize => PTCSolver_initialize
    procedure :: step => PTCSolver_step
    procedure :: solve => PTCSolver_solve
    procedure :: check_convergence => PTCSolver_check_convergence

    procedure :: set_verify_timestep => PTCSolver_set_verify_timestep
    procedure :: set_compute_timestep => PTCSolver_set_compute_timestep
  end type PTCSolver

contains

  !> Initialize solver state, allocate work arrays, and register user callbacks.
  !!
  !! Configures dense or banded Jacobian storage, sets PETSc-like defaults,
  !! and optionally applies user-provided tolerances and stepping controls.
  subroutine PTCSolver_initialize(self, x0, f, jacobian_type, dt0, jac, kl, ku, fatol, frtol, dt_increment, dt_max, increment_dt_from_initial_dt, max_reject, max_steps, weighted_rtol, weighted_atol)
    class(PTCSolver), intent(inout) :: self  !! Solver object to initialize.
    real(wp), intent(in) :: x0(:)  !! Initial state guess.
    procedure(rhs_fcn) :: f  !! User residual callback.
    integer, intent(in) :: jacobian_type  !! Jacobian mode (`PTC_JAC_DENSE` or `PTC_JAC_BAND`).
    real(wp), intent(in) :: dt0  !! Initial pseudo-time step.
    procedure(jac_fcn) :: jac  !! User Jacobian callback.
    integer, intent(in), optional :: kl  !! Number of sub-diagonals for banded Jacobian mode.
    integer, intent(in), optional :: ku  !! Number of super-diagonals for banded Jacobian mode.
    integer, intent(in), optional :: max_reject  !! Maximum rejections allowed per `step()` call.
    integer, intent(in), optional :: max_steps  !! Maximum accepted steps allowed in `solve()`.
    real(wp), intent(in), optional :: fatol  !! Absolute residual-norm convergence tolerance.
    real(wp), intent(in), optional :: frtol  !! Relative residual-norm convergence tolerance.
    real(wp), intent(in), optional :: dt_increment  !! Default timestep growth factor.
    real(wp), intent(in), optional :: dt_max  !! Maximum allowed timestep (non-positive means no cap).
    logical, intent(in), optional :: increment_dt_from_initial_dt  !! Optional switch for initial-reference dt adaptation.
    real(wp), intent(in), optional :: weighted_rtol  !! Relative tolerance in WRMS weighting (`rtol` term).
    real(wp), intent(in), optional :: weighted_atol(:)  !! Per-component absolute tolerances in WRMS weighting.

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
    if (present(weighted_rtol) .or. present(weighted_atol)) then
      if (.not. present(weighted_rtol) .or. .not. present(weighted_atol)) then
        self%reason = PTC_DIVERGED_INVALID_INPUT
        return
      end if
      if (weighted_rtol <= 0.0_wp .or. size(weighted_atol) /= self%neq) then
        self%reason = PTC_DIVERGED_INVALID_INPUT
        return
      end if
      if (any(weighted_atol <= 0.0_wp)) then
        self%reason = PTC_DIVERGED_INVALID_INPUT
        return
      end if
      self%use_weighted_norm = .true.
      self%weighted_rtol = weighted_rtol
      allocate(self%weighted_atol(self%neq))
      self%weighted_atol = weighted_atol
    end if

    self%fnorm = -1.0_wp
    self%fnorm_initial = -1.0_wp
    self%fnorm_previous = -1.0_wp
    self%fnorm_l2 = -1.0_wp
    self%fnorm_wrms = -1.0_wp

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

  !> Register a callback to verify/possibly reject each candidate pseudo-step.
  subroutine PTCSolver_set_verify_timestep(self, verify)
    class(PTCSolver), intent(inout) :: self  !! Solver object to update.
    procedure(verify_step_fcn) :: verify  !! User callback for accept/reject decisions.

    self%verify => verify
  end subroutine PTCSolver_set_verify_timestep

  !> Register a callback that overrides default pseudo-time-step adaptation.
  subroutine PTCSolver_set_compute_timestep(self, compute_dt)
    class(PTCSolver), intent(inout) :: self  !! Solver object to update.
    procedure(timestep_fcn) :: compute_dt  !! User callback for computing next timestep.

    self%compute_dt => compute_dt
  end subroutine PTCSolver_set_compute_timestep

  !> Advance the solver by one accepted pseudo-step (with internal retries).
  !!
  !! Performs linearized PTC update, optional verify callback, default/custom
  !! timestep adaptation, and convergence checks.
  subroutine PTCSolver_step(self)
    class(PTCSolver), intent(inout) :: self  !! Solver object advanced by one accepted pseudo-step.

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

  !> Repeatedly call `step()` until convergence or a terminal failure reason.
  subroutine PTCSolver_solve(self)
    class(PTCSolver), intent(inout) :: self  !! Solver object advanced until termination.

    if (.not. self%initialized) then
      self%reason = PTC_DIVERGED_NOT_INITIALIZED
      return
    end if

    do while (self%reason == PTC_REASON_NONE)
      if (self%steps >= self%max_steps) then
        self%reason = PTC_DIVERGED_MAX_STEPS
        exit
      end if
      call self%step()
    end do
  end subroutine PTCSolver_solve

  !> Evaluate stopping criteria based on absolute and relative residual norms.
  subroutine PTCSolver_check_convergence(self)
    class(PTCSolver), intent(inout) :: self  !! Solver object whose residual norms are tested.

    if (self%use_weighted_norm) then
      if (self%fnorm < 1.0_wp) then
        self%reason = PTC_CONVERGED_PSEUDO_FATOL
        return
      end if
    else
      if (self%fnorm < self%fatol) then
        self%reason = PTC_CONVERGED_PSEUDO_FATOL
        return
      end if
    end if

    if (self%fnorm_initial > 0.0_wp) then
      if ((self%fnorm / self%fnorm_initial) < self%frtol) then
        self%reason = PTC_CONVERGED_PSEUDO_FRTOL
        return
      end if
    end if

    self%reason = PTC_REASON_NONE
  end subroutine PTCSolver_check_convergence

  !> Compute residual vector and its active metric norm at a given state.
  subroutine PTCSolver_compute_residual(self, x, fvec, fnorm, ierr)
    class(PTCSolver), intent(inout) :: self  !! Solver object providing residual callback.
    real(wp), intent(in) :: x(:)  !! State at which to evaluate residual.
    real(wp), intent(out) :: fvec(:)  !! Residual vector `f(x)`.
    real(wp), intent(out) :: fnorm  !! Residual metric norm used by timestep/adaptation logic.
    integer, intent(out) :: ierr  !! Callback status (`0` success, nonzero failure).

    call self%f(x, fvec, ierr)
    if (ierr /= 0) then
      fnorm = -1.0_wp
      return
    end if

    self%fnorm_l2 = norm2(fvec)
    if (self%use_weighted_norm) then
      self%fnorm_wrms = PTCSolver_compute_wrms_norm(self, x, fvec)
      fnorm = self%fnorm_wrms
    else
      self%fnorm_wrms = -1.0_wp
      fnorm = self%fnorm_l2
    end if
  end subroutine PTCSolver_compute_residual

  !> Compute WRMS norm for a residual vector using CVODE-style weights.
  function PTCSolver_compute_wrms_norm(self, x, rvec) result(wrms_norm)
    class(PTCSolver), intent(in) :: self  !! Solver object providing WRMS weighting parameters.
    real(wp), intent(in) :: x(:)  !! Current state vector.
    real(wp), intent(in) :: rvec(:)  !! Residual vector.
    real(wp) :: wrms_norm  !! Weighted root-mean-square norm.

    real(wp) :: scale(size(x))

    scale = self%weighted_atol + self%weighted_rtol * abs(x)
    wrms_norm = sqrt(sum((rvec / scale) ** 2) / real(size(x), wp))
  end function PTCSolver_compute_wrms_norm

  !> Perform one linearized pseudo-transient correction solve and state update.
  !!
  !! Solves `(I/dt - J) s = f(x)` using dense (`dgesv`) or banded (`dgbsv`)
  !! LAPACK routines, then updates `x <- x + s`.
  subroutine PTCSolver_take_newton_update(self, ierr)
    class(PTCSolver), intent(inout) :: self  !! Solver object updated in place.
    integer, intent(out) :: ierr  !! Status (`0` success, positive rejectable failure, negative fatal failure).

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

  !> Compute the next pseudo-time step from residual norms.
  !!
  !! Uses either user callback `compute_dt` or the default TSPSEUDO-style
  !! residual-ratio formula with optional cap `dt_max`.
  subroutine PTCSolver_compute_next_dt(self, next_dt, ierr)
    class(PTCSolver), intent(inout) :: self  !! Solver object providing residual history and controls.
    real(wp), intent(out) :: next_dt  !! Computed timestep for next accepted step.
    integer, intent(out) :: ierr  !! Status (`0` success, nonzero failure).

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

  !> Reject current attempt, rollback state, update `dt`, and count rejection.
  subroutine PTCSolver_reject_step(self, rejections, new_dt)
    class(PTCSolver), intent(inout) :: self  !! Solver object to rollback and update.
    integer, intent(inout) :: rejections  !! Rejection counter for current `step()` call.
    real(wp), intent(in) :: new_dt  !! Timestep to use after rejection.

    self%x = self%x_old
    self%dt = max(new_dt, tiny(1.0_wp))
    self%rejects_total = self%rejects_total + 1

    rejections = rejections + 1
    if (rejections > self%max_reject) then
      self%reason = PTC_DIVERGED_STEP_REJECTED
    end if
  end subroutine PTCSolver_reject_step

  !> Deallocate solver work arrays and clear callback pointers.
  subroutine reset_storage(self)
    class(PTCSolver), intent(inout) :: self  !! Solver object whose allocations/pointers are cleared.

    if (allocated(self%x)) deallocate(self%x)
    if (allocated(self%x_old)) deallocate(self%x_old)
    if (allocated(self%fvec)) deallocate(self%fvec)
    if (allocated(self%step_vec)) deallocate(self%step_vec)
    if (allocated(self%rhs_mat)) deallocate(self%rhs_mat)
    if (allocated(self%jac_mat)) deallocate(self%jac_mat)
    if (allocated(self%a_dense)) deallocate(self%a_dense)
    if (allocated(self%a_band)) deallocate(self%a_band)
    if (allocated(self%weighted_atol)) deallocate(self%weighted_atol)
    if (allocated(self%ipiv)) deallocate(self%ipiv)

    self%initialized = .false.
    self%use_weighted_norm = .false.
    self%weighted_rtol = 0.0_wp
    self%f => null()
    self%jac => null()
    self%verify => null()
    self%compute_dt => null()
  end subroutine reset_storage

end module pseudo_transient
