program robertson_ptc_demo
  use pseudo_transient, only: PTCSolver, wp, PTC_JAC_DENSE, PTC_JAC_BAND, PTC_REASON_NONE
  implicit none

  call test()

contains

  subroutine test()
    real(wp), parameter :: y0(3) = [1.0_wp, 0.0_wp, 0.0_wp]
    real(wp) :: x_dense(3), x_band(3)
    real(wp) :: fnorm_dense, fnorm_band
    real(wp) :: mass_drift_dense, mass_drift_band, diff_norm
    integer :: reason_dense, reason_band
    integer :: steps_dense, steps_band
    integer :: rejects_dense, rejects_band

    call run_case("dense", PTC_JAC_DENSE, y0, x_dense, fnorm_dense, reason_dense, steps_dense, rejects_dense)
    call run_case("band", PTC_JAC_BAND, y0, x_band, fnorm_band, reason_band, steps_band, rejects_band)

    mass_drift_dense = sum(x_dense) - sum(y0)
    mass_drift_band = sum(x_band) - sum(y0)
    diff_norm = norm2(x_dense - x_band)

    write (*, '(a)') ''
    write (*, '(a)') 'Dense vs band comparison'
    write (*, '(a,1x,es12.4)') '  ||x_dense - x_band||_2 =', diff_norm
    write (*, '(a,1x,es12.4)') '  Dense mass drift       =', mass_drift_dense
    write (*, '(a,1x,es12.4)') '  Band mass drift        =', mass_drift_band
    write (*, '(a,1x,l1)') '  Both converged         =', (reason_dense > PTC_REASON_NONE .and. reason_band > PTC_REASON_NONE)
  end subroutine test

  subroutine run_case(label, jacobian_type, y0, x_final, fnorm_final, reason, steps, rejects)
    character(len=*), intent(in) :: label
    integer, intent(in) :: jacobian_type
    real(wp), intent(in) :: y0(:)
    real(wp), intent(out) :: x_final(:)
    real(wp), intent(out) :: fnorm_final
    integer, intent(out) :: reason
    integer, intent(out) :: steps
    integer, intent(out) :: rejects

    type(PTCSolver) :: solver
    real(wp), parameter :: dt0 = 1.0e-8_wp
    real(wp) :: relnorm

    call solver%initialize(y0, robertson_rhs, jacobian_type, dt0, robertson_jac, &
      kl=1, ku=2, max_steps=200000, dt_max=1.0e20_wp)

    write (*, '(a,a)') 'Stepping case: ', trim(label)
    write (*, '(a)') '  step                  dt            abs_norm            rel_norm'
    do while (solver%reason == PTC_REASON_NONE)
      call solver%step()
      if (solver%steps > 0 .and. solver%fnorm_initial > 0.0_wp) then
        relnorm = solver%fnorm / solver%fnorm_initial
        write (*, '(2x,i6,3(2x,es18.10))') solver%steps, solver%dt, solver%fnorm, relnorm
      end if
    end do

    x_final = solver%x
    fnorm_final = solver%fnorm
    reason = solver%reason
    steps = solver%steps
    rejects = solver%rejects_total

    write (*, '(a)') ''
    write (*, '(a,a)') 'Case: ', trim(label)
    write (*, '(a,1x,i0,1x,a)') '  reason  =', reason, trim(reason_name(reason))
    write (*, '(a,1x,i0)') '  steps   =', steps
    write (*, '(a,1x,i0)') '  rejects =', rejects
    write (*, '(a,1x,es12.4)') '  fnorm   =', fnorm_final
    write (*, '(a,1x,es12.4,1x,es12.4,1x,es12.4)') '  x       =', x_final(1), x_final(2), x_final(3)
    write (*, '(a,1x,es12.4)') '  sum(x)-1=', sum(x_final) - 1.0_wp
  end subroutine run_case

  subroutine robertson_rhs(u, udot, ierr)
    real(wp), intent(in) :: u(:)
    real(wp), intent(out) :: udot(:)
    integer, intent(out) :: ierr

    if (size(u) /= 3 .or. size(udot) /= 3) then
      ierr = 1
      return
    end if

    udot(1) = -0.04_wp * u(1) + 1.0e4_wp * u(2) * u(3)
    udot(2) = 0.04_wp * u(1) - 1.0e4_wp * u(2) * u(3) - 3.0e7_wp * u(2) * u(2)
    udot(3) = 3.0e7_wp * u(2) * u(2)

    ierr = 0
  end subroutine robertson_rhs

  subroutine robertson_jac(u, jac, ierr)
    real(wp), intent(in) :: u(:)
    real(wp), intent(out) :: jac(:, :)
    integer, intent(out) :: ierr

    if (size(u) /= 3 .or. size(jac, 2) /= 3) then
      ierr = 1
      return
    end if

    jac = 0.0_wp

    if (size(jac, 1) == 3) then
      jac(1, 1) = -0.04_wp
      jac(1, 2) = 1.0e4_wp * u(3)
      jac(1, 3) = 1.0e4_wp * u(2)

      jac(2, 1) = 0.04_wp
      jac(2, 2) = -1.0e4_wp * u(3) - 6.0e7_wp * u(2)
      jac(2, 3) = -1.0e4_wp * u(2)

      jac(3, 2) = 6.0e7_wp * u(2)

    else if (size(jac, 1) == 4) then
      ! Compact banded storage for n=3, kl=1, ku=2:
      !   jac(ku+1+i-j, j) = J(i,j)
      jac(3, 1) = -0.04_wp
      jac(4, 1) = 0.04_wp
      jac(2, 2) = 1.0e4_wp * u(3)
      jac(3, 2) = -1.0e4_wp * u(3) - 6.0e7_wp * u(2)
      jac(4, 2) = 6.0e7_wp * u(2)
      jac(1, 3) = 1.0e4_wp * u(2)
      jac(2, 3) = -1.0e4_wp * u(2)

    else
      ierr = 1
      return
    end if

    ierr = 0
  end subroutine robertson_jac

  function reason_name(reason) result(name)
    use pseudo_transient, only: PTC_CONVERGED_PSEUDO_FATOL, PTC_CONVERGED_PSEUDO_FRTOL, &
      PTC_DIVERGED_STEP_REJECTED, PTC_DIVERGED_CALLBACK_FATAL, PTC_DIVERGED_NOT_INITIALIZED, &
      PTC_DIVERGED_INVALID_INPUT, PTC_DIVERGED_MAX_STEPS, PTC_REASON_NONE
    integer, intent(in) :: reason
    character(len=40) :: name

    select case (reason)
    case (PTC_REASON_NONE)
      name = 'PTC_REASON_NONE'
    case (PTC_CONVERGED_PSEUDO_FATOL)
      name = 'PTC_CONVERGED_PSEUDO_FATOL'
    case (PTC_CONVERGED_PSEUDO_FRTOL)
      name = 'PTC_CONVERGED_PSEUDO_FRTOL'
    case (PTC_DIVERGED_STEP_REJECTED)
      name = 'PTC_DIVERGED_STEP_REJECTED'
    case (PTC_DIVERGED_CALLBACK_FATAL)
      name = 'PTC_DIVERGED_CALLBACK_FATAL'
    case (PTC_DIVERGED_NOT_INITIALIZED)
      name = 'PTC_DIVERGED_NOT_INITIALIZED'
    case (PTC_DIVERGED_INVALID_INPUT)
      name = 'PTC_DIVERGED_INVALID_INPUT'
    case (PTC_DIVERGED_MAX_STEPS)
      name = 'PTC_DIVERGED_MAX_STEPS'
    case default
      name = 'PTC_REASON_UNKNOWN'
    end select
  end function reason_name

end program robertson_ptc_demo
