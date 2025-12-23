! Module for procedures used for convergence testing.
! Standard: Fortran 2018
! Preprocessor: none
! Author: Ben Trettel (<http://trettel.us/>)
! Project: [flt](https://github.com/btrettel/flt)
! License: [GPLv3](https://www.gnu.org/licenses/gpl-3.0.en.html)

module convergence

use fmad, only: ad
use prec, only: WP
implicit none
private

public :: dnorm, convergence_test, logspace

interface dnorm
    ! <https://numpy.org/doc/stable/reference/generated/numpy.linalg.norm.html>
    ! <https://help.imsl.com/fortran/fnlmath/current/NORM.htm>
    ! <https://gcc.gnu.org/onlinedocs/gfortran/NORM2.html>
    
    ! <https://www.petercheng.me/blog/discrete-L2-norm>
    ! <https://scicomp.stackexchange.com/questions/5095/discrete-lp-norms-for-non-uniform-grid>
    ! roy_review_2005 eq. 7: the discrete norm takes into account `n` so that it's comparable to a point-wise value.
    ! <https://youtu.be/v9OnNeYH4Ok?t=2216>
    
    module procedure dnorm_real_rank_1
    module procedure dnorm_ad_rank_1
end interface dnorm

contains

subroutine convergence_test(n_arr, solver_ne, p_expected, message, tests, p_tol, p_d_tol)
    use, intrinsic :: iso_fortran_env, only: ERROR_UNIT
    use checks, only: assert, assert_dimension
    use unittest, only: test_results_type
    use fmad, only: log
    
    integer, intent(in)                     :: n_arr(:) ! array of number of grid cells, time steps, Monte Carlo samples, etc.
    real(WP), intent(in)                    :: p_expected(:) ! expected order of convergence
    character(len=*), intent(in)            :: message
    type(test_results_type), intent(in out) :: tests
    real(WP), intent(in), optional          :: p_tol(:)
    real(WP), intent(in), optional          :: p_d_tol(:)
    
    interface
        subroutine solver_ne(n, ne, ne_d)
            use fmad, only: ad
            use prec, only: WP
            
            integer, intent(in)                :: n          ! number of grid cells, time steps, Monte Carlo samples, etc.
            type(ad), intent(out), allocatable :: ne(:)      ! numerical error for value(s) (n_var)
            real(WP), intent(out), allocatable :: ne_d(:, :) ! numerical error for derivatives (n_var, n_d)
            
            ! `n_var` is the number of variables to return numerical error for.
            ! `n_d` is the size of the gradient vector for AD.
            
            ! This is not `pure` to make debugging easier.
            
            ! Exact or manufactured solutions are called in this function.
            
            ! Numerical error (usually discretization error) is calculated in here.
            ! A norm can be used or a local metric can be used.
            
            ! Instead pass out `ne` and calculate `ne_d` in `convergence_test`?
            ! Start as-is, later figure out how to refactor to simplify.
            
            ! Previously, additional tests could be added to be used with `tests`.
            ! However, I decided that feature wasn't going to be used often and wasn't worth adding here.
            ! I could use normal assertions.
        end subroutine solver_ne
    end interface
    
    integer               :: i_n, n_n     ! index of `n_arr` and size of `n_arr`
    integer               :: i_var, n_var ! index for dependent variables and number of dependent variables
    integer               :: i_d, n_d     ! index for derivatives and number of derivatives
    integer               :: n_failures
    type(ad), allocatable :: ne_i_n(:), ne(:, :) ! ne(n_n, n_var)
    real(WP), allocatable :: ne_d_i_n(:, :), ne_d(:, :, :) ! ne_d(n_n, n_var, n_d)
    type(ad), allocatable :: p_v(:)
    real(WP), allocatable :: p_tol_(:), p_d_tol_(:), p_d(:, :)
    character(len=32)     :: i_var_string, i_d_string
    logical               :: stdout
    
    if (present(p_tol)) then
        p_tol_ = p_tol
    else
        allocate(p_tol_(size(p_expected)))
        p_tol_ = 0.05_WP
    end if
    
    if (present(p_d_tol)) then
        p_d_tol_ = p_d_tol
    else
        allocate(p_d_tol_(size(p_expected)))
        p_d_tol_ = 0.05_WP
    end if
    
    call assert_dimension(p_tol_, p_expected)
    call assert_dimension(p_d_tol_, p_expected)
    
    call assert(size(n_arr) >= 2,       "convergence (convergence_test): n_arr should have at least 2 elements to calculate p")
    call assert(all(n_arr > 0),         "convergence (convergence_test): n can not be zero or negative")
    call assert(len(message) > 0,       "convergence (convergence_test): message can not be empty")
    call assert(all(p_tol_ > 0.0_WP),   "convergence (convergence_test): p_tol is too small")
    call assert(all(p_d_tol_ > 0.0_WP), "convergence (convergence_test): p_d_tol is too small")
    
    n_n = size(n_arr)
    
    ! Suppress test messages while printing the table.
    stdout       = tests%stdout
    tests%stdout = .false.
    n_failures   = tests%n_failures
    
    print "(2a)", message, ":"
    print "(3a6, 2a14)", "n", "var #", "v/d", "ne", "p"
    ! MAYBE: Run convergence tests in parallel later?
    do i_n = 1, n_n ! SERIAL
        call solver_ne(n_arr(i_n), ne_i_n, ne_d_i_n)
        
        if (i_n == 1) then
            n_var = size(ne_i_n)
            allocate(ne(n_n, n_var))
            allocate(p_v(n_var))
            n_d = size(ne_i_n(1)%d)
            allocate(p_d(n_var, n_d))
            allocate(ne_d(n_n, n_var, n_d))
            call assert_dimension(p_v, p_expected)
            
            do i_var = 1, n_var ! SERIAL
                call assert_numerical_error(ne_i_n(i_var)%v, i_var, 0)
                ne(i_n, i_var) = ne_i_n(i_var)
                
                if (i_var == 1) then
                    print "(2i6, a6, es14.5)", n_arr(i_n), i_var, "v", ne(i_n, i_var)%v
                else
                    print "(a6, i6, a6, es14.5)", "", i_var, "v", ne(i_n, i_var)%v
                end if
                do i_d = 1, n_d ! SERIAL
                    call assert_numerical_error(ne_d_i_n(i_var, i_d), i_var, i_d)
                    ne_d(i_n, i_var, i_d) = ne_d_i_n(i_var, i_d)
                    print "(a6, 2i6, es14.5)", "", i_var, i_d, ne_d(i_n, i_var, i_d)
                end do
            end do
        else
            do i_var = 1, n_var ! SERIAL
                call assert_numerical_error(ne_i_n(i_var)%v, i_var, 0)
                ne(i_n, i_var) = ne_i_n(i_var)
                
                ! order of accuracy; see roy_review_2005 eq. 6 or 8
                p_v(i_var) = log(ne(i_n, i_var) / ne(i_n - 1, i_var)) &
                                / log(real(n_arr(i_n - 1), WP) / real(n_arr(i_n), WP))
                
                if (i_var == 1) then
                    print "(2i6, a6, es14.5, f14.6)", n_arr(i_n), i_var, "v", ne(i_n, i_var)%v, p_v(i_var)%v
                else
                    print "(a6, i6, a6, es14.5, f14.6)", "", i_var, "v", ne(i_n, i_var)%v, p_v(i_var)%v
                end if
                
                do i_d = 1, n_d ! SERIAL
                    call assert_numerical_error(ne_d_i_n(i_var, i_d), i_var, i_d)
                    ne_d(i_n, i_var, i_d) = ne_d_i_n(i_var, i_d)
                    
                    ! order of accuracy; see roy_review_2005 eq. 6 or 8
                    p_d(i_var, i_d) = log(ne_d(i_n, i_var, i_d) / ne_d(i_n - 1, i_var, i_d)) &
                                            / log(real(n_arr(i_n - 1), WP) / real(n_arr(i_n), WP))
                    print "(a6, 2i6, es14.5, f14.6)", "", i_var, i_d, ne_d(i_n, i_var, i_d), p_d(i_var, i_d)
                end do
            end do
        end if
    end do
    
    ! Re-enable test failure messages.
    if (tests%n_failures > n_failures) then
        write(unit=ERROR_UNIT, fmt="(a)") "One or more test failures were suppressed during convergence test. Check the test log."
    end if
    tests%stdout = stdout
    
    ! Check that the orders of accuracy are as expected.
    do i_var = 1, n_var ! SERIAL
        call assert(p_expected(i_var) > 0.0_WP, "convergence (convergence_test): p_expected is zero or negative, " &
                                                // "which probably isn't desired")
        
        write(unit=i_var_string, fmt="(i0)") i_var
        call tests%real_eq(p_v(i_var)%v, p_expected(i_var), message // ", p_v=expected, var=" // trim(i_var_string), &
                                abs_tol=p_tol_(i_var))
        
        do i_d = 1, n_d ! SERIAL
            write(unit=i_d_string, fmt="(i0)") i_d
            call tests%real_eq(p_v(i_var)%d(i_d), 0.0_WP, message // ", p_v%d(" // trim(i_d_string) &
                                    // ")=0, var=" // trim(i_var_string), abs_tol=p_d_tol_(i_var))
            call tests%real_eq(p_d(i_var, i_d), p_expected(i_var), message // ", p_d(" // trim(i_d_string) &
                                    // ")=expected, var=" // trim(i_var_string), abs_tol=p_d_tol_(i_var))
        end do
    end do
end subroutine convergence_test

subroutine assert_numerical_error(ne, i_var, i_d)
    use checks, only: assert, TOL_FACTOR
    
    real(WP), intent(in) :: ne
    integer, intent(in)  :: i_var, i_d
    
    character(len=32) :: ne_string, arg_string
    
    write(unit=ne_string, fmt="(es14.5)") ne
    
    if (i_d == 0) then ! `i_d = 0` is used for the value here.
        write(unit=arg_string, fmt="(a, i0, a, i0, a)") "ne(i_var=", i_var, ", i_d=", i_d, ")"
    else
        write(unit=arg_string, fmt="(a, i0, a, i0, a)") "ne_d(i_var=", i_var, ", i_d=", i_d, ")"
    end if
    
    call assert(ne > TOL_FACTOR * spacing(0.0_WP), &
                    "convergence (convergence_test): " // trim(adjustl(arg_string)) // "=" // trim(adjustl(ne_string)) &
                        // ", but it must be > 0. " &
                        // "If one or more variables are expected to be exact, test that separately with real_eq.")
end subroutine assert_numerical_error

pure function dnorm_real_rank_1(x, ord, lower, upper)
    use checks, only: assert
    
    real(WP), intent(in)          :: x(:)
    integer, intent(in), optional :: ord, lower, upper
    
    ! `lower` and `upper` are used for cases where not all indices are to be summed over.
    ! For example, if you have ghost cells, those cells have fictitious data that should not be summed over.
    
    real(WP) :: dnorm_real_rank_1
    
    integer :: ord_, lower_, upper_, i, n
    
    if (present(ord)) then
        ord_ = ord
    else
        ord_ = 2
    end if
    
    if (present(lower)) then
        lower_ = lower
    else
        lower_ = 1
    end if
    
    if (present(upper)) then
        upper_ = upper
    else
        upper_ = size(x)
    end if
    
    call assert(lower_ >= lbound(x, dim=1), "convergence (dnorm_ad_rank_1): lower index bound must be >= lbound")
    call assert(upper_ <= ubound(x, dim=1), "convergence (dnorm_ad_rank_1): upper index bound must be >= ubound")
    call assert(lower_ < upper_, "convergence (dnorm_ad_rank_1): lower index bound must be above upper index bound")
    
    dnorm_real_rank_1 = 0.0_WP
    if (ord_ == huge(1)) then
        ! $l_\infty$ norm
        
        do i = lower_, upper_ ! SERIAL
            dnorm_real_rank_1 = max(dnorm_real_rank_1, abs(x(i)))
        end do
    else
        n = 0
        do i = lower_, upper_ ! SERIAL
            n = n + 1
            dnorm_real_rank_1 = dnorm_real_rank_1 + abs(x(i))**ord_
        end do
        dnorm_real_rank_1 = (dnorm_real_rank_1 / real(n, WP))**(1.0_WP/real(ord_, WP))
    end if
    
    call assert(dnorm_real_rank_1 >= 0.0_WP, "convergence (dnorm_real_rank_1): negative norm?")
end function dnorm_real_rank_1

pure function dnorm_ad_rank_1(x, ord, lower, upper)
    use fmad, only: max, abs
    use checks, only: assert
    
    type(ad), intent(in)          :: x(:)
    integer, intent(in), optional :: ord, lower, upper
    
    type(ad) :: dnorm_ad_rank_1
    
    integer :: ord_, lower_, upper_, i, n
    
    if (present(ord)) then
        ord_ = ord
    else
        ord_ = 2
    end if
    
    if (present(lower)) then
        lower_ = lower
    else
        lower_ = 1
    end if
    
    if (present(upper)) then
        upper_ = upper
    else
        upper_ = size(x)
    end if
    
    call assert(lower_ >= lbound(x, dim=1), "convergence (dnorm_ad_rank_1): lower index bound must be >= lbound")
    call assert(upper_ <= ubound(x, dim=1), "convergence (dnorm_ad_rank_1): upper index bound must be >= ubound")
    call assert(lower_ < upper_, "convergence (dnorm_ad_rank_1): lower index bound must be above upper index bound")
    
    call dnorm_ad_rank_1%init_const(0.0_WP, size(x(1)%d))
    if (ord_ == huge(1)) then
        ! $l_\infty$ norm
        
        do i = lower_, upper_ ! SERIAL
            dnorm_ad_rank_1 = max(dnorm_ad_rank_1, abs(x(i)))
        end do
    else
        n = 0
        do i = lower_, upper_ ! SERIAL
            n = n + 1
            dnorm_ad_rank_1 = dnorm_ad_rank_1 + abs(x(i))**ord_
        end do
        dnorm_ad_rank_1 = (dnorm_ad_rank_1 / real(n, WP))**(1.0_WP/real(ord_, WP))
    end if
    
    call assert(dnorm_ad_rank_1%v >= 0.0_WP, "convergence (dnorm_ad_rank_1): negative norm?")
end function dnorm_ad_rank_1

pure function logspace(loglower, logupper, n)
    ! <https://numpy.org/doc/stable/reference/generated/numpy.logspace.html>
    ! <https://www.mathworks.com/help/matlab/ref/logspace.html>
    
    use checks, only: assert, is_close
    
    real(WP), intent(in)          :: loglower, logupper
    integer, intent(in), optional :: n
    
    integer, allocatable :: logspace(:)
    
    integer  :: i, n_
    real(WP) :: logdelta
    
    if (present(n)) then
        n_ = n
    else
        n_ = 50
    end if
    
    allocate(logspace(n_))
    
    call assert(loglower < logupper, "convergence (logspace): loglower must be less than logupper")
    call assert(n_ >= 2, "convergence (logspace): n must be greater than or equal to 2")
    
    logdelta = (logupper - loglower) / real(n_ - 1, WP)
    do concurrent (i = 1:n_)
        logspace(i) = nint(10.0_WP**(loglower + logdelta * real(i - 1, WP)))
    end do
end function logspace

end module convergence
