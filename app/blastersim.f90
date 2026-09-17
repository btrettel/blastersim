! BlasterSim: a spring and pneumatic blaster simulator
! Standard: Fortran 2018
! Preprocessor: none
! Author: Ben Trettel (<http://trettel.us/>)
! Project: [BlasterSim](https://github.com/btrettel/blastersim)
! License: [AGPLv3](https://www.gnu.org/licenses/agpl-3.0.en.html)

program blastersim

use, intrinsic :: iso_fortran_env, only: IOSTAT_END, ERROR_UNIT, OUTPUT_UNIT
use prec, only: CL, WP
use cli, only: get_input_file_name_from_cli
use io, only: I_BARREL, I_BARREL_ATM, read_pneumatic_namelist, read_springer_namelist
use cva, only: run_config_type, cv_system_type, run_status_type, T_STOP_DEFAULT, run, MAX_ITERS_TIME_LOOP, &
                    SUCCESS_RC, TIMEOUT_RUN_RC, NEGATIVE_CV_M_TOTAL_RUN_RC, NEGATIVE_CV_TEMP_RUN_RC, &
                    MASS_TOLERANCE_RUN_RC, ENERGY_TOLERANCE_RUN_RC, MASS_DERIV_TOLERANCE_RUN_RC, &
                    ENERGY_DERIV_TOLERANCE_RUN_RC, IDEAL_EOS_RUN_RC, MIRROR_X_TOLERANCE_RUN_RC, &
                    NEGATIVE_CV_X_RUN_RC, MAX_ITERS_TIME_LOOP_RUN_RC, MAX_ITERS_GET_SYS_AT_X_RUN_RC, &
                    RK_STAGE_NEGATIVE_MASS_RC, RK_STAGE_NEGATIVE_ENERGY_RC, MAX_VELOCITY_EXCEEDED_RC
use stopcodes, only: EX_OK, EX_USAGE, EX_SOFTWARE
use rev, only: TAG, REVISION_DATE, MODIFIED
use checks, only: assert, is_close
use build, only: FUZZ
use units
implicit none

character(len=CL)                 :: input_file, extra, modified_string
type(run_config_type)             :: config
type(cv_system_type), allocatable :: sys_start, sys_end
integer                           :: rc, out_unit
type(run_status_type)             :: status
type(si_length)                   :: l_end, l_travel
real(WP)                          :: f, sum_g

extra = "<http://trettel.us/blastersim/>" // new_line("a") // "Written by Ben Trettel."

call get_input_file_name_from_cli("blastersim", input_file, extra=extra)

if (MODIFIED) then
   modified_string = ", modified"
else
   modified_string = "" 
end if

write(unit=*, fmt="(a)") "BlasterSim " // TAG // " (" // REVISION_DATE // trim(modified_string) // ")"
write(unit=*, fmt="(a)") "Running simulation..."

nml_blk: block
    call read_pneumatic_namelist(trim(input_file), sys_start, config, rc)
    if (rc == 0) then
        exit nml_blk
    else if (rc /= IOSTAT_END) then
        stop EX_USAGE, quiet=.true.
    end if
    
    call read_springer_namelist(trim(input_file), sys_start, config, rc)
    if (rc == 0) then
        exit nml_blk
    else if (rc /= IOSTAT_END) then
        stop EX_USAGE, quiet=.true.
    end if
    
    call assert(rc == IOSTAT_END, "blastersim: rc == IOSTAT_END violated")
    write(unit=ERROR_UNIT, fmt="(a)") "ERROR: Empty input file? No pneumatic or springer namelists detected. " // &
        "If your input file does have a pneumatic or springer namelist, make sure there is an empty line after the " // &
        "final /. Some Fortran compilers require a return character after the slash for a namelist to be read properly."
    stop EX_USAGE, quiet=.true.
end block nml_blk

call run(config, sys_start, sys_end, status)

call post_run_checks(sys_start, sys_end, rc)

!tripwire$ begin 724E768C Update `if (rc_read /= 0) then` sections of io.f90 to account for different total `sum_g` here.
! `sum_g` there needs to strictly be higher than `sum_g` here to encourage going through input validation.
if (FUZZ) then
    ! Write out data used in feedback-based fuzzing.
    
    ! More time steps indicates more opportunities for things to go wrong, so incentivize that.
    ! Scale it so that it's not huge.
    f = -real(status%i, WP)/real(MAX_ITERS_TIME_LOOP, WP)
    
    if (allocated(status%data)) then
        call assert(sum(status%data) >= 0.0_WP, &
                        "blastersim: all elements of status%data should be non-negative to enable fuzz testing")
        
        f = f - sum(status%data)
    end if
    
    if (status%rc < SUCCESS_RC) then
        ! If successful, no constraints are violated.
        sum_g = 0.0_WP
    else
        ! If not successful, set a constraint to incentivize the projectile leaving the barrel.
        ! With purely random testing, the vast majority of cases do not leave the barrel.
        ! So I want to test more cases that leave the barrel.
        
        ! One part of the constraint is whether the pressure is high enough to cause the projectile to move at all.
        sum_g = max(0.0_WP, (max(sys_end%cv(I_BARREL)%p_fs%v%v, sys_end%cv(I_BARREL)%p_fd%v%v) &
                                - (sys_end%cv(I_BARREL)%p_peak%v%v - sys_end%cv(I_BARREL_ATM)%p_const%v%v)) &
                                    / max(sys_end%cv(I_BARREL)%p_fs%v%v, sys_end%cv(I_BARREL)%p_fd%v%v))
        
        ! The other part of the constraint is how far the projectile moves down the barrel.
        l_travel = sys_start%cv(I_BARREL)%x_stop - sys_start%cv(I_BARREL)%x
        l_end    = sys_end%cv(I_BARREL)%x        - sys_start%cv(I_BARREL)%x
        sum_g    = sum_g + (l_travel%v%v - l_end%v%v)/l_travel%v%v
        
        call assert(.not. is_close(l_travel%v%v, 0.0_WP), "blastersim: l_travel must be /= 0")
        
        call assert(sum_g >= 0.0_WP, "blastersim: sum_g >= 0 violated", &
                        print_real=[sum_g, l_end%v%v, l_travel%v%v])
        
        ! Commented out as until I fix the friction to not have any backwards motion, `sum_g` can go above 1.
        !call assert(sum_g <= 1.0_WP, "blastersim: sum_g <= 1 violated", &
                        !print_real=[sum_g, l_end%v%v, l_travel%v%v])
    end if
    
    open(newunit=out_unit, action="write", status="replace", position="rewind", &
            file=trim(input_file) // ".out")
    write(unit=out_unit, fmt="(es24.17, 1x, es24.17)") f, sum_g
    close(unit=out_unit)
end if
!tripwire$ end

!tripwire$ begin 36CD5B66 Update `\secref{return-codes}` of usage.tex.
if (status%rc < SUCCESS_RC) then
    write(unit=OUTPUT_UNIT, fmt="(a)") "SUCCESS!"
    write(unit=OUTPUT_UNIT, fmt="(a, f0.2, a)") "muzzle velocity: ", sys_end%cv(I_BARREL)%x_dot%v%v, " m/s"
    
    stop EX_OK, quiet=.true.
else
    write(unit=ERROR_UNIT, fmt="(a, i0)") "ERROR: return code ", status%rc
    select case (status%rc)
        case (TIMEOUT_RUN_RC)
            write(unit=ERROR_UNIT, fmt="(a, f3.1, a)") "Projectile did not leave barrel after ", &
                                                            T_STOP_DEFAULT, " seconds."
            call refer_to_docs()
            stop EX_USAGE, quiet=.true.
        case (NEGATIVE_CV_M_TOTAL_RUN_RC, NEGATIVE_CV_TEMP_RUN_RC)
            write(unit=ERROR_UNIT, fmt="(2a)") "Negative mass or temperature of control volume. ", &
                                                "This is a bug that should be reported."
            call refer_to_docs()
            stop EX_SOFTWARE, quiet=.true.
        case (MASS_TOLERANCE_RUN_RC, ENERGY_TOLERANCE_RUN_RC, MASS_DERIV_TOLERANCE_RUN_RC, &
                    ENERGY_DERIV_TOLERANCE_RUN_RC)
            write(unit=ERROR_UNIT, fmt="(2a)") "Mass or energy tolerance exceeded. ", &
                                                "Decrease dt by a factor of 10 and report a bug if that doesn't help."
            call refer_to_docs()
            stop EX_USAGE, quiet=.true.
        case (IDEAL_EOS_RUN_RC)
            write(unit=ERROR_UNIT, fmt="(a)") "Critical pressure exceeded. The ideal gas law is inaccurate here. ", &
                    "BlasterSim can not handle pressures this high at the moment."
            call refer_to_docs()
            stop EX_USAGE, quiet=.true.
        case (MIRROR_X_TOLERANCE_RUN_RC)
            write(unit=ERROR_UNIT, fmt="(2a)") "Plunger position desynchronization. ", &
                                                "This is a bug that should be reported."
            call refer_to_docs()
            stop EX_SOFTWARE, quiet=.true.
        case (NEGATIVE_CV_X_RUN_RC)
            write(unit=ERROR_UNIT, fmt="(2a)") "The plunger has moved past its stopping point. ", &
                                                "This is a bug that should be reported."
            call refer_to_docs()
            stop EX_SOFTWARE, quiet=.true.
        case (MAX_ITERS_TIME_LOOP_RUN_RC, MAX_ITERS_GET_SYS_AT_X_RUN_RC)
            write(unit=ERROR_UNIT, fmt="(2a)") "Maximum number of iterations reached. ", &
                                                "This might be a bug or possibly the projectile will not leave the barrel."
            call refer_to_docs()
            stop EX_SOFTWARE, quiet=.true.
        case (RK_STAGE_NEGATIVE_MASS_RC, RK_STAGE_NEGATIVE_ENERGY_RC)
            write(unit=ERROR_UNIT, fmt="(2a)") "Negative mass of a gas species or energy detected during a Runge-Kutta stage. ", &
                    "Check whether d_e is too large, or possibly if dt is too large."
            call refer_to_docs()
            stop EX_USAGE, quiet=.true.
        case (MAX_VELOCITY_EXCEEDED_RC)
            write(unit=ERROR_UNIT, fmt="(a)") "Muzzle velocity exceeded what is physically possible. ", &
                                                "This is a bug that should be reported."
            write(unit=ERROR_UNIT, fmt="(a, f0.2, a)") "muzzle velocity: ", sys_end%cv(I_BARREL)%x_dot%v%v, " m/s"
            call refer_to_docs()
            stop EX_SOFTWARE, quiet=.true.
        case default
            write(unit=ERROR_UNIT, fmt="(a)") "Unknown error. This is a bug that should be reported."
            stop EX_SOFTWARE, quiet=.true.
    end select
end if
!tripwire$ end

contains

subroutine refer_to_docs()
    write(unit=ERROR_UNIT, fmt="(a)") "Refer to BlasterSim User's Guide for possibly more information."
    write(unit=ERROR_UNIT, fmt="(a)") "<http://trettel.us/blastersim/docs/return-codes.html>"
end subroutine refer_to_docs

pure subroutine post_run_checks(sys_start, sys_end, rc)
    use io, only: I_SOURCE
    
    type(cv_system_type), allocatable, intent(in) :: sys_start, sys_end
    integer, intent(in out)                       :: rc
    
    type(unitless)    :: y(size(sys_start%cv(1)%m_k))
    type(si_velocity) :: v_escape
    
    ! The reason this check isn't in `check_sys` is that it shouldn't apply to every CV, only the barrel.
    
    y = sys_start%cv(I_SOURCE)%y()
    
    ! corner_theory_1950 p. 364, eq. 74
    ! seigel_theory_1965 p. 19, eq. 11-8
    v_escape = 2.0_WP * sqrt(sys_start%cv(I_SOURCE)%gamma(y) * sys_start%cv(I_SOURCE)%r() * sys_start%cv(I_SOURCE)%temp()) &
                            / (sys_start%cv(I_SOURCE)%gamma(y) - 1.0_WP)
    
    if (rc < SUCCESS_RC) then
        if (sys_end%cv(I_BARREL)%x_dot > v_escape) then
            rc = MAX_VELOCITY_EXCEEDED_RC
        end if
    end if
end subroutine post_run_checks

end program blastersim
