! BlasterSim: a spring and pneumatic blaster simulator
! Standard: Fortran 2018
! Preprocessor: none
! Author: Ben Trettel (<http://trettel.us/>)
! Project: [BlasterSim](https://github.com/btrettel/blastersim)
! License: [AGPLv3](https://www.gnu.org/licenses/agpl-3.0.en.html)

program blastersim

use, intrinsic :: iso_fortran_env, only: IOSTAT_END, ERROR_UNIT, OUTPUT_UNIT
use prec, only: CL
use cli, only: get_input_file_name_from_cli
use io, only: I_BARREL, read_pneumatic_namelist, read_springer_namelist
use cva, only: run_config_type, cv_system_type, run_status_type, T_STOP_DEFAULT, run, &
                    SUCCESS_RC, TIMEOUT_RUN_RC, NEGATIVE_CV_M_TOTAL_RUN_RC, NEGATIVE_CV_TEMP_RUN_RC, &
                    MASS_TOLERANCE_RUN_RC, ENERGY_TOLERANCE_RUN_RC, MASS_DERIV_TOLERANCE_RUN_RC, &
                    ENERGY_DERIV_TOLERANCE_RUN_RC, IDEAL_EOS_RUN_RC, MIRROR_X_TOLERANCE_RUN_RC, &
                    NEGATIVE_CV_X_RUN_RC, MAX_ITERS_TIME_LOOP_RUN_RC, MAX_ITERS_GET_SYS_AT_X_RUN_RC, &
                    RK_STAGE_NEGATIVE_MASS_RC, RK_STAGE_NEGATIVE_ENERGY_RC
use stopcodes, only: EX_OK, EX_USAGE, EX_SOFTWARE
use rev, only: TAG, REVISION_DATE, MODIFIED
use checks, only: assert
implicit none

character(len=CL)                 :: input_file, extra, modified_string
type(run_config_type)             :: config
type(cv_system_type), allocatable :: sys_start, sys_end
integer                           :: rc
type(run_status_type)             :: status

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

!tripwire$ begin 4DA5CFBA Update `\secref{run-time-checks}` of verval.tex.
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
                                                "This is a bug that should be reported."
            call refer_to_docs()
            stop EX_SOFTWARE, quiet=.true.
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
            write(unit=ERROR_UNIT, fmt="(2a)") "Maximum number of iterations exceeded. ", &
                                                "This is a bug that should be reported."
            call refer_to_docs()
            stop EX_SOFTWARE, quiet=.true.
        case (RK_STAGE_NEGATIVE_MASS_RC, RK_STAGE_NEGATIVE_ENERGY_RC)
            write(unit=ERROR_UNIT, fmt="(2a)") "Negative mass of a gas species or energy detected during a Runge-Kutta stage. ", &
                    "Check whether d_e is too large, or possibly if dt is too large."
            call refer_to_docs()
            stop EX_USAGE, quiet=.true.
        case default
            write(unit=ERROR_UNIT, fmt="(a)") "Unknown error. This is a bug that should be reported."
            stop EX_SOFTWARE, quiet=.true.
    end select
end if
!tripwire$ end

contains

subroutine refer_to_docs()
    write(unit=ERROR_UNIT, fmt="(a)") "Refer to BlasterSim User's Guide for possibly more information."
    write(unit=ERROR_UNIT, fmt="(a)") "<http://trettel.us/blastersim/docs/verification.html#run-time-checks>"
end subroutine refer_to_docs

end program blastersim
