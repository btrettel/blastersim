! tests for I/O procedures
! Standard: Fortran 2018
! Preprocessor: none
! Author: Ben Trettel (<http://trettel.us/>)
! Project: [BlasterSim](https://github.com/btrettel/blastersim)
! License: [AGPLv3](https://www.gnu.org/licenses/agpl-3.0.en.html)

program test_io

use prec, only: WP
use units
use unittest, only: test_results_type
implicit none

type(test_results_type) :: tests

call tests%start_tests("io.nml")

! TODO: test `write_latex_engineering`
! TODO: test `create_barrel`

call test_read_pneumatic_namelist_non_default(tests)
call test_read_pneumatic_namelist_default(tests)

call test_read_springer_namelist_non_default(tests)
call test_read_springer_namelist_default(tests)

call test_calculate_eta_pneumatic(tests)
call test_calculate_eta_springer(tests)

call tests%end_tests()

contains

subroutine test_read_pneumatic_namelist_non_default(tests)
    use prec, only: CL, PI
    use cva, only: cv_system_type, run_config_type
    use port, only: path_join
    use io, only: I_BARREL, read_pneumatic_namelist
    
    type(test_results_type), intent(in out) :: tests
    
    character(len=CL)                 :: input_file, path_array(2)
    type(cv_system_type), allocatable :: sys
    type(run_config_type)             :: config
    integer                           :: rc
    
    integer, parameter :: I_CHAMBER = 2, I_BARREL_ATM  = 3
    
    type(si_length)   :: d_barrel, d_e, d_chamber, l_travel
    type(si_volume)   :: vol_dead
    type(si_velocity) :: actual_v_muzzle
    integer           :: actual_rc
    
    path_array(1) = "test"
    path_array(2) = "pneumatic-non-default.nml"
    input_file = path_join(path_array)
    
    call read_pneumatic_namelist(input_file, sys, config, rc, actual_v_muzzle_=actual_v_muzzle, actual_rc_=actual_rc)
    
    call tests%integer_eq(rc, 0, "read_pneumatic_namelist, non-default values, rc")
    
    if (rc == 0) then
        call tests%real_eq(config%dt%v%v, 2.0e-7_WP, "read_pneumatic_namelist, non-default values, dt")
        call tests%logical_false(config%csv_output, "read_pneumatic_namelist, non-default values, csv_output")
        call tests%integer_eq(config%csv_frequency, 1, "read_pneumatic_namelist, non-default values, csv_frequency")
        
        d_barrel = sqrt((4.0_WP/PI)*sys%cv(I_BARREL)%csa)
        call tests%real_eq(d_barrel%v%v, 13.0e-3_WP, "read_pneumatic_namelist, non-default values, d_barrel")
        
        d_e = sqrt((4.0_WP/PI)*sys%con(I_CHAMBER, I_BARREL)%a_e)
        call tests%real_eq(d_e%v%v, 3.0e-3_WP, "read_pneumatic_namelist, non-default values, d_e")
        
        d_chamber = sqrt((4.0_WP/PI)*sys%cv(I_CHAMBER)%csa)
        call tests%real_eq(d_chamber%v%v, 12.5e-3_WP, "read_pneumatic_namelist, non-default values, d_chamber")
        
        l_travel = sys%cv(I_BARREL)%x_stop - sys%cv(I_BARREL)%x
        call tests%real_eq(l_travel%v%v, 40.0e-2_WP, "read_pneumatic_namelist, non-default values, l_travel")
        
        call tests%real_eq(1.0_WP/sys%cv(I_BARREL)%rm_p%v%v, 1.0e-3_WP, "read_pneumatic_namelist, non-default values, m_proj")
        call tests%real_eq(sys%cv(I_BARREL)%p_fd%v%v, 3.0e3_WP, "read_pneumatic_namelist, non-default values, p_fd_proj")
        call tests%real_eq(sys%cv(I_BARREL)%p_fs%v%v, 7.0e3_WP, "read_pneumatic_namelist, non-default values, p_fs_proj")
        
        vol_dead = sys%cv(I_BARREL)%x*sys%cv(I_BARREL)%csa
        call tests%real_eq(vol_dead%v%v, 1.0e-5_WP, "read_pneumatic_namelist, non-default values, vol_dead")
        
        call tests%character_eq(config%id, "pneumatic-non-default", "read_pneumatic_namelist, non-default values, id")
        
        call tests%real_eq(sys%con(I_CHAMBER, I_BARREL)%b%v%v, 0.6_WP, "read_pneumatic_namelist, non-default values, b")
        call tests%real_eq(sys%cv(I_BARREL_ATM)%p_const%v%v, 102000.0_WP, "read_pneumatic_namelist, non-default values, p_atm")
        call tests%real_eq(sys%cv(I_BARREL_ATM)%temp_const%v%v, 295.0_WP, &
                                    "read_pneumatic_namelist, non-default values, temp_atm")
        call tests%real_eq(actual_v_muzzle%v%v, 70.0_WP, "read_pneumatic_namelist, non-default values, actual_v_muzzle")
        call tests%integer_eq(actual_rc, 2, "read_pneumatic_namelist, non-default values, actual_rc")
    end if
end subroutine test_read_pneumatic_namelist_non_default

subroutine test_read_pneumatic_namelist_default(tests)
    use prec, only: CL, PI
    use cva, only: DT_DEFAULT, CSV_FREQUENCY_DEFAULT, X_GT_X_STOP_RUN_RC, cv_system_type, run_config_type
    use port, only: path_join
    use io, only: I_BARREL, read_pneumatic_namelist
    use gasdata, only: P_ATM, TEMP_ATM
    
    type(test_results_type), intent(in out) :: tests
    
    character(len=CL)                 :: input_file, path_array(2)
    type(cv_system_type), allocatable :: sys
    type(run_config_type)             :: config
    integer                           :: rc
    
    integer, parameter :: I_CHAMBER = 2, I_BARREL_ATM  = 3
    
    type(si_length)   :: d_barrel, d_e, d_chamber, l_travel
    type(si_volume)   :: vol_dead
    type(si_velocity) :: actual_v_muzzle
    integer           :: actual_rc
    
    path_array(1) = "test"
    path_array(2) = "pneumatic-default.nml"
    input_file = path_join(path_array)
    
    call read_pneumatic_namelist(input_file, sys, config, rc, actual_v_muzzle_=actual_v_muzzle, actual_rc_=actual_rc)
    
    call tests%integer_eq(rc, 0, "read_pneumatic_namelist, default values, rc")
    
    if (rc == 0) then
        call tests%real_eq(config%dt%v%v, DT_DEFAULT, "read_pneumatic_namelist, default values, dt")
        call tests%logical_true(config%csv_output, "read_pneumatic_namelist, default values, csv_output")
        call tests%integer_eq(config%csv_frequency, CSV_FREQUENCY_DEFAULT, &
                                "read_pneumatic_namelist, default values, csv_frequency")
        
        d_barrel = sqrt((4.0_WP/PI)*sys%cv(I_BARREL)%csa)
        call tests%real_eq(d_barrel%v%v, 13.0e-3_WP, "read_pneumatic_namelist, default values, d_barrel")
        
        d_e = sqrt((4.0_WP/PI)*sys%con(I_CHAMBER, I_BARREL)%a_e)
        call tests%real_eq(d_e%v%v, 3.0e-3_WP, "read_pneumatic_namelist, default values, d_e")
        
        d_chamber = sqrt((4.0_WP/PI)*sys%cv(I_CHAMBER)%csa)
        call tests%real_eq(d_chamber%v%v, 12.5e-3_WP, "read_pneumatic_namelist, default values, d_chamber")
        
        l_travel = sys%cv(I_BARREL)%x_stop - sys%cv(I_BARREL)%x
        call tests%real_eq(l_travel%v%v, 40.0e-2_WP, "read_pneumatic_namelist, default values, l_travel")
        
        call tests%real_eq(1.0_WP/sys%cv(I_BARREL)%rm_p%v%v, 1.0e-3_WP, "read_pneumatic_namelist, default values, m_proj")
        call tests%real_eq(sys%cv(I_BARREL)%p_fd%v%v, 3.0e3_WP, "read_pneumatic_namelist, default values, p_fd_proj")
        call tests%real_eq(sys%cv(I_BARREL)%p_fs%v%v, 7.0e3_WP, "read_pneumatic_namelist, default values, p_fs_proj")
        
        vol_dead = sys%cv(I_BARREL)%x*sys%cv(I_BARREL)%csa
        call tests%real_eq(vol_dead%v%v, 1.0e-5_WP, "read_pneumatic_namelist, default values, vol_dead")
        
        call tests%character_eq(config%id, "pneumatic-default", "read_pneumatic_namelist, default values, id")
        
        call tests%real_eq(sys%con(I_CHAMBER, I_BARREL)%b%v%v, 0.5_WP, "read_pneumatic_namelist, default values, b")
        call tests%real_eq(sys%cv(I_BARREL_ATM)%p_const%v%v, P_ATM, "read_pneumatic_namelist, default values, p_atm")
        call tests%real_eq(sys%cv(I_BARREL_ATM)%temp_const%v%v, TEMP_ATM, "read_pneumatic_namelist, default values, temp_atm")
        call tests%real_eq(actual_v_muzzle%v%v, 0.0_WP, "read_pneumatic_namelist, default values, actual_v_muzzle")
        call tests%integer_eq(actual_rc, X_GT_X_STOP_RUN_RC, "read_pneumatic_namelist, default values, actual_rc")
    end if
end subroutine test_read_pneumatic_namelist_default

subroutine test_read_springer_namelist_non_default(tests)
    use prec, only: CL, PI
    use cva, only: cv_system_type, run_config_type
    use port, only: path_join
    use io, only: I_BARREL, read_springer_namelist
    
    type(test_results_type), intent(in out) :: tests
    
    character(len=CL)                 :: input_file, path_array(2)
    type(cv_system_type), allocatable :: sys
    type(run_config_type)             :: config
    integer                           :: rc
    
    integer, parameter :: I_PLUNGER = 2, I_BARREL_ATM  = 3, I_PLUNGER_ATM = 4
    
    type(si_length)   :: d_barrel, d_e, d_plunger, l_travel
    type(si_volume)   :: vol_dead
    type(si_velocity) :: actual_v_muzzle
    integer           :: actual_rc
    
    path_array(1) = "test"
    path_array(2) = "springer-non-default.nml"
    input_file = path_join(path_array)
    
    call read_springer_namelist(input_file, sys, config, rc, actual_v_muzzle_=actual_v_muzzle, actual_rc_=actual_rc)
    
    call tests%integer_eq(rc, 0, "read_springer_namelist, non-default values, rc")
    
    if (rc == 0) then
        call tests%real_eq(config%dt%v%v, 1.0e-7_WP, "read_springer_namelist, non-default values, dt")
        call tests%logical_false(config%csv_output, "read_springer_namelist, non-default values, csv_output")
        call tests%integer_eq(config%csv_frequency, 1, "read_springer_namelist, non-default values, csv_frequency")
        
        d_barrel = sqrt((4.0_WP/PI)*sys%cv(I_BARREL)%csa)
        call tests%real_eq(d_barrel%v%v, 13.0e-3_WP, "read_springer_namelist, non-default values, d_barrel")
        
        d_e = sqrt((4.0_WP/PI)*sys%con(I_PLUNGER, I_BARREL)%a_e)
        call tests%real_eq(d_e%v%v, 5.0e-3_WP, "read_springer_namelist, non-default values, d_e")
        
        d_plunger = sqrt((4.0_WP/PI)*sys%cv(I_PLUNGER)%csa)
        call tests%real_eq(d_plunger%v%v, 30.0e-3_WP, "read_springer_namelist, non-default values, d_plunger")
        
        call tests%real_eq(sys%cv(I_PLUNGER)%k%v%v, 500.0_WP, "read_springer_namelist, non-default values, k")
        call tests%real_eq(sys%cv(I_PLUNGER)%x%v%v - sys%cv(I_PLUNGER)%x_min%v%v, 20.0e-2_WP, &
                                "read_springer_namelist, non-default values, l_draw")
        call tests%real_eq(sys%cv(I_PLUNGER)%delta_pre%v%v, 5.0e-2_WP, "read_springer_namelist, non-default values, delta_pre")
        call tests%real_eq(sys%cv(I_PLUNGER)%cor%v%v, 0.5_WP, "read_springer_namelist, non-default values, cor")
        
        l_travel = sys%cv(I_BARREL)%x_stop - sys%cv(I_BARREL)%x
        call tests%real_eq(l_travel%v%v, 50.0e-2_WP, "read_springer_namelist, non-default values, l_travel")
        
        call tests%real_eq(1.0_WP/sys%cv(I_PLUNGER)%rm_p%v%v, 30.0e-3_WP, &
                            "read_springer_namelist, non-default values, m_plunger")
        call tests%real_eq(1.0_WP/sys%cv(I_BARREL)%rm_p%v%v, 1.1e-3_WP, "read_springer_namelist, non-default values, m_proj")
        call tests%real_eq(sys%cv(I_PLUNGER)%m_spring%v%v, 10.0e-3_WP, "read_springer_namelist, non-default values, m_spring")
        call tests%real_eq(sys%cv(I_PLUNGER)%p_fd%v%v, 1.0e3_WP, "read_springer_namelist, non-default values, p_fd_plunger")
        call tests%real_eq(sys%cv(I_BARREL)%p_fd%v%v, 2.0e3_WP, "read_springer_namelist, non-default values, p_fd_proj")
        call tests%real_eq(sys%cv(I_PLUNGER)%p_fs%v%v, 4.0e3_WP, "read_springer_namelist, non-default values, p_fs_plunger")
        call tests%real_eq(sys%cv(I_BARREL)%p_fs%v%v, 5.0e3_WP, "read_springer_namelist, non-default values, p_fs_proj")
        
        vol_dead = sys%cv(I_BARREL)%x*sys%cv(I_BARREL)%csa
        call tests%real_eq(vol_dead%v%v, 10.0e-6_WP, "read_springer_namelist, non-default values, vol_dead (barrel)")
        
        vol_dead = sys%cv(I_PLUNGER)%x_min*sys%cv(I_PLUNGER)%csa
        call tests%real_eq(vol_dead%v%v, 10.0e-6_WP, "read_springer_namelist, non-default values, vol_dead (plunger)")
        
        call tests%character_eq(config%id, "springer-non-default", "read_springer_namelist, non-default values, id")
        
        call tests%real_eq(sys%con(I_PLUNGER, I_BARREL)%b%v%v, 0.528_WP, "read_springer_namelist, non-default values, b")
        call tests%real_eq(sys%cv(I_PLUNGER_ATM)%p_const%v%v, 100000.0_WP, &
                                    "read_springer_namelist, non-default values, p_atm (1)")
        call tests%real_eq(sys%cv(I_BARREL_ATM)%p_const%v%v, 100000.0_WP, &
                                    "read_springer_namelist, non-default values, p_atm (2)")
        call tests%real_eq(sys%cv(I_PLUNGER_ATM)%temp_const%v%v, 310.0_WP, &
                            "read_springer_namelist, non-default values, temp_atm (1)")
        call tests%real_eq(sys%cv(I_BARREL_ATM)%temp_const%v%v, 310.0_WP, &
                            "read_springer_namelist, non-default values, temp_atm (2)")
        call tests%real_eq(actual_v_muzzle%v%v, 60.0_WP, "read_springer_namelist, non-default values, actual_v_muzzle")
        call tests%integer_eq(actual_rc, 1, "read_springer_namelist, non-default values, actual_rc")
    end if
end subroutine test_read_springer_namelist_non_default

subroutine test_read_springer_namelist_default(tests)
    use prec, only: CL, PI
    use cva, only: DT_DEFAULT, CSV_FREQUENCY_DEFAULT, X_GT_X_STOP_RUN_RC, cv_system_type, run_config_type
    use port, only: path_join
    use io, only: I_BARREL, read_springer_namelist
    use gasdata, only: P_ATM, TEMP_ATM
    
    type(test_results_type), intent(in out) :: tests
    
    character(len=CL)                 :: input_file, path_array(2)
    type(cv_system_type), allocatable :: sys
    type(run_config_type)             :: config
    integer                           :: rc
    
    integer, parameter :: I_PLUNGER = 2, I_BARREL_ATM  = 3, I_PLUNGER_ATM = 4
    
    type(si_length)   :: d_barrel, d_e, d_plunger, l_travel
    type(si_volume)   :: vol_dead
    type(si_velocity) :: actual_v_muzzle
    integer           :: actual_rc
    
    path_array(1) = "test"
    path_array(2) = "springer-default.nml"
    input_file = path_join(path_array)
    
    call read_springer_namelist(input_file, sys, config, rc, actual_v_muzzle_=actual_v_muzzle, actual_rc_=actual_rc)
    
    call tests%integer_eq(rc, 0, "read_springer_namelist, default values, rc")
    
    if (rc == 0) then
        call tests%real_eq(config%dt%v%v, DT_DEFAULT, "read_springer_namelist, default values, dt")
        call tests%logical_true(config%csv_output, "read_springer_namelist, default values, csv_output")
        call tests%integer_eq(config%csv_frequency, CSV_FREQUENCY_DEFAULT, &
                                "read_springer_namelist, default values, csv_frequency")
        
        d_barrel = sqrt((4.0_WP/PI)*sys%cv(I_BARREL)%csa)
        call tests%real_eq(d_barrel%v%v, 13.0e-3_WP, "read_springer_namelist, default values, d_barrel")
        
        d_e = sqrt((4.0_WP/PI)*sys%con(I_PLUNGER, I_BARREL)%a_e)
        call tests%real_eq(d_e%v%v, 5.0e-3_WP, "read_springer_namelist, default values, d_e")
        
        d_plunger = sqrt((4.0_WP/PI)*sys%cv(I_PLUNGER)%csa)
        call tests%real_eq(d_plunger%v%v, 30.0e-3_WP, "read_springer_namelist, default values, d_plunger")
        
        call tests%real_eq(sys%cv(I_PLUNGER)%k%v%v, 500.0_WP, "read_springer_namelist, default values, k")
        call tests%real_eq(sys%cv(I_PLUNGER)%x%v%v - sys%cv(I_PLUNGER)%x_min%v%v, 20.0e-2_WP, &
                                    "read_springer_namelist, default values, l_draw")
        call tests%real_eq(sys%cv(I_PLUNGER)%delta_pre%v%v, 5.0e-2_WP, "read_springer_namelist, default values, delta_pre")
        call tests%real_eq(sys%cv(I_PLUNGER)%cor%v%v, 0.0_WP, "read_springer_namelist, default values, cor")
        
        l_travel = sys%cv(I_BARREL)%x_stop - sys%cv(I_BARREL)%x
        call tests%real_eq(l_travel%v%v, 50.0e-2_WP, "read_springer_namelist, default values, d_plunger")
        
        call tests%real_eq(1.0_WP/sys%cv(I_PLUNGER)%rm_p%v%v, 30.0e-3_WP, "read_springer_namelist, default values, m_plunger")
        call tests%real_eq(1.0_WP/sys%cv(I_BARREL)%rm_p%v%v, 1.1e-3_WP, "read_springer_namelist, default values, m_proj")
        call tests%real_eq(sys%cv(I_PLUNGER)%m_spring%v%v, 10.0e-3_WP, "read_springer_namelist, default values, m_spring")
        call tests%real_eq(sys%cv(I_PLUNGER)%p_fd%v%v, 1.0e3_WP, "read_springer_namelist, default values, p_fd_plunger")
        call tests%real_eq(sys%cv(I_BARREL)%p_fd%v%v, 2.0e3_WP, "read_springer_namelist, default values, p_fd_proj")
        call tests%real_eq(sys%cv(I_PLUNGER)%p_fs%v%v, 4.0e3_WP, "read_springer_namelist, default values, p_fs_plunger")
        call tests%real_eq(sys%cv(I_BARREL)%p_fs%v%v, 5.0e3_WP, "read_springer_namelist, default values, p_fs_proj")
        
        vol_dead = sys%cv(I_BARREL)%x*sys%cv(I_BARREL)%csa
        call tests%real_eq(vol_dead%v%v, 10.0e-6_WP, "read_springer_namelist, default values, vol_dead (barrel)")
        
        vol_dead = sys%cv(I_PLUNGER)%x_min*sys%cv(I_PLUNGER)%csa
        call tests%real_eq(vol_dead%v%v, 10.0e-6_WP, "read_springer_namelist, default values, vol_dead (plunger)")
        
        call tests%character_eq(config%id, "springer-default", "read_springer_namelist, default values, id")
        
        call tests%real_eq(sys%con(I_PLUNGER, I_BARREL)%b%v%v, 0.5_WP, "read_springer_namelist, default values, b")
        call tests%real_eq(sys%cv(I_PLUNGER_ATM)%p_const%v%v, P_ATM, "read_springer_namelist, default values, p_atm (1)")
        call tests%real_eq(sys%cv(I_BARREL_ATM)%p_const%v%v, P_ATM, "read_springer_namelist, default values, p_atm (2)")
        call tests%real_eq(sys%cv(I_PLUNGER_ATM)%temp_const%v%v, TEMP_ATM, &
                                    "read_springer_namelist, default values, temp_atm (1)")
        call tests%real_eq(sys%cv(I_BARREL_ATM)%temp_const%v%v, TEMP_ATM, &
                                    "read_springer_namelist, default values, temp_atm (2)")
        call tests%real_eq(actual_v_muzzle%v%v, 0.0_WP, "read_springer_namelist, default values, actual_v_muzzle")
        call tests%integer_eq(actual_rc, X_GT_X_STOP_RUN_RC, "read_springer_namelist, default values, actual_rc")
    end if
end subroutine test_read_springer_namelist_default

subroutine test_calculate_eta_pneumatic(tests)
    use cva, only: run_config_type, cv_system_type, run_status_type, run
    use io, only: I_BARREL, PNEUMATIC_MODE, read_pneumatic_namelist, calculate_eta
    use port, only: path_join
    use prec, only: CL
    
    type(test_results_type), intent(in out) :: tests
    
    character(len=CL)     :: path_array(2), input_file
    type(run_config_type) :: config
    integer               :: rc
    type(si_velocity)     :: actual_v_muzzle
    type(unitless)        :: eta
    type(cv_system_type), allocatable :: sys_start, sys_end
    
    path_array(1) = "examples"
    path_array(2) = "pneumatic-2010-08-07-70-psi.nml"
    input_file    = path_join(path_array)
    call read_pneumatic_namelist(input_file, sys_start, config, rc, actual_v_muzzle_=actual_v_muzzle)
    
    sys_end = sys_start
    sys_end%cv(I_BARREL)%x_dot = actual_v_muzzle
    
    call calculate_eta(sys_start, sys_end, PNEUMATIC_MODE, eta)
    
    ! This number was calculated by an old spreadsheet using rounded `v_muzzle` and `vol_chamber` values matching the input file.
    call tests%real_eq(eta%v%v, 0.16049644048116973_WP, "calculate_eta, pneumatic", abs_tol=1.0e-8_WP)
end subroutine test_calculate_eta_pneumatic

subroutine test_calculate_eta_springer(tests)
    use cva, only: run_config_type, cv_system_type, run_status_type, run
    use io, only: I_BARREL, SPRINGER_MODE, read_springer_namelist, calculate_eta
    use port, only: path_join
    use prec, only: CL
    
    type(test_results_type), intent(in out) :: tests
    
    character(len=CL)     :: path_array(2), input_file
    type(run_config_type) :: config
    integer               :: rc
    type(si_velocity)     :: actual_v_muzzle
    type(unitless)        :: eta
    type(cv_system_type), allocatable :: sys_start, sys_end
    
    path_array(1) = "examples"
    path_array(2) = "radioactive-ranger-mega-sniper-k11.nml"
    input_file    = path_join(path_array)
    call read_springer_namelist(input_file, sys_start, config, rc, actual_v_muzzle_=actual_v_muzzle)
    
    sys_end = sys_start
    sys_end%cv(I_BARREL)%x_dot = actual_v_muzzle
    
    call calculate_eta(sys_start, sys_end, SPRINGER_MODE, eta)
    
    ! <https://discord.com/channels/825852031239061545/1177799131268382902/1178501509722427414>
    call tests%real_eq(eta%v%v, 0.5152_WP, "calculate_eta, springer", abs_tol=1.0e-3_WP)
end subroutine test_calculate_eta_springer

end program test_io
