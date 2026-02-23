! validation tests (experimental data)
! Standard: Fortran 2018
! Preprocessor: none
! Author: Ben Trettel (<http://trettel.us/>)
! Project: [BlasterSim](https://github.com/btrettel/blastersim)
! License: [AGPLv3](https://www.gnu.org/licenses/agpl-3.0.en.html)

program test_validation

use units
use unittest, only: test_results_type
implicit none

type(test_results_type) :: tests

call tests%start_tests("validation.nml")

call pneumatic_2010_08_07(tests)

call tests%end_tests()

contains

subroutine run_pneumatic_get_mv(input_file, v_muzzle_predicted, rc_predicted)
    use cva, only: run_config_type, cv_system_type, run_status_type, run
    use io, only: I_BARREL, read_pneumatic_namelist
    
    character(len=*), intent(in)   :: input_file
    type(si_velocity), intent(out) :: v_muzzle_predicted
    integer, intent(out)           :: rc_predicted
    
    type(cv_system_type), allocatable :: sys_start, sys_end
    type(run_config_type)             :: config
    integer                           :: rc
    type(run_status_type)             :: status
    
    call read_pneumatic_namelist(input_file, sys_start, config, rc)
    call run(config, sys_start, sys_end, status)
    
    v_muzzle_predicted = sys_end%cv(I_BARREL)%x_dot
    rc_predicted       = status%rc
end subroutine run_pneumatic_get_mv

!!!!!!!!!!!!!!!!!!!!
! 2010-08-07 tests !
!!!!!!!!!!!!!!!!!!!!

! Why start here? I have validation data I collected back in the day. I wasn't able to model this set of experiments properly in the past (with my theoretical model, not computational), so I should first figure out how to model this properly before making too many assumptions about what the remainder of the model will loo like.

! /home/ben/svn/old/ballistics/data/2010-08-07 muzzle velocity test.gnumeric
! Possibly the test gun: /home/ben/photos/14/DSCN0070.JPG

! I believe the QEV was a Clippard JEV-F2M2.
! I recall doing the tests in the garage in Jefferson. It was hot outside.

! <https://web.archive.org/web/20220927075508/http://btrettel.nerfers.com/forums/topic20.html>:
! > In order to validate my simplified model of the internal ballistics of a pneumatic gun, I did some tests in 2010. I first determined what pressure was necessary for the projectile to leave the barrel. I then fired 15 shots at 30 psi, 40 psi, 50 psi, 60 psi, and 70 psi. I measured the muzzle velocity in ft/s with my chronograph.
! > These tests forced me to conclude that my model was too simple as it did not capture the dynamics of the gun as well as I had hoped. At the moment I'm planning a more complicated model suited better for spud guns, but also applicable for Nerf guns.

! TODO: Check handwritten notes for more information. I can't find anything in the scans I have.

subroutine pneumatic_2010_08_07(tests)
    use convert, only: fps_const
    use prec, only: CL, WP
    use port, only: path_join
    use cva, only: TIMEOUT_RUN_RC, SUCCESS_RUN_RC
    
    type(test_results_type), intent(in out) :: tests
    
    character(len=CL) :: path_array(2), input_file
    type(si_velocity) :: v_muzzle_predicted, v_exp
    integer           :: rc_predicted
    
    path_array(1) = "examples"
    
    path_array(2) = "pneumatic-2010-08-07-25-psi.nml"
    input_file    = path_join(path_array)
    call run_pneumatic_get_mv(input_file, v_muzzle_predicted, rc_predicted)
    call tests%integer_eq(rc_predicted, TIMEOUT_RUN_RC, &
                            "pneumatic, 2010-08-07 experiments, 25 psi, status%rc (projectile did not exit)")
    
    path_array(2) = "pneumatic-2010-08-07-30-psi.nml"
    input_file    = path_join(path_array)
    call run_pneumatic_get_mv(input_file, v_muzzle_predicted, rc_predicted)
    call tests%integer_eq(rc_predicted, SUCCESS_RUN_RC, &
                            "pneumatic, 2010-08-07 experiments, 30 psi, status%rc (projectile exited)")
    call tests%real_eq(v_muzzle_predicted%v%v, 23.80285338515869_WP, &
                        "pneumatic, 2010-08-07 experiments, 30 psi, muzzle velocity (characterization)")
    v_exp = fps_const(80.465_WP, 0)
    call tests%real_eq(v_muzzle_predicted%v%v, v_exp%v%v, &
                        "pneumatic, 2010-08-07 experiments, 30 psi, muzzle velocity (validation)", abs_tol=1.0_WP)
    
    path_array(2) = "pneumatic-2010-08-07-40-psi.nml"
    input_file    = path_join(path_array)
    call run_pneumatic_get_mv(input_file, v_muzzle_predicted, rc_predicted)
    call tests%integer_eq(rc_predicted, SUCCESS_RUN_RC, &
                            "pneumatic, 2010-08-07 experiments, 40 psi, status%rc (projectile exited)")
    call tests%real_eq(v_muzzle_predicted%v%v, 35.761901429472438_WP, &
                        "pneumatic, 2010-08-07 experiments, 40 psi, muzzle velocity (characterization)")
    v_exp = fps_const(120.353_WP, 0)
    call tests%real_eq(v_muzzle_predicted%v%v, v_exp%v%v, &
                        "pneumatic, 2010-08-07 experiments, 40 psi, muzzle velocity (validation)", abs_tol=1.0_WP)
    
    path_array(2) = "pneumatic-2010-08-07-50-psi.nml"
    input_file    = path_join(path_array)
    call run_pneumatic_get_mv(input_file, v_muzzle_predicted, rc_predicted)
    call tests%integer_eq(rc_predicted, SUCCESS_RUN_RC, &
                            "pneumatic, 2010-08-07 experiments, 50 psi, status%rc (projectile exited)")
    call tests%real_eq(v_muzzle_predicted%v%v, 45.359461217896239_WP, &
                        "pneumatic, 2010-08-07 experiments, 50 psi, muzzle velocity (characterization)")
    v_exp = fps_const(145.664_WP, 0)
    call tests%real_eq(v_muzzle_predicted%v%v, v_exp%v%v, &
                        "pneumatic, 2010-08-07 experiments, 50 psi, muzzle velocity (validation)", abs_tol=1.0_WP)
    
    path_array(2) = "pneumatic-2010-08-07-60-psi.nml"
    input_file    = path_join(path_array)
    call run_pneumatic_get_mv(input_file, v_muzzle_predicted, rc_predicted)
    call tests%integer_eq(rc_predicted, SUCCESS_RUN_RC, &
                            "pneumatic, 2010-08-07 experiments, 60 psi, status%rc (projectile exited)")
    call tests%real_eq(v_muzzle_predicted%v%v, 49.27168246923017_WP, &
                        "pneumatic, 2010-08-07 experiments, 60 psi, muzzle velocity (characterization)")
    v_exp = fps_const(158.567_WP, 0)
    call tests%real_eq(v_muzzle_predicted%v%v, v_exp%v%v, &
                        "pneumatic, 2010-08-07 experiments, 60 psi, muzzle velocity (validation)", abs_tol=1.0_WP)
    
    path_array(2) = "pneumatic-2010-08-07-70-psi.nml"
    input_file    = path_join(path_array)
    call run_pneumatic_get_mv(input_file, v_muzzle_predicted, rc_predicted)
    call tests%integer_eq(rc_predicted, SUCCESS_RUN_RC, &
                            "pneumatic, 2010-08-07 experiments, 70 psi, status%rc (projectile exited)")
    call tests%real_eq(v_muzzle_predicted%v%v, 59.13020469015979_WP, &
                        "pneumatic, 2010-08-07 experiments, 70 psi, muzzle velocity (characterization)")
    v_exp = fps_const(190.987_WP, 0)
    call tests%real_eq(v_muzzle_predicted%v%v, v_exp%v%v, &
                        "pneumatic, 2010-08-07 experiments, 70 psi, muzzle velocity (validation)", abs_tol=1.0_WP)
end subroutine pneumatic_2010_08_07

end program test_validation
