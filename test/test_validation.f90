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

call pneumatic_validation(tests)

call tests%end_tests()

contains

!subroutine v_muzzle_predicted_vs_observed(v_muzzle_predicted, actual_v_muzzle, actual_v_muzzle_stdev)
!    type(si_velocity), intent(in) :: v_muzzle_predicted(:), actual_v_muzzle(:), actual_v_muzzle_stdev(:)
    
    
!end subroutine v_muzzle_predicted_vs_observed

subroutine run_pneumatic_get_mv(input_file, v_muzzle_predicted, rc_predicted, actual_v_muzzle, actual_v_muzzle_stdev, &
                                actual_v_muzzle_n, actual_rc)
    use cva, only: run_config_type, cv_system_type, run_status_type, run
    use io, only: I_BARREL, read_pneumatic_namelist
    
    character(len=*), intent(in)   :: input_file
    type(si_velocity), intent(out) :: v_muzzle_predicted
    integer, intent(out)           :: rc_predicted
    type(si_velocity), intent(out) :: actual_v_muzzle, actual_v_muzzle_stdev
    integer, intent(out)           :: actual_v_muzzle_n, actual_rc
    
    type(cv_system_type), allocatable :: sys_start, sys_end
    type(run_config_type)             :: config
    integer                           :: rc
    type(run_status_type)             :: status
    
    call read_pneumatic_namelist(input_file, sys_start, config, rc, actual_v_muzzle_=actual_v_muzzle, &
                                    actual_v_muzzle_stdev_=actual_v_muzzle_stdev, actual_v_muzzle_n_=actual_v_muzzle_n, &
                                    actual_rc_=actual_rc)
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

subroutine pneumatic_validation(tests)
    ! For now, I'll compare the estimated mean against the confidence interval of the mean for the experiment.
    ! Later, when I have UQ, I'll do something more sophisticated like looking at the K-L divergence.
    ! This current setup will fail if *any* data point is off.
    ! Given the lack of data and me calibrating `d_e` in each input file specifically rather than having one `d_e` value
    ! for each valve, this is a weaker form of validation. It's more of a conceptual validation, looking at whether BlasterSim
    ! can in principle reproduce the data.
    ! Later, overall validation will be determined by a global statistic like R^2 (though not R^2 itself as that's not the best).
    
    use prec, only: CL, WP
    use port, only: path_join
    
    type(test_results_type), intent(in out) :: tests
    
    character(len=CL) :: INPUT_FILES(6) = ["pneumatic-2010-08-07-25-psi.nml", &
                                           "pneumatic-2010-08-07-30-psi.nml", &
                                           "pneumatic-2010-08-07-40-psi.nml", &
                                           "pneumatic-2010-08-07-50-psi.nml", &
                                           "pneumatic-2010-08-07-60-psi.nml", &
                                           "pneumatic-2010-08-07-70-psi.nml"]
    real(WP), parameter :: Z_HALF_ALPHA = 1.96_WP
    
    character(len=CL)                               :: path_array(2), input_file
    type(si_velocity), dimension(size(INPUT_FILES)) :: v_muzzle_predicted, actual_v_muzzle, actual_v_muzzle_stdev
    integer                                         :: i, rc_predicted, actual_v_muzzle_n, actual_rc
    
    path_array(1) = "examples"
    
    do i = 1, size(INPUT_FILES)
        path_array(2) = trim(INPUT_FILES(i))
        input_file    = path_join(path_array)
        call run_pneumatic_get_mv(input_file, v_muzzle_predicted(i), rc_predicted, actual_v_muzzle(i), &
                                    actual_v_muzzle_stdev(i), actual_v_muzzle_n, actual_rc)
        
        call tests%integer_eq(rc_predicted, actual_rc, trim(INPUT_FILES(i)) // ", status%rc")
        if (actual_rc == 0) call tests%real_eq(v_muzzle_predicted(i)%v%v, actual_v_muzzle(i)%v%v, &
                                    trim(INPUT_FILES(i)) // ", muzzle velocity (validation)", &
                                    abs_tol=Z_HALF_ALPHA*actual_v_muzzle_stdev(i)%v%v/sqrt(real(actual_v_muzzle_n, WP)))
    end do
end subroutine pneumatic_validation

end program test_validation
