! validation tests (experimental data)
! Standard: Fortran 2018
! Preprocessor: none
! Author: Ben Trettel (<http://trettel.us/>)
! Project: [BlasterSim](https://github.com/btrettel/blastersim)
! License: [AGPLv3](https://www.gnu.org/licenses/agpl-3.0.en.html)

program test_validation

use units
use prec, only: WP
use unittest, only: test_results_type
implicit none

real(WP), parameter :: Z_HALF_ALPHA = 1.96_WP

type(test_results_type) :: tests

call tests%start_tests("validation.nml")

call pneumatic_validation(tests)

call tests%end_tests()

contains

subroutine predicted_v_muzzle_vs_observed(output_basename, predicted_v_muzzle, actual_v_muzzle, actual_v_muzzle_stdev, &
                                            actual_v_muzzle_n)
    character(len=*), intent(in)  :: output_basename
    type(si_velocity), intent(in) :: predicted_v_muzzle(:), actual_v_muzzle(:), actual_v_muzzle_stdev(:)
    integer, intent(in)           :: actual_v_muzzle_n(:)
    
    integer           :: gp_unit, i_data, n_data
    type(si_velocity) :: maximum_muzzle_velocity, muzzle_velocity_error
    
    n_data = size(predicted_v_muzzle)
    call maximum_muzzle_velocity%v%init_const(0.0_WP, size(predicted_v_muzzle(1)%v%d))
    
    open(newunit=gp_unit, action="write", status="replace", position="rewind", file=output_basename//".gp", delim="quote")
    write(unit=gp_unit, fmt="(a)") '# auto-generated'
    
    write(unit=gp_unit, fmt="(a)") '$data << EOD'
    write(unit=gp_unit, fmt="(a)") '# predicted_v_muzzle actual_v_muzzle actual_v_muzzle_error'
    do i_data = 1, n_data
        muzzle_velocity_error = Z_HALF_ALPHA*actual_v_muzzle_stdev(i_data)/sqrt(real(actual_v_muzzle_n(i_data), WP))
        write(unit=gp_unit, fmt="(g0, a, g0, a, g0)") predicted_v_muzzle(i_data)%v%v, " ", actual_v_muzzle(i_data)%v%v, " ", &
                                                        muzzle_velocity_error%v%v
        maximum_muzzle_velocity = max(maximum_muzzle_velocity, predicted_v_muzzle(i_data))
        maximum_muzzle_velocity = max(maximum_muzzle_velocity, actual_v_muzzle(i_data) + muzzle_velocity_error)
    end do
    write(unit=gp_unit, fmt="(a)") 'EOD'
    
    maximum_muzzle_velocity%v%v = real(ceiling(maximum_muzzle_velocity%v%v / 10.0_WP), WP)*10.0_WP
    
    write(unit=gp_unit, fmt="(a)") '$perfect << EOD'
    write(unit=gp_unit, fmt="(a)") '0 0'
    write(unit=gp_unit, fmt="(g0, a, g0)") maximum_muzzle_velocity%v%v, " ", maximum_muzzle_velocity%v%v
    
    write(unit=gp_unit, fmt="(a)") 'EOD'
    write(unit=gp_unit, fmt="(a)") 'set term tikz'
    write(unit=gp_unit, fmt="(3a)") 'set output "', output_basename, '.tikz"'
    write(unit=gp_unit, fmt="(a, g0, a)") 'set xrange [0:', maximum_muzzle_velocity%v%v, ']'
    write(unit=gp_unit, fmt="(a, g0, a)") 'set yrange [0:', maximum_muzzle_velocity%v%v, ']'
    write(unit=gp_unit, fmt="(a)") 'set xlabel "predicted average muzzle velocity (m/s)"'
    write(unit=gp_unit, fmt="(a)") 'set ylabel "actual average muzzle velocity (m/s)" offset -1,0'
    write(unit=gp_unit, fmt="(a)") 'set grid'
    write(unit=gp_unit, fmt="(a)") 'unset key'
    write(unit=gp_unit, fmt="(2a)") 'plot $perfect with lines linecolor "black" dashtype ".", ', &
                                        '$data with yerrorbars linecolor "black"'
    
    ! Why use `pngcairo` and not `png`? `png` doesn't support dashed lines.
    write(unit=gp_unit, fmt="(a)") 'set term pngcairo'
    write(unit=gp_unit, fmt="(3a)") 'set output "', output_basename, '.png"'
    write(unit=gp_unit, fmt="(a)") 'set xlabel "predicted average muzzle velocity (m/s)"'
    write(unit=gp_unit, fmt="(a)") 'set ylabel "actual average muzzle velocity (m/s)"'
    write(unit=gp_unit, fmt="(a)") 'set grid'
    write(unit=gp_unit, fmt="(a)") 'unset key'
    write(unit=gp_unit, fmt="(2a)") 'plot $perfect with lines linecolor "black" dashtype ". ", ', &
                                        '$data with yerrorbars linecolor "black"'
    close(gp_unit)
end subroutine predicted_v_muzzle_vs_observed

subroutine run_pneumatic_get_mv(input_file, predicted_v_muzzle, rc_predicted, actual_v_muzzle, actual_v_muzzle_stdev, &
                                actual_v_muzzle_n, actual_rc)
    use cva, only: run_config_type, cv_system_type, run_status_type, run
    use io, only: I_BARREL, read_pneumatic_namelist
    
    character(len=*), intent(in)   :: input_file
    type(si_velocity), intent(out) :: predicted_v_muzzle
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
    
    if (status%rc == 0) then
        predicted_v_muzzle = sys_end%cv(I_BARREL)%x_dot
    else
        call predicted_v_muzzle%v%init_const(0.0_WP, size(sys_end%cv(I_BARREL)%x_dot%v%d))
    end if
    rc_predicted = status%rc
end subroutine run_pneumatic_get_mv

subroutine pneumatic_validation(tests)
    ! For now, I'll compare the estimated mean against the confidence interval of the mean for the experiment.
    ! Later, when I have UQ, I'll do something more sophisticated like looking at the K-L divergence.
    ! This current setup will fail if *any* data point is off.
    ! Given the lack of data and me calibrating `d_e` in each input file specifically rather than having one `d_e` value
    ! for each valve, this is a weaker form of validation. It's more of a conceptual validation, looking at whether BlasterSim
    ! can in principle reproduce the data.
    ! Later, overall validation will be determined by a global statistic like R^2 (though not R^2 itself as that's not the best).
    
    use prec, only: CL
    use port, only: path_join
    
    type(test_results_type), intent(in out) :: tests
    
    character(len=CL) :: INPUT_FILES(6) = ["pneumatic-2010-08-07-25-psi.nml", &
                                           "pneumatic-2010-08-07-30-psi.nml", &
                                           "pneumatic-2010-08-07-40-psi.nml", &
                                           "pneumatic-2010-08-07-50-psi.nml", &
                                           "pneumatic-2010-08-07-60-psi.nml", &
                                           "pneumatic-2010-08-07-70-psi.nml"]
    
    character(len=CL)                               :: path_array(2), input_file
    type(si_velocity), dimension(size(INPUT_FILES)) :: predicted_v_muzzle, actual_v_muzzle, actual_v_muzzle_stdev
    integer                                         :: i_data, rc_predicted, actual_v_muzzle_n(size(INPUT_FILES)), actual_rc
    
    path_array(1) = "examples"
    
    do i_data = 1, size(INPUT_FILES)
        path_array(2) = trim(INPUT_FILES(i_data))
        input_file    = path_join(path_array)
        call run_pneumatic_get_mv(input_file, predicted_v_muzzle(i_data), rc_predicted, actual_v_muzzle(i_data), &
                                    actual_v_muzzle_stdev(i_data), actual_v_muzzle_n(i_data), actual_rc)
        
        call tests%integer_eq(rc_predicted, actual_rc, trim(INPUT_FILES(i_data)) // ", status%rc")
        if (actual_rc == 0) call tests%real_eq(predicted_v_muzzle(i_data)%v%v, actual_v_muzzle(i_data)%v%v, &
                                trim(INPUT_FILES(i_data)) // ", muzzle velocity (validation)", &
                                abs_tol=Z_HALF_ALPHA*actual_v_muzzle_stdev(i_data)%v%v/sqrt(real(actual_v_muzzle_n(i_data), WP)))
    end do
    
    call predicted_v_muzzle_vs_observed("pneumatic-validation", predicted_v_muzzle, actual_v_muzzle, actual_v_muzzle_stdev, &
                                            actual_v_muzzle_n)
end subroutine pneumatic_validation

end program test_validation
