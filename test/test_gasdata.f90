! tests for gas thermochemical data
! Standard: Fortran 2018
! Preprocessor: none
! Author: Ben Trettel (<http://trettel.us/>)
! Project: [BlasterSim](https://github.com/btrettel/blastersim)
! License: [AGPLv3](https://www.gnu.org/licenses/agpl-3.0.en.html)

program test_gasdata

use prec, only: WP
use units
use unittest, only: test_results_type
implicit none

type(test_results_type) :: tests

call tests%start_tests("gasdata.nml")

call test_u_h(tests)
call test_c_p_c_v(tests)

call test_p_v_h2o(tests)

call tests%end_tests()

contains

subroutine test_u_h(tests)
    use gasdata, only: DRY_AIR
    
    type(test_results_type), intent(in out) :: tests

    type(si_temperature)     :: temp
    type(si_specific_energy) :: u, h
    
    call temp%v%init_const(300.0_WP, 0)
    
    u = DRY_AIR%u(temp)
    h = DRY_AIR%h(temp)
    
    ! Data from moran_fundamentals_2008 table A-22.
    call tests%real_eq(u%v%v, DRY_AIR%u_0, "u (gas), 300 K")
    call tests%real_eq(h%v%v, DRY_AIR%h_0, "h (gas), 300 K")
    
    call temp%v%init_const(400.0_WP, 0)
    
    u = DRY_AIR%u(temp)
    h = DRY_AIR%h(temp)
    
    ! Data from moran_fundamentals_2008 table A-22.
    call tests%real_eq(u%v%v, 286.16e3_WP, "u (gas), 400 K", abs_tol=1.0e3_WP)
    call tests%real_eq(h%v%v, 400.98e3_WP, "h (gas), 400 K", abs_tol=1.0e3_WP)
    
    call temp%v%init_const(200.0_WP, 0)
    
    u = DRY_AIR%u(temp)
    h = DRY_AIR%h(temp)
    
    ! Data from moran_fundamentals_2008 table A-22.
    call tests%real_eq(u%v%v, 142.56e3_WP, "u (gas), 200 K", abs_tol=1.0e3_WP)
    call tests%real_eq(h%v%v, 199.97e3_WP, "h (gas), 200 K", abs_tol=1.0e3_WP)
end subroutine test_u_h

subroutine test_c_p_c_v(tests)
    use gasdata, only: DRY_AIR, AR
    
    type(test_results_type), intent(in out) :: tests

    type(si_specific_heat) :: c_p, c_v
    
    c_p = DRY_AIR%c_p(0)
    c_v = DRY_AIR%c_v(0)
    
    ! Data from moran_fundamentals_2008 table A-20.
    call tests%real_eq(c_p%v%v, 1.005e3_WP, "c_p (DRY_AIR), 300 K", abs_tol=1.0_WP)
    call tests%real_eq(c_v%v%v, 0.718e3_WP, "c_v (DRY_AIR), 300 K", abs_tol=1.0_WP)
    
    c_p = AR%c_p(0)
    c_v = AR%c_v(0)
    
    ! See `AR` in cva.f90 for the source of these values.
    ! The NIST data seems to be less consistent with ideal gas equations than moran_fundamentals_2008's data.
    call tests%real_eq(c_p%v%v, 0.52154e3_WP, "c_p (AR), 300 K", abs_tol=3.0_WP)
    call tests%real_eq(c_v%v%v, 0.31239e3_WP, "c_v (AR), 300 K", abs_tol=2.0_WP)
end subroutine test_c_p_c_v

subroutine test_p_v_h2o(tests)
    use convert, only: TEMP_C_TO_K
    use gasdata, only: p_v_h2o
    
    type(test_results_type), intent(in out) :: tests
    
    type(si_temperature) :: temp
    type(si_pressure)    :: p_v
    
    ! Data from moran_fundamentals_2008 table A-2.
    ! This won't match up exactly with the regression equation as it's a different data source.
    
    call temp%v%init_const(TEMP_C_TO_K + 0.01_WP, 0)
    p_v = p_v_h2o(temp)
    call tests%real_eq(p_v%v%v, 0.00611e5_WP, "water vapor pressure at 0.01 C", abs_tol=30.0_WP)
    
    call temp%v%init_const(TEMP_C_TO_K + 20.0_WP, 0)
    p_v = p_v_h2o(temp)
    call tests%real_eq(p_v%v%v, 0.02339e5_WP, "water vapor pressure at 20 C", abs_tol=10.0_WP)
    
    call temp%v%init_const(TEMP_C_TO_K + 50.0_WP, 0)
    p_v = p_v_h2o(temp)
    call tests%real_eq(p_v%v%v, 0.1235e5_WP, "water vapor pressure at 50 C", abs_tol=200.0_WP)
end subroutine test_p_v_h2o

end program test_gasdata
