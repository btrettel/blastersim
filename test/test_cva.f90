! tests for the port module
! Standard: Fortran 2018
! Preprocessor: none
! Author: Ben Trettel (<http://trettel.us/>)
! Project: [BlasterSim](https://github.com/btrettel/blastersim)
! License: [AGPLv3](https://www.gnu.org/licenses/agpl-3.0.en.html)

program test_cva

use unittest, only: test_results_type
implicit none

type(test_results_type) :: tests

call tests%start_tests("cva.nml")

call test_p_eos(tests)
call test_p_f(tests)

call tests%end_tests()

contains

subroutine test_p_eos(tests)
    use units, only: si_mass_density => unit_m30_p10_p00_p00, &
                     si_temperature  => unit_p00_p00_p00_p10, &
                     si_pressure     => unit_m10_p10_m20_p00
    use prec, only: WP
    use cva, only: P_ATM, T_ATM, RHO_ATM, p_eos
    
    type(test_results_type), intent(in out) :: tests

    type(si_mass_density) :: rho
    type(si_temperature)  :: temp
    type(si_pressure)     :: p
    
    call rho%v%init_const(RHO_ATM, 0)
    call temp%v%init_const(T_ATM, 0)
    
    p = p_eos(rho, temp)
    
    call tests%real_eq(p%v%v, P_ATM, "p_eos, atmospheric", abs_tol=1.0_WP)
end subroutine test_p_eos

subroutine test_p_f(tests)
    use units, only: si_pressure => unit_m10_p10_m20_p00
    use cva, only: cv_type
    
    type(test_results_type), intent(in out) :: tests

    type(cv_type)     :: cv
    type(si_pressure) :: p_f0, p_f
    
    call cv%x%v%init_const(0.0_WP, 0)
    call cv%p_fs%v%init_const(0.1e5_WP, 0)
    call cv%p_fd%v%init_const(0.05e5_WP, 0)
    
    call cv%x_dot%v%init_const(-10.0_WP, 0)
    call p_f0%v%init_const(-1.0e5_WP, 0)
    p_f = cv%p_f(p_f0)
    call tests%real_eq(p_f%v%v, cv%p_fs, "p_f, x_dot < 0, p_f0 < 0")
end subroutine test_p_f

end program test_cva
