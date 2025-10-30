- adjacency matrix
- RHS functions for RK4
- `cv_system` type
    - `type(cv_type), allocatable  :: cvs(:)`
    - `type(con_type), allocatable :: cons(:, :)`

***

- gnuplot subroutine for testing and generic output
- water vapor
    - moran_thermodynamics_2008 pp. 666--667: can use $h_g(T)$, $u_g(T)$, or could use Table A-22 (but recognize that Table A-22 has a different datum than the steam tables and can not be used when liquid water is present, which is irrelevant here)
    - Have function to construct `y` given relative humidity.
    - How to handle water condensing out of the air is not clear at the moment, but I suppose I can ignore that to start. All I should have to do to handle water vapor is add the right `gas_type`.
    - This is more complicated than I originally thought: Relative humidity from weather data is given at atmospheric conditions, not pressurized conditions. This means that as the air is compressed, the relative humidity decreases because the vapor pressure increases. So I need to model the compression process. Take ambient air at given relative humidity, and increase the pressure.
- Property test to compare BlasterSim derivatives against numerical derivatives of BlasterSim input.
    - test_fmad.f90: `test_num_deriv`
- Write program to generate namelist reader code, AD `d` indices, and fuzz testing bounds. (Part of FLT but for BlasterSim.)
    - Namelist/ad generator allows for disabling AD for some inputs. Perhaps a second AD namelist where variables can be set to .true. if you want AD. A UQ namelist could set the standard deviation. Combining all of this, you could make a UQ namelist where not defining something in the UQ namelist disables AD for that variable. Read the UQ namelist first, then the value namelist.
