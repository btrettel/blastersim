### v0.1.0

- Test `m_s` in `d_xdot_d_t`.
- Make CSV file includes interpolated stop time.
- transonic corrections in the barrel (Corner eq. 123)
    - Make this the default but optional if desired for testing.
    - Input validation at first to not use this with RK EOS (if that's added first)
    - Assert `CONST_EOS` or `IDEAL_EOS` with `p < p_c` to satisfy requirements of version in Corner's book.
- Order-of-accuracy test for single control volume
    - adiabatic test without atmospheric pressure and friction
        - Can test the following: `x_dot`, `e`, `e_f`, `p`, `temp`
        - gas `m` is constant
        - Make all input variables have derivatives to test the derivatives too.
        - Need actually constant friction pressure for testing.
    - constant pressure test with atmospheric pressure and friction
        - seigel_theory_1965 eq. 3-2 is a simpler version of this (just subtract atmospheric and friction pressures from the pressure to factor those in)
        - Can test the following: `x`, `x_dot`, `e_f`
- `check_sys`: test each `status%rc` code
    - Add a comment for each check to more easily figure out where each check is done.
    - All state variables don't exceed certain amounts.
        - Also check derivatives in stability checks.
    - `x` > 0
        - Could simply abort if there is piston bounce during shot rather than handle it. Piston bounce is undesirable anyway as it would reduce accuracy of the projectile.
        - Make test case where piston will impact end ($x = 0$) to see what happens and where piston impact detection needs to be added. Inside a RK stage? If the check is in `check_sys`, before getting to that point, an assertion might trigger or BlasterSim might go haywire.
    - Add a length invariant: `x + x_mirror` = const.
    - Add check for whether transonic corrections are needed with the speed of sound.
    - Add check for validity of lumped parameter approximation
        - Does this require the pressure gradient? That's basically how the Biot number works.
        - Some sort of CFL number? $u \, \Delta t / x$? What would $u$ be if there are multiple connections? It doesn't make sense to use this as whether the lumped parameter approximation is valid or not does not depend on $\Delta t$.        - Ask Gemini for ideas on what to do here.
- `write_csv_row`
    - Make methods to get gas kinetic energy and internal energy, use in CSV output
        - Test cases for temperature and internal energy with non-zero gas velocity.
    - Print internal energy and gas kinetic energy in CSV output
    - Write each component of energy and total energy so that someone can create their own energy balance plots.
- Input file reader generator
    - Input validation:
        - Any diameter is too large or too small to not only make sure that it's physically possible, but also that they use the correct units. Perhaps allow the latter to be disabled with `suggestions = .false.`.
        - p < p_c until RK EOS added
- documentation
    - quick start
    - Put all the drawings and equations on paper first.
    - drawings of pneumatic and springer guns with lengths labeled
    - explanation of governing equations
    - list of all inputs
    - FAQ:
        - TODO: Collect comments making these points.
        - Accuracy of the simulation
        - Answer complaints about units: Writing `2.0e-2` is nearly as easy as writing 2.0 cm
            - <https://www.reddit.com/r/nerfhomemades/comments/1p16zi8/interest_in_springer_physics_simulator/npoj6xa/>
            - <https://discord.com/channels/727038380054937610/1172390267890958366/1431436991752572968>
        - Too complex to simulate.
            - <https://www.reddit.com/r/Nerf/comments/1ordwlk/optimising_for_fps_designing_springer_blasters/nnpwewz/>
            - <https://www.reddit.com/r/Nerf/comments/dmgxeb/science_and_math_of_nerf_formulae_and_how_to/f50mx2c/>
            - <https://discord.com/channels/146386512873783296/146680423173455872/1389713797173743788>
        - Simulation isn't worthwhile because it takes less time to figure things out experimentally. Not true from my perspective. In practice, the alternative to simulation has been speculation. People spend a huge amount of time working on things that a simulation could show is not plausible. That's a waste of time. If anything, simulation saves time by reducing the amount of experiments that need to be done. Simulation and experimentation are complementary.
        - How do I calculate optimal barrel length?
            - Comment on some common formulas, give better approximate formulas based on adiabatic process relations, and say how to do it in BlasterSim.
            - <https://www.reddit.com/r/Nerf/comments/1k428am/ideal_barrel_length_math/>
            - <https://www.reddit.com/r/Nerf/comments/186ckcn/formula_for_optimal_barrel_length/>
            - <http://btrettel.nerfers.com/archives/54>
    - API
        - `i_cv_mirror = 0` disables mirror CVs; use for constant volume chambers
    - Write Fortran code to output gas data table to put in documentation.
        - Can I use similar code generation to get other important values from the code?
- Create subroutines in io.f90 to create different types of CVs. Use these in the tests.
- At termination, print:
    - If a success, say so.
    - If a failure, say so.
    - Muzzle velocity including units.
    - efficiency
    - whether under- or over-barreled
- Upload Windows BlasterSim to malware scanner to check. 
- Check that Windows BlasterSim works in Wine to make sure it doesn't require extra libraries.
- Valve opening time, valve poppet model using pressures from CVs
    - Find what you saved on valve opening profiles.
- Time step estimate?
    - <https://www.spudfiles.com/viewtopic.php?p=391877#p391877>: > So I try to pick a time step intelligently. I first make a very crude guesstimate of muzzle energy. That gives me a (crude estimate of) muzzle velocity. I then assume constant acceleration and determine how long it would take a projectile to clear the muzzle.

***

- Plunger head motion bounds (lower and upper) (plunger impact)
    - Lower is not necessarily zero.
    - Use forcing to set x_min and x_max? Make how far the force extends out depend on `dt`. You'd have to track energy lost to this forcing for the energy balance and also the estimate of plunger impact energy.
        - This seems to be the easiest way to do it. It would not require changes to handle zero volume as the piston should always stop before volume decreases to zero, unlike the sudden impact approach.
        - Do I need "impact" force in both directions to reach an equilibrium position? I guess not as overdamped systems exist. This would then allow me to not have any impact force as the piston moves away. Also, the piston could oscillate back and forth even if there is no damping on the back direction as there still would be force on the back direction.
    - Would suddenly stopping the piston when it goes past the stop cause problems with the derivatives? Yes. After impact, $\dot{x} = 0$ and $x = x_\text{min}$. But what do I pick for the derivatives of $\dot{x}$ and $x$? Unlikely. If $t$ is constant, impact might not have occurred if one of the differentiable variables were different. So $\dot{x}$ and $x$ should have specific values that I don't know.
    - Related: Make BlasterSim handle zero volume CVs. For zero volume, use pressure on other side of connection? What should I do if there are multiple connections? Minimum connected pressure?
        - Changing the $p$ equation of state at low volume would introduce a complexity. I'd need to feed in an alternative pressure value to use from the linked control volumes. Picking the minimum pressure of the linked control volumes would satisfy the requirement that nothing can flow to a higher pressure. If I simply required that $x$ be finite, that would lead to unphysical oscillations in pressure as the first time step could flow too much mass into the small CV, leading to high pressure, leading to backflow, etc. What do I pick for the switching point? It would make sense to compare the volumes of the different connected control volumes. If the current control volume is much smaller than a connected control volume, then it could get the pressure from the connected control volume.
        - It would be best to not change the $m$ equation as that would violate mass conservation.
        - Link $x$ and $m$ so that as $x \rightarrow 0$, $m \rightarrow 0$ so $p$ will remain finite?
        - A correction approach might work: If in the next time step, $x$ is anticipated to be zero or below, change $\dot{m}_{i\,\rightarrow\,j}$ so that $m$ goes to zero and change $x$ to zero. I suppose it would be reasonable to scale $\dot{m}_{i\,\rightarrow\,j}$ when there are multiple non-zero entries. Redo the entire time step? If this is done, how can $x$ increase above zero again? Mass entering while $x$ is zero would trigger the mass to be adjusted, which wouldn't work for inflows. I guess I need logic so that for inflows, $x$ is adjusted, and for outflows $\dot{m}_{i\,\rightarrow\,j}$ is adjusted.
    - Test cases for piston impact:
        - Doesn't overshoot when approaching.
        - No effect on piston trajectory when moving away from stopping point.
        - Make sure derivatives are correct after impact. I'm not sure what they should be, though.
- Test `sys_interp`.
- Check that derivatives are correct in special cases where something is set to zero with no derivatives. Check for "TODO: Not sure the derivatives of this should be zero."
- tests for `test_const`, `p_eos`, `temp_cv`, and others for `MIRROR_CV_TYPE`
- Tests for io.f90
- Test `d_x_d_t`, `d_xdot_d_t`, and `d_e_f_d_t` for `CONST_EOS` and `MIRROR_CV_TYPE`.
- Functions to calculate input PE (spring and adiabatic compression) for efficiency calculation.
- function to calculate efficiency
- gnuplot subroutine for testing and generic output
- water vapor
    - moran_thermodynamics_2008 pp. 666--667: can use $h_g(T)$, $u_g(T)$, or could use Table A-22 (but recognize that Table A-22 has a different datum than the steam tables and can not be used when liquid water is present, which is irrelevant here)
    - Have function to construct `y` given relative humidity.
    - How to handle water condensing out of the air is not clear at the moment, but I suppose I can ignore that to start. All I should have to do to handle water vapor is add the right `gas_type`.
    - This is more complicated than I originally thought: Relative humidity from weather data is given at atmospheric conditions, not pressurized conditions. This means that as the air is compressed, the relative humidity decreases because the vapor pressure increases. So I need to model the compression process. Take ambient air at given relative humidity, and increase the pressure.
- Property test to compare BlasterSim derivatives against numerical derivatives of BlasterSim input.
    - test_fmad.f90: `test_num_deriv`
- Readd `smooth_min` assertions including new one from Wikipedia including some extra gap for floating point error
- Add elevation angle
- pressure gradient
    - Is the pressure gradient necessary? The multiple control volumes will provide a pressure gradient of sorts.
    - After adding the pressure gradient, connections will need to know where they are connected.
    - Model pressure drop over long tubes from the friction factor.
- exterior ballistics
- Error messages:
    - When mass or temperature goes negative, suggest that perhaps the effective area is too large.
- Make subroutine to fit `sys%con(:, :)%a_e`, `sys%con(:, :)%b`, `sys%cv(:)%p_fs`, `sys%cv(:)%p_fd` for all `con_types` and `cv_types`.
- Test if correct `sys` is output for `run` (old or new)
- Data to collect from the open literature:
    - sudden contraction data for $A_\text{e}$ and $b$ would be useful for springers
    - Use loss coefficients to estimate effective area?
    - Make regression for loss coefficient considering contraction ratio and entrance radius of curvature?
- Couple BlasterSim with some sort of geometric analysis. Optimal transfer port geometry (considering dead space), optimal notch geometry to minimize weight, etc.
- Have ability to run multiple cycles and terminate when mass gets low in any particular CV.
- If I use a constant pressure/temperature CV to model a HPA or CO2 tank, then I'll still need a way to estimate the real gas internal energy and enthalpy. Going all the way with a better equation of state and thermodynamic properties might not be much more complex. I could make each control volume use a different EOS if I want to avoid iterations associated with a different EOS.
- Arbitrary displacement vs. force curves, using cubic splines for smoothness
- Why is there a small amount of backwards motion of the projectile in `test_conservation`? Is it due to friction?
- Post-projectile-exit analysis to wait until plunger impact, but still interpolate to get muzzle velocity.
- Constraints
    - plunger impact energy
    - to prevent projectile damage: maximum acceleration or maximum projectile pressure difference
- `check_sys`
    - Time step criteria based on flow rate to empty CV? This wouldn't work right if the CV should empty as might be the case for springers. I could still check $\Delta m/m$ for each CV.
    - Try "Lipschitz constant estimate" suggested by Gemini.
    - Message for check_sys error: `CRITICAL_ERROR_MESSAGE = "Please report this input file to the GitHub. https://github.com/btrettel/blastersim/issues"`
- parallel optimization
    - need parallel RNG
- Stopping criteria based on acceleration to find optimal barrel length
- Check entropy conservation.

***

Questions to think about:

- Why doesn't `u_cv` just get `u` from `e`?
- Does plunger impact before the dart exits the barrel cause inaccuracy?
- Why is the derivative of total energy with respect to `x_dot` unstable in the one CV exact solution test?

***

Atean Armory data:

- Post drawing showing interpretation of numbers in efficiency board asking for confirmation.
- Other variables to ask for:
    - plunger mass (everything that moves, not including the spring)
    - spring mass
    - dart mass
    - barrel diameter
    - barrel length
    - dead volume in spring chamber
    - dead volume in barrel
    - whether the barrel can function as a blowgun with these darts

***

Notes:

<https://discord.com/channels/825852031239061545/825852073382772758/1441916424942649354>: > if i'm not mistaken, they're the length in mm of the spring at rest in the blaster (so including precomp if there is any) and the length of the spring when primed [L1 and L2 respectively]
