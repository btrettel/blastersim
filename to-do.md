### v0.1.0

- Add BlasterSim URL to cli.f90.
- `make dist`
- Try exact solution for springers
- Try multiple CV exact solution. One constant pressure chamber, one barrel?
- docs:
    - In BlasterSim code, refer to docs in procedures as appropriate. Exact solution test, governing equations, run time checks.
    - Also refer to code in docs and link to GitHub.
    - Look at Zotero "BibTeX quality report" lines
        - Export to BibLaTeX for even more checks?
    - chktex custom regex option
        - Example: /etc/chktexrc
        - `UserWarn`
        - `UserWarnRegex`
    - BlasterSim output for a working case
    - Put springer and pneumatic governing equations in usage chapter
    - Test building docs on Windows.
    - Thanks appendix
        - Andrew Trettel for macOS binary
    - Development appendix/chapter
    - Add BlasterSim usage in documentation
    - In V&V chapter, discuss running all the tests with `make check` (or jom or NMAKE)
    - Add index to docs.
        - <https://www.overleaf.com/learn/latex/Indices>
        - <https://en.wikibooks.org/wiki/LaTeX/Indexing>
    - Use [Intel Fortran's namelist terminology](https://www.intel.com/content/www/us/en/docs/fortran-compiler/developer-guide-reference/2024-2/namelist.html): namelist group and variable.
    - Example springer input file in documentation (use listings package?)
    - Document friction model to understand what's there are present
    - Quick-start guide for Windows
    - Code documentation
        - <https://fortranwiki.org/fortran/show/Automatic+documentation>
            - Not listed: Sphinx-Fortran
        - LaTeX output: ROBODoc, Doxygen
- BlasterSim output
    - Check acceleration to know if under-barreled or over-barreled.
    - Explain the meaning of the return code if there is an error.
    - efficiency
- Make static friction force actually cancel out properly and not approximately, or at the very least prevent the backwards motion
    - New requirement: $p_\text{f0} \leq p_\text{fe}$ to avoid backwards motion.
- Optimal barrel length mode where the barrel length is not specified and BlasterSim stops where acceleration is zero.
    - It would be important to stop the backwards motion before adding this, otherwise BlasterSim will stop at the wrong time.
- Make going on level deeper (`%v`) optional in geninput when using genunits.
- Add option for `*_stdev` variables to geninput.
- Pneumatic mode.
- Get documentation done before sensitivity analysis.
    - Get sensitivity analysis done before UQ.
    - Get UQ done before optimization so that all optimization is robust for simplicity (no need to have both non-robust and robust optimization set up).

***

- Check that $\dot{x}_0$ derivative is now good with exact solution
- Make friction plot for debugging. Try typical case and also `p_fs = p_fd` to help debug what's going on with that. Why does `p_f` go so much higher than `p_fs`/`p_fd` in that case?
- Test `m_spring` in `d_xdot_d_t` and `m_p_ke`.
- transonic corrections in the barrel (corner_theory_1950 eq. 123)
    - Make this the default but optional if desired for testing.
    - Input validation at first to not use this with RK EOS (if that's added first)
    - Assert `CONST_EOS` or `IDEAL_EOS` with `p < p_c` to satisfy requirements of version in Corner's book.
- pressure gradient
    - Is the pressure gradient necessary? The multiple control volumes will provide a pressure gradient of sorts.
        - <https://www.spudfiles.com/viewtopic.php?f=5&t=26154>
    - After adding the pressure gradient, connections will need to know where they are connected.
    - $p_\text{l}$, $p_\text{r}$
    - Related:
        - Model pressure drop over long tubes from the friction factor.
        - gas kinetic energy
            - Need integrated gas kinetic energy equation.
    - projectile base pressure is used in force calculation
    - take spatial average of equation of state to relate average pressure to other average thermodynamic functions
    - <https://apps.dtic.mil/sti/citations/ADA222590>
    - also consider that the other side of the chamber may have a non-zero velocity (unlike normal Lagrange solution)
    - One issue with the pressure gradient is that the assumed geometry might not be the actual geometry. The pressure gradient solution might be okay for the barrel and plunger tube where the assumed geometry matches the actual geometry, but other flow paths could deviate dramatically.
- Include $\alpha$ in the energy equations.
- Input file reader generator
    - Input validation:
        - Any diameter is too large or too small to not only make sure that it's physically possible, but also that they use the correct units. Perhaps allow the latter to be disabled with `suggestions = .false.`.
        - `p < p_c` until RK EOS added
        - Return error if ambient temperature is too low. Likely they gave the temperature in C or F.
        - `a_e` is not likely larger than the barrel diameter.
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
- Print some derived results initially.
    - Plunger volume
- Determine terminology to use
    - $y_0$: draw or draw length
    - "core" for pressure chamber?
    - GGDT confuses people, so have clearer documentation and names:
        - <https://www.spudfiles.com/viewtopic.php?f=26&t=27224>
- Upload Windows BlasterSim to malware scanner to check. 
- Check that Windows BlasterSim works in Wine to make sure it doesn't require extra libraries.
- Valve opening time, valve poppet model using pressures from CVs
    - Find what you saved on valve opening profiles/valve characteristic curves.
    - Have ability to model valve internals by getting pressures from other control volumes? Then you could model the movement of poppets and whatnot. You'll need some way to handle the "valve profile" or whatever it's called: relationship between poppet location and flow cross-sectional area.
- Time step estimate?
    - <https://www.spudfiles.com/viewtopic.php?p=391877#p391877>: > So I try to pick a time step intelligently. I first make a very crude guesstimate of muzzle energy. That gives me a (crude estimate of) muzzle velocity. I then assume constant acceleration and determine how long it would take a projectile to clear the muzzle.
    - Minimum of multiple time scales?

***

- Use linters including fortitude
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
- Check that derivatives are correct in special cases where something is set to zero with no derivatives. Check for "TODO: Not sure the derivatives of this should be zero."
- tests for `test_const`, `p_eos`, `temp_cv`, and others for `MIRROR_CV_TYPE`
- Tests for io.f90
- Test `d_x_d_t`, `d_x_dot_d_t`, and `d_e_f_d_t` for `CONST_EOS` and `MIRROR_CV_TYPE`.
- Functions to calculate input PE (spring and adiabatic compression) for efficiency calculation.
- function to calculate efficiency
- gnuplot subroutine for testing and generic output
- water vapor
    - moran_thermodynamics_2008 pp. 666--667: can use $h_g(T)$, $u_g(T)$, or could use Table A-22 (but recognize that Table A-22 has a different datum than the steam tables and can not be used when liquid water is present, which is irrelevant here)
    - Have function to construct `y` given relative humidity.
    - How to handle water condensing out of the air is not clear at the moment, but I suppose I can ignore that to start. All I should have to do to handle water vapor is add the right `gas_type`.
    - This is more complicated than I originally thought: Relative humidity from weather data is given at atmospheric conditions, not pressurized conditions. This means that as the air is compressed, the relative humidity decreases because the vapor pressure increases. So I need to model the compression process. Take ambient air at given relative humidity, and increase the pressure.
    - other humidity effects
        - <https://discord.com/channels/825852031239061545/1011329333949907105/1318103337568174121>
        - The inconsistency could be from the dart material or other materials and not the gas.
        - <https://discord.com/channels/1196986370573467668/1247341528561487892/1322273804029923358>: > The weirdest shit happens in the northeast from temperature and humidity fluctuations
- Property test to compare BlasterSim derivatives against numerical derivatives of BlasterSim input.
    - test_fmad.f90: `test_num_deriv`
- Readd `smooth_min` assertions including new one from Wikipedia including some extra gap for floating point error
- Add elevation angle
- Maybe: Add energy loss term for valves?
- Maybe: Quadratic `u` and `h` can be easily solved explicitly
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
- elastic tubing and other non-linear springs
- Post-projectile-exit analysis to wait until plunger impact, but still interpolate to get muzzle velocity.
- `check_sys`
    - Time step criteria based on flow rate to empty CV? This wouldn't work right if the CV should empty as might be the case for springers. I could still check $\Delta m/m$ for each CV.
    - Try "Lipschitz constant estimate" suggested by Gemini.
    - Message for check_sys error: `CRITICAL_ERROR_MESSAGE = "Please report this input file to the GitHub. https://github.com/btrettel/blastersim/issues"`
- optimization
    - common optimization activities: optimal barrel length
    - objective functions: maximize muzzle velocity, minimize input energy, minimize input gas mass (to maximize number of shots per tank in a simpler way as minimum input gas mass assumes all gas in the tank can be extracted)
    - optimize GA parameters to most quickly optimize an example pneumatic blaster
    - parallel optimization
        - need parallel RNG
    - Constraints
        - plunger impact energy
        - to prevent projectile damage: maximum acceleration or maximum projectile pressure difference
        - Make BlasterSim have a generic constraint that can handle any output. `output_constraint` namelist group?
        - Kinetic energy density, muzzle velocity (target, upper limit, lower limit)
        - min/max dart mass
        - dart head decapitation
        - pneumatics might want to use less gas mass per shot
        - spring compression
            - <https://discord.com/channels/727038380054937610/1172390267890958366/1466253503151476877>
- Stopping criteria based on acceleration to find optimal barrel length
- Check entropy conservation.
- Order-of-accuracy tests
    - single control volume constant pressure test with atmospheric pressure and friction
        - seigel_theory_1965 eq. 3-2 is a simpler version of this (just subtract atmospheric and friction pressures from the pressure to factor those in)
        - Can test the following: `x`, `x_dot`, `e_f`
- Fitting model coefficients to data
    - Have a way to put muzzle velocity measurements in the input file as an array so there's only one file with everything.
- Isometric icon and logo for BlasterSim? Check Super Soaker icons you have.
    - Sell BlasterSim stickers to put on blasters that were designed using it? Getting a good logo for this is key.
- Sensitivity analysis for muzzle velocity and optimal barrel length. For barrel length in particular it can be useful to show (if true) that the length doesn't depend much on the plunger mass, spring stiffness.
- Determine when dart heads will be blown off and include that in BlasterSim. Add as a constraint too.
- spring fatigue life, add as a constraint too
- time integration continues after projectile leaves barrel to get the little extra bit of acceleration there
- Make BlasterSim able to handle light gas guns (kinda) by making each side of a piston potentially different diameters. Force is transmitted between two control volumes.
- Allow for modeling of porting with $x$-dependent connections.
- Market research: Ask which types of blasters people want from a simulator and which features people want.
    - <https://discord.com/channels/825852031239061545/825852033898774543/1257084202407428158>
- springers
    - Model the air behind the plunger too as it might be pulling a vacuum.
        - <https://discord.com/channels/825852031239061545/825852033898774543/1219837257653944422>
    - Data:
        - Plunger position: <https://discord.com/channels/727038380054937610/1172390267890958366/1177752285703573504>
        - pressure traces if available
- In `check_sys`, use something with less cancellation error? Pick different points for the derivative calculation to avoid catastrophic cancellation? ash_optimal_1981 eq. 2 won't be the best as it would require 4 function evaluations per iteration.
- Daniel Beaver validation cases:
    - <http://nerfhaven.com/forums/topic/21832-experimental-methods-for-determining-and-predicting-blaster-power/?p=307341>
    - <http://www.danielbeaver.net/storage/projects/nerf/SpringerTesting/>
    - Doesn't provide enough data to make a complete test case.
- Test `sys_interp`.
    - Comparison with `x_stop`. Is a test needed given the assertion?
    - Test with something that can be solved exactly by RK4. Then `x_dot` can be known exactly.
- `check_sys`
    - test each `status%rc` code
    - All state variables don't exceed certain amounts.
        - Also check derivatives in stability checks.
    - `x` > 0
        - Could simply abort if there is piston bounce during shot rather than handle it. Piston bounce is undesirable anyway as it would reduce accuracy of the projectile.
        - Make test case where piston will impact end ($x = 0$) to see what happens and where piston impact detection needs to be added. Inside a RK stage? If the check is in `check_sys`, before getting to that point, an assertion might trigger or BlasterSim might go haywire.
    - Add check for whether transonic corrections are needed with the speed of sound.
    - Add check for validity of lumped parameter approximation
        - Does this require the pressure gradient? That's basically how the Biot number works.
        - Ask Gemini for ideas on what to do here.
- `write_csv_row`
    - Add Python code to test reading the output file.
        - `rc` column is -1 until changing to something different on last row.
        - Same number of columns for each row.
    - Make methods to get gas kinetic energy and internal energy, use in CSV output
        - Test cases for temperature and internal energy with non-zero gas velocity.
    - Print internal energy and gas kinetic energy in CSV output
- Check /home/ben/svn/old/ballistics/text/ for more ideas.
- Add spring $k$ calculator
    - Then you could optimize the length of the spring to cut to.
- Make `a_e` and `b` different for reverse flow in the springer case.
- Allow for negative precompression. This would require changing the spring force law. The spring typically would not apply a restoring force (opposite direction) as it's not firmly attached to the plunger.
- Try hevea, tex4ht, and pandoc for HTML version of the docs.
    - [HEVEA](https://hevea.inria.fr/)
    - [latex2html](https://www.latex2html.org/)
    - [lwarp](https://ctan.org/pkg/lwarp?lang=en)
        - <https://github.com/DominikPeters/pgf-tikz-html-manual>
    - pandoc
        - <https://www.danwjoyce.com/data-blog/2018/2/20/latex-to-html-via-pandoc>
        - <https://tex.stackexchange.com/questions/431719/how-to-use-pandoc-to-derive-output-from-latex-and-tikz-to-a-docx-file>
    - [TeX4ht](https://www.tug.org/tex4ht/)
- Test `d_m_k_d_t` with multiple gases.
- promotion
    - Post on:
        - r/nerf
        - r/nerfhomemades
        - r/hpanerf
        - Discord
            - Atean Armory #the-workbench
            - Kelly Industries #springer-science
            - Roboman Automation #showoff
            - Nerfcord #showoff
            - Private server
        - SpudFiles
        - Title: BlasterSim: simulate your blasters, calculate optimal barrel length, and more
    - Distribution:
        - trettel.us
        - <https://blasterdownloads.com/are-you-foam-dart-blaster-designer/>

***

Questions to think about:

- Why doesn't `u_cv` just get `u` from `e`?
- Does plunger impact before the dart exits the barrel cause inaccuracy?
- Problems maybe with the friction model:
    - Why is the derivative of total energy with respect to initial `x_dot` unstable in `test_one_cv`?
        - Does using constant `p_f` solve the `x_dot` derivative instability problem? `p_f` would change dramatically as `x_dot` changes for `x_dot == 0`.
    - Why is there a small amount of backwards motion of the projectile in `test_conservation`? Is it due to friction?
- How can I handle projectiles on the outside of the barrel?
    - Initial barrel volume is the dead volume. CSA is the area the axial force is acting on.
- How can I handle having a rod inside the barrel?
    - If I use CSA instead of diameter, use the actual CSA and not the (outer) diameter.
- How can I handle dead volume in the dart?
    - This just increases the dead volume for the barrel.
    - Could have normal dead volume and a separate dead volume for the projectile.

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

corner_theory_1950 p. 58: > so the $R$ associated with the $A$ is in heat units, namely, calories per mole per degree.

***

Outline of planned advanced mode input file:

    &bsim
    stdout = 
    /

    &barrel
    ! cv_id = 1 ! implicit
    projectile_mass = 1.0e-3 ! kg
    /

    &pressure_chamber
    ! cv_id = 2 ! implicit
    p_initial = 5.0e5 ! Pa
    /

    &connection

    effective_area = 1.0e-3 ! m2
    b = 0.5
    /

    &external
    ! cv_id = 1 ! implicit, to get projectile mass
    C_d = 0.67
    /

inspiration for UI: SPICE, <https://en.wikipedia.org/wiki/Netlist>

Simplified modes:

- `springer` namelist
- `pneumatic` namelist

***

GUI ideas:

- Try WebAssembly to see what the limitations are. Having something that people can use without much effort and on their phones is important.
    - Might need `bind(c)` subset to interface with webpage.
    - I could write the core so that it will compile with LFortran as that seems to be the easiest way to get WebAssembly.
    - LFortran
        - I'd need to avoid custom derived type operators as I believe lfortran doesn't support those as of this writing (2024-06-29). Perhaps `interface` operators are okay? Then I could use genunits.
        - <https://lfortran.org/blog/2024/05/fortran-on-web-using-lfortran/>
            - <https://fortran-lang.discourse.group/t/fortran-on-web-using-lfortran/7957>
            - <https://github.com/lfortran/mnist-classifier-blas-wasm/>
        - <https://fortran-lang.discourse.group/t/flang-wasm-compiler/7589/8>
        - <https://github.com/lfortran/Fortran-On-Web>
    - Flang
        - <https://gws.phd/posts/fortran_wasm/>
            - <https://news.ycombinator.com/item?id=39944275>
            - <https://hackaday.com/2024/04/08/fortran-and-webassembly-bringing-zippy-linear-algebra-to-nodejs-browsers/>
        - <https://fortran-lang.discourse.group/t/flang-wasm-compiler/7589>
        - <https://niconiconi.neocities.org/tech-notes/fortran-in-webassembly-and-field-solver/>
- Look into running with a GUI using `iso_c_binding` (or otherwise). This should be started with early as it likely will limit the Fortran code in some way. Or perhaps not, if some sort of conversion subroutines will be needed.
    - <https://docs.python.org/3/library/tk.html>
- Order inputs by sensitivity. Do a sensitivity study first to know the order.
- Get some ideas about how to make a good UI from Engine Simulator.
    - <https://www.youtube.com/watch?v=ndUyNJEk4ng>
    - <https://github.com/Engine-Simulator/engine-sim-community-edition>

***

- Make official website similar to some of the old 90s-style software websites you have links to saved. LaTeX documentation is available on the website per LaTeXML.
    - Make old Nerf ballistics PDF file redirect to BlasterSim's docs when BlasterSim is released.
    - Examples:
        - <https://flatassembler.net/>
            - <https://board.flatassembler.net/index.php>
        - <https://systemd.io/>
        - <https://zeta.asie.pl/>
        - <https://webalizer.net/>
        - <https://ice-wm.org/>
            - <https://web.archive.org/web/20080302175539/http://www.icewm.org/>
            - <https://web.archive.org/web/20020325091254/http://icewm.sourceforge.net/>
            - <https://web.archive.org/web/20020328083952/http://www.icewm.com/>
        - <https://unhaut.epizy.com/psxsdk/>
        - <https://chenthread.asie.pl/fromage/>
        - <https://plasma-gate.weizmann.ac.il/Grace/>
        - <https://www.nongnu.org/chktex/>
        - <http://aspell.net/>
        - Dark mode?
            - <https://speeddemosarchive.com/>
            - <https://hardforum.com/>
- BlasterSim tutorial on YouTube, with slides showing each control volume and how they are connected to what's in the input file.
    - This could also be good marketing.
