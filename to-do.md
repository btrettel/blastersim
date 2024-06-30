# to-do

## Before minimum viable product

- Market research: Ask which types of blasters people want from a simulator and which features people want.
    - <https://discord.com/channels/825852031239061545/825852033898774543/1257084202407428158>

## Minimum viable product

- Connected control volumes
    - Two types of connections:
        - pistons
            - The motion of the "piston" will be calculated. The projectile will be a "piston".
            - The overall length, area, and volume of each piston-connected control volume must be the same. And piston-connected control volumes must only refer to each other.
        - valves
    - Connections represented with adjacency matrices.
        - <https://en.wikipedia.org/wiki/Adjacency_matrix>
    - `state_type` contains an array of `dual_control_volume_type`s and the adjacency matrix for connections between control volumes.
    - Start with an input file where you specify the control volumes specifically, one per `control_volume` namelist, then later add modes which automatically create the control volumes.
        - Make `control_volume_type` and use a separate control volume for each side of a piston.
        - Control volume 1: atmosphere
        - Control volume 2: barrel
        - Control volume 3: pressure chamber
- MC data type
    - Why part of the MVP? Because retrofitting the code to add this later will take too much time, and it's simple to add at the beginning.
    - Why Monte Carlo? BlasterSim should be fast, so computation time is not a problem. Uncertainty also is likely to be large, and Monte Carlo would handle that better than FOSM.
    - Write MC module for Fortran.
- Lagrange pressure gradient
    - projectile base pressure is used in force calculation
    - take spatial average of equation of state to relate average pressure to other average thermodynamic functions
    - <https://apps.dtic.mil/sti/citations/ADA222590>
    - also consider that the other side of the chamber may have a non-zero velocity
- ISO valve model with laminar modification
- units.f90
- Documentation
- High order RK scheme
- Simplicity:
    - Interior ballistics only to start.
    - ideal gas
    - constant specific heat
    - one species (air) by default
- Tests
    - mass, momentum, and energy conservation (Using a time integration scheme which preserves these would be nice, but is not necessary for a MVP.)
- Checks:
    - all velocities are much less than speed of sound
- Air gun
    - Why start here? I have validation data I collected back in the day.
        - /home/ben/svn/old/ballistics/data/2010-08-07 muzzle velocity test.gnumeric
- springers
    - Model the air behind the plunger too as it might be pulling a vacuum.
        - <https://discord.com/channels/825852031239061545/825852033898774543/1219837257653944422>
    - Data:
        - <http://nerfhaven.com/forums/topic/21832-experimental-methods-for-determining-and-predicting-blaster-power/?p=307341>
            - <http://www.danielbeaver.net/storage/projects/nerf/SpringerTesting/>
- `make release` compiles for more generic architecture (not `-march=native` in gfortran, etc.) for portability.
- Use GA when simulation is fast. Combine with MC for robust optimization. Also have fitting model coefficients to data, with uncertainties if possible.
- optimization (optimal barrel length, optimal dart mass, optimal 
    - objective functions: maximize muzzle velocity, minimize input energy, maximize number of shots per tank
    - constraints: max KED, max dart mass, max muzzle velocity

### Files

- bsim.f90

## 

- Make Wasm version for a simple GUI. Might need `bind(c)` subset to interface with webpage.
    - I could write the core so that it will compile with LFortran as that seems to be the easiest way to get WebAssembly.
    - This is not a MVP task as I'd need to avoid custom derived type operators as I believe lfortran doesn't support those as of this writing (2024-06-29).
    - <https://fortran-lang.discourse.group/t/flang-wasm-compiler/7589/8>
- Interior Ballistics
- Transitional Ballistics
- Exterior Ballistics
    - Have some model for "fishtailing", which would require 3 DOF (2D trajectory, 1 degree of freedom for angle from velocity vector)
- Terminal Ballistics (Safety)
- The BlasterSim user's guide might be a good replacement for your old ballistics notes, to teach people about Nerf engineering.
- time integration continues after projectile leaves barrel to get the little extra bit of acceleration there
- input file specifies mean, 95% error bar, and units
- input file
    - `interior` namelist
    - `transitional` namelist
    - `exterior` namelist
    - `terminal` namelist
    - other
        - `optimization` namelist
        - `calibration` namelist
- Better validation data: pressure traces
- GGDT confuses people, so have clearer documentation and names:
    - <https://www.spudfiles.com/viewtopic.php?f=26&t=27224>
- Get some ideas about how to make a good UI from Engine Simulator.
    - <https://www.youtube.com/watch?v=ndUyNJEk4ng>
    - <https://github.com/Engine-Simulator/engine-sim-community-edition>
- <http://www.zdspb.com/tech/index.html>
    - Make sure that you can model paintball guns easily too to get that part of the market.
