## 0.3.0

- In namelist input files, to improve understanding by users, the unit of volume were changed from m3 to mL the unit of length was changed from m to mm.
- Simulation progress and events are now printed on standard output for the user to better understand how the simulation is progressing.
- Made BlasterSim more robust to floating point error in plunger impact.
- Adaptive time stepping was added.
- Certain error messages now include instructions for the user to fix the error.
- Dead space in springers is now split between the plunger tube and barrel to eliminate numerical issues as the plunger approaches the end of the plunger tube.
- Plunger impact and bounce are handled.
- The documentation now covers all the topics I wanted to cover at a minimum.
- Fuzz testing is now supported via nmlfuzz.
    - BlasterSim was made overall more robust through issues found with fuzz testing.
- Added ability to disable CSV output.
