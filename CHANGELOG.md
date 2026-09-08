## 0.3.0

- Made BlasterSim more robust to floating point error in plunger impact.
- Adaptive time stepping was added.
- Certain error messages now include instructions for the user to fix the error.
- Dead space in springers is now split between the plunger tube and barrel to eliminate numerical issues as the plunger approaches the end of the plunger tube.
- Plunger impact and bounce are handled.
- The documentation now covers all the topics I wanted to cover at a minimum.
- Fuzz testing is now supported via nmlfuzz.
    - BlasterSim was made overall more robust through issues found with fuzz testing.
- Added ability to disable CSV output.
