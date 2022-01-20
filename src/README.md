**Notice:**

The source codes provided in this folder are versions under development and obsolete prototypes. 

In any case they are not stable and are not ready for conducting comparative numerical experiments. 

Currently the source codes are documented mainly in french.

***

**Source codes:**

- `GMmain.jl` version under development (January 2022)
- others | `testCones.jl` code for displaying the cones
- others | `testDirections`code for computing the directions
- obsolete | vRefactoring2021 | `gravityMachineV4.jl` starting point of the refactoring
- obsolete | vConferences2021 | `SPA_GM.jl` very first version used for talks presented in 2021

***

**Configuration required:**

Packages Julia to add:

- JuMP
- GLPK
- PyPlot

MIP Solver required:

- GLPK (to install it with Homebrew on macOS: [https://brew.sh](https://brew.sh) and [https://formulae.brew.sh/formula/glpk](https://formulae.brew.sh/formula/glpk) )

***

**Running GravityMachine:**

- get a local copy of this repository
- open a Julia REPL
- set the current directory to `src`
- type into the REPL `include("GMmain.jl")`

