**Notice:**

The source codes provided in this folder are versions under development and obsolete prototypes. 

In any case they are not stable and are not ready for conducting comparative numerical experiments. 

Currently the source codes are documented mainly in french.

***

**Source codes:**

- `GM.jl` version under development (January 2022)
- `testcones.jl` code for displaying the cones
- `gravityMachineV3.jl` starting point of the refactoring; obsolete
- `SPA_GM.jl` version of July 2021 (used for talks presented in 2021); obsolete 

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
- type into the REPL `include("GM.jl")`

