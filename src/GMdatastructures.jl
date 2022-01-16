# types ------------------------------------------------------------------------

# type corresponding to a solution
mutable struct tSolution{T}
    x :: Vector{T}                # vector variables x (1..n)
    y :: Vector{T}                # vector outcomes  y (1..p)
end

# type corresponding to a point generator
mutable struct tGenerateur
    sRel :: tSolution{Float64}    # initial relaxed solution
    sInt :: tSolution{Int64}      # integer solution
    sPrj :: tSolution{Float64}    # projected solution
    sFea :: Bool                  # indicate if sInt is feasible or not
end

# type of a point (x,y) in the objective space
mutable struct tPoint
    x :: Float64
    y :: Float64
end

# Global variables -------------------------------------------------------------

# listes de points pour les affichages graphiques
xLf1  = (Float64)[]; yLf1  = (Float64)[] # liste des points (x,y) relaches
xLf2  = (Float64)[]; yLf2  = (Float64)[] # liste des points (x,y) relaches
xL    = (Float64)[]; yL    = (Float64)[] # liste des points (x,y) relaches
XInt  = (Int64)[];   YInt  = (Int64)[]   # liste des points (x,y) entiers
XProj = (Float64)[]; YProj = (Float64)[] # liste des points (x,y) projetes
XFeas = (Int64)[];   YFeas = (Int64)[]   # liste des points (x,y) admissibles
XPert = (Int64)[];   YPert = (Int64)[]   # liste des points (x,y) perturbes


# ==============================================================================
# Initialisation structure donnees contenant tous les generateurs

function allocateDatastructure(nbgen::Int64, nbvar::Int64, nbobj::Int64)

    verbose ? println("\n  → Allocation memoire pour ",nbgen," generateurs\n") : nothing

    vg = Vector{tGenerateur}(undef, nbgen)
    for k = 1:nbgen
        vg[k] = tGenerateur( tSolution{Float64}(zeros(Float64,nbvar),zeros(Float64,nbobj)),
                              tSolution{Int64}(zeros(Int64,nbvar),zeros(Int64,nbobj)),
                              tSolution{Float64}(zeros(Float64,nbvar),zeros(Float64,nbobj)),
                              false
                            )
    end
    return vg
end