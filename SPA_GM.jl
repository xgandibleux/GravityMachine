using JuMP, GLPK, PyPlot, Printf , vOptGeneric
include("crea_read.jl") # permet de créer les fichiers : données bi-obj, solution Y_N, solutions relachées
include("xp_num.jl")
# https://github.com/fuofuo-gg/TER

#const GUROBI_ENV = Gurobi.Env()

# ==============================================================================
#= Retourne un booléen indiquant si un point se trouve dans un secteur défini dans
  le sens de rotation trigonométrique (repère X de gauche à droite, Y du haut vers
  le bas).
  https://www.stashofcode.fr/presence-dun-point-dans-un-secteur-angulaire/#more-328
  M    Point dont la position est à tester (point resultant a tester)
  O    Point sommet du secteur (point generateur)
  A    Point de départ du secteur (point adjacent inferieur)
  B    Point d'arrivée du secteur (point adjacent superieur)
  sortie : Booléen indiquant si le point est dans le secteur ou non.

  Exemple :

  B=point(2.0,1.0)
  O=point(2.5,2.5)
  A=point(5.0,5.0)

  M=point(5.0,4.0)
  InSector(M, O, A, B)
=#

mutable struct point
    x::Float64
    y::Float64
end

function InSector(M, O, A, B)
    cpAB = (A.x - O.x) * (B.y - O.y) - (A.y - O.y) * (B.x - O.x)
    cpAM = (A.x - O.x) * (M.y - O.y) - (A.y - O.y) * (M.x - O.x)
    cpBM = (B.x - O.x) * (M.y - O.y) - (B.y - O.y) * (M.x - O.x)

    if (cpAB > 0)
        if ((cpAM >= 0) && (cpBM <= 0))
            return true
        else
            return false
        end
    else
        if (!((cpAM < 0) && (cpBM > 0)))
            return true
        else
            return false
        end
    end
end

# ----------------------------------------------------------------
# Model solutions optimales SPA bi obj
function SPA_vopt(nbvar::Int, nbcontraintes::Int, L::Array{Int,2}, c1::Array{Int,1}, c2::Array{Int,1})
    model = vModel(with_optimizer(GLPK.Optimizer))
    JuMP.set_silent(model)

    @variable(model, x[1:nbvar], Bin)

    @constraint(model, [i=1:nbcontraintes],(sum((x[j]*L[i,j]) for j in 1:nbvar)) == 1)

    @addobjective(model, Min, sum(c1[i]*x[i] for i in 1:nbvar))
    @addobjective(model, Min, sum(c2[i]*x[i] for i in 1:nbvar))


    vSolve(model, method=:epsilon, step = 0.5)

    Y_N = getY_N(model)

    return Y_N
end

# ==============================================================================
# Model solution relaxation linéaire sur un objectif
function relaxLinXG(nbvar::Int, nbcontraintes::Int, L::Array{Int,2}, c1::Array{Int,1}, c2::Array{Int,1}, vobj, delta, obj)

    model = Model(with_optimizer(GLPK.Optimizer))
#    model = Model(with_optimizer(Gurobi.Optimizer, GUROBI_ENV))
#    set_optimizer_attributes(model, "OutputFlag" => 0)

    @variable(model, 0.0 <= x[1:nbvar] <= 1.0 )
    @constraint(model, [i=1:nbcontraintes],(sum((x[j]*L[i,j]) for j in 1:nbvar)) == 1)
    if obj == 1
        @objective(model, Min, sum((c1[i])*x[i] for i in 1:nbvar))
        @constraint(model, sum((c2[i])*x[i] for i in 1:nbvar) <= vobj-delta)
    else
        @objective(model, Min, sum((c2[i])*x[i] for i in 1:nbvar))
        @constraint(model, sum((c1[i])*x[i] for i in 1:nbvar) <= vobj-delta)
    end
    optimize!(model)
    return objective_value(model), value.(x)
end

# ==============================================================================
# Model solution relaxation linéaire sur un objectif vertical
function relaxLinXG4(nbvar::Int, nbcontraintes::Int, L::Array{Int,2}, c1::Array{Int,1}, c2::Array{Int,1}, epsilon, obj)

    model = Model(with_optimizer(GLPK.Optimizer))
#    model = Model(with_optimizer(Gurobi.Optimizer, GUROBI_ENV))
#    set_optimizer_attributes(model, "OutputFlag" => 0)

    @variable(model, 0.0 <= x[1:nbvar] <= 1.0 )
    @constraint(model, [i=1:nbcontraintes],(sum((x[j]*L[i,j]) for j in 1:nbvar)) == 1)
    if obj == 1
        @objective(model, Min, sum((c1[i])*x[i] for i in 1:nbvar))
        @constraint(model, sum((c2[i])*x[i] for i in 1:nbvar) <= epsilon)
    else
        @objective(model, Min, sum((c2[i])*x[i] for i in 1:nbvar))
        @constraint(model, sum((c1[i])*x[i] for i in 1:nbvar) <= epsilon)
    end
    optimize!(model)
    return objective_value(model), value.(x)
end

# ==============================================================================
# Elabore 2 ensembles d'indices selon que xTilde[i] vaut 0 ou 1
function splitXG(xTilde)

   indices0 = (Int64)[]
   indices1 = (Int64)[]

   for i=1:length(xTilde)
       if xTilde[i] == 0
           push!(indices0,i)
       else
           push!(indices1,i)
       end
    end

   return indices0, indices1
end

#xTilde = [0,0,1,1,1,0,1]

# ==============================================================================
# Projete xTilde sur le polyedre X
function Δ2(L::Array{Int,2}, xTilde::Array{Int,1})

    nbcontraintes = size(L,1)
    nbvar = size(L,2)
    idxTilde0, idxTilde1 = splitXG(xTilde)

    proj=Model(GLPK.Optimizer)
    @variable(proj, 0.0 <= x[1:length(xTilde)] <= 1.0 )
    @objective(proj, Min, sum(x[i] for i in idxTilde0) + sum((1-x[i]) for i in idxTilde1) )
    @constraint(proj, [i=1:nbcontraintes],(sum((x[j]*L[i,j]) for j in 1:nbvar)) == 1)
    optimize!(proj)
    return objective_value(proj), value.(x)
end

# ==============================================================================
# test si une solution est admissible en verifiant si sa relaxation lineaire
# conduit a une solution entiere
function estAdmissible(x)
    admissible = true
    i=1
    while admissible && i<=length(x)
        if x[i] > 0.001 && x[i] < 0.999
            return false
        end
        i+=1
    end
    return admissible
end

# ==============================================================================
# calcule la performance d'une solution sur les 2 objectifs
function evaluerSolution(x, c1, c2)
    z1 = 0.0; z2 = 0.0
    for i in 1:length(x)
        z1 += x[i] * c1[i]
        z2 += x[i] * c2[i]
    end
    return round(z1, digits=2), round(z2, digits=2)
end

# ==============================================================================
# The gravity machine (Man of Steel) -> to terraform the world
#   e-contrainte avec z1 comme fonction a minimiser
#   arrondi qui evite de conduire vers un point domine par le point issu de la RL

# structure corresponding to one point generator
mutable struct tGenerateur
    xRlx :: Vector{Float64} # relaxed solution x
    yRlx :: Vector{Float64} # point corresponding to the relaxed solution x
    #
    xInt :: Vector{Int64}   # integer solution x
    yInt :: Vector{Int64}   # point corresponding to the integer solution x
    #
    xPrj :: Vector{Float64} # projected solution x
    yPrj :: Vector{Float64} # point corresponding to the projected solution x
    #
    xFea :: Bool            # indicate if a solution x is feasible or not
end

mutable struct tEnsBornant
    xEBD :: Vector{Float64} # solution x member of the dual bound set
    yEBD :: Vector{Float64} # performance y member of the dual bound set
end

mutable struct tPoint
    x :: Vector{Float64} # vector solution x (1..n)
    y :: Vector{Float64} # vector performances y (1..p)
end

# ============================================================================
#= On effectue une vérification si tous les fichiers de donnees ont bien ete
 crees avant d'executer gravity machine
=#

function verification(s::String)
    if !isfile("SPA_Databio/bio"*s*".txt")
        println("Création du fichier bi-objectif : ")
        dirname(pwd())
        SPA_creation("$s")
        println("Fichier bi-objectif cree : ")
    end
    nbvar, nbcontraintes, L, c1, c2 = SPA_loaddemandebi("SPA_Databio/bio"*s*".txt")
    if !isfile("Archive/Y_N_"*s*".txt")
        println("Création du fichier Y_N_ : ")
        start = time()
        Y_N = SPA_vopt(nbvar, nbcontraintes, L, c1, c2) #solver vopt
        temps = time() - start
        archiver_Y_N(s, Y_N, temps)
        println("Fichier Y_N_ créé : ")
    end

end

# ----------------------------------------------------------------
# Nettoyage de la valeur de vg[j].xRlx et calcul du point bi-objectif
function calcNoInteger(vg::Array{tGenerateur,1}, nbvar::Int, j::Int)
    for i in 1:nbvar
        if round(vg[j].xRlx[i], digits=3) == 0.0
            vg[j].xRlx[i] = 0.0
            vg[j].xInt[i] = 0
        elseif round(vg[j].xRlx[i], digits=3) == 1.0
            vg[j].xRlx[i] = 1.0
            vg[j].xInt[i] = 1
        else
            vg[j].xRlx[i] = round(vg[j].xRlx[i], digits=3)
            vg[j].xInt[i] = -1 # valeur sentinelle quand non entier
        end
    end
    return vg
end

# ----------------------------------------------------------------
# Selectionne les points pour constituer le cone d'interet centre sur le point generateur O
function selectPoints(j::Int, tailleSampling::Int, eb::Array{tEnsBornant,1})
    O=point(eb[j].yEBD[1], eb[j].yEBD[2])

    #choix de B
    flag = j
    while flag>1 && (eb[flag-1].yEBD[1] > eb[j].yEBD[1])
        flag -= 1
    end
    if flag == 1
        B=point(eb[1].yEBD[1], eb[1].yEBD[2]+1.0)
    else
        B=point(eb[flag-1].yEBD[1], eb[flag-1].yEBD[2])
    end

    #choix de A
    flag = j
    while flag<tailleSampling && (eb[flag+1].yEBD[1] < eb[j].yEBD[1])
        flag += 1
    end
    if flag == tailleSampling
        A=point(eb[end].yEBD[1]+1.0, eb[end].yEBD[2])
    else
        A=point(eb[flag+1].yEBD[1], eb[flag+1].yEBD[2])
    end
    return B, O, A
end

# ----------------------------------------------------------------
# Operations pour arrondir les valeurs non-entieres d'un générateur
function roundOneCone(j::Int, nbvar::Int, vg::Array{tGenerateur,1}, c1::Array{Int64,1}, c2::Array{Int64,1}, O::point, A::point, B::point)
    vg[j].yInt[1] = 0
    vg[j].yInt[2] = 0
    z1 = vg[j].yRlx[1]
    z2 = vg[j].yRlx[2]

    nbVarNonEntiere=0

    for i in 1:nbvar
        if vg[j].xInt[i] == -1
            # la variable i est non entiere
            nbVarNonEntiere += 1
            M = point( z1 - (vg[j].xRlx[i] * c1[i]) , z2 - (vg[j].xRlx[i] * c2[i]))
            #print(M)
            if !InSector(M, O, A, B)
                # le point M obtenu est hors du cone => variable i fixee a 1
                vg[j].xInt[i] = 1
                z1 = z1 - (vg[j].xRlx[i] * c1[i]) + c1[i]
                z2 = z2 - (vg[j].xRlx[i] * c2[i]) + c2[i]
            else
                # le point M obtenu est dans le cone => variable i fixee a 0
                vg[j].xInt[i] = 0
                z1 = z1 - (vg[j].xRlx[i] * c1[i])
                z2 = z2 - (vg[j].xRlx[i] * c2[i])
            end
        end
        vg[j].yInt[1] += vg[j].xInt[i] * c1[i]
        vg[j].yInt[2] += vg[j].xInt[i] * c2[i]
    end

    return vg, nbVarNonEntiere
end

# ----------------------------------------------------------------
# Operations pour arrondir les valeurs non-entieres de la projection
function roundTwoCones(j::Int, nbvar::Int, vg::Array{tGenerateur,1}, c1::Array{Int64,1}, c2::Array{Int64,1}, O::point, A::point, B::point, OO::point, param_round2::Int)
    vg[j].yInt[1] = 0
    vg[j].yInt[2] = 0
    Xapx=(Float64)[]; Yapx=(Float64)[]
    z1 = vg[j].yRlx[1]
    z2 = vg[j].yRlx[2]
    nbVarNonEntiere=0

    if param_round2 == 1
        for i in 1:nbvar
            if vg[j].xInt[i] == -1
                M = point( z1 - vg[j].xRlx[i] * c1[i] , z2 - vg[j].xRlx[i] * c2[i])
                push!(Xapx,M.x)
                push!(Yapx,M.y)
                # la variable i est non entiere
                nbVarNonEntiere += 1
                if !InSector(M, O, A, B)
                    # le point M obtenu est hors du cone => variable i fixee a 1
                    vg[j].xInt[i] = 1
                    z1 = z1 - (vg[j].xRlx[i] * c1[i]) + c1[i]
                    z2 = z2 - (vg[j].xRlx[i] * c2[i]) + c2[i]
                else
                    if !InSector(M, OO, B, A)
                        MM = point( z1 - vg[j].xRlx[i] * c1[i] + c1[i] , z2 - vg[j].xRlx[i] * c2[i] + c2[i])
                        if !InSector(MM, OO, B, A)
                            # le point MM obtenu est hors du cone => choix de 0 ou 1 avec la plus courte distance à O et OO
                            distm = sqrt((O.x-M.x)^2+(O.y-M.y)^2)+sqrt((OO.x-M.x)^2+(OO.y-M.y)^2)
                            dismm = sqrt((O.x-MM.x)^2+(O.y-MM.y)^2)+sqrt((OO.x-MM.x)^2+(OO.y-MM.y)^2)
                            if distm < dismm
                                vg[j].xInt[i] = 0
                                z1 = z1 - (vg[j].xRlx[i] * c1[i])
                                z2 = z2 - (vg[j].xRlx[i] * c2[i])
                            else
                                vg[j].xInt[i] = 1
                                z1 = z1 - (vg[j].xRlx[i] * c1[i]) + c1[i]
                                z2 = z2 - (vg[j].xRlx[i] * c2[i]) + c2[i]
                            end
                        else
                            vg[j].xInt[i] = 1
                            z1 = z1 - (vg[j].xRlx[i] * c1[i]) + c1[i]
                            z2 = z2 - (vg[j].xRlx[i] * c2[i]) + c2[i]
                        end
                        MM = point( z1 - (vg[j].xRlx[i] * c1[i]) + c1[i] , z2 - (vg[j].xRlx[i] * c2[i]) + c2[i])
                    else
                        # le point M obtenu est dans du cone => variable i fixee a 0
                        vg[j].xInt[i] = 0
                        #print(" ",i)
                        z1 = z1 - (vg[j].xRlx[i] * c1[i])
                        z2 = z2 - (vg[j].xRlx[i] * c2[i])
                    end
                end
            end
            vg[j].yInt[1] += vg[j].xInt[i] * c1[i]
            vg[j].yInt[2] += vg[j].xInt[i] * c2[i]
        end
    elseif param_round2 == 2
        i = 1
        while i <= nbvar
            if vg[j].xInt[i] == -1
                nbVarNonEntiere += 1
                k = i + 1
                while k <= nbvar
                    if vg[j].xInt[k] == -1
                        nbVarNonEntiere += 1
                        MOO = point( z1 - vg[j].xRlx[i] * c1[i] - vg[j].xRlx[k] * c1[k] , z2 - vg[j].xRlx[i] * c2[i] - vg[j].xRlx[k] * c2[k])
                        MOI = point( z1 - vg[j].xRlx[i] * c1[i] - vg[j].xRlx[k] * c1[k] + c1[k] , z2 - vg[j].xRlx[i] * c2[i] - vg[j].xRlx[k] * c2[k] + c2[k])
                        MIO = point( z1 - vg[j].xRlx[i] * c1[i] + c1[i] - vg[j].xRlx[k] * c1[k] , z2 - vg[j].xRlx[i] * c2[i] + c2[i] - vg[j].xRlx[k] * c2[k])
                        MII = point( z1 - vg[j].xRlx[i] * c1[i] + c1[i] - vg[j].xRlx[k] * c1[k] + c1[k] , z2 - vg[j].xRlx[i] * c2[i] + c2[i] - vg[j].xRlx[k] * c2[k] + c2[k])
                        if InSector(MOO, O, A, B) && InSector(MOO, OO, B, A)
                            vg[j].xInt[i] = 0
                            z1 = z1 - (vg[j].xRlx[i] * c1[i])
                            z2 = z2 - (vg[j].xRlx[i] * c2[i])
                            vg[j].xInt[k] = 0
                            z1 = z1 - (vg[j].xRlx[k] * c1[k])
                            z2 = z2 - (vg[j].xRlx[k] * c2[k])
                        elseif InSector(MOI, O, A, B) && InSector(MOI, OO, B, A)
                            vg[j].xInt[i] = 0
                            z1 = z1 - (vg[j].xRlx[i] * c1[i])
                            z2 = z2 - (vg[j].xRlx[i] * c2[i])
                            vg[j].xInt[k] = 1
                            z1 = z1 - (vg[j].xRlx[k] * c1[k]) + c1[k]
                            z2 = z2 - (vg[j].xRlx[k] * c2[k]) + c2[k]
                        elseif InSector(MIO, O, A, B) && InSector(MIO, OO, B, A)
                            vg[j].xInt[i] = 1
                            z1 = z1 - (vg[j].xRlx[i] * c1[i]) + c1[i]
                            z2 = z2 - (vg[j].xRlx[i] * c2[i]) + c2[i]
                            vg[j].xInt[k] = 0
                            z1 = z1 - (vg[j].xRlx[k] * c1[k])
                            z2 = z2 - (vg[j].xRlx[k] * c2[k])
                        elseif InSector(MII, O, A, B) && InSector(MII, OO, B, A)
                            vg[j].xInt[i] = 1
                            z1 = z1 - (vg[j].xRlx[i] * c1[i]) + c1[i]
                            z2 = z2 - (vg[j].xRlx[i] * c2[i]) + c2[i]
                            vg[j].xInt[k] = 1
                            z1 = z1 - (vg[j].xRlx[k] * c1[k]) + c1[k]
                            z2 = z2 - (vg[j].xRlx[k] * c2[k]) + c2[k]
                        else
                            distMOO = sqrt((O.x-MOO.x)^2+(O.y-MOO.y)^2)+sqrt((OO.x-MOO.x)^2+(OO.y-MOO.y)^2)
                            distMOI = sqrt((O.x-MOI.x)^2+(O.y-MOI.y)^2)+sqrt((OO.x-MOI.x)^2+(OO.y-MOI.y)^2)
                            distMIO = sqrt((O.x-MIO.x)^2+(O.y-MIO.y)^2)+sqrt((OO.x-MIO.x)^2+(OO.y-MIO.y)^2)
                            distMII = sqrt((O.x-MII.x)^2+(O.y-MII.y)^2)+sqrt((OO.x-MII.x)^2+(OO.y-MII.y)^2)
                            if InSector(MOO, O, A, B) && distMOO <= distMOI && distMOO <= distMIO && distMOO <= distMII
                                vg[j].xInt[i] = 0
                                z1 = z1 - (vg[j].xRlx[i] * c1[i])
                                z2 = z2 - (vg[j].xRlx[i] * c2[i])
                                vg[j].xInt[k] = 0
                                z1 = z1 - (vg[j].xRlx[k] * c1[k])
                                z2 = z2 - (vg[j].xRlx[k] * c2[k])
                            elseif InSector(MOI, O, A, B) && distMOI <= distMIO && distMOI <= distMII
                                vg[j].xInt[i] = 0
                                z1 = z1 - (vg[j].xRlx[i] * c1[i])
                                z2 = z2 - (vg[j].xRlx[i] * c2[i])
                                vg[j].xInt[k] = 1
                                z1 = z1 - (vg[j].xRlx[k] * c1[k]) + c1[k]
                                z2 = z2 - (vg[j].xRlx[k] * c2[k]) + c2[k]
                            elseif InSector(MIO, O, A, B) && distMIO <= distMII
                                vg[j].xInt[i] = 1
                                z1 = z1 - (vg[j].xRlx[i] * c1[i]) + c1[i]
                                z2 = z2 - (vg[j].xRlx[i] * c2[i]) + c2[i]
                                vg[j].xInt[k] = 0
                                z1 = z1 - (vg[j].xRlx[k] * c1[k])
                                z2 = z2 - (vg[j].xRlx[k] * c2[k])
                            else
                                vg[j].xInt[i] = 1
                                z1 = z1 - (vg[j].xRlx[i] * c1[i]) + c1[i]
                                z2 = z2 - (vg[j].xRlx[i] * c2[i]) + c2[i]
                                vg[j].xInt[k] = 1
                                z1 = z1 - (vg[j].xRlx[k] * c1[k]) + c1[k]
                                z2 = z2 - (vg[j].xRlx[k] * c2[k]) + c2[k]
                            end
                        end
                        k += nbvar + 1
                    end
                    k += 1
                end
                if k == nbvar + 1
                    nbVarNonEntiere += 1
                    MO = point( z1 - vg[j].xRlx[i] * c1[i] , z2 - vg[j].xRlx[i] * c2[i])
                    MI = point( z1 - vg[j].xRlx[i] * c1[i] + c1[i] , z2 - vg[j].xRlx[i] * c2[i] + c2[i])
                    if InSector(MO, O, A, B) && InSector(MO, OO, B, A)
                        vg[j].xInt[i] = 0
                        z1 = z1 - (vg[j].xRlx[i] * c1[i])
                        z2 = z2 - (vg[j].xRlx[i] * c2[i])
                    elseif InSector(MI, O, A, B) && InSector(MI, OO, B, A)
                        vg[j].xInt[i] = 1
                        z1 = z1 - (vg[j].xRlx[i] * c1[i]) + c1[i]
                        z2 = z2 - (vg[j].xRlx[i] * c2[i]) + c2[i]
                    else
                        distMO = sqrt((O.x-MO.x)^2+(O.y-MO.y)^2)+sqrt((OO.x-MO.x)^2+(OO.y-MO.y)^2)
                        distMI = sqrt((O.x-MI.x)^2+(O.y-MI.y)^2)+sqrt((OO.x-MI.x)^2+(OO.y-MI.y)^2)
                        if InSector(MO, O, A, B) && distMO < distMI
                            vg[j].xInt[i] = 0
                            z1 = z1 - (vg[j].xRlx[i] * c1[i])
                            z2 = z2 - (vg[j].xRlx[i] * c2[i])
                        else
                            vg[j].xInt[i] = 1
                            z1 = z1 - (vg[j].xRlx[i] * c1[i]) + c1[i]
                            z2 = z2 - (vg[j].xRlx[i] * c2[i]) + c2[i]
                        end
                    end
                end
            end
            vg[j].yInt[1] += vg[j].xInt[i] * c1[i]
            vg[j].yInt[2] += vg[j].xInt[i] * c2[i]
            i += 1
        end
    elseif param_round2 == 3
        for i in 1:nbvar
            if vg[j].xInt[i] == -1
                M = point( z1 - vg[j].xRlx[i] * c1[i] , z2 - vg[j].xRlx[i] * c2[i])
                MM = point( z1 - vg[j].xRlx[i] * c1[i] + c1[i] , z2 - vg[j].xRlx[i] * c2[i] + c2[i])
                if !InSector(M, O, A, B)
                    # le point M obtenu est hors du cone => variable i fixee a 1
                    vg[j].xInt[i] = 1
                    z1 = z1 - (vg[j].xRlx[i] * c1[i]) + c1[i]
                    z2 = z2 - (vg[j].xRlx[i] * c2[i]) + c2[i]
                else
                    distm = sqrt((O.x-M.x)^2+(O.y-M.y)^2)+sqrt((OO.x-M.x)^2+(OO.y-M.y)^2)
                    distmm = sqrt((O.x-MM.x)^2+(O.y-MM.y)^2)+sqrt((OO.x-MM.x)^2+(OO.y-MM.y)^2)
                    if distm < distmm
                        vg[j].xInt[i] = 0
                        z1 = z1 - (vg[j].xRlx[i] * c1[i])
                        z2 = z2 - (vg[j].xRlx[i] * c2[i])
                    else
                        vg[j].xInt[i] = 1
                        z1 = z1 - (vg[j].xRlx[i] * c1[i]) + c1[i]
                        z2 = z2 - (vg[j].xRlx[i] * c2[i]) + c2[i]
                    end
                end
            end
            vg[j].yInt[1] += vg[j].xInt[i] * c1[i]
            vg[j].yInt[2] += vg[j].xInt[i] * c2[i]
        end
    elseif param_round2 == 4
        i = 1
        while i <= nbvar
            if vg[j].xInt[i] == -1
                nbVarNonEntiere += 1
                k = i + 1
                while k <= nbvar
                    if vg[j].xInt[k] == -1
                        nbVarNonEntiere += 1
                        MOO = point( z1 - vg[j].xRlx[i] * c1[i] - vg[j].xRlx[k] * c1[k] , z2 - vg[j].xRlx[i] * c2[i] - vg[j].xRlx[k] * c2[k])
                        MOI = point( z1 - vg[j].xRlx[i] * c1[i] - vg[j].xRlx[k] * c1[k] + c1[k] , z2 - vg[j].xRlx[i] * c2[i] - vg[j].xRlx[k] * c2[k] + c2[k])
                        MIO = point( z1 - vg[j].xRlx[i] * c1[i] + c1[i] - vg[j].xRlx[k] * c1[k] , z2 - vg[j].xRlx[i] * c2[i] + c2[i] - vg[j].xRlx[k] * c2[k])
                        MII = point( z1 - vg[j].xRlx[i] * c1[i] + c1[i] - vg[j].xRlx[k] * c1[k] + c1[k] , z2 - vg[j].xRlx[i] * c2[i] + c2[i] - vg[j].xRlx[k] * c2[k] + c2[k])
                        distMOO = sqrt((O.x-MOO.x)^2+(O.y-MOO.y)^2)+sqrt((OO.x-MOO.x)^2+(OO.y-MOO.y)^2)
                        distMOI = sqrt((O.x-MOI.x)^2+(O.y-MOI.y)^2)+sqrt((OO.x-MOI.x)^2+(OO.y-MOI.y)^2)
                        distMIO = sqrt((O.x-MIO.x)^2+(O.y-MIO.y)^2)+sqrt((OO.x-MIO.x)^2+(OO.y-MIO.y)^2)
                        distMII = sqrt((O.x-MII.x)^2+(O.y-MII.y)^2)+sqrt((OO.x-MII.x)^2+(OO.y-MII.y)^2)
                        if InSector(MOO, O, A, B) && distMOO <= distMOI && distMOO <= distMIO && distMOO <= distMII
                            vg[j].xInt[i] = 0
                            z1 = z1 - (vg[j].xRlx[i] * c1[i])
                            z2 = z2 - (vg[j].xRlx[i] * c2[i])
                            vg[j].xInt[k] = 0
                            z1 = z1 - (vg[j].xRlx[k] * c1[k])
                            z2 = z2 - (vg[j].xRlx[k] * c2[k])
                        elseif InSector(MOI, O, A, B) && distMOI <= distMIO && distMOI <= distMII
                            vg[j].xInt[i] = 0
                            z1 = z1 - (vg[j].xRlx[i] * c1[i])
                            z2 = z2 - (vg[j].xRlx[i] * c2[i])
                            vg[j].xInt[k] = 1
                            z1 = z1 - (vg[j].xRlx[k] * c1[k]) + c1[k]
                            z2 = z2 - (vg[j].xRlx[k] * c2[k]) + c2[k]
                        elseif InSector(MIO, O, A, B) && distMIO <= distMII
                            vg[j].xInt[i] = 1
                            z1 = z1 - (vg[j].xRlx[i] * c1[i]) + c1[i]
                            z2 = z2 - (vg[j].xRlx[i] * c2[i]) + c2[i]
                            vg[j].xInt[k] = 0
                            z1 = z1 - (vg[j].xRlx[k] * c1[k])
                            z2 = z2 - (vg[j].xRlx[k] * c2[k])
                        else
                            vg[j].xInt[i] = 1
                            z1 = z1 - (vg[j].xRlx[i] * c1[i]) + c1[i]
                            z2 = z2 - (vg[j].xRlx[i] * c2[i]) + c2[i]
                            vg[j].xInt[k] = 1
                            z1 = z1 - (vg[j].xRlx[k] * c1[k]) + c1[k]
                            z2 = z2 - (vg[j].xRlx[k] * c2[k]) + c2[k]
                        end
                        k = nbvar + 1
                    end
                    k += 1
                end
                if k == nbvar + 1
                    nbVarNonEntiere += 1
                    M = point( z1 - vg[j].xRlx[i] * c1[i] , z2 - vg[j].xRlx[i] * c2[i])
                    MM = point( z1 - vg[j].xRlx[i] * c1[i] + c1[i] , z2 - vg[j].xRlx[i] * c2[i] + c2[i])
                    distm = sqrt((O.x-M.x)^2+(O.y-M.y)^2)+sqrt((OO.x-M.x)^2+(OO.y-M.y)^2)
                    distmm = sqrt((O.x-MM.x)^2+(O.y-MM.y)^2)+sqrt((OO.x-MM.x)^2+(OO.y-MM.y)^2)
                    if InSector(M, O, A, B) && distm < distmm
                        vg[j].xInt[i] = 0
                        z1 = z1 - (vg[j].xRlx[i] * c1[i])
                        z2 = z2 - (vg[j].xRlx[i] * c2[i])
                    else
                        vg[j].xInt[i] = 1
                        z1 = z1 - (vg[j].xRlx[i] * c1[i]) + c1[i]
                        z2 = z2 - (vg[j].xRlx[i] * c2[i]) + c2[i]
                    end
                end
            end
            vg[j].yInt[1] += vg[j].xInt[i] * c1[i]
            vg[j].yInt[2] += vg[j].xInt[i] * c2[i]
            i += 1
        end
    elseif param_round2 == 5
        i = 1
        while i <= nbvar
            if vg[j].xInt[i] == -1
                nbVarNonEntiere += 1
                k = i + 1
                while k <= nbvar
                    if vg[j].xInt[k] == -1
                        nbVarNonEntiere += 1
                        l = k + 1
                        while l <= nbvar
                            if vg[j].xInt[l] == -1
                                nbVarNonEntiere += 1
                                MOOO = point( z1 - vg[j].xRlx[i] * c1[i] - vg[j].xRlx[k] * c1[k] - vg[j].xRlx[l] * c1[l],
                                              z2 - vg[j].xRlx[i] * c2[i] - vg[j].xRlx[k] * c2[k] - vg[j].xRlx[l] * c2[l])
                                MOOI = point( z1 - vg[j].xRlx[i] * c1[i] - vg[j].xRlx[k] * c1[k] - vg[j].xRlx[l] * c1[l] + c1[l],
                                              z2 - vg[j].xRlx[i] * c2[i] - vg[j].xRlx[k] * c2[k] - vg[j].xRlx[l] * c2[l] + c2[l])
                                MOIO = point( z1 - vg[j].xRlx[i] * c1[i] - vg[j].xRlx[k] * c1[k] - vg[j].xRlx[l] * c1[l] + c1[k],
                                              z2 - vg[j].xRlx[i] * c2[i] - vg[j].xRlx[k] * c2[k] - vg[j].xRlx[l] * c2[l] + c2[k])
                                MOII = point( z1 - vg[j].xRlx[i] * c1[i] - vg[j].xRlx[k] * c1[k] - vg[j].xRlx[l] * c1[l] + c1[k] + c1[l],
                                              z2 - vg[j].xRlx[i] * c2[i] - vg[j].xRlx[k] * c2[k] - vg[j].xRlx[l] * c2[l] + c2[k] + c2[l])
                                MIOO = point( z1 - vg[j].xRlx[i] * c1[i] - vg[j].xRlx[k] * c1[k] - vg[j].xRlx[l] * c1[l] + c1[i],
                                              z2 - vg[j].xRlx[i] * c2[i] - vg[j].xRlx[k] * c2[k] - vg[j].xRlx[l] * c2[l] + c2[i])
                                MIOI = point( z1 - vg[j].xRlx[i] * c1[i] - vg[j].xRlx[k] * c1[k] - vg[j].xRlx[l] * c1[l] + c1[i] + c1[l],
                                              z2 - vg[j].xRlx[i] * c2[i] - vg[j].xRlx[k] * c2[k] - vg[j].xRlx[l] * c2[l] + c2[i] + c2[l])
                                MIIO = point( z1 - vg[j].xRlx[i] * c1[i] - vg[j].xRlx[k] * c1[k] - vg[j].xRlx[l] * c1[l] + c1[i] + c1[k],
                                              z2 - vg[j].xRlx[i] * c2[i] - vg[j].xRlx[k] * c2[k] - vg[j].xRlx[l] * c2[l] + c2[i] + c2[k])
                                MIII = point( z1 - vg[j].xRlx[i] * c1[i] - vg[j].xRlx[k] * c1[k] - vg[j].xRlx[l] * c1[l] + c1[i] + c1[k] + c1[l],
                                              z2 - vg[j].xRlx[i] * c2[i] - vg[j].xRlx[k] * c2[k] - vg[j].xRlx[l] * c2[l] + c2[i] + c2[k] + c2[l])
                                distMOOO = sqrt((O.x-MOOO.x)^2+(O.y-MOOO.y)^2)+sqrt((OO.x-MOOO.x)^2+(OO.y-MOOO.y)^2)
                                distMOOI = sqrt((O.x-MOOI.x)^2+(O.y-MOOI.y)^2)+sqrt((OO.x-MOOI.x)^2+(OO.y-MOOI.y)^2)
                                distMOIO = sqrt((O.x-MOIO.x)^2+(O.y-MOIO.y)^2)+sqrt((OO.x-MOIO.x)^2+(OO.y-MOIO.y)^2)
                                distMOII = sqrt((O.x-MOII.x)^2+(O.y-MOII.y)^2)+sqrt((OO.x-MOII.x)^2+(OO.y-MOII.y)^2)
                                distMIOO = sqrt((O.x-MIOO.x)^2+(O.y-MIOO.y)^2)+sqrt((OO.x-MIOO.x)^2+(OO.y-MIOO.y)^2)
                                distMIOI = sqrt((O.x-MIOI.x)^2+(O.y-MIOI.y)^2)+sqrt((OO.x-MIOI.x)^2+(OO.y-MIOI.y)^2)
                                distMIIO = sqrt((O.x-MIIO.x)^2+(O.y-MIIO.y)^2)+sqrt((OO.x-MIIO.x)^2+(OO.y-MIIO.y)^2)
                                distMIII = sqrt((O.x-MIII.x)^2+(O.y-MIII.y)^2)+sqrt((OO.x-MIII.x)^2+(OO.y-MIII.y)^2)
                                if InSector(MOOO, O, A, B) && distMOOO <= distMOOI && distMOOO <= distMOIO && distMOOO <= distMOII && distMOOO <= distMIOO && distMOOO <= distMIOI && distMOOO <= distMIIO && distMOOO <= distMIII
                                    vg[j].xInt[i] = 0
                                    z1 = z1 - (vg[j].xRlx[i] * c1[i])
                                    z2 = z2 - (vg[j].xRlx[i] * c2[i])
                                    vg[j].xInt[k] = 0
                                    z1 = z1 - (vg[j].xRlx[k] * c1[k])
                                    z2 = z2 - (vg[j].xRlx[k] * c2[k])
                                    vg[j].xInt[l] = 0
                                    z1 = z1 - (vg[j].xRlx[l] * c1[l])
                                    z2 = z2 - (vg[j].xRlx[l] * c2[l])
                                elseif InSector(MOOI, O, A, B) && distMOOI <= distMOIO && distMOOI <= distMOII && distMOOI <= distMIOO && distMOOI <= distMIOI && distMOOI <= distMIIO && distMOOI <= distMIII
                                    vg[j].xInt[i] = 0
                                    z1 = z1 - (vg[j].xRlx[i] * c1[i])
                                    z2 = z2 - (vg[j].xRlx[i] * c2[i])
                                    vg[j].xInt[k] = 0
                                    z1 = z1 - (vg[j].xRlx[k] * c1[k])
                                    z2 = z2 - (vg[j].xRlx[k] * c2[k])
                                    vg[j].xInt[l] = 1
                                    z1 = z1 - (vg[j].xRlx[l] * c1[l]) + c1[l]
                                    z2 = z2 - (vg[j].xRlx[l] * c2[l]) + c2[l]
                                elseif InSector(MOIO, O, A, B) && distMOIO <= distMOII && distMOIO <= distMIOO && distMOIO <= distMIOI && distMOIO <= distMIIO && distMOIO <= distMIII
                                    vg[j].xInt[i] = 0
                                    z1 = z1 - (vg[j].xRlx[i] * c1[i])
                                    z2 = z2 - (vg[j].xRlx[i] * c2[i])
                                    vg[j].xInt[k] = 1
                                    z1 = z1 - (vg[j].xRlx[k] * c1[k]) + c1[k]
                                    z2 = z2 - (vg[j].xRlx[k] * c2[k]) + c2[k]
                                    vg[j].xInt[l] = 0
                                    z1 = z1 - (vg[j].xRlx[l] * c1[l])
                                    z2 = z2 - (vg[j].xRlx[l] * c2[l])
                                elseif InSector(MOII, O, A, B) && distMOII <= distMIOO && distMOII <= distMIOI && distMOII <= distMIIO && distMOII <= distMIII
                                    vg[j].xInt[i] = 0
                                    z1 = z1 - (vg[j].xRlx[i] * c1[i])
                                    z2 = z2 - (vg[j].xRlx[i] * c2[i])
                                    vg[j].xInt[k] = 1
                                    z1 = z1 - (vg[j].xRlx[k] * c1[k]) + c1[k]
                                    z2 = z2 - (vg[j].xRlx[k] * c2[k]) + c2[k]
                                    vg[j].xInt[l] = 1
                                    z1 = z1 - (vg[j].xRlx[l] * c1[l]) + c1[l]
                                    z2 = z2 - (vg[j].xRlx[l] * c2[l]) + c2[l]
                                elseif InSector(MIOO, O, A, B) && distMIOO <= distMIOI && distMIOO <= distMIIO && distMIOO <= distMIII
                                    vg[j].xInt[i] = 1
                                    z1 = z1 - (vg[j].xRlx[i] * c1[i]) + c1[i]
                                    z2 = z2 - (vg[j].xRlx[i] * c2[i]) + c2[i]
                                    vg[j].xInt[k] = 0
                                    z1 = z1 - (vg[j].xRlx[k] * c1[k])
                                    z2 = z2 - (vg[j].xRlx[k] * c2[k])
                                    vg[j].xInt[l] = 0
                                    z1 = z1 - (vg[j].xRlx[l] * c1[l])
                                    z2 = z2 - (vg[j].xRlx[l] * c2[l])
                                elseif InSector(MIOI, O, A, B) && distMIOI <= distMIIO && distMIOI <= distMIII
                                    vg[j].xInt[i] = 1
                                    z1 = z1 - (vg[j].xRlx[i] * c1[i]) + c1[i]
                                    z2 = z2 - (vg[j].xRlx[i] * c2[i]) + c2[i]
                                    vg[j].xInt[k] = 0
                                    z1 = z1 - (vg[j].xRlx[k] * c1[k])
                                    z2 = z2 - (vg[j].xRlx[k] * c2[k])
                                    vg[j].xInt[l] = 1
                                    z1 = z1 - (vg[j].xRlx[l] * c1[l]) + c1[l]
                                    z2 = z2 - (vg[j].xRlx[l] * c2[l]) + c2[l]
                                elseif InSector(MIIO, O, A, B) && distMIIO <= distMIII
                                    vg[j].xInt[i] = 1
                                    z1 = z1 - (vg[j].xRlx[i] * c1[i]) + c1[i]
                                    z2 = z2 - (vg[j].xRlx[i] * c2[i]) + c2[i]
                                    vg[j].xInt[k] = 1
                                    z1 = z1 - (vg[j].xRlx[k] * c1[k]) + c1[k]
                                    z2 = z2 - (vg[j].xRlx[k] * c2[k]) + c2[k]
                                    vg[j].xInt[l] = 0
                                    z1 = z1 - (vg[j].xRlx[l] * c1[l])
                                    z2 = z2 - (vg[j].xRlx[l] * c2[l])
                                else
                                    vg[j].xInt[i] = 1
                                    z1 = z1 - (vg[j].xRlx[i] * c1[i]) + c1[i]
                                    z2 = z2 - (vg[j].xRlx[i] * c2[i]) + c2[i]
                                    vg[j].xInt[k] = 1
                                    z1 = z1 - (vg[j].xRlx[k] * c1[k]) + c1[k]
                                    z2 = z2 - (vg[j].xRlx[k] * c2[k]) + c2[k]
                                    vg[j].xInt[l] = 1
                                    z1 = z1 - (vg[j].xRlx[l] * c1[l]) + c1[l]
                                    z2 = z2 - (vg[j].xRlx[l] * c2[l]) + c2[l]
                                end
                                l = nbvar + 1
                                k = nbvar + 1
                            end
                            l += 1
                        end
                        if l == nbvar + 1
                            MOO = point( z1 - vg[j].xRlx[i] * c1[i] - vg[j].xRlx[k] * c1[k] , z2 - vg[j].xRlx[i] * c2[i] - vg[j].xRlx[k] * c2[k])
                            MOI = point( z1 - vg[j].xRlx[i] * c1[i] - vg[j].xRlx[k] * c1[k] + c1[k] , z2 - vg[j].xRlx[i] * c2[i] - vg[j].xRlx[k] * c2[k] + c2[k])
                            MIO = point( z1 - vg[j].xRlx[i] * c1[i] + c1[i] - vg[j].xRlx[k] * c1[k] , z2 - vg[j].xRlx[i] * c2[i] + c2[i] - vg[j].xRlx[k] * c2[k])
                            MII = point( z1 - vg[j].xRlx[i] * c1[i] + c1[i] - vg[j].xRlx[k] * c1[k] + c1[k] , z2 - vg[j].xRlx[i] * c2[i] + c2[i] - vg[j].xRlx[k] * c2[k] + c2[k])
                            distMOO = sqrt((O.x-MOO.x)^2+(O.y-MOO.y)^2)+sqrt((OO.x-MOO.x)^2+(OO.y-MOO.y)^2)
                            distMOI = sqrt((O.x-MOI.x)^2+(O.y-MOI.y)^2)+sqrt((OO.x-MOI.x)^2+(OO.y-MOI.y)^2)
                            distMIO = sqrt((O.x-MIO.x)^2+(O.y-MIO.y)^2)+sqrt((OO.x-MIO.x)^2+(OO.y-MIO.y)^2)
                            distMII = sqrt((O.x-MII.x)^2+(O.y-MII.y)^2)+sqrt((OO.x-MII.x)^2+(OO.y-MII.y)^2)
                            if InSector(MOO, O, A, B) && distMOO <= distMOI && distMOO <= distMIO && distMOO <= distMII
                                vg[j].xInt[i] = 0
                                z1 = z1 - (vg[j].xRlx[i] * c1[i])
                                z2 = z2 - (vg[j].xRlx[i] * c2[i])
                                vg[j].xInt[k] = 0
                                z1 = z1 - (vg[j].xRlx[k] * c1[k])
                                z2 = z2 - (vg[j].xRlx[k] * c2[k])
                            elseif InSector(MOI, O, A, B) && distMOI <= distMIO && distMOI <= distMII
                                vg[j].xInt[i] = 0
                                z1 = z1 - (vg[j].xRlx[i] * c1[i])
                                z2 = z2 - (vg[j].xRlx[i] * c2[i])
                                vg[j].xInt[k] = 1
                                z1 = z1 - (vg[j].xRlx[k] * c1[k]) + c1[k]
                                z2 = z2 - (vg[j].xRlx[k] * c2[k]) + c2[k]
                            elseif InSector(MIO, O, A, B) && distMIO <= distMII
                                vg[j].xInt[i] = 1
                                z1 = z1 - (vg[j].xRlx[i] * c1[i]) + c1[i]
                                z2 = z2 - (vg[j].xRlx[i] * c2[i]) + c2[i]
                                vg[j].xInt[k] = 0
                                z1 = z1 - (vg[j].xRlx[k] * c1[k])
                                z2 = z2 - (vg[j].xRlx[k] * c2[k])
                            else
                                vg[j].xInt[i] = 1
                                z1 = z1 - (vg[j].xRlx[i] * c1[i]) + c1[i]
                                z2 = z2 - (vg[j].xRlx[i] * c2[i]) + c2[i]
                                vg[j].xInt[k] = 1
                                z1 = z1 - (vg[j].xRlx[k] * c1[k]) + c1[k]
                                z2 = z2 - (vg[j].xRlx[k] * c2[k]) + c2[k]
                            end
                            k = nbvar + 1
                        end
                    end
                    k += 1
                end
                if k == nbvar + 1
                    M = point( z1 - vg[j].xRlx[i] * c1[i] , z2 - vg[j].xRlx[i] * c2[i])
                    MM = point( z1 - vg[j].xRlx[i] * c1[i] + c1[i] , z2 - vg[j].xRlx[i] * c2[i] + c2[i])
                    distm = sqrt((O.x-M.x)^2+(O.y-M.y)^2)+sqrt((OO.x-M.x)^2+(OO.y-M.y)^2)
                    distmm = sqrt((O.x-MM.x)^2+(O.y-MM.y)^2)+sqrt((OO.x-MM.x)^2+(OO.y-MM.y)^2)
                    if InSector(M, O, A, B) && distm < distmm
                        vg[j].xInt[i] = 0
                        z1 = z1 - (vg[j].xRlx[i] * c1[i])
                        z2 = z2 - (vg[j].xRlx[i] * c2[i])
                    else
                        vg[j].xInt[i] = 1
                        z1 = z1 - (vg[j].xRlx[i] * c1[i]) + c1[i]
                        z2 = z2 - (vg[j].xRlx[i] * c2[i]) + c2[i]
                    end
                end
            end
            vg[j].yInt[1] += vg[j].xInt[i] * c1[i]
            vg[j].yInt[2] += vg[j].xInt[i] * c2[i]
            i += 1
        end
    end

    #scatter(Xapx,Yapx, color="black", marker=".")
    return vg, nbVarNonEntiere
end

# ----------------------------------------------------------------
# Cherche l'indice de la variable |x*-x~| la plus grande
function hightest(xstar::Array{Float64, 1}, xtild::Array{Int, 1}, supr::Array{Int, 1}, nbvar::Int)
    indice = 0
    max = 0
    for i in 1:nbvar
        if i in supr
            #rien
        elseif abs(xstar[i] - xtild[i]) > max
            max = abs(xstar[i] - xtild[i])
            indice = i
        end
    end
    return indice
end

# ----------------------------------------------------------------
# Sélectionne les points pour former l'ensemble bornant primal
function kung(XFeas, YFeas)
    S = []
    for i=1:length(XFeas)
        push!(S, (XFeas[i] , YFeas[i]) )
    end
    sort!(S, by = x -> x[1])
    SN=[] ; push!(SN, S[1]) ; minYFeas = S[1][2]
    for i=2:length(XFeas)
        if S[i][2] < minYFeas
            push!(SN, S[i]) ; minYFeas = S[i][2]
        end
    end
    return SN
end

function LS(nbvar::Int, nbcontraintes::Int, L::Array{Int,2}, c1::Array{Int,1}, c2::Array{Int,1}, nb_de_un::Array{Int,1}, L_LS::Array{Int,2})
    origine = []
    desti = []
    i = 1
    while i <= nbvar-1
        j = i+1
        while j <= nbvar
            if nb_de_un[i] == nb_de_un[j]
                k = 1
            else
                k = 0
                i = j-1
                j = nbvar+1
            end
            while k > 0
                if L_LS[k,i] != L_LS[k,j]
                    k = 0
                    i = j-1
                    j = nbvar+1
                else
                    k += 1
                end
                if k > 0 && L_LS[k,i] == 0
                    if c1[i] + c2[i] > c1[j] + c2[j]
                        append!(origine, i)
                        append!(desti, j)
                        i += 1
                    else
                        append!(origine, j)
                        append!(desti, i)
                        i += 1
                    end
                    k = 0
                end
            end
            j += 1
        end
        i += 1
    end
    return origine, desti
end

function main(fname::String, tailleSampling::Int64, terraform::Int64, param_round2 = 5, param_LS = 1)

    #@printf("Running the gravity machine...\n\n")
    @assert tailleSampling>=3 "Erreur : Au moins 3 sont requis"

    verification(fname)

    # chargement de l'instance numerique ---------------------------------------
    nbvar, nbcontraintes, L, c1, c2, nb_de_un, L_LS = SPA_loaddemandebi("SPA_Databio/bio"*fname*".txt") # instance numerique de SPA
    nbobj = 2
    nbgen = tailleSampling

    start = time()

    # allocation de memoire pour la structure de donnees -----------------------
    vg=Vector{tGenerateur}(undef,nbgen)
    for j=1:nbgen
        vg[j]=tGenerateur( zeros(Float64,nbvar), zeros(Float64,nbobj),
                           zeros(Int64,nbvar),   zeros(Int64,nbobj),
                           zeros(Float64,nbvar), zeros(Float64,nbobj),
                           false
                         )
    end

    # allocation de memoire pour les ensembles bornants ------------------------
    eb=Vector{tEnsBornant}(undef,nbgen)
    for j=1:nbgen
        eb[j]=tEnsBornant( zeros(Float64,nbvar), zeros(Float64,nbobj)
                         )
    end

    # listes de points pour les affichages graphiques --------------------------
    XEBD=(Float64)[]; YEBD=(Float64)[]    # liste des points (x,y) relaches
    XInt=(Int64)[];  YInt=(Int64)[]       # liste des points (x,y) entiers
    XProj=(Float64)[]; YProj=(Float64)[]  # liste des points (x,y) projetes
    XFeas=(Int64)[]; YFeas=(Int64)[]      # liste des points (x,y) admissibles

    # --------------------------------------------------------------------------
    # --------------------------------------------------------------------------
    #@printf("1) calcule les etendues de valeurs sur les 2 objectifs\n\n")

    # calcule la valeur optimale relachee de f1 seule et le point (z1,z2) correspondant
    f1RL, xf1RL = relaxLinXG(nbvar, nbcontraintes, L, c1, c2, typemax(Int), 0, 1) # opt fct 1
    z1RL1, z2RL1 = evaluerSolution(xf1RL, c1, c2)

    # calcule la valeur optimale relachee de f2 seule et le point (z1,z2) correspondant
    f2RL, xf2RL = relaxLinXG(nbvar, nbcontraintes, L, c1, c2, typemax(Int), 0, 2) # opt fct 2
    z1RL2, z2RL2 = evaluerSolution(xf2RL, c1, c2)

    #@printf("  f1_min=%8.2f ↔ f1_max=%8.2f (Δ=%.2f) \n",z1RL1, z1RL2, z1RL2-z1RL1)
    #@printf("  f2_min=%8.2f ↔ f2_max=%8.2f (Δ=%.2f) \n\n",z2RL2, z2RL1, z2RL1-z2RL2)

    if param_LS == 1
        origine, desti = LS(nbvar, nbcontraintes, L, c1, c2, nb_de_un, L_LS)
    end

    # --------------------------------------------------------------------------
    # --------------------------------------------------------------------------
    #@printf("2) calcule les generateurs par e-contrainte avec z1 la fonction a minimiser\n\n")

    pasSampleV = 1.6*(z2RL1 - z2RL2) / ((tailleSampling-1)) # pas vertical de l'echantillonage
    pasSampleH = 1.6*(z1RL2 - z1RL1) / ((tailleSampling-1)) # pas horizontal de l'echantillo

    for j = 1:tailleSampling
        if j<floor(Int, tailleSampling/2)+1
            #@printf("  %2d : ϵ = %8.2f  ", j, z2RL1 - (j-1) * pasSampleV)

            # calcul d'une solution epsilon-contrainte -----------------------------
            fRL, eb[j].xEBD = relaxLinXG4(nbvar, nbcontraintes, L, c1, c2, z2RL1 - (j-1) * pasSampleV, 1)
            #@printf("fRL = %8.2f  ",round(fRL, digits=2))

        else
            #@printf("  %2d : ϵ = %8.2f  ", j, z1RL2 - (2 * floor(Int, tailleSampling/2) - j) * pasSampleH)

            # calcul d'une solution epsilon-contrainte -----------------------------
            fRL, eb[j].xEBD = relaxLinXG4(nbvar, nbcontraintes, L, c1, c2, z1RL2 - (2 * floor(Int, tailleSampling/2) - j) * pasSampleH, 2)
            #@printf("fRL = %8.2f  ",round(fRL, digits=2))
        end

        # nettoyage de la valeur de eb[j].xEBD et calcul du point bi-objectif --
        vg = calcNoInteger(vg, nbvar, j)

        eb[j].yEBD[1], eb[j].yEBD[2] = evaluerSolution(eb[j].xEBD, c1, c2)
        #@printf("[ %8.2f , %8.2f ] ", eb[j].yEBD[1], eb[j].yEBD[2])

        # copie de l'ensemble bornant dual dans la stru de donnees iterative ---
        vg[j].xRlx = deepcopy(eb[j].xEBD)
        vg[j].yRlx[1] = eb[j].yEBD[1]
        vg[j].yRlx[2] = eb[j].yEBD[2]
        push!(XEBD,eb[j].yEBD[1])
        push!(YEBD,eb[j].yEBD[2])

        # test d'admissibilite et marquage de la solution le cas echeant -------
        if estAdmissible(eb[j].xEBD)
            for i = 1:nbvar
                vg[j].xInt[i] = convert(Int, round(eb[j].xEBD[i]))
            end
            vg[j].yInt[1] = eb[j].yEBD[1]
            vg[j].yInt[2] = eb[j].yEBD[2]
            vg[j].xFea = true
            push!(XFeas,eb[j].yEBD[1])
            push!(YFeas,eb[j].yEBD[2])

            #ajout de la recherche locale
            if param_LS == 1
                for i in 1:length(origine)
                    if eb[j].xEBD[origine[i]] == 1

                        eb[j].xEBD[origine[i]] = 0
                        eb[j].xEBD[desti[i]] = 1

                        vg[j].yInt[1] = vg[j].yInt[1] - c1[origine[i]] + c1[desti[i]]
                        vg[j].yInt[2] = vg[j].yInt[2] - c2[origine[i]] + c2[desti[i]]
                    end
                end
                push!(XFeas, vg[j].yInt[1])
                push!(YFeas, vg[j].yInt[2])
            end

            #@printf(" Admissible \n")
        else
            vg[j].xFea = false
            #@printf(" x \n")
        end

    end
    #println("")

    # ==========================================================================
    #@printf("3) terraformation en cours...\n\n")

    for j in 1:tailleSampling
        if vg[j].xFea == false
            # Nettoyage de la valeur de vg[j].xRlx et calcul du point bi-objectif
            vg = calcNoInteger(vg, nbvar, j)

            vg[j].yRlx[1],vg[j].yRlx[2] = evaluerSolution(vg[j].xRlx, c1, c2)
            #@printf("  %2d : [ %8.2f , %8.2f ] ", j, vg[j].yRlx[1], vg[j].yRlx[2])

            # ----------------------------------------------------------------
            # Selectionne les points pour constituer le cone d'interet centre sur le point generateur O
            B, O, A = selectPoints(j, tailleSampling, eb)

            # ----------------------------------------------------------------
            # Operations pour arrondir les valeurs non-entieres de la solution relachee
            vg, nbVarNonEntiere = roundOneCone(j, nbvar, vg, c1, c2, O, A, B)

            #@printf("  → #round : %4d → [ %5d , %5d ] ", nbVarNonEntiere, vg[j].yInt[1], vg[j].yInt[2])
            push!(XInt,vg[j].yInt[1])
            push!(YInt,vg[j].yInt[2])

            trial = 0
        end

        while vg[j].xFea == false && trial < terraform
            #print(" -->")
            # ----------------------------------------------------------------
            # Projete la solution entiere sur le polytope X avec norme-L1
            fPrj, vg[j].xPrj = Δ2(L, vg[j].xInt)

            # Nettoyage de la valeur de vg[j].xPrj et calcul du point bi-objectif
            vg = calcNoInteger(vg, nbvar, j)

            vg[j].yPrj[1], vg[j].yPrj[2] = evaluerSolution(vg[j].xPrj, c1, c2)

            #@printf(" [ %8.2f , %8.2f ] ", vg[j].yPrj[1], vg[j].yPrj[2])
            push!(XProj, vg[j].yPrj[1])
            push!(YProj, vg[j].yPrj[2])

            # ----------------------------------------------------------------
            # Teste si la projection est admissible
            if estAdmissible(vg[j].xPrj)
                vg[j].xFea = true
                vg[j].yInt[1] = vg[j].yPrj[1]
                vg[j].yInt[2] = vg[j].yPrj[2]
                push!(XFeas, vg[j].yInt[1])
                push!(YFeas, vg[j].yInt[2])

                #ajout de la recherche locale
                if param_LS == 1
                    for i in 1:length(origine)
                        if vg[j].xPrj[origine[i]] == 1

                            vg[j].xPrj[origine[i]] = 0
                            vg[j].xPrj[desti[i]] = 1

                            vg[j].yInt[1] = vg[j].yInt[1] - c1[origine[i]] + c1[desti[i]]
                            vg[j].yInt[2] = vg[j].yInt[2] - c2[origine[i]] + c2[desti[i]]
                        end
                    end
                    push!(XFeas, vg[j].yInt[1])
                    push!(YFeas, vg[j].yInt[2])
                end
                #@printf(" Admissible \n")
            else
                vg[j].xFea = false
                #@printf(" x \n")
            end
            if vg[j].xFea == false && trial < (terraform - 1)

                # prepare pour l'iteration suivante
                OO = point(vg[j].yPrj[1], vg[j].yPrj[2])

                # ----------------------------------------------------------------
                # Operations pour arrondir les valeurs non-entieres de la solution relachee
                vg, nbVarNonEntiere = roundTwoCones(j, nbvar, vg, c1, c2, O, A, B, OO, param_round2)

                # ----------------------------------------------------------------
                # Gestion de 1-cycle
                if vg[j].yInt[1] == XInt[end] && vg[j].yInt[2] == YInt[end] && trial > 0
                    T = 30
                    supr = Int64[]
                    p = rand(T/2:3*T/4)
                    i = 1
                    while i <= p
                        indice = hightest(vg[j].xPrj, vg[j].xInt, supr, nbvar)
                        if indice != 0
                            append!(supr, indice)
                            if vg[j].xInt[indice] == 0
                                M = point( vg[j].yInt[1] + c1[indice] , vg[j].yInt[2] + c2[indice])
                                if InSector(M, O, A, B)
                                    if InSector(M, OO, B, A)
                                        vg[j].xInt[indice] = 1
                                        vg[j].yInt[1] += c1[indice]
                                        vg[j].yInt[2] += c2[indice]
                                    end
                                end
                            elseif vg[j].xInt[indice] == 1
                                M = point( vg[j].yInt[1] - c1[indice] , vg[j].yInt[2] - c2[indice])
                                if InSector(M, O, A, B)
                                    if InSector(M, OO, B, A)
                                        vg[j].xInt[indice] = 0
                                        vg[j].yInt[1] -= c1[indice]
                                        vg[j].yInt[2] -= c2[indice]
                                    end
                                end
                            end
                        i += 1
                        else
                            i = p+1
                        end
                    end
                end
                #@printf("  → #round : %4d → [ %5d , %5d ] ", nbVarNonEntiere, vg[j].yInt[1], vg[j].yInt[2])
                push!(XInt,vg[j].yInt[1])
                push!(YInt,vg[j].yInt[2])
            end
            trial += 1
            if length(XInt) > 1 && XInt[end] == XInt[end-1] && YInt[end] == YInt[end-1]
                trial = terraform
            end
        end
        #println("")
    end

    temps = time() - start

    # ==========================================================================
    #@printf("\n3) Edition des resultats \n\n")

    # Donne les points relaches initiaux ---------------------------------------
    #scatter(XEBD,YEBD,color="blue", marker="x", label = "images relaxations lineaires")
    #@show XEBD
    #@show YEBD

    # Donne les points entiers -------------------------------------------------
    #scatter(XInt,YInt,color="orange", marker="s", label = " images arrondis rela lin")
    #@show XInt
    #@show YInt

    # Donne les points apres projection Δ(x,x̃) ---------------------------------
    #scatter(XProj,YProj, color="red", marker="x", label = "images projections")
    #@show XProj
    #@show YProj

    # Donne les points admissibles ---------------------------------------------
    #scatter(XFeas,YFeas, color="green", marker="o", label = "images sol admissibles")
    #@show XFeas
    #@show YFeas
    XSN = []
    YSN = []
    XSNR = []
    YSNR = []
    if length(XFeas) > 0
        SN = kung(XFeas, YFeas)
        push!(XSN, SN[1][1])
        push!(YSN, 1.1 * maximum(YInt))
        for i in 1:length(SN)-1
            push!(XSN, SN[i][1])
            push!(YSN, SN[i][2])
            push!(XSN, SN[i+1][1])
            push!(YSN, SN[i][2])
            push!(XSNR, SN[i][1])
            push!(YSNR, SN[i][2])
        end
        push!(XSN, SN[end][1])
        push!(YSN, SN[end][2])
        push!(XSN, 1.1 * maximum(XInt))
        push!(YSN, SN[end][2])
        push!(XSNR, SN[end][1])
        push!(YSNR, SN[end][2])
        #plot(XSN, YSN, color="green", marker="o", label = "ensemble bornant primal")
        #scatter(XSNR, YSNR, color="green", label = "images sol admissibles non dominées", s = 150, alpha = 0.3)
    end

    #Donne les solutions optimales du problème ---------------------------------
    X_Y_N, Y_Y_N = lecture_Y_N(fname::String)
    SM = kung(X_Y_N, Y_Y_N)
    XSM = []
    YSM = []
    XSMR = []
    YSMR = []
    for i in 1:length(SM)-1
        push!(XSM, SM[i][1])
        push!(YSM, SM[i][2])
        push!(XSM, SM[i+1][1])
        push!(YSM, SM[i][2])
        push!(XSMR, SM[i][1])
        push!(YSMR, SM[i][2])
    end
    push!(XSM, SM[end][1])
    push!(YSM, SM[end][2])
    push!(XSMR, SM[end][1])
    push!(YSMR, SM[end][2])
    #plot(XSM, YSM, color="black", label = "ensemble solutions optimales")
    #scatter(XSM, YSM, color="black", marker="+", label = "images points non dominés")
    #@show X_Y_N
    #@show Y_Y_N

    print(" & ", length(XSNR) )

    print(" & ", round(temps, digits = 2))

    println(" & ", round(area(XSM, YSM, XSN, YSN)*100, digits = 2), " \\\\")


    #legend(bbox_to_anchor=[1,1], loc=0, borderaxespad=0, fontsize = "x-small", title = "Legende")
    #PyPlot.title("f-p bi-obj conique 1 round: $fname")

end

#=
fname :
    nom de l'instance sans le .txt
tailleSampling :
    nombre de générateurs
terraform :
    nombre d'itération par générateur
Param_round2 :
    1. arrondi 1 variabale par 1 varaible avec 2 cones
    2. arrondi 2 varaibles par 2 variables avec 2 cones
    3. arrondi 1 varaible par 1 varaible avec 1 cone et 1 direction
    4. arrondi 2 varaibles par 2 varaibles avec 1 cone et 1 direction
    5. arrondi 3 varaibles par 3 varaibles avec 1 cone et 1 direction
=#

#println("Entrer nom de l'instance sans le .txt : ")
#s = chomp(readline())
#print("$s")


main("sppnw19", 30, 5)

#=
println("")
println("\\hline")
print("sppaa03 ")
main("sppaa03", 30, 5)
println("\\hline")
print("sppaa05 ")
main("sppaa05", 30, 5)
println("\\hline")
print("sppnw01 ")
main("sppnw01", 30, 5)
println("\\hline")
for i in 3:9
    b = string(i)
    s = "sppnw0"*b
    print("$s ")
    main("$s", 30, 5)
    println("\\hline")
end

for i in 19:43
    b = string(i)
    s = "sppnw"*b
    print("$s ")
    main("$s", 30, 5)
    println("\\hline")
end
=#
nothing
