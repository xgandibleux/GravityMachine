# ==============================================================================
# The gravity machine (Man of Steel) -> to terraform the world

println("""\nAlgorithme "Gravity machine" --------------------------------\n""")

const verbose = true
const graphic = true

println("-) Active les packages requis\n")
using JuMP, GLPK, PyPlot, Printf, Random
verbose ? println("  Fait \n") : nothing


# ==============================================================================

include("GMdatastructures.jl") # types, datastructures and global variables specially defined for GM
include("GMparsers.jl")        # parsers of instances and non-dominated points
include("GMjumpModels.jl")     # JuMP models for computing optima
include("GMmopPrimitives.jl")  # usuals algorithms in multiobjective optimization
include("GMperturbation.jl")   # routines dealing with the perturbation of a solution when a cycle is detected


# ==============================================================================
# Ajout d'une solution relachee initiale a un generateur

function ajouterX0!(vg::Vector{tGenerateur}, k::Int64, s::tSolution{Float64})

    vg[k].sRel = deepcopy(s) # met en place le generateur \bar{x}^k
    vg[k].sPrj = deepcopy(s) # le generateur est la premiere projection \bar{x}^{k,0}
    return nothing
end


# ==============================================================================
# Ajout d'une solution entiere (arrondie ou perturbee) a un generateur

function ajouterXtilde!(vg::Vector{tGenerateur}, k::Int64, x::Vector{Int64}, y::Vector{Int64})

    vg[k].sInt.x = copy(x)
    vg[k].sInt.y = copy(y)
    return nothing
end


# ==============================================================================
# Ajout d'une solution fractionnaire (projetee) a un generateur

function ajouterXbar!(vg::Vector{tGenerateur}, k::Int64, x::Vector{Float64}, y::Vector{Float64})

    vg[k].sPrj.x = copy(x)
    vg[k].sPrj.y = copy(y)
    return nothing
end


# ==============================================================================
# Elabore 2 ensembles d'indices selon que xTilde[i] vaut 0 ou 1

function split01(xTilde::Array{Int,1})

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


# ==============================================================================
# Projete xTilde sur le polyedre X

function Δ2SPA(A::Array{Int,2}, xTilde::Array{Int,1})

    nbctr = size(A,1)
    nbvar = size(A,2)
    idxTilde0, idxTilde1 = split01(xTilde)

    proj = Model(GLPK.Optimizer)
    @variable(proj, 0.0 <= x[1:length(xTilde)] <= 1.0 )
    @objective(proj, Min, sum(x[i] for i in idxTilde0) + sum((1-x[i]) for i in idxTilde1) )
    @constraint(proj, [i=1:nbctr],(sum((x[j]*A[i,j]) for j in 1:nbvar)) == 1)
    optimize!(proj)
    return objective_value(proj), value.(x)
end

# ==============================================================================
# Projete xTilde sur le polyedre X

function Δ2SPAbis(A::Array{Int,2}, xTilde::Array{Int,1}, 
                  c1::Array{Int,1}, c2::Array{Int,1}, k::Int64, λ1::Vector{Float64}, λ2::Vector{Float64})

    nbctr = size(A,1)
    nbvar = size(A,2)
    idxTilde0, idxTilde1 = split01(xTilde)

#    cλ = 0.5.*c1 + 0.5.*c2
    proj = Model(GLPK.Optimizer)
    @variable(proj, 0.0 <= x[1:length(xTilde)] <= 1.0 )
    @objective(proj, Min, sum(λ1[k]*x[i] for i in idxTilde0) + sum(λ2[k]*(1-x[i]) for i in idxTilde1) )
#    @objective(proj, Min, sum(cλ[i]*x[i] for i in idxTilde0) + sum(cλ[i]*(1-x[i]) for i in idxTilde1) )
    @constraint(proj, [i=1:nbctr],(sum((x[j]*A[i,j]) for j in 1:nbvar)) == 1)
    optimize!(proj)
    return objective_value(proj), value.(x)
end


# ==============================================================================
# test si une solution est admissible en verifiant si sa relaxation lineaire
# conduit a une solution entiere

function estAdmissible(x::Vector{Float64})

    admissible = true
    i=1
    while admissible && i<=length(x)
        if round(x[i], digits=3)!=0.0 && round(x[i], digits=3)!=1.0
            admissible = false
        end
        i+=1
    end
    return admissible
end


# ==============================================================================
# calcule la performance z d'une solution x sur les 2 objectifs

function evaluerSolution(x::Vector{Float64}, c1::Array{Int,1}, c2::Array{Int,1})

    z1 = 0.0; z2 = 0.0
    for i in 1:length(x)
        z1 += x[i] * c1[i]
        z2 += x[i] * c2[i]
    end
    return round(z1, digits=2), round(z2, digits=2)
end


# ==============================================================================
# nettoyage des valeurs des variables d'une solution x relachee sur [0,1]

function nettoyageSolution!(x::Vector{Float64})

    nbvar=length(x)
    for i in 1:nbvar
        if     round(x[i], digits=3) == 0.0
                   x[i] = 0.0
        elseif round(x[i], digits=3) == 1.0
                   x[i] = 1.0
        else
                   x[i] = round(x[i], digits=3)
        end
    end
end


# ==============================================================================
# Calcul des generateurs avec 2 ϵ-contraintes alternees jusqu'a leur rencontre

function calculGenerateurs(A::Array{Int,2}, c1::Array{Int,1}, c2::Array{Int,1}, 
                           tailleSampling::Int64, 
                           minf1RL::Float64, maxf2RL::Float64, maxf1RL::Float64, minf2RL::Float64,
                           d::tListDisplay)

    nbctr = size(A,1)
    nbvar = size(A,2)

    L1 = (tSolution{Float64})[]
    L2 = (tSolution{Float64})[]

    # Premiere generation de points (avec z1 la fonction a minimiser) ----------
    pasSample2 = (maxf2RL - minf2RL) / (tailleSampling-1) # pas de l'echantillonage sur z2
    j1 = 1

    # Seconde generation de points (avec z2 la fonction a minimiser) -----------
    pasSample1 = (maxf1RL-minf1RL) / (tailleSampling-1) # pas de l'echantillonage sur z1
    j2 = 1

    alternance = 1
    maxf2RLlimite = maxf2RL
    minf2RLlimite = minf2RL

    while maxf2RLlimite > minf2RLlimite
        alternance += 1
        if alternance % 2 == 0

            # minimise f1 avec ϵ-contrainte sur f2 -----------------------------
            verbose ? @printf("  z1 %2d : ϵ = %8.2f  ", j1, maxf2RL - (j1-1) * pasSample2) : nothing # echantillonage sur z2

            # calcul d'une solution epsilon-contrainte
            f1RL, xf1RL = computeLinearRelax2SPA(nbvar, nbctr, A, c1, c2, maxf2RL - (j1-1) * pasSample2, 1)

            # reconditionne les valeurs 0 et 1 et arrondi les autres valeurs
            nettoyageSolution!(xf1RL)
            verbose ? @printf("fRL = %8.2f  ",round(f1RL, digits=2)) : nothing

            # recalcule la solution au regard des 2 objectifs
            z1f1RLcourant, z2f1RLcourant = evaluerSolution(xf1RL, c1, c2)
            verbose ? @printf("[ %8.2f , %8.2f ] \n", z1f1RLcourant, z2f1RLcourant) : nothing

            # maj la valeur limite sur l'objectif 2 pour la solution courante
            maxf2RLlimite = z2f1RLcourant
            if maxf2RLlimite > minf2RLlimite
                push!(L1, (tSolution{Float64})(xf1RL,[z1f1RLcourant, z2f1RLcourant]))
                push!(d.xL,z1f1RLcourant);push!(d.yL,z2f1RLcourant)
                push!(d.xLf1,z1f1RLcourant);push!(d.yLf1,z2f1RLcourant)
            end
            j1 = j1+1

        else

            # minimise f2 avec ϵ-contrainte sur f1 -----------------------------
            verbose ? @printf("  z2 %2d : ϵ = %8.2f  ", j2, maxf1RL - (j2-1) * pasSample1) : nothing # echantillonage sur z1

            # calcul d'une solution epsilon-contrainte
            f2RL, xf2RL = computeLinearRelax2SPA(nbvar, nbctr, A, c1, c2, maxf1RL - (j2-1) * pasSample1, 2)

            # reconditionne les valeurs 0 et 1 et arrondi les autres valeurs
            nettoyageSolution!(xf2RL)
            verbose ? @printf("fRL = %8.2f  ",round(f2RL, digits=2)) : nothing

            # recalcule la solution au regard des 2 objectifs
            z1f2RLcourant, z2f2RLcourant = evaluerSolution(xf2RL, c1, c2)
            verbose ? @printf("[ %8.2f , %8.2f ] \n", z1f2RLcourant, z2f2RLcourant) : nothing

            # maj la valeur limite sur l'objectif 2 pour la solution courante
            minf2RLlimite = z2f2RLcourant
            if maxf2RLlimite > minf2RLlimite
                push!(L2, (tSolution{Float64})(xf2RL,[z1f2RLcourant, z2f2RLcourant]))
                push!(d.xL,z1f2RLcourant);push!(d.yL,z2f2RLcourant)
                push!(d.xLf2,z1f2RLcourant);push!(d.yLf2,z2f2RLcourant)
            end
            j2 = j2+1
        end
    end
    return  length(L1)+length(L2), vcat(L1,reverse(L2))
end


# ==============================================================================
# predicat : verifie si une solution entiere est realisable
function isFeasible(vg::Vector{tGenerateur}, k::Int64)
    #verbose && vg[k].sFea == true ? println("   feasible") : nothing
    return (vg[k].sFea == true)
end


# ==============================================================================
# predicat : verifie si le nombre d'essai maximum a ete tente
function isFinished(trial::Int64, maxTrial::Int64)
#    verbose && trial > maxTrial ? println("   maxTrial") : nothing
    return (trial > maxTrial)
end


# ==============================================================================
# predicat : verifie si le budget de calcul maximum a ete consomme
function isTimeout(temps, maxTime)
#    verbose && time()- temps > maxTime ? println("   maxTime") : nothing
    return (time()- temps > maxTime)
end


# ==============================================================================
# elabore pC le pointeur du cone ouvert vers L

function elaborePointConeOuvertversL(vg::Vector{tGenerateur}, k::Int64, pB::tPoint, pA::tPoint)

    # recupere les coordonnees du point projete
    pC=tPoint(vg[k].sPrj.y[1], vg[k].sPrj.y[2])

    # etablit le point nadir pN au depart des points pA et pB adjacents au generateur k
    pN = tPoint( pA.x , pB.y )

#    print("Coordonnees du cone 2 : ")
#    @show pC, pN

    # retient pN si pC domine pN (afin d'ouvrir le cone)
    if pC.x < pN.x  &&  pC.y < pN.y
        # remplace pC par pN
        pC=tPoint( pA.x , pB.y )
    end

    return pC
end


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
  inSector(M, O, A, B)
=#

function inSector(M, O, A, B)

    cpAB = (A.y - O.y) * (B.x - O.x) - (A.x - O.x) * (B.y - O.y)
    cpAM = (A.y - O.y) * (M.x - O.x) - (A.x - O.x) * (M.y - O.y)
    cpBM = (B.y - O.y) * (M.x - O.x) - (B.x - O.x) * (M.y - O.y)

    if (cpAB > 0)
        if ((cpAM > 0) && (cpBM < 0))
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

function inCone(pOrg, pDeb, pFin, pCur)
    # pOrg : point origine du cone (la ou il est pointe)
    # pDeb : point depart du cone (point du rayon [pOrg,pDeb])
    # pFin : point final du cone (point du rayon [pOrg,pFin])
    # pCur : point courant a tester
    # retourne VRAI si pCur est dans le cone pDeb-pFin-pOrg, FAUX sinon

    cp_pDeb_pFin = (pDeb.x - pOrg.x) * (pFin.y - pOrg.y) - (pDeb.y - pOrg.y) * (pFin.x - pOrg.x)
    cp_pDeb_pCur = (pDeb.x - pOrg.x) * (pCur.y - pOrg.y) - (pDeb.y - pOrg.y) * (pCur.x - pOrg.x)
    cp_pFin_pCur = (pFin.x - pOrg.x) * (pCur.y - pOrg.y) - (pFin.y - pOrg.y) * (pCur.x - pOrg.x)

    if (cp_pDeb_pFin > 0)
        if ((cp_pDeb_pCur >= 0) && (cp_pFin_pCur <= 0))
            return true
        else
            return false
        end
    else
        if (!((cp_pDeb_pCur < 0) && (cp_pFin_pCur > 0)))
            return true
        else
            return false
        end
    end
end

function inCone1VersZ(pOrg, pDeb, pFin, pCur)
    return inCone(pOrg, pDeb, pFin, pCur)
end

function inCone2Vers0(pOrg, pDeb, pFin, pCur)
    return !inCone(pOrg, pDeb, pFin, pCur)
end


# ==============================================================================
# Selectionne les points pour le cone pointe sur le generateur k (pCour) et ouvert vers Y
function selectionPoints(vg::Vector{tGenerateur}, k::Int64)
    nbgen = size(vg,1)
    if k==1
        # premier generateur (point predecesseur fictif)
        pPrec = tPoint(vg[k].sRel.y[1], vg[k].sRel.y[2]+1.0)
        pCour = tPoint(vg[k].sRel.y[1], vg[k].sRel.y[2])
        pSuiv = tPoint(vg[k+1].sRel.y[1], vg[k+1].sRel.y[2])
    elseif k==nbgen
        # dernier generateur (point suivant fictif)
        pPrec = tPoint(vg[k-1].sRel.y[1], vg[k-1].sRel.y[2])
        pCour = tPoint(vg[k].sRel.y[1], vg[k].sRel.y[2])
        pSuiv = tPoint(vg[k].sRel.y[1]+1.0, vg[k].sRel.y[2])
    else
        # generateur non extreme
        pPrec = tPoint(vg[k-1].sRel.y[1], vg[k-1].sRel.y[2])
        pCour = tPoint(vg[k].sRel.y[1], vg[k].sRel.y[2])
        pSuiv = tPoint(vg[k+1].sRel.y[1], vg[k+1].sRel.y[2])
    end
#    print("Coordonnees du cone 1 : ")
#    @show pPrec, pCour, pSuiv
    return pPrec, pCour, pSuiv
end


# ==============================================================================
# arrondi la solution correspondant au generateur (pas d'historique donc)
# version avec cone inferieur seulement
function roundingSolution!(vg::Vector{tGenerateur}, k::Int64, c1::Array{Int,1}, c2::Array{Int,1}, d::tListDisplay)

    nbvar = length(vg[k].sInt.x)
    nbgen = size(vg,1)

    # --------------------------------------------------------------------------
    # identifie les variables fractionnaires et marquage par valeur sentinelle -1
    for i in 1:nbvar
        if  isapprox(vg[k].sPrj.x[i] , 0.0, atol=1e-3)
            vg[k].sInt.x[i] = 0
        elseif isapprox(vg[k].sPrj.x[i] , 1.0, atol=1e-3)
            vg[k].sInt.x[i] = 1
        else
            vg[k].sInt.x[i] = -1
        end
    end

    # --------------------------------------------------------------------------
    # Selectionne les points pour le cone pointe sur le generateur k (pCour) et ouvert vers Y
    pPrec, pCour, pSuiv = selectionPoints(vg, k)

    # --------------------------------------------------------------------------
    # elabore le point pointeur du cone ouvert vers L
    pC = elaborePointConeOuvertversL(vg, k, pPrec, pSuiv)

    # --------------------------------------------------------------------------
    # Arrondi les valeurs non-entieres d'une solution fractionnaire
    # 1) applique un arrondi qui fixe 1 variable (la premiere venue) par iteration

    vg[k].sInt.y[1] = 0
    vg[k].sInt.y[2] = 0

    # calcule en differentiel par rapport a l'image de la solution fractionnaire
    z1 = vg[k].sPrj.y[1]
    z2 = vg[k].sPrj.y[2]
    pM = tPoint( z1 , z2 )

    print("Depart a arrondir : ")
    @show pPrec, pCour, pSuiv, pM
    if length(vg[k].sPrj.x) ≤ 20 @show vg[k].sPrj.x end
    if length(vg[k].sInt.x) ≤ 20 @show vg[k].sInt.x end

    nbVarNonEntiere = 0
    for i in 1:nbvar
        if vg[k].sInt.x[i] == -1

            # la variable i est non entiere
            nbVarNonEntiere += 1

            # defalque la contribution fractionnaire de la variable aux objectifs => pose x[i]=0
            pM = tPoint( z1 - vg[k].sPrj.x[i] * c1[i] , z2 - vg[k].sPrj.x[i] * c2[i])

            if inSector(pM, pCour, pSuiv, pPrec)
                # le point pM obtenu est hors du cone => valide x[i]=1
                vg[k].sInt.x[i] = 1
                z1 = z1 - (vg[k].sPrj.x[i] * c1[i]) + c1[i]
                z2 = z2 - (vg[k].sPrj.x[i] * c2[i]) + c2[i]
                @printf("  x[%2d]=1 |  : [ %12.5f , %12.5f ] \n",i, z1, z2)
            else
                # le point pM obtenu est dans le cone => valide x[i]=0
                vg[k].sInt.x[i] = 0
                z1 = z1 - (vg[k].sPrj.x[i] * c1[i])
                z2 = z2 - (vg[k].sPrj.x[i] * c2[i])
                @printf("  x[%2d]=0 |  : [ %12.5f , %12.5f ] \n",i, z1, z2)
            end
        end

        # calcule la performance de la solution entiere
        vg[k].sInt.y[1] += vg[k].sInt.x[i] * c1[i]
        vg[k].sInt.y[2] += vg[k].sInt.x[i] * c2[i]
    end

    if length(vg[k].sInt.x) ≤ 20 @show vg[k].sInt.x end
    @printf("→ #round : %4d → [ %5d , %5d ] ", nbVarNonEntiere, vg[k].sInt.y[1], vg[k].sInt.y[2])
    push!(d.XInt,vg[k].sInt.y[1])
    push!(d.YInt,vg[k].sInt.y[2])

end


# ==============================================================================
# arrondi la solution correspondant au generateur (pas d'historique donc)
# version avec cone inferieur et superieur
function roundingSolutionnew24!(vg::Vector{tGenerateur}, k::Int64, c1::Array{Int,1}, c2::Array{Int,1}, d::tListDisplay)

    nbvar = length(vg[k].sInt.x)
    nbgen = size(vg,1)

    # --------------------------------------------------------------------------
    # identifie les variables fractionnaires et marquage par valeur sentinelle -1
    for i in 1:nbvar
        if  isapprox(vg[k].sPrj.x[i] , 0.0, atol=1e-3)
            vg[k].sInt.x[i] = 0
        elseif isapprox(vg[k].sPrj.x[i] , 1.0, atol=1e-3)
            vg[k].sInt.x[i] = 1
        else
            vg[k].sInt.x[i] = -1
        end
    end

    # --------------------------------------------------------------------------
    # Selectionne les points pour le cone pointe sur le generateur k (pCour) et ouvert vers Y
    pPrec, pCour, pSuiv = selectionPoints(vg, k)

    # --------------------------------------------------------------------------
    # elabore le point pointeur du cone ouvert vers L
    pC = elaborePointConeOuvertversL(vg, k, pPrec, pSuiv)

    # --------------------------------------------------------------------------
    # Arrondi les valeurs non-entieres d'une solution fractionnaire
    # 1) applique un arrondi qui fixe 1 variable (la premiere venue) par iteration

    vg[k].sInt.y[1] = 0
    vg[k].sInt.y[2] = 0

    # calcule en differentiel par rapport a l'image de la solution fractionnaire
    z1 = vg[k].sPrj.y[1]
    z2 = vg[k].sPrj.y[2]
    pM = tPoint( z1 , z2 )

#    print("Depart a arrondir : ")
#    @show pPrec, pCour, pSuiv, pM
#    if length(vg[k].sPrj.x) ≤ 20 @show vg[k].sPrj.x end
#    if length(vg[k].sInt.x) ≤ 20 @show vg[k].sInt.x end

    verbose ? @printf("  %2dR : [ %8.2f , %8.2f ] ", k, z1, z2) : nothing

    nbVarNonEntiere = 0
    for i in 1:nbvar
        if vg[k].sInt.x[i] == -1

            # la variable i est non entiere
            nbVarNonEntiere += 1

            # defalque la contribution fractionnaire de la variable aux objectifs => pose x[i]=0
            pM0 = tPoint( z1 - vg[k].sPrj.x[i] * c1[i] , z2 - vg[k].sPrj.x[i] * c2[i])
            pM1 = tPoint( z1 - vg[k].sPrj.x[i] * c1[i]+c1[i] , z2 - vg[k].sPrj.x[i] * c2[i]+c2[i])

            if !inCone1VersZ(pCour, pSuiv, pPrec, pM0)
                # le point pM0 obtenu est hors du cone inferieur => valide x[i]=1
                vg[k].sInt.x[i] = 1
                z1 = z1 - (vg[k].sPrj.x[i] * c1[i]) + c1[i]
                z2 = z2 - (vg[k].sPrj.x[i] * c2[i]) + c2[i]
#                @printf("  x[%2d]=1 |  : [ %12.5f , %12.5f ] \n",i, z1, z2)
            elseif inCone2Vers0(pC, pSuiv, pPrec, pM0)
                # le point pM0 obtenu est dans le cone => valide x[i]=0
                vg[k].sInt.x[i] = 0
                z1 = z1 - (vg[k].sPrj.x[i] * c1[i])
                z2 = z2 - (vg[k].sPrj.x[i] * c2[i])
#                @printf("  x[%2d]=0 |  : [ %12.5f , %12.5f ] \n",i, z1, z2)
            elseif inCone2Vers0(pC, pSuiv, pPrec, pM1)
                # le point pM1 obtenu est dans le cone => valide x[i]=1
                vg[k].sInt.x[i] = 1
                z1 = z1 - (vg[k].sPrj.x[i] * c1[i]) + c1[i]
                z2 = z2 - (vg[k].sPrj.x[i] * c2[i]) + c2[i]
#                @printf("  x[%2d]=0 |  : [ %12.5f , %12.5f ] \n",i, z1, z2)
            else
#                println("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
                # choisir entre pM0 et pM1 selon la distance a la droite [O;OO]
                O = tPoint(vg[k].sRel.y[1], vg[k].sRel.y[2]) # generateur
                OO = tPoint(vg[k].sPrj.y[1], vg[k].sPrj.y[2]) # projection
                # les points pM0 et pM1 obtenus sont hors du cone => choix de 0 ou 1 avec la plus courte distance à O et OO
                distM0 = sqrt((O.x-pM0.x)^2+(O.y-pM0.y)^2)+sqrt((OO.x-pM0.x)^2+(OO.y-pM0.y)^2)
                distM1 = sqrt((O.x-pM1.x)^2+(O.y-pM1.y)^2)+sqrt((OO.x-pM1.x)^2+(OO.y-pM1.y)^2)
                if distM0 < distM1
                    vg[k].sInt.x[i] = 0
                    z1 = z1 - (vg[k].sPrj.x[i] * c1[i])
                    z2 = z2 - (vg[k].sPrj.x[i] * c2[i])
                else
                    vg[k].sInt.x[i] = 1
                    z1 = z1 - (vg[k].sPrj.x[i] * c1[i]) + c1[i]
                    z2 = z2 - (vg[k].sPrj.x[i] * c2[i]) + c2[i]
                end
            end
        end

        # calcule la performance de la solution entiere
        vg[k].sInt.y[1] += vg[k].sInt.x[i] * c1[i]
        vg[k].sInt.y[2] += vg[k].sInt.x[i] * c2[i]
    end

    if length(vg[k].sInt.x) ≤ 20 @show vg[k].sInt.x end
    @printf("→ #round : %4d → [ %5d , %5d ] ", nbVarNonEntiere, vg[k].sInt.y[1], vg[k].sInt.y[2])
    push!(d.XInt,vg[k].sInt.y[1])
    push!(d.YInt,vg[k].sInt.y[2])

end


# ==============================================================================
# arrondi la solution correspondant au generateur (pas d'historique donc)
# version avec voisinage et selection d'un voisin selon distance L1 avec generateur
function roundingSolutionNew23!(vg::Vector{tGenerateur}, k::Int64, c1::Array{Int,1}, c2::Array{Int,1}, d::tListDisplay)

    nbvar = length(vg[k].sInt.x)
    nbgen = size(vg,1)
    lstIdxFrac =(Int64)[]

    vg[k].sInt.y[1] = 0
    vg[k].sInt.y[2] = 0

    # --------------------------------------------------------------------------
    # identifie les variables fractionnaires et marquage par valeur sentinelle -1
    for i in 1:nbvar
        if  isapprox(vg[k].sPrj.x[i] , 0.0, atol=1e-3)
            vg[k].sInt.x[i] = 0
        elseif isapprox(vg[k].sPrj.x[i] , 1.0, atol=1e-3)
            vg[k].sInt.x[i] = 1
            vg[k].sInt.y[1] += c1[i] # Comptabilise la contribution des var a 1
            vg[k].sInt.y[2] += c2[i]
        else
            vg[k].sInt.x[i] = -1
            push!(lstIdxFrac,i)
        end
    end

    # --------------------------------------------------------------------------
    # Selectionne les points pour le cone pointe sur le generateur k (pCour) et ouvert vers Y
    pPrec, pCour, pSuiv = selectionPoints(vg, k)

    # --------------------------------------------------------------------------
    # elabore le point pointeur du cone ouvert vers L
    pC = elaborePointConeOuvertversL(vg, k, pPrec, pSuiv)

    # --------------------------------------------------------------------------
    # Arrondi les valeurs non-entieres d'une solution fractionnaire
    # 1) applique un arrondi qui fixe 1 variable (la premiere venue) par iteration

    # calcule en differentiel par rapport a l'image de la solution fractionnaire
    z1 = vg[k].sPrj.y[1]
    z2 = vg[k].sPrj.y[2]
    pM = tPoint( z1 , z2 )

#    if length(vg[k].sPrj.x) ≤ 20 @show vg[k].sPrj.x end
#    if length(vg[k].sInt.x) ≤ 20 @show vg[k].sInt.x end

    nbVarNonEntiere = length(lstIdxFrac)
    # traite toutes les variables non entieres

#=
    @printf("\n\n")
    @show c1
    @show c2
    @printf("VARIABLES FIXEES : \n")
    @show vg[k].sPrj.x
    @show vg[k].sInt.x
    @show vg[k].sInt.y
    @show pM
    @show lstIdxFrac
=#
    verbose ? @printf("  %2dR : [ %8.2f , %8.2f ] ", k, z1, z2) : nothing

    # traitement de toutes les variables fractionnaires une a une
    for i= 1:nbVarNonEntiere
#        @printf("\n")

        distMin = Inf; ind = -1; val = -1
        # voisinage a 1 variable sur toutes les variables
        for iFrac in [i for i=1:nbvar if vg[k].sInt.x[i] == -1]

            # Evalue x[iFrac]=0
            # defalque la contribution fractionnaire de la variable aux objectifs => pose x[iFrac]=0
            pM = tPoint( z1 - c1[iFrac] * vg[k].sPrj.x[iFrac] , z2 - c2[iFrac] * vg[k].sPrj.x[iFrac])
            if inCone1VersZ(pCour, pSuiv, pPrec, pM)
                distL1 = abs(vg[k].sRel.y[1]-pM.x) + abs(vg[k].sRel.y[2]-pM.y)
                if distL1 < distMin
                    distMin = distL1; ind = iFrac; val = 0
                end
            end

            # Evalue x[iFrac]=1
            # ajoute la contribution entiere de la variable aux objectifs => pose x[iFrac]=1
            pM = tPoint( z1 - c1[iFrac] * vg[k].sPrj.x[iFrac] + c1[iFrac], z2 - c2[iFrac] * vg[k].sPrj.x[iFrac] + c2[iFrac])
            distL1 = abs(vg[k].sRel.y[1]-pM.x) + abs(vg[k].sRel.y[2]-pM.y)
            if distL1 < distMin
                distMin = distL1; ind = iFrac; val = 1
            end
        end

        # applique le voisin le plus interessant trouve
        vg[k].sInt.x[ind] = val
        vg[k].sInt.y[1] += vg[k].sInt.x[ind] * c1[ind] * val
        vg[k].sInt.y[2] += vg[k].sInt.x[ind] * c2[ind] * val
        z1 = z1 - (vg[k].sPrj.x[ind] * c1[ind]) + c1[ind] * val
        z2 = z2 - (vg[k].sPrj.x[ind] * c2[ind]) + c2[ind] * val
    #    @printf("  x[%2d]=%2d |  : [ %12.5f , %12.5f ] \n",ind, val, z1, z2)

    end

#    if length(vg[k].sInt.x) ≤ 20 @show vg[k].sInt.x end
    @printf("→ #round : %4d → [ %5d , %5d ] ", nbVarNonEntiere, vg[k].sInt.y[1], vg[k].sInt.y[2])
    push!(d.XInt,vg[k].sInt.y[1])
    push!(d.YInt,vg[k].sInt.y[2])

end


# ==============================================================================
# projecte la solution entiere correspondant au generateur k et test d'admissibilite
function projectingSolution!(vg::Vector{tGenerateur}, k::Int64, 
                             A::Array{Int,2}, c1::Array{Int,1}, c2::Array{Int,1}, 
                             λ1::Vector{Float64}, λ2::Vector{Float64},
                             d::tListDisplay)

    # --------------------------------------------------------------------------
    # Projete la solution entiere sur le polytope X avec norme-L1

#    fPrj, vg[k].sPrj.x = Δ2SPA(A,vg[k].sInt.x)
    fPrj, vg[k].sPrj.x = Δ2SPAbis(A,vg[k].sInt.x,c1,c2,k,λ1,λ2)

    # Nettoyage de la valeur de vg[k].sPrj.x et calcul du point bi-objectif
    # reconditionne les valeurs 0 et 1 et arrondi les autres valeurs
    nettoyageSolution!(vg[k].sPrj.x)
#    verbose ? @printf("  %2dP : fPrj = %8.2f  ",k, round(fPrj, digits=2)) : nothing

    # recalcule la solution au regard des 2 objectifs
    vg[k].sPrj.y[1], vg[k].sPrj.y[2] = evaluerSolution(vg[k].sPrj.x, c1, c2)
    verbose ? @printf("  %2dP : [ %8.2f , %8.2f ] ",k, vg[k].sPrj.y[1], vg[k].sPrj.y[2]) : nothing

    push!(d.XProj, vg[k].sPrj.y[1])
    push!(d.YProj, vg[k].sPrj.y[2])

    # ----------------------------------------------------------------
    # Teste si la projection est admissible

    if estAdmissible(vg[k].sPrj.x)
        vg[k].sInt.x = deepcopy(vg[k].sPrj.x)
        vg[k].sInt.y[1] = vg[k].sPrj.y[1]
        vg[k].sInt.y[2] = vg[k].sPrj.y[2]
        vg[k].sFea = true
        push!(d.XFeas, vg[k].sPrj.y[1])
        push!(d.YFeas, vg[k].sPrj.y[2])
        @printf("→ Admissible "); print("                       ")
    else
        vg[k].sFea = false
        @printf("→ x          "); print("                       ")
        # prepare pour l'iteration suivante
#        vg[k].xRlx = deepcopy(vg[k].sPrj.x) !!!!!!!!!!!!!
    end

end


# ==============================================================================
# Calcule la direction d'interet du nadir vers le milieu de segment reliant deux points generateurs
function calculerDirections(L::Vector{tSolution{Float64}}, vg::Vector{tGenerateur})
   # function calculerDirections(L, vg::Vector{tGenerateur})    

    nbgen = size(vg,1)
    for k in 2:nbgen

        n1 = L[end].y[1]
        n2 = L[1].y[2]

        x1,y1 = vg[k-1].sRel.y[1], vg[k-1].sRel.y[2]
        x2,y2 = vg[k].sRel.y[1], vg[k].sRel.y[2]
        xm=(x1+x2)/2.0
        ym=(y1+y2)/2.0
        Δx = abs(n1-xm)
        Δy = abs(n2-ym)
        λ1 =  1 - Δx / (Δx+Δy)
        λ2 =  1 - Δy / (Δx+Δy)
        @printf("  x1= %7.2f   y1= %7.2f \n",x1,y1)
        @printf("  x2= %7.2f   y2= %7.2f \n",x2,y2)
        @printf("  Δx= %7.2f    Δy= %7.2f \n",Δx,Δy)
        @printf("  λ1= %6.5f    λ2= %6.5f \n",λ1,λ2)
        plot(n1, n2, xm, ym, linestyle="-", color="blue", marker="+")
        annotate("",
                 xy=[xm;ym],# Arrow tip
                 xytext=[n1;n2], # Text offset from tip
                 arrowprops=Dict("arrowstyle"=>"->"))
        println("")
    end

end


# ==============================================================================
# Calcule la direction d'interet du nadir vers un point generateur
function calculerDirections2(L::Vector{tSolution{Float64}}, vg::Vector{tGenerateur})
    #function calculerDirections2(L, vg::Vector{tGenerateur})

    nbgen = size(vg,1)
    λ1=Vector{Float64}(undef, nbgen)
    λ2=Vector{Float64}(undef, nbgen)
    for k in 1:nbgen

        n1 = L[end].y[1]
        n2 = L[1].y[2]

        xm=vg[k].sRel.y[1]
        ym=vg[k].sRel.y[2]
        Δx = abs(n1-xm)
        Δy = abs(n2-ym)
        λ1[k] =  1 - Δx / (Δx+Δy)
        λ2[k] =  1 - Δy / (Δx+Δy)
        @printf("  k= %3d   ",k)
        @printf("  xm= %7.2f   ym= %7.2f ",xm,ym)
        @printf("  Δx= %8.2f    Δy= %8.2f ",Δx,Δy)
        @printf("  λ1= %6.5f    λ2= %6.5f \n",λ1[k],λ2[k])
        plot(n1, n2, xm, ym, linestyle="-", color="blue", marker="+")
        annotate("",
                 xy=[xm;ym],# Arrow tip
                 xytext=[n1;n2], # Text offset from tip
                 arrowprops=Dict("arrowstyle"=>"->"))
        #println("")
    end
    return λ1, λ2
end


# ==============================================================================
# point d'entree principal

function GM( fname::String,
             tailleSampling::Int64,
             maxTrial::Int64,
             maxTime::Int64
           )

    @assert tailleSampling>=3 "Erreur : Au moins 3 sont requis"  

    @printf("0) instance et parametres \n\n")
    verbose ? println("  instance = $fname | tailleSampling = $tailleSampling | maxTrial = $maxTrial | maxTime = $maxTime\n\n") : nothing

    # chargement de l'instance numerique ---------------------------------------
    c1, c2, A = loadInstance2SPA(fname) # instance numerique de SPA
    nbctr = size(A,1)
    nbvar = size(A,2)
    nbobj = 2

    # structure pour les points qui apparaitront dans l'affichage graphique
    d = tListDisplay([],[], [],[], [],[], [],[], [],[], [],[], [],[])

    # --------------------------------------------------------------------------
    # --------------------------------------------------------------------------
    @printf("1) calcule les etendues de valeurs sur les 2 objectifs\n\n")

    # calcule la valeur optimale relachee de f1 seule et le point (z1,z2) correspondant
    f1RL, xf1RL = computeLinearRelax2SPA(nbvar, nbctr, A, c1, c2, typemax(Int), 1) # opt fct 1
    minf1RL, maxf2RL = evaluerSolution(xf1RL, c1, c2)

    # calcule la valeur optimale relachee de f2 seule et le point (z1,z2) correspondant
    f2RL, xf2RL = computeLinearRelax2SPA(nbvar, nbctr, A, c1, c2, typemax(Int), 2) # opt fct 2
    maxf1RL, minf2RL = evaluerSolution(xf2RL, c1, c2)

    verbose ? @printf("  f1_min=%8.2f ↔ f1_max=%8.2f (Δ=%.2f) \n",minf1RL, maxf1RL, maxf1RL-minf1RL) : nothing
    verbose ? @printf("  f2_min=%8.2f ↔ f2_max=%8.2f (Δ=%.2f) \n\n",minf2RL, maxf2RL, maxf2RL-minf2RL) : nothing


    # --------------------------------------------------------------------------
    # --------------------------------------------------------------------------
    @printf("2) calcule les generateurs par e-contrainte alternant minimiser z1 et z2\n\n")

    nbgen, L = calculGenerateurs(A, c1, c2, tailleSampling, minf1RL, maxf2RL, maxf1RL, minf2RL, d)

    # --------------------------------------------------------------------------
    # --------------------------------------------------------------------------
    # allocation de memoire pour la structure de donnees -----------------------

    vg = allocateDatastructure(nbgen, nbvar, nbobj)

    # --------------------------------------------------------------------------
    # --------------------------------------------------------------------------
    @printf("3) place L dans structure et verifie l'admissibilite de chaque generateur\n\n")

    for k=1:nbgen

        verbose ? @printf("  %2d  : [ %8.2f , %8.2f ] ", k, L[k].y[1], L[k].y[2]) : nothing

        # copie de l'ensemble bornant inferieur dans la stru de donnees iterative ---
        ajouterX0!(vg, k, L[k])

        # test d'admissibilite et marquage de la solution le cas echeant -------
        if estAdmissible(vg[k].sRel.x)
            ajouterXtilde!(vg, k, convert.(Int, vg[k].sRel.x), convert.(Int, L[k].y))
            vg[k].sFea   = true
            push!(d.XFeas,vg[k].sInt.y[1])
            push!(d.YFeas,vg[k].sInt.y[2])
            verbose ? @printf("→ Admissible \n") : nothing
        else
            vg[k].sFea   = false
            verbose ? @printf("→ x          \n") : nothing
        end

    end
    verbose ? println("") : nothing

    # --------------------------------------------------------------------------
    # --------------------------------------------------------------------------
    # Sortie graphique

    figure("Gravity Machine",figsize=(6.5,5))
    #xlim(25000,45000)
    #ylim(20000,40000)
    xlabel(L"z^1(x)")
    ylabel(L"z^2(x)")
    PyPlot.title("Cone | 1 rounding | 2-$fname")

    # --------------------------------------------------------------------------
    # --------------------------------------------------------------------------
    # calcule les directions (λ1,λ2) pour chaque generateur a utiliser lors des projections
    λ1,λ2 = calculerDirections2(L,vg)

    # ==========================================================================

    @printf("4) terraformation generateur par generateur \n\n")

    for k in [i for i in 1:nbgen if !isFeasible(vg,i)]
        temps = time()
        trial = 0
        H =(Vector{Int64})[]

#perturbSolution30!(vg,k,c1,c2,d)

        # rounding solution : met a jour sInt dans vg --------------------------
        #roundingSolution!(vg,k,c1,c2,d)  # un cone
        #roundingSolutionnew24!(vg,k,c1,c2,d) # deux cones
        roundingSolutionNew23!(vg,k,c1,c2,d) # un cone et LS sur generateur

        push!(H,[vg[k].sInt.y[1],vg[k].sInt.y[2]])
        println("   t=",trial,"  |  Tps=", round(time()- temps, digits=4))

        while !(t1=isFeasible(vg,k)) && !(t2=isFinished(trial, maxTrial)) && !(t3=isTimeout(temps, maxTime))

            trial+=1

            # projecting solution : met a jour sPrj, sInt, sFea dans vg --------
            projectingSolution!(vg,k,A,c1,c2,λ1,λ2,d)
            println("   t=",trial,"  |  Tps=", round(time()- temps, digits=4))

            if !isFeasible(vg,k)

                # rounding solution : met a jour sInt dans vg --------------------------
                #roundingSolution!(vg,k,c1,c2,d)
                #roundingSolutionnew24!(vg,k,c1,c2,d)
                roundingSolutionNew23!(vg,k,c1,c2,d)
                println("   t=",trial,"  |  Tps=", round(time()- temps, digits=4))

                # test detection cycle sur solutions entieres ------------------
                cycle = [vg[k].sInt.y[1],vg[k].sInt.y[2]] in H
                if (cycle == true)
                    println("CYCLE!!!!!!!!!!!!!!!")
                    # perturb solution
                    perturbSolution30!(vg,k,c1,c2,d)
                end
                push!(H,[vg[k].sInt.y[1],vg[k].sInt.y[2]])

            end
        end
        if t1
            println("   feasible \n")
        elseif t2
            println("   maxTrial \n")
        elseif t3
            println("   maxTime \n")
        end


    end

    println("");

    # ==========================================================================

    @printf("5) Extraction des resultats\n\n")


    for k=1:nbgen
        verbose ? @printf("  %2d  : [ %8.2f , %8.2f ] ", k, vg[k].sInt.y[1],vg[k].sInt.y[2]) : nothing
        # test d'admissibilite et marquage de la solution le cas echeant -------
        if vg[k].sFea
            verbose ? @printf("→ Admissible \n") : nothing
        else
            verbose ? @printf("→ x          \n") : nothing
        end
    end

    # allocation de memoire pour les ensembles bornants ------------------------
    U = Vector{tSolution{Int64}}(undef,nbgen)
    for j = 1:nbgen
        U[j] = tSolution{Int64}(zeros(Int64,nbvar),zeros(Int64,nbobj))
    end
    #--> TODO : stocker l'EBP dans U proprement


    # ==========================================================================
    @printf("6) Edition des resultats \n\n")

#    figure("Gravity Machine",figsize=(6.5,5))
    #xlim(25000,45000)
    #ylim(20000,40000)
#    xlabel(L"z^1(x)")
#    ylabel(L"z^2(x)")
    # Donne les points relaches initiaux ---------------------------------------
#    scatter(d.xLf1,d.yLf1,color="blue", marker="x")
#    scatter(d.xLf2,d.yLf2,color="red", marker="+")
    graphic ? scatter(d.xL,d.yL,color="blue", marker="x", label = L"y \in L") : nothing

    # Donne les points entiers -------------------------------------------------
    graphic ? scatter(d.XInt,d.YInt,color="orange", marker="s", label = L"y"*" rounded") : nothing
#    @show d.XInt
#    @show d.YInt

    # Donne les points apres projection Δ(x,x̃) ---------------------------------
    graphic ? scatter(d.XProj,d.YProj, color="red", marker="x", label = L"y"*" projected") : nothing
#    @show d.XProj
#    @show d.YProj

    # Donne les points admissibles ---------------------------------------------
    graphic ? scatter(d.XFeas,d.YFeas, color="green", marker="o", label = L"y \in F") : nothing
#    @show d.XFeas
#    @show d.YFeas

    # Donne l'ensemble bornant primal obtenu + la frontiere correspondante -----
    #--> TODO : stocker l'EBP dans U proprement  
    X_EBP_frontiere, Y_EBP_frontiere, X_EBP, Y_EBP = ExtractEBP(d.XFeas, d.YFeas)
    plot(X_EBP_frontiere, Y_EBP_frontiere, color="green", markersize=3.0, marker="x")
    scatter(X_EBP, Y_EBP, color="green", s = 150, alpha = 0.3, label = L"y \in U")  
    @show X_EBP
    @show Y_EBP 

    # Donne les points qui ont fait l'objet d'une perturbation -----------------
     scatter(d.XPert,d.YPert, color="magenta", marker="s", label ="pertub")

    # Donne les points non-domines exacts de cette instance --------------------
     XN,YN = loadNDPoints2SPA(fname)
     plot(XN, YN, color="black", linewidth=0.75, marker="+", markersize=1.0, linestyle=":", label = L"y \in Y_N")
     scatter(XN, YN, color="black", marker="+")
#    @show X_Y_N
#    @show Y_Y_N

    # Affiche le cadre avec les legendes des differents traces -----------------
    legend(bbox_to_anchor=[1,1], loc=0, borderaxespad=0, fontsize = "x-small")
    #PyPlot.title("Cone | 1 rounding | 2-$fname")

    #--> TODO : integrer la mesure de qualite de EBP

end

# ==============================================================================

#testidee()
#@time GM("sppaa02.txt", 6, 20, 20)
#@time GM("sppnw03.txt", 6, 20, 20) #pb glpk
#@time GM("sppnw10.txt", 6, 20, 20)
@time GM("didactic5.txt", 5, 5, 10)
#@time GM("sppnw29.txt", 6, 30, 20)
nothing
