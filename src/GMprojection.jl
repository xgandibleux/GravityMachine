# ==============================================================================
# Projete xTilde sur le polyedre X du SPA avec norme-L1
# version FP 2005

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
# Projete xTilde sur le polyedre X du SPA avec norme-L1
# version avec somme ponderee donnant la direction vers le generateur k

function Δ2SPAbis(A::Array{Int,2}, xTilde::Array{Int,1}, 
                  c1::Array{Int,1}, c2::Array{Int,1}, k::Int64, λ1::Vector{Float64}, λ2::Vector{Float64})

    nbctr = size(A,1)
    nbvar = size(A,2)
    idxTilde0, idxTilde1 = split01(xTilde)

    cλ = λ1[k].*c1 + λ2[k].*c2
    proj = Model(GLPK.Optimizer)
    @variable(proj, 0.0 <= x[1:length(xTilde)] <= 1.0 )
#    @objective(proj, Min, sum(λ1[k]*x[i] for i in idxTilde0) + sum(λ2[k]*(1-x[i]) for i in idxTilde1) )
    @objective(proj, Min, sum(cλ[i]*x[i] for i in idxTilde0) + sum(cλ[i]*(1-x[i]) for i in idxTilde1) )
    @constraint(proj, [i=1:nbctr],(sum((x[j]*A[i,j]) for j in 1:nbvar)) == 1)
    optimize!(proj)
    return objective_value(proj), value.(x)
end

# ==============================================================================
# Projete xTilde sur le polyedre X
# version avec somme ponderee donnant la direction vers le milieu de segment reliant deux points generateurs

function Δ2SPAbis2(A::Array{Int,2}, vg::Vector{tGenerateur},
    c1::Array{Int,1}, c2::Array{Int,1}, k::Int64, λ1::Vector{Float64}, λ2::Vector{Float64})

    xTilde = vg[k].sInt.x
    nbctr = size(A, 1)
    nbvar = size(A, 2)
    idxTilde0, idxTilde1 = split01(xTilde)
    
    if k == length(vg)
        k -= 1
    end

    proj = Model(GLPK.Optimizer)
    @variable(proj, 0.0 <= x[1:length(xTilde)] <= 1.0)
#    @objective(proj, Min, sum(λ1[k] * x[i] for i in idxTilde0) + sum(λ2[k] * (1 - x[i]) for i in idxTilde1))
    cλ = λ1[k]*c1 + λ2[k]*c2
    @objective(proj, Min, sum(cλ[i]*x[i] for i in idxTilde0) + sum(cλ[i]*(1-x[i]) for i in idxTilde1) )
    @constraint(proj, [i = 1:nbctr], (sum((x[j] * A[i, j]) for j in 1:nbvar)) == 1)
    optimize!(proj)
    return objective_value(proj), value.(x)
end

# ==============================================================================
# projecte la solution entiere correspondant au generateur k et test d'admissibilite
function projectingSolution!(vg::Vector{tGenerateur}, k::Int64, 
                             A::Array{Int,2}, c1::Array{Int,1}, c2::Array{Int,1}, 
                             λ1::Vector{Float64}, λ2::Vector{Float64},
                             d::tListDisplay)

    # --------------------------------------------------------------------------
    # Projete la solution entiere sur le polytope X 

#    fPrj, vg[k].sPrj.x = Δ2SPA(A,vg[k].sInt.x)
    fPrj, vg[k].sPrj.x = Δ2SPAbis(A,vg[k].sInt.x,c1,c2,k,λ1,λ2)
#    fPrj, vg[k].sPrj.x = Δ2SPAbis2(A,vg,c1,c2,k,λ1,λ2)

    # Nettoyage de la valeur de vg[k].sPrj.x et calcul du point bi-objectif
    # reconditionne les valeurs 0 et 1 et arrondi les autres valeurs
    nettoyageSolution!(vg[k].sPrj.x)
#    verbose ? @printf("  %2dP : fPrj = %8.2f  ",k, round(fPrj, digits=2)) : nothing

    # recalcule la solution au regard des 2 objectifs
    vg[k].sPrj.y[1], vg[k].sPrj.y[2] = evaluerSolution(vg[k].sPrj.x, c1, c2)
    verbose ? @printf("  %2dP : [ %8.2f , %8.2f ] ",k, vg[k].sPrj.y[1], vg[k].sPrj.y[2]) : nothing

    # archive le point obtenu pour les besoins d'affichage
    if generateurVisualise == -1 
        # archivage pour tous les generateurs
        push!(d.XProj, vg[k].sPrj.y[1])
        push!(d.YProj, vg[k].sPrj.y[2])
    elseif generateurVisualise == k
        # archivage seulement pour le generateur k
        push!(d.XProj, vg[k].sPrj.y[1])
        push!(d.YProj, vg[k].sPrj.y[2])
    end            

    # ----------------------------------------------------------------
    # Teste si la projection est admissible

    if estAdmissible(vg[k].sPrj.x)

        # sauvegarde de la solution entiere admissible obtenue
        vg[k].sInt.x = deepcopy(vg[k].sPrj.x)
        vg[k].sInt.y[1] = vg[k].sPrj.y[1]
        vg[k].sInt.y[2] = vg[k].sPrj.y[2]
        vg[k].sFea = true
        @printf("→ Admissible "); print("                       ")

        # archive le point obtenu pour les besoins d'affichage
        if generateurVisualise == -1 
            # archivage pour tous les generateurs
            push!(d.XFeas, vg[k].sPrj.y[1])
            push!(d.YFeas, vg[k].sPrj.y[2])
        elseif generateurVisualise == k
            # archivage seulement pour le generateur k
            push!(d.XFeas, vg[k].sPrj.y[1])
            push!(d.YFeas, vg[k].sPrj.y[2])
        end  

    else

        vg[k].sFea = false
        @printf("→ x          "); print("                       ")
        # prepare pour l'iteration suivante
#        vg[k].xRlx = deepcopy(vg[k].sPrj.x) !!!!!!!!!!!!!

    end

end
