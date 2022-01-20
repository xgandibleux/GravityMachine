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