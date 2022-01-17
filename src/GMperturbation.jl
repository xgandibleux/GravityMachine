# ==============================================================================
# applique une perturbation sur la solution entiere faisant l'objet d'un cycle

function perturbSolution!(vg::Vector{tGenerateur}, k::Int64, c1::Array{Int,1}, c2::Array{Int,1}, d::tListDisplay)
    # nombre de variables maximum a considerer
    T = 20
    # nombre effectif de variables a flipper
    TT = rand(T/2:3*T/2)

    # liste des candidats (valeur, indice) et tri decroissant
    nbvar = length(vg[k].sInt.x)
    candidats=[( abs( vg[k].sPrj.x[i] - vg[k].sInt.x[i] ) , i ) for i=1:nbvar]
    sort!(candidats, rev=true, by = x -> x[1])
#    sort!(candidats,  by = x -> x[1])

#@show vg[k].sPrj.x
#@show vg[k].sInt.x
#@show candidats
#@show TT
#@show nbvar

    i = 1
    while (i<= nbvar) && (i<=TT)
        j=candidats[i][2]
#        @show candidats[i][2]
        if vg[k].sInt.x[j] == 0
            vg[k].sInt.x[j] = 1
            vg[k].sInt.y[1] = vg[k].sInt.y[1] + c1[j]
            vg[k].sInt.y[2] = vg[k].sInt.y[2] + c2[j]
        else
            vg[k].sInt.x[j] = 0
            vg[k].sInt.y[1] = vg[k].sInt.y[1] - c1[j]
            vg[k].sInt.y[2] = vg[k].sInt.y[2] - c2[j]
        end
        i+=1
    end
    @printf("  %2dC : [ %8.2f , %8.2f ] \n", k, vg[k].sInt.y[1], vg[k].sInt.y[2])
    push!(d.XPert,vg[k].sInt.y[1]); push!(d.YPert,vg[k].sInt.y[2])
#    @show vg[k].sInt.x
    return nothing
end


# ==============================================================================
# applique une perturbation sur la solution entiere faisant l'objet d'un cycle

function perturbSolution30!(vg::Vector{tGenerateur}, k::Int64, c1::Array{Int,1}, c2::Array{Int,1}, d::tListDisplay)

    # liste des candidats (valeur, indice) et tri decroissant
    nbvar = length(vg[k].sInt.x)
    idxTilde0, idxTilde1 = split01(vg[k].sInt.x)

#    candidats=[( abs( vg[k].sPrj.x[i] - vg[k].sInt.x[i] ) , i ) for i=1:nbvar if vg[k].sPrj.x[i]>0 && vg[k].sPrj.x[i]<1]
    candidats=[( vg[k].sPrj.x[i] , i ) for i=1:nbvar if vg[k].sPrj.x[i]>0 && vg[k].sPrj.x[i]<1]
#    sort!(candidats, rev=true, by = x -> x[1])


#@show vg[k].sPrj.x
#@show vg[k].sInt.x
#@show candidats
#@show nbvar

seq = randperm(length(candidats)) # melange les candidats afin d'avoir une composante variee
etat = vg[k].sInt.x[ candidats[seq[1]][2] ] # etat 0 ou 1 de la premiere variable candidate
for i = 1:length(candidats)
    j=candidats[seq[i]][2]
    #@show candidats[seq[i]][2]
    if etat == 0
        if vg[k].sInt.x[j] == 0
            vg[k].sInt.x[j] = 1
            vg[k].sInt.y[1] = vg[k].sInt.y[1] + c1[j]
            vg[k].sInt.y[2] = vg[k].sInt.y[2] + c2[j]
        end
    else
        if vg[k].sInt.x[j] == 1
            vg[k].sInt.x[j] = 0
            vg[k].sInt.y[1] = vg[k].sInt.y[1] - c1[j]
            vg[k].sInt.y[2] = vg[k].sInt.y[2] - c2[j]
        end
    end
    etat=(etat+1)%2
end

    @printf("  %2dC : [ %8.2f , %8.2f ] \n", k, vg[k].sInt.y[1], vg[k].sInt.y[2])
    push!(d.XPert,vg[k].sInt.y[1]); push!(d.YPert,vg[k].sInt.y[2])
#    @show vg[k].sInt.x

    return nothing
end
