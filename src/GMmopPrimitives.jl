# ==============================================================================
# Algorithme de Kung (extrait S_N d'un ensemble statique de points S de IR^2)

function getNonDominatedPoints(XFeas, YFeas)
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


# ==============================================================================
# Extraction de l'ensemble bornant primal + sa frontiere de l'ensemble des points realisables

function ExtractEBP(XFeas, YFeas)
    X_EBP_frontiere = (Int64)[];   Y_EBP_frontiere = (Int64)[]
    X_EBP = (Int64)[] ;            Y_EBP = (Int64)[]

    if length(XFeas) > 0
        SN = getNonDominatedPoints(XFeas, YFeas)

        push!(X_EBP_frontiere, SN[1][1]) ;   
        push!(Y_EBP_frontiere, ceil(Int64, 1.1 * maximum(YFeas)))

        for i in 1:length(SN)-1
            push!(X_EBP_frontiere, SN[i][1]);     push!(Y_EBP_frontiere, SN[i][2])
            push!(X_EBP_frontiere, SN[i+1][1]);   push!(Y_EBP_frontiere, SN[i][2])
            push!(X_EBP, SN[i][1]);               push!(Y_EBP, SN[i][2])
        end
        push!(X_EBP_frontiere, SN[end][1]);       push!(Y_EBP_frontiere, SN[end][2])

        push!(X_EBP_frontiere, ceil(Int64, 1.1 * maximum(XFeas)))
        push!(Y_EBP_frontiere, SN[end][2])

        push!(X_EBP, SN[end][1]);   push!(Y_EBP, SN[end][2])
    end
    return X_EBP_frontiere, Y_EBP_frontiere,   X_EBP, Y_EBP
end