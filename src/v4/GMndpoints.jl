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