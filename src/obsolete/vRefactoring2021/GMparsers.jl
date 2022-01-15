# ==============================================================================
# Parseur lisant une instance de probleme de partitionnement d'ensembles (SPA) bi-objectif

function loadInstance2SPA(fname::String)           
    f = open(fname)
    nbcontraintes, nbvar = parse.(Int, split(readline(f))) # nombre de contraintes , nombre de variables
    L = zeros(Int, nbcontraintes, nbvar)    # matrice des contraintes
    c1 = zeros(Int, nbvar)                  # vecteur des couts
    c2 = zeros(Int, nbvar)                  # deuxième vecteur des couts
    nb = zeros(Int, nbvar)
    for i in 1:nbvar
        flag = 1
        for valeur in split(readline(f))
            if flag == 1
                c1[i] = parse(Int, valeur)
                flag +=1
            elseif flag == 2
                c2[i] = parse(Int, valeur)
                flag +=1
            elseif flag == 3
                nb[i] = parse(Int, valeur)
                flag +=1
            else
                j = parse(Int, valeur)
                L[j,i] = 1
            end
        end
    end
    close(f)
    return nbvar, nbcontraintes, L, c1, c2
end