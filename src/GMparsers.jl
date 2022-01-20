# ==============================================================================
# Parseur lisant une instance de probleme de partitionnement d'ensembles bi-objectif (2-SPA)

function loadInstance2SPA(fname::String)

    f = open("../SPA/instances/"*"bio"*fname)
    nbctr, nbvar = parse.(Int, split(readline(f))) # nombre de contraintes , nombre de variables
    A = zeros(Int, nbctr, nbvar)                   # matrice des contraintes
    c1 = zeros(Int, nbvar)                         # vecteur des couts
    c2 = zeros(Int, nbvar)                         # deuxième vecteur des couts
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
                A[j,i] = 1
            end
        end
    end
    close(f)
    return c1, c2, A
end


# ==============================================================================
# Parseur lisant un fichier de Y_N (dans dossier SPA/Y) pour le 2-SPA traite

function loadNDPoints2SPA(fname::String)
    fname = "../SPA/Y/Y_N_"*fname
    f = open(fname)
    temps = parse.(Float64, readline(f))
    len = parse.(Int, readline(f))
    x = Vector{Float64}(undef, len)
    y = Vector{Float64}(undef, len)
    for i in 1:len
        x[i], y[i] =  parse.(Float64, split(readline(f)))
    end
    close(f)
    return x, y
end