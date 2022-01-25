# Create an instance of bi-objective set partitionning problem (2-SPA)


# ---- Parser reading an instance of SPA (format of instances compliant with OR-library)
function loadInstanceSPA(fname::String)

    f = open(fname)
    nbctr, nbvar = parse.(Int, split(readline(f)))
    L = zeros(Int, nbctr, nbvar)       # matrix of constraints
    c = zeros(Int, nbvar)              # vector of costs
    nb = zeros(Int, nbvar)
    for i in 1:nbvar
        flag = 1
        for valeur in split(readline(f))
            if flag == 1
                c[i] = parse(Int, valeur)
                flag +=1
            elseif flag == 2
                nb[i] = parse(Int, valeur)
                flag +=1
            else
                j = parse(Int, valeur)
                L[j,i] = 1
            end
        end
    end
    close(f)
    return nbvar, nbctr, L, c, nb
end


# ---- Create and save the 2-SPA instance (format of instances compliant with vOptLib)
function CreateSave2SPA(fname::String)

    nbvar, nbctr, L, c, nb = loadInstanceSPA(fname)
    cc = copy(c)
    way = "bio"*fname
    println(way)
    open(way, "w") do f
        write(f, "$nbctr $nbvar \n")
        for i in 1:nbvar
            randi = rand(1:length(cc))
            write(f, "$(c[i]) $(cc[randi])")
            write(f, " $(nb[i])")
            splice!(cc, randi)
            for j in 1:nbctr
                if L[j,i] == 1
                    write(f, " $j")
                end
            end
            write(f, "\n")
        end
    end
end


# ---- Entry point
function main(fname::String)

    print("Create an instance of bi-objective set partitionning problem (2-SPA) \n\n")
    CreateSave2SPA(fname)

end

# ---- Example

#main("didacticSPA.txt")
