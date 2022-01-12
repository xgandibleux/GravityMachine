using JuMP, GLPK, PyPlot, Printf #, vOptGeneric
# https://github.com/fuofuo-gg/TER

#const GUROBI_ENV = Gurobi.Env()

# ==============================================================================
function loaddemandebi(fname::String)           #parseur lisant le fichier text entrant
    f = open(fname)
    nbcontraintes, nbvar = parse.(Int, split(readline(f))) # nombre de contraintes , nombre de variables
    L = zeros(Int, nbcontraintes, nbvar)   #matrice des contraintes
    c1 = zeros(Int, nbvar)                  #vecteur des couts
    c2 = zeros(Int, nbvar)                  #deuxième vecteur des couts
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

mutable struct point
    x::Float64
    y::Float64
end

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
# Model solution relaxation linéaire sur un objectif
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
        if round(x[i], digits=3)!=0.0 && round(x[i], digits=3)!=1.0
            admissible = false
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



function mainXG4(fname::String, tailleSampling::Int64, terraform::Int64)

    @printf("Running the gravity machine...\n\n")
    @assert tailleSampling>=3 "Erreur : Au moins 3 sont requis"

    # chargement de l'instance numerique ---------------------------------------
    nbvar, nbcontraintes, L, c1, c2 = loaddemandebi(fname) # instance numerique de SPA
    nbobj = 2
    nbgen = tailleSampling

    # allocation de memoire pour la structure de donnees -----------------------
    vg=Vector{tGenerateur}(undef,2*nbgen-2)
    for j=1:2*nbgen-2
        vg[j]=tGenerateur( zeros(Float64,nbvar), zeros(Float64,nbobj),
                           zeros(Int64,nbvar),   zeros(Int64,nbobj),
                           zeros(Float64,nbvar), zeros(Float64,nbobj),
                           false
                         )
    end

    # allocation de memoire pour les ensembles bornants ------------------------
    eb=Vector{tEnsBornant}(undef,2*nbgen-2)
    for j=1:2*nbgen-2
        eb[j]=tEnsBornant( zeros(Float64,nbvar), zeros(Float64,nbobj)
                         )
    end
    ebtmp1=Vector{tEnsBornant}(undef,nbgen)
    ebtmp2=Vector{tEnsBornant}(undef,nbgen)
    for j=1:nbgen
        ebtmp1[j]=tEnsBornant( zeros(Float64,nbvar), zeros(Float64,nbobj)
                             )
        ebtmp2[j]=tEnsBornant( zeros(Float64,nbvar), zeros(Float64,nbobj)
                             )
    end

    # listes de points pour les affichages graphiques --------------------------
    XEBD=(Float64)[]; YEBD=(Float64)[]    # liste des points (x,y) relaches
    XInt=(Int64)[];  YInt=(Int64)[]       # liste des points (x,y) entiers
    XProj=(Float64)[]; YProj=(Float64)[]  # liste des points (x,y) projetes
    XFeas=(Int64)[]; YFeas=(Int64)[]      # liste des points (x,y) admissibles


    # --------------------------------------------------------------------------
    # --------------------------------------------------------------------------
    @printf("1) calcule les etendues de valeurs sur les 2 objectifs\n\n")

    # calcule la valeur optimale relachee de f1 seule et le point (z1,z2) correspondant
    f1RL, xf1RL = relaxLinXG(nbvar, nbcontraintes, L, c1, c2, typemax(Int), 0, 1) # opt fct 1
    z1RL1, z2RL1 = evaluerSolution(xf1RL, c1, c2)

    # calcule la valeur optimale relachee de f2 seule et le point (z1,z2) correspondant
    f2RL, xf2RL = relaxLinXG(nbvar, nbcontraintes, L, c1, c2, typemax(Int), 0, 2) # opt fct 2
    z1RL2, z2RL2 = evaluerSolution(xf2RL, c1, c2)

    @printf("  f1_min=%8.2f ↔ f1_max=%8.2f (Δ=%.2f) \n",z1RL1, z1RL2, z1RL2-z1RL1)
    @printf("  f2_min=%8.2f ↔ f2_max=%8.2f (Δ=%.2f) \n\n",z2RL2, z2RL1, z2RL1-z2RL2)


    # --------------------------------------------------------------------------
    # --------------------------------------------------------------------------
    @printf("2) calcule les generateurs par e-contrainte avec z1 et ensuite z2 a minimiser\n\n")

    # Premiere generation de points (avec z1 la fonction a minimiser) ----------
    pasSample = (z2RL1 - z2RL2) / (tailleSampling-1) # pas de l'echantillonage sur z2
    for j=1:tailleSampling
        @printf("  %2d : ϵ = %8.2f  ", j, z2RL1 - (j-1) * pasSample) # echantillonage sur z2

        # calcul d'une solution epsilon-contrainte -----------------------------
        fRL, ebtmp1[j].xEBD = relaxLinXG4(nbvar, nbcontraintes, L, c1, c2, z2RL1 - (j-1) * pasSample, 1)
        @printf("fRL = %8.2f  ",round(fRL, digits=2))

        # nettoyage de la valeur de ebtmp1[j].xEBD et calcul du point bi-objectif --
        for i in 1:nbvar
            if round(ebtmp1[j].xEBD[i], digits=3) == 0.0
                ebtmp1[j].xEBD[i] = 0.0
            elseif round(ebtmp1[j].xEBD[i], digits=3) == 1.0
                ebtmp1[j].xEBD[i] = 1.0
            else
                ebtmp1[j].xEBD[i] = round(ebtmp1[j].xEBD[i], digits=3)
            end
        end
        ebtmp1[j].yEBD[1], ebtmp1[j].yEBD[2] = evaluerSolution(ebtmp1[j].xEBD, c1, c2)
        @printf("[ %8.2f , %8.2f ] \n", ebtmp1[j].yEBD[1], ebtmp1[j].yEBD[2])
    end

    println("")

    # Seconde generation de points (avec z2 la fonction a minimiser) ----------
    pasSample = (z1RL2-z1RL1) / (tailleSampling-1) # pas de l'echantillonage sur z1
    for j in reverse(2:tailleSampling-1)
        @printf("  %2d : ϵ = %8.2f  ", tailleSampling-j, z1RL2 - (j-1) * pasSample)  # echantillonage sur z1

        # calcul d'une solution epsilon-contrainte -----------------------------
        fRL, ebtmp2[tailleSampling-j].xEBD = relaxLinXG4(nbvar, nbcontraintes, L, c1, c2, z1RL2 - (j-1) * pasSample, 2)
        @printf("fRL = %8.2f  ",round(fRL, digits=2))

        # nettoyage de la valeur de ebtmp2[j].xEBD et calcul du point bi-objectif --
        for i in 1:nbvar
            if round(ebtmp2[tailleSampling-j].xEBD[i], digits=3) == 0.0
                ebtmp2[tailleSampling-j].xEBD[i] = 0.0
            elseif round(ebtmp2[tailleSampling-j].xEBD[i], digits=3) == 1.0
                ebtmp2[tailleSampling-j].xEBD[i] = 1.0
            else
                ebtmp2[tailleSampling-j].xEBD[i] = round(ebtmp2[tailleSampling-j].xEBD[i], digits=3)
            end
        end
        ebtmp2[tailleSampling-j].yEBD[1], ebtmp2[tailleSampling-j].yEBD[2] = evaluerSolution(ebtmp2[tailleSampling-j].xEBD, c1, c2)
        @printf("[ %8.2f , %8.2f ] \n", ebtmp2[tailleSampling-j].yEBD[1], ebtmp2[tailleSampling-j].yEBD[2])
    end

    # --------------------------------------------------------------------------
    # --------------------------------------------------------------------------
    @printf("\n 3) fusion ordonnee des points generateur pour constituer l'EBDual\n\n")

    # fusion de ebtmp1 et ebtmp2 dans eb ---------------------------------------
    j  = 1
    j1 = 1
    j2 = 1

    # fusion ordonnee
    while j1 ≤ tailleSampling && j2 ≤ tailleSampling-2
        if ebtmp1[j1].yEBD[1] ≤ ebtmp2[j2].yEBD[1]
            @printf("   1 → [ %8.2f , %8.2f ] \n", ebtmp1[j1].yEBD[1], ebtmp1[j1].yEBD[2])
            eb[j] = deepcopy(ebtmp1[j1])
            j  += 1
            j1 += 1
        else
            @printf("   2 → [ %8.2f , %8.2f ] \n", ebtmp2[j2].yEBD[1], ebtmp2[j2].yEBD[2])
            eb[j] = deepcopy(ebtmp2[j2])
            j  += 1
            j2 += 1
        end
    end
    # queue restante de ebtmp1
    while j1 ≤ tailleSampling
        @printf("   1 → [ %8.2f , %8.2f ] \n", ebtmp1[j1].yEBD[1], ebtmp1[j1].yEBD[2])
        eb[j] = deepcopy(ebtmp1[j1])
        j  += 1
        j1 += 1
    end
    # queue restante de ebtmp2 !! partie normalement jamais executee
    while j2 ≤ tailleSampling-2
        @printf("   2 → [ %8.2f , %8.2f ] \n", ebtmp2[j2].yEBD[1], ebtmp2[j2].yEBD[2])
        eb[j] = deepcopy(ebtmp2[j2])
        j  += 1
        j2 += 1
    end

    # --------------------------------------------------------------------------
    # --------------------------------------------------------------------------
    @printf("\n 4) copie EBDual et verifie l'admissibilite\n\n")

    for j=1:2*tailleSampling-2

        @printf("  %2d : [ %8.2f , %8.2f ] ", j, eb[j].yEBD[1], eb[j].yEBD[2])

        # copie de l'ensemble bornant dual dans la stru de donnees iterative ---
        vg[j].xRlx = deepcopy(eb[j].xEBD)
        vg[j].yRlx[1] = eb[j].yEBD[1]
        vg[j].yRlx[2] = eb[j].yEBD[2]
        push!(XEBD,eb[j].yEBD[1])
        push!(YEBD,eb[j].yEBD[2])

        # test d'admissibilite et marquage de la solution le cas echeant -------
        if estAdmissible(eb[j].xEBD)
            for i = 1:nbvar
                vg[j].xInt[i] = convert(Int, eb[j].xEBD[i])
            end
            vg[j].yInt[1] = eb[j].yEBD[1]
            vg[j].yInt[2] = eb[j].yEBD[2]
            vg[j].xFea = true
            push!(XFeas,eb[j].yEBD[1])
            push!(YFeas,eb[j].yEBD[2])
            @printf(" Admissible \n")
        else
            vg[j].xFea = false
            @printf(" x \n")
        end

    end
    println("")

#        @assert false "fin"

    n1 = eb[end].yEBD[1]
    n2 = eb[1].yEBD[2]

    # ==========================================================================
    # ==========================================================================

    @printf("5) terraformation en cours...\n\n")

    for trial=1:terraform
      println("--- ",trial," ---")

    #  Lst=[]

      for j=1:2*tailleSampling-2
          if vg[j].xFea == false

              # ----------------------------------------------------------------
              # Nettoyage de la valeur de vg[j].xRlx et calcul du point bi-objectif
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
              vg[j].yRlx[1],vg[j].yRlx[2] = evaluerSolution(vg[j].xRlx, c1, c2)
              @printf("  %2d : [ %8.2f , %8.2f ] ", j, vg[j].yRlx[1], vg[j].yRlx[2])


              # ----------------------------------------------------------------
              # Selectionne les points pour constituer le cone d'interet centre sur le point generateur O
              if j==1
                  # premier point (predecesseur fictif)
                  B=point(eb[1].yEBD[1], eb[1].yEBD[2]+1.0)
                  O=point(eb[1].yEBD[1], eb[1].yEBD[2])
                  A=point(eb[2].yEBD[1], eb[2].yEBD[2])
              elseif j==2*tailleSampling-2
                  # dernier points (successeur fictif)
                  B=point(eb[end-1].yEBD[1], eb[end-1].yEBD[2])
                  O=point(eb[end].yEBD[1], eb[end].yEBD[2])
                  A=point(eb[end].yEBD[1]+1.0, eb[end].yEBD[2])
              else
                  B=point(eb[j-1].yEBD[1], eb[j-1].yEBD[2])
                  O=point(eb[j].yEBD[1], eb[j].yEBD[2])
                  A=point(eb[j+1].yEBD[1], eb[j+1].yEBD[2])
              end

              # ----------------------------------------------------------------
              # Selectionne les points pour constituer le cone d'interet centre sur le point projete P
              P=point(vg[j].yPrj[1], vg[j].yPrj[2])
              @show P
              corner = point( A.x , B.y )
              @show corner
              if P.x < corner.x && P.y < corner.y
                  # replace P par corner
                  P=point( A.x , B.y )
              end

              # ----------------------------------------------------------------
              # Calcule la direction d'interet pour le point generateur

              x1,y1 = B.x, B.y
              x2,y2 = A.x, A.y
              xm=(x1+x2)/2.0
              ym=(y1+y2)/2.0
              Δx = abs(n1-xm)
              Δy = abs(n2-ym)
              λ1 =  1 - Δx / (Δx+Δy)
              λ2 =  1 - Δy / (Δx+Δy)
              @show Δx
              @show Δy
              println("directions deduites globalement")
              @show λ1
              @show λ2
              plot(n1, n2, xm, ym, linestyle="-", color="blue", marker="+")
              annotate("",
                       xy=[xm;ym],# Arrow tip
                       xytext=[n1;n2], # Text offset from tip
                       arrowprops=Dict("arrowstyle"=>"->"))

              # ----------------------------------------------------------------
              # Operations pour arrondir les valeurs non-entieres de la solution relachee
              vg[j].yInt[1] = 0
              vg[j].yInt[2] = 0
              z1 = vg[j].yRlx[1]
              z2 = vg[j].yRlx[2]

              @show B
              @show P
              @show A
              M = point( z1 , z2 )
              @show M

              @show inSector(M, P, A, B) # true -> cote (0,0) du cone
              #@assert false "brol"

              nbVarNonEntiere=0
              for i in 1:nbvar
                  if vg[j].xInt[i] == -1
                    #  if j==10
                    #      push!(Lst,M)
                    #  end

                      # la variable i est non entiere
                      nbVarNonEntiere += 1
                      M = point( z1 - vg[j].xRlx[i] * c1[i] , z2 - vg[j].xRlx[i] * c2[i])

                      if inSector(M, O, A, B)
                          # le point M obtenu est hors du cone => variable i fixee a 1
                          vg[j].xInt[i] = 1
                          #print(" ",i)
                          z1 = z1 - (vg[j].xRlx[i] * c1[i]) + c1[i]
                          z2 = z2 - (vg[j].xRlx[i] * c2[i]) + c2[i]
                      else
                          # le point M obtenu est dans le cone => variable i fixee a 0
                          vg[j].xInt[i] = 0
                          #print(" ",i)
                          z1 = z1 - (vg[j].xRlx[i] * c1[i])
                          z2 = z2 - (vg[j].xRlx[i] * c2[i])
                      end
                  end
                  vg[j].yInt[1] += vg[j].xInt[i] * c1[i]
                  vg[j].yInt[2] += vg[j].xInt[i] * c2[i]
              end

              @printf("→ #round : %4d → [ %5d , %5d ] ", nbVarNonEntiere, vg[j].yInt[1], vg[j].yInt[2])
              push!(XInt,vg[j].yInt[1])
              push!(YInt,vg[j].yInt[2])

    #          if j==10
    #          for i=1:length(Lst)
    #              scatter(Lst[i].x,Lst[i].y, color="black", marker=".")
    #          end
    #      end


              # ----------------------------------------------------------------
              # Projete la solution entiere sur le polytope X avec norme-L1

              fPrj, vg[j].xPrj = Δ2(L, vg[j].xInt)

              # Nettoyage de la valeur de vg[j].xPrj et calcul du point bi-objectif
              for i in 1:nbvar
                  if round(vg[j].xPrj[i], digits=3) == 0.0
                      vg[j].xPrj[i] = 0.0
                  elseif round(vg[j].xPrj[i], digits=3) == 1.0
                      vg[j].xPrj[i] = 1.0
                  else
                      vg[j].xPrj[i] = round(vg[j].xPrj[i], digits=3)
                  end
              end

              vg[j].yPrj[1], vg[j].yPrj[2] = evaluerSolution(vg[j].xPrj, c1, c2)

              @printf(" [ %8.2f , %8.2f ] ", vg[j].yPrj[1], vg[j].yPrj[2])
              push!(XProj, vg[j].yPrj[1])
              push!(YProj, vg[j].yPrj[2])


              # ----------------------------------------------------------------
              # Teste si la projection est admissible

              if estAdmissible(vg[j].xPrj)
                  vg[j].xFea = true
                  vg[j].yInt[1] = vg[j].yPrj[1]
                  vg[j].yInt[2] = vg[j].yPrj[2]
                  push!(XFeas, vg[j].yPrj[1])
                  push!(YFeas, vg[j].yPrj[2])
                  @printf(" Admissible \n")
              else
                  vg[j].xFea = false
                  @printf(" x \n")
                  # prepare pour l'iteration suivante
                  vg[j].xRlx = deepcopy(vg[j].xPrj)
              end

          end
      end
    end # for pump




    # ==========================================================================
    @printf("\n6) Edition des resultats \n\n")

    # Donne les points relaches initiaux ---------------------------------------
    scatter(XEBD,YEBD,color="blue", marker="x")
    @show XEBD
    @show YEBD

    # Donne les points entiers -------------------------------------------------
    scatter(XInt,YInt,color="orange", marker="s", label="int")
    @show XInt
    @show YInt

    # Donne les points apres projection Δ(x,x̃) ---------------------------------
    scatter(XProj,YProj, color="red", marker="x")
    @show XProj
    @show YProj

    # Donne les points admissibles ---------------------------------------------
    scatter(XFeas,YFeas, color="green", marker="o")
    @show XFeas
    @show YFeas

end

using PyPlot
using LinearAlgebra, Printf
function testidee()
    c1 = rand(0:100,20)
    c2 = rand(0:100,20)
    x  = rand(0.0:1.0,20)
    for i=1:5
        x[rand(1:20)]=rand()
    end
    Lfrac=[]
    for i=1:20
        if x[i]!=0.0 && x[i]!=1.0
            push!(Lfrac,i)
        end
    end
    @show x
    @show Lfrac

    # --------------------------------------------------------------------------
    P=[]
    push!(P,[rand(1:200) rand(800:1000) ])
    while P[end][1]<900 && P[end][2]>100
        push!(P,[rand(P[end][1]+1:min(1000,P[end][1]+500)) rand(max(1,P[end][2]-500):P[end][2]-1)])
    end
    # ---
    fig = figure("test idees",figsize=(8,8))
    xlim(0,1500)
    ylim(0,1500)
    for i=1:length(P)
        scatter(P[i][1],P[i][2], color="black")
    end
    i=rand(1:length(P)-1)
    a=P[i]
    x1,y1 = a
    b=P[i+1]
    x2,y2 = b
    Δx = abs(x1-x2)
    Δy = abs(y1-y2)
    λ1 =  1 - Δx / (Δx+Δy)
    λ2 =  1 - Δy / (Δx+Δy)
    @show P
    @show a
    @show b
    @show Δx
    @show Δy
    println("directions deduites localement")
    @show λ1
    @show λ2
    xlim(0,1000)
    ylim(0,1000)
    scatter(x1, y1, color="red", marker="x")
    scatter(x2, y2, color="red", marker="x")
    zλ = λ1 * dot(c1,x) + λ2 * dot(c2,x)
    #----
    n1=P[end][1]
    n2=P[1][2]
    xm=(x1+x2)/2.0
    ym=(y1+y2)/2.0
    Δx = abs(n1-xm)
    Δy = abs(n2-ym)
    λ1 =  1 - Δx / (Δx+Δy)
    λ2 =  1 - Δy / (Δx+Δy)
    @show Δx
    @show Δy
    println("directions deduites globalement")
    @show λ1
    @show λ2
    plot(n1, n2, xm, ym, linestyle="-", color="blue", marker="+")
    annotate("",
	xy=[xm;ym],# Arrow tip
	xytext=[n1;n2], # Text offset from tip
	arrowprops=Dict("arrowstyle"=>"->"))

    # --------------------------------------------------------------------------
    for f in Lfrac  @printf("%5d ",f) end ; @printf("\n")
    for f in Lfrac  @printf("%5d ",c1[f]) end; @printf("\n")
    for f in Lfrac  @printf("%5d ",c2[f]) end; @printf("\n")
    for f in Lfrac  @printf("%5.2f ",λ1*c1[f]+λ2*c2[f]) end; @printf("\n")
    for f in Lfrac  @printf("%5.2f ",x[f]) end; @printf("\n")

end

#testidee()
#mainXG4("biosppaa02.txt", 6, 5)
mainXG4("../SPA/instances/biosppaa02.txt", 6, 5)
nothing
