# Code qui determine un jeu de poids soit local a 2 points, soit global a tous les points

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

testidee()