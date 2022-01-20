using PyPlot

mutable struct point
    x::Float64
    y::Float64
end


function inCone(pOrg, pDeb, pFin, pCur)
    # pOrg : point origine du cone (la ou il est pointe)
    # pDeb : point depart du cone (point du rayon [pOrg,pDeb])
    # pFin : point final du cone (point du rayon [pOrg,pFin])
    # pCur : point courant a tester
    # retourne VRAI si pCur est dans le cone pDeb-pFin-pOrg, FAUX sinon

    cp_pDeb_pFin = (pDeb.x - pOrg.x) * (pFin.y - pOrg.y) - (pDeb.y - pOrg.y) * (pFin.x - pOrg.x)
    cp_pDeb_pCur = (pDeb.x - pOrg.x) * (pCur.y - pOrg.y) - (pDeb.y - pOrg.y) * (pCur.x - pOrg.x)
    cp_pFin_pCur = (pFin.x - pOrg.x) * (pCur.y - pOrg.y) - (pFin.y - pOrg.y) * (pCur.x - pOrg.x)

    if (cp_pDeb_pFin > 0)
        if ((cp_pDeb_pCur >= 0) && (cp_pFin_pCur <= 0))
            return true
        else
            return false
        end
    else
        if (!((cp_pDeb_pCur < 0) && (cp_pFin_pCur > 0)))
            return true
        else
            return false
        end
    end
end

function inCone1VersZ(pOrg, pDeb, pFin, pCur)
    return inCone(pOrg, pDeb, pFin, pCur)
end

function inCone2Vers0(pOrg, pDeb, pFin, pCur)
    return !inCone(pOrg, pDeb, pFin, pCur)
end

xlim(0,10)
ylim(0,10)
B=point(2.0,5.0); scatter(2.0,5.0,marker=".",color="black")
O=point(2.5,2.5); scatter(2.5,2.5,marker="x",color="red")
A=point(5.0,2.0); scatter(5.0,2.0,marker=".",color="black")
#N=point(5.0,5.0); scatter(5.0,5.0,marker="+")

x=round(rand()*8+2;digits=2); y=round(rand()*8+2;digits=2)
N=point(x,y); scatter(x,y,marker="x",color="blue")

x=round(rand()*10;digits=2); y=round(rand()*10;digits=2)
M=point(x,y); scatter(x,y,marker="o",color="magenta")

lx1=[O.x,B.x]; ly1=[O.y,B.y]; plot(lx1,ly1,color="red")
lx1=[O.x,0]; ly1=[O.y,(B.y-O.y)/(B.x-O.x)*(0-O.x)+O.y]; plot(lx1,ly1,color="red",linestyle="dotted")
lx2=[O.x,A.x]; ly2=[O.y,A.y]; plot(lx2,ly2,color="red")
lx2=[O.x,10]; ly2=[O.y,(A.y-O.y)/(A.x-O.x)*(10-O.x)+O.y]; plot(lx2,ly2,color="red",linestyle="dotted")

lx3=[N.x,B.x]; ly3=[N.y,B.y]; plot(lx3,ly3,color="blue",linestyle="dashed")
lx4=[N.x,A.x]; ly4=[N.y,A.y]; plot(lx4,ly4,color="blue",linestyle="dashed")

@show inCone1VersZ(O, A, B, M)
@show inCone2Vers0(N, A, B, M)
nothing
