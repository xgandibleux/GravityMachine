function qualityMeasure(xN::Array{Int64,1}, yN::Array{Int64,1},
                        xU::Array{Int64,1}, yU::Array{Int64,1}
                       )

    # ---- Set the reference point
    xRef = xN[end] + 1
    yRef = yN[1] + 1

    # ---- Measure the area corresponding to the set of non-dominated points Y_N
    areaN = 0
    for i = 1:length(xN)-1
        areaN = areaN + (xN[i+1]-xN[i]) * (yRef - yN[i])
    end
    areaN = areaN + (xRef-xN[end]) * (yRef - yN[end])

    # ---- Measure the (restricted) area corresponding to the bound set U
    areaU = 0
    for i = 1:length(xU)-1
        areaU = areaU + ( min(xRef,xU[i+1]) - min(xRef,xU[i]) ) * (yRef - min(yRef,yU[i]) )
    end
    areaU = areaU + ( xRef - min(xRef,xU[end]) ) * (yRef - min(yRef,yU[end]) )

    # ---- Computer the ratio
    result = areaU/areaN
    return result

end
