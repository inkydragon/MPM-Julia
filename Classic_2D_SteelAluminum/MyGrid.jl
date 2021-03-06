module moduleGrid

using Printf
import ..moduleMaterialPoint #sina, do not use include here, since you have already included the module in Main.jl

# fTime = 0.0

function index2DTo1D(i::Int64, j::Int64, nColumns::Int64, nRows::Int64)
    index = nColumns*(j-1) + i

    if(index > nRows*nColumns || index < 1)
        @printf("Index out of bounds: i, j: %d, %d \n", i, j)
    end

    return(Int64(index))
end

mutable struct mpmGridPoint
    v2Fixed::Vector{Bool}
    fMass::Float64
    v2Position::Vector{Float64}
    v2Velocity::Vector{Float64}
    v2Momentum::Vector{Float64}
    v2Force::Vector{Float64}

    function mpmGridPoint()
        new(
            [false; false],
            0.0,
            zeros(2),
            zeros(2),
            zeros(2),
            zeros(2)
        )
    end
end

mutable struct mpmGrid   #grid container
    v2Length_Grid::Vector{Float64}
    v2Nodes::Vector{Int64}
    iNodes::Int64
    v2Length_Cell::Vector{Float64}

    GridPoints::Vector{mpmGridPoint}

    function mpmGrid(fGL_x, fGL_y, iN_x, iN_y)

        v2CL = zeros(2)
        v2CL[1] = fGL_x / Float64(iN_x - 1.0)
        v2CL[2] = fGL_y / Float64(iN_y - 1.0)

        thisGridPoint = Vector{mpmGridPoint}(undef, iN_x * iN_y)
        for j=1:1:iN_y
            for i=1:1:iN_x
                x = (i-1) * v2CL[1]
                y = (j-1) * v2CL[2]
                index = index2DTo1D(i, j, iN_x, iN_y)

                bFixed_x = false
                bFixed_y = false
                if(i == 1 || i == iN_x)
                    bFixed_x = true
                    bFixed_y = true
                end
                if(j == 1 || j == iN_y)
                    bFixed_x = true
                    bFixed_y = true
                end

                thisGridPoint[index] = mpmGridPoint() # starts with every member equal to zero

                thisGridPoint[index].v2Fixed = [bFixed_x; bFixed_y]

                thisGridPoint[index].v2Position = [x; y]
            end
        end

        new([fGL_x; fGL_y], [iN_x; iN_y], iN_x*iN_y, v2CL, thisGridPoint)
    end
end

function getAdjacentGridPoints(thisMaterialPoint::moduleMaterialPoint.mpmMaterialPoint_2D_Classic, thisGrid::mpmGrid)
    thisAdjacentGridPoints = Vector{Int64}(undef, 0)

    v2Coordinate = thisMaterialPoint.v2Centroid

    fLength_Cell_x = thisGrid.v2Length_Cell[1]
    fLength_Cell_y = thisGrid.v2Length_Cell[2]

    iBottomLeft_i    = (floor(v2Coordinate[1] / fLength_Cell_x) + 1.0)
    iBottomLeft_j    = (floor(v2Coordinate[2] / fLength_Cell_y) + 1.0)

    if(iBottomLeft_j < 1 || iBottomLeft_j > thisGrid.v2Nodes[2])
        @printf("Index out of bounds: j: %d \n", iBottomLeft_j)
        @printf("v2Coordinate[2]: %e \n", v2Coordinate[2])
    end

    for i = iBottomLeft_i:1:iBottomLeft_i+1
        for j = iBottomLeft_j:1:iBottomLeft_j+1
            iIndex = index2DTo1D(Int64(i), Int64(j), thisGrid.v2Nodes[1], thisGrid.v2Nodes[2])

            push!(thisAdjacentGridPoints, iIndex)
        end
    end

    return((thisAdjacentGridPoints))
end

end # module moduleGrid