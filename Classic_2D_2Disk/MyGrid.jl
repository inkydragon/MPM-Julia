module moduleGrid

using ..moduleMaterialPoint #sina, do not use include here, since you have already included the module in Main.jl
export mpmGridPoint, mpmGrid,
    getAdjacentGridPoints


mutable struct mpmGridPoint
    v2Fixed     :: Vector{Bool}
    fMass       :: Float64
    v2Position  :: Vector{Float64}
    v2Momentum  :: Vector{Float64}
    v2Force     :: Vector{Float64}

    function mpmGridPoint()
        new(
            [false; false],
            0.0,
            zeros(2),
            zeros(2),
            zeros(2)
        )
    end
end

function index2DTo1D(i::Int64, j::Int64, nColumns::Int64, nRows::Int64)
    index = nColumns*(j-1) + i

    if(index > nRows*nColumns || index < 1)
        @warn("Index out of bounds: i, j: ", i, j)
    end

    index
end

"grid container"
mutable struct mpmGrid
    "length in x and y dirs"
    v2Length_Grid   :: Vector{Float64}
    v2Nodes         :: Vector{Int64}
    "number of grid nodes"
    iNodes          :: Int64

    "size of each cell/element"
    v2Length_Cell   :: Vector{Float64}
    "inverse of size of each cell/element"
    v2Length_CellI  :: Vector{Float64}

    "array of all grid points"
    GridPoints      :: Vector{mpmGridPoint}

    """
    constructor, GL_x is length of the grid in x dir
    iN_x: number of nodes in x dir
    """
    function mpmGrid(fGL_x, fGL_y, iN_x, iN_y)
        v2CL   = zeros(2)
        v2CLI  = zeros(2)
        v2CL[1]   = fGL_x / (iN_x - 1.0)
        v2CL[2]   = fGL_y / (iN_y - 1.0)
        v2CLI[1]  = 1.0 / v2CL[1]
        v2CLI[2]  = 1.0 / v2CL[2]

        thisGridPoint = Vector{mpmGridPoint}(undef, iN_x * iN_y)
        for j = 1:iN_y, i = 1:iN_x
            x = (i-1) * v2CL[1]
            y = (j-1) * v2CL[2]
            index = index2DTo1D(i, j, iN_x, iN_y)

            # starts with every member equal to zero
            thisGridPoint[index] = mpmGridPoint()
            thisGridPoint[index].v2Fixed = [false; false]
            thisGridPoint[index].v2Position = [x; y]
        end

        new([fGL_x; fGL_y], [iN_x; iN_y], iN_x*iN_y, v2CL, v2CLI, thisGridPoint)
    end
end

function getAdjacentGridPoints(thisMaterialPoint::mpmMaterialPoint_2D_Classic, thisGrid::mpmGrid)
    v2Coordinate   = thisMaterialPoint.v2Centroid
    fLength_Cell_x = thisGrid.v2Length_CellI[1]
    fLength_Cell_y = thisGrid.v2Length_CellI[2]

    iBottomLeft_i  = floor(Int64, v2Coordinate[1] * fLength_Cell_x + 1.0)
    iBottomLeft_j  = floor(Int64, v2Coordinate[2] * fLength_Cell_y + 1.0)

    if(iBottomLeft_j < 1 || iBottomLeft_j > thisGrid.v2Nodes[2])
        @warn("Index out of bounds: j: ", iBottomLeft_j)
        @warn("v2Coordinate[2]: ", v2Coordinate[2])
    end

    iIndex = index2DTo1D(iBottomLeft_i, iBottomLeft_j,   thisGrid.v2Nodes[1], thisGrid.v2Nodes[2])
    jIndex = index2DTo1D(iBottomLeft_i, iBottomLeft_j+1, thisGrid.v2Nodes[1], thisGrid.v2Nodes[2])

    [ iIndex, jIndex, iIndex+1, jIndex+1 ]
end

end # module moduleGrid