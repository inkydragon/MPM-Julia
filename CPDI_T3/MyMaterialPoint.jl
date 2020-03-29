module moduleMaterialPoint

using LinearAlgebra
using ..moduleMath #sina, do not use include here, since you have already included the module in Main.jl

export mpmMaterialPoint,
    createMaterialDomain_Rectangle,
    createMaterialDomain_Circle

mutable struct mpmMaterialPoint   #material point container
    fMass::Float64
    fVolumeInitial::Float64
    fVolume::Float64
    v2Length::Vector2D

    fElasticModulus::Float64
    fPoissonRatio::Float64

    v2Position::Vector2D
    v2PositionIncrement::Vector2D
    v2Velocity::Vector2D
    v2Momentum::Vector2D
    v2ExternalForce::Vector2D
    v2Restraint::Vector2D    # 0.0=no restraint, 1.0=fully restrained

    mCorner::Array{Float64,2}
    mCorner_Increment::Array{Float64,2}

    mRadial1::Vector{Float64}
    mRadial2::Vector{Float64}

    mDeformationGradient::Array{Float64,2}
    mDeformationGradientIncrement::Array{Float64,2}

    v3Strain::Vector3D
    v3Stress::Vector3D
    v3StrainIncrement::Vector3D
    v3StressIncrement::Vector3D

    function mpmMaterialPoint()
        new(
            1.0, 1.0, 1.0,
            Vector2D(0.0, 0.0),
            1.0, 0.3,
            Vector2D(0.0, 0.0),
            Vector2D(0.0, 0.0),
            Vector2D(0.0, 0.0),
            Vector2D(0.0, 0.0),
            Vector2D(0.0, 0.0),
            Vector2D(0.0, 0.0),
            zeros(0,2),
            zeros(0,2),
            ones(2),
            ones(2),
            Matrix{Float64}(I, 2, 2),
            Matrix{Float64}(I, 2, 2),
            Vector3D(0.0, 0.0, 0.0),
            Vector3D(0.0, 0.0, 0.0),
            Vector3D(0.0, 0.0, 0.0),
            Vector3D(0.0, 0.0, 0.0)
        )
    end

    function mpmMaterialPoint(
        fM::Float64,
        fV::Float64,
        fEM::Float64,
        fPR::Float64,
        v2P::Vector2D,
        v2V::Vector2D,
        v2M::Vector2D,
        v2ExternalForce::Vector2D,
        v3Strain::Vector3D,
        v3Stress::Vector3D
    )
        new(fM, fV, fV,
            Vector2D(0.0, 0.0),
            fEM, fPR,
            v2P,
            Vector2D(0.0, 0.0),
            v2V,
            v2M,
            v2ExternalForce,
            Vector2D(0.0, 0.0),
            zeros(0,2),
            zeros(0,2),
            ones(2),
            ones(2),
            Matrix{Float64}(I, 2, 2),
            Matrix{Float64}(I, 2, 2),
            Vector3D(0.0, 0.0, 0.0),
            Vector3D(0.0, 0.0, 0.0),
            Vector3D(0.0, 0.0, 0.0),
            Vector3D(0.0, 0.0, 0.0))
    end
end

function createMaterialDomain_Rectangle(sParticleShape::String, fCenter_x::Float64, fCenter_y::Float64, fWidth::Float64, fHeight::Float64, fOffset::Float64)
    thisMaterialDomain = Vector{mpmMaterialPoint}(undef, 0)

    if(sParticleShape == "triangle")
        # left triangles
        for fBaseCorner_y = -0.5*fHeight:fOffset:+0.5*fHeight
            for fBaseCorner_x = -0.5*fWidth:fOffset:+0.5*fWidth
                v2Corner = Array{Float64, 2}(undef, 3,2)
                v2Corner[1,1] = fBaseCorner_x + 0.0
                v2Corner[1,2] = fBaseCorner_y + 0.0

                v2Corner[2,1] = fBaseCorner_x + 0.5*fOffset
                v2Corner[2,2] = fBaseCorner_y + 0.5*fOffset

                v2Corner[3,1] = fBaseCorner_x + 0.0
                v2Corner[3,2] = fBaseCorner_y + 1.0*fOffset

                v2Centroid = [0.0 0.0]

                fWeight = 1.0 / size(v2Corner,1)
                for iIndex_Corner = 1:1:size(v2Corner,1)
                    v2Centroid[1,1] += fWeight * v2Corner[iIndex_Corner,1]
                    v2Centroid[1,2] += fWeight * v2Corner[iIndex_Corner,2]
                end

                if(-0.5*fHeight < v2Centroid[1,2] && v2Centroid[1,2] < +0.5fHeight)
                    if(-0.5*fWidth < v2Centroid[1,1] && v2Centroid[1,1] < +0.5fWidth)
                        thisMaterialPoint = mpmMaterialPoint()
                        thisMaterialPoint.v2Position.fx = fCenter_x + v2Centroid[1,1]
                        thisMaterialPoint.v2Position.fy = fCenter_y + v2Centroid[1,2]

                        for iIndex_Corner = 1:1:size(v2Corner,1) # 4 corners for rectangle
                            newCorner = [0.0 0.0]
                            newCorner[1,1] = fCenter_x + v2Corner[iIndex_Corner,1]
                            newCorner[1,2] = fCenter_y + v2Corner[iIndex_Corner,2]
                            thisMaterialPoint.mCorner = vcat(thisMaterialPoint.mCorner, newCorner)

                            thisMaterialPoint.mCorner_Increment = vcat(thisMaterialPoint.mCorner_Increment, [0.0 0.0])
                        end

                        push!(thisMaterialDomain, thisMaterialPoint)
                    end
                end
            end
        end

        # bottom triangles
        for fBaseCorner_y = -0.5*fHeight:fOffset:+0.5*fHeight
            for fBaseCorner_x = -0.5*fWidth:fOffset:+0.5*fWidth
                v2Corner = Array{Float64, 2}(undef, 3,2)
                v2Corner[1,1] = fBaseCorner_x + 0.0
                v2Corner[1,2] = fBaseCorner_y + 0.0

                v2Corner[2,1] = fBaseCorner_x + 1.0*fOffset
                v2Corner[2,2] = fBaseCorner_y + 0.0

                v2Corner[3,1] = fBaseCorner_x + 0.5*fOffset
                v2Corner[3,2] = fBaseCorner_y + 0.5*fOffset

                v2Centroid = [0.0 0.0]

                fWeight = 1.0 / size(v2Corner,1)
                for iIndex_Corner = 1:1:size(v2Corner,1)
                    v2Centroid[1,1] += fWeight * v2Corner[iIndex_Corner,1]
                    v2Centroid[1,2] += fWeight * v2Corner[iIndex_Corner,2]
                end

                if(-0.5*fHeight < v2Centroid[1,2] && v2Centroid[1,2] < +0.5fHeight)
                    if(-0.5*fWidth < v2Centroid[1,1] && v2Centroid[1,1] < +0.5fWidth)
                        thisMaterialPoint = mpmMaterialPoint()
                        thisMaterialPoint.v2Position.fx = fCenter_x + v2Centroid[1,1]
                        thisMaterialPoint.v2Position.fy = fCenter_y + v2Centroid[1,2]

                        for iIndex_Corner = 1:1:size(v2Corner,1) # 4 corners for rectangle
                            newCorner = [0.0 0.0]
                            newCorner[1,1] = fCenter_x + v2Corner[iIndex_Corner,1]
                            newCorner[1,2] = fCenter_y + v2Corner[iIndex_Corner,2]
                            thisMaterialPoint.mCorner = vcat(thisMaterialPoint.mCorner, newCorner)

                            thisMaterialPoint.mCorner_Increment = vcat(thisMaterialPoint.mCorner_Increment, [0.0 0.0])
                        end

                        push!(thisMaterialDomain, thisMaterialPoint)
                    end
                end
            end
        end

        # right triangles
        for fBaseCorner_y = -0.5*fHeight:fOffset:+0.5*fHeight
            for fBaseCorner_x = -0.5*fWidth:fOffset:+0.5*fWidth
                v2Corner = Array{Float64, 2}(undef, 3,2)
                v2Corner[1,1] = fBaseCorner_x + 1.0*fOffset
                v2Corner[1,2] = fBaseCorner_y + 0.0

                v2Corner[2,1] = fBaseCorner_x + 1.0*fOffset
                v2Corner[2,2] = fBaseCorner_y + 1.0*fOffset

                v2Corner[3,1] = fBaseCorner_x + 0.5*fOffset
                v2Corner[3,2] = fBaseCorner_y + 0.5*fOffset

                v2Centroid = [0.0 0.0]

                fWeight = 1.0 / size(v2Corner,1)
                for iIndex_Corner = 1:1:size(v2Corner,1)
                    v2Centroid[1,1] += fWeight * v2Corner[iIndex_Corner,1]
                    v2Centroid[1,2] += fWeight * v2Corner[iIndex_Corner,2]
                end

                if(-0.5*fHeight < v2Centroid[1,2] && v2Centroid[1,2] < +0.5fHeight)
                    if(-0.5*fWidth < v2Centroid[1,1] && v2Centroid[1,1] < +0.5fWidth)
                        thisMaterialPoint = mpmMaterialPoint()
                        thisMaterialPoint.v2Position.fx = fCenter_x + v2Centroid[1,1]
                        thisMaterialPoint.v2Position.fy = fCenter_y + v2Centroid[1,2]

                        for iIndex_Corner = 1:1:size(v2Corner,1) # 4 corners for rectangle
                            newCorner = [0.0 0.0]
                            newCorner[1,1] = fCenter_x + v2Corner[iIndex_Corner,1]
                            newCorner[1,2] = fCenter_y + v2Corner[iIndex_Corner,2]
                            thisMaterialPoint.mCorner = vcat(thisMaterialPoint.mCorner, newCorner)

                            thisMaterialPoint.mCorner_Increment = vcat(thisMaterialPoint.mCorner_Increment, [0.0 0.0])
                        end

                        push!(thisMaterialDomain, thisMaterialPoint)
                    end
                end
            end
        end

        # top triangles
        for fBaseCorner_y = -0.5*fHeight:fOffset:+0.5*fHeight
            for fBaseCorner_x = -0.5*fWidth:fOffset:+0.5*fWidth
                v2Corner = Array{Float64, 2}(undef, 3,2)
                v2Corner[1,1] = fBaseCorner_x + 1.0*fOffset
                v2Corner[1,2] = fBaseCorner_y + 1.0*fOffset

                v2Corner[2,1] = fBaseCorner_x + 0.0
                v2Corner[2,2] = fBaseCorner_y + 1.0*fOffset

                v2Corner[3,1] = fBaseCorner_x + 0.5*fOffset
                v2Corner[3,2] = fBaseCorner_y + 0.5*fOffset

                v2Centroid = [0.0 0.0]

                fWeight = 1.0 / size(v2Corner,1)
                for iIndex_Corner = 1:1:size(v2Corner,1)
                    v2Centroid[1,1] += fWeight * v2Corner[iIndex_Corner,1]
                    v2Centroid[1,2] += fWeight * v2Corner[iIndex_Corner,2]
                end

                if(-0.5*fHeight < v2Centroid[1,2] && v2Centroid[1,2] < +0.5fHeight)
                    if(-0.5*fWidth < v2Centroid[1,1] && v2Centroid[1,1] < +0.5fWidth)
                        thisMaterialPoint = mpmMaterialPoint()
                        thisMaterialPoint.v2Position.fx = fCenter_x + v2Centroid[1,1]
                        thisMaterialPoint.v2Position.fy = fCenter_y + v2Centroid[1,2]

                        for iIndex_Corner = 1:1:size(v2Corner,1) # 4 corners for rectangle
                            newCorner = [0.0 0.0]
                            newCorner[1,1] = fCenter_x + v2Corner[iIndex_Corner,1]
                            newCorner[1,2] = fCenter_y + v2Corner[iIndex_Corner,2]
                            thisMaterialPoint.mCorner = vcat(thisMaterialPoint.mCorner, newCorner)

                            thisMaterialPoint.mCorner_Increment = vcat(thisMaterialPoint.mCorner_Increment, [0.0 0.0])
                        end

                        push!(thisMaterialDomain, thisMaterialPoint)
                    end
                end
            end
        end
    end

    if(sParticleShape == "rectangle")
        for fBaseCorner_y = -0.5*fHeight:fOffset:+0.5*fHeight
            for fBaseCorner_x = -0.5*fWidth:fOffset:+0.5*fWidth
                v2Corner = Array{Float64, 2}(undef, 4,2)
                v2Corner[1,1] = fBaseCorner_x + 0.0
                v2Corner[1,2] = fBaseCorner_y + 0.0

                v2Corner[2,1] = fBaseCorner_x + fOffset
                v2Corner[2,2] = fBaseCorner_y + 0.0

                v2Corner[3,1] = fBaseCorner_x + fOffset
                v2Corner[3,2] = fBaseCorner_y + fOffset

                v2Corner[4,1] = fBaseCorner_x + 0.0
                v2Corner[4,2] = fBaseCorner_y + fOffset

                v2Centroid = [0.0 0.0]

                for iIndex_Corner = 1:1:size(v2Corner,1)
                    v2Centroid[1,1] += 0.25 * v2Corner[iIndex_Corner,1]
                    v2Centroid[1,2] += 0.25 * v2Corner[iIndex_Corner,2]
                end

                if(-0.5*fHeight < v2Centroid[1,2] && v2Centroid[1,2] < +0.5fHeight)
                    if(-0.5*fWidth < v2Centroid[1,1] && v2Centroid[1,1] < +0.5fWidth)
                        thisMaterialPoint = mpmMaterialPoint()
                        thisMaterialPoint.v2Position.fx = fCenter_x + v2Centroid[1,1]
                        thisMaterialPoint.v2Position.fy = fCenter_y + v2Centroid[1,2]

                        for iIndex_Corner = 1:1:size(v2Corner,1)
                            newCorner = [0.0 0.0]
                            newCorner[1,1] = fCenter_x + v2Corner[iIndex_Corner,1]
                            newCorner[1,2] = fCenter_y + v2Corner[iIndex_Corner,2]
                            thisMaterialPoint.mCorner = vcat(thisMaterialPoint.mCorner, newCorner)

                            thisMaterialPoint.mCorner_Increment = vcat(thisMaterialPoint.mCorner_Increment, [0.0 0.0])
                        end

                        push!(thisMaterialDomain, thisMaterialPoint)
                    end
                end
            end
        end
    end

    return(thisMaterialDomain)
end

function createMaterialDomain_Rectangle(fCenter_x::Float64, fCenter_y::Float64, fWidth::Float64, fHeight::Float64, fOffset::Float64)
    thisMaterialDomain = Vector{mpmMaterialPoint}(undef, 0)

    fWidth = floor(fWidth/fOffset) * fOffset    #just in case width is not a multiple of offset
    fHeight = floor(fHeight/fOffset) * fOffset    #just in case height is not a multiple of offset

    for fy in -0.5*fHeight+0.5*fOffset:fOffset:+0.5*fHeight-0.5*fOffset
        for fx in -0.5*fWidth+0.5*fOffset:fOffset:+0.5*fWidth-0.5*fOffset
            thisMaterialPoint = mpmMaterialPoint()
            thisMaterialPoint.v2Position.fx = fCenter_x + fx
            thisMaterialPoint.v2Position.fy = fCenter_y + fy
            push!(thisMaterialDomain, thisMaterialPoint)
        end
    end

    return(thisMaterialDomain)
end

function createMaterialDomain_Circle(fCenter_x::Float64, fCenter_y::Float64, fRadius::Float64, fOffset::Float64)
    thisMaterialDomain = Vector{mpmMaterialPoint}(undef, 0)

    fRadius = floor(fRadius/fOffset) * fOffset    #just in case radius is not a multiple of offset

    for fy in -fRadius+0.5*fOffset:fOffset:+fRadius-0.5*fOffset
        for fx in -fRadius+0.5*fOffset:fOffset:+fRadius-0.5*fOffset
            if(fx^2 + fy^2 < fRadius^2)
                thisMaterialPoint = mpmMaterialPoint()
                thisMaterialPoint.v2Position.fx = fCenter_x + fx
                thisMaterialPoint.v2Position.fy = fCenter_y + fy
                push!(thisMaterialDomain, thisMaterialPoint)
            end
        end
    end

    return(thisMaterialDomain)
end

end # module moduleMaterialPoint