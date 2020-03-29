module moduleMaterialPoint

using LinearAlgebra
import ..moduleMath #sina, do not use include here, since you have already included the module in Main.jl

mutable struct mpmMaterialPoint   #material point container
    fMass::Float64
    fVolumeInitial::Float64
    fVolume::Float64
    v2Length::moduleMath.Vector2D

    fElasticModulus::Float64
    fPoissonRatio::Float64

    v2Position::moduleMath.Vector2D
    v2PositionIncrement::moduleMath.Vector2D
    v2Velocity::moduleMath.Vector2D
    v2Momentum::moduleMath.Vector2D
    v2ExternalForce::moduleMath.Vector2D
    v2Restraint::moduleMath.Vector2D	# 0.0=no restraint, 1.0=fully restrained

    mCorner::Array{Float64,2}
    mCorner_Increment::Array{Float64,2}

    mRadial1::Vector{Float64}
    mRadial2::Vector{Float64}

    mDeformationGradient::Array{Float64, 2}
    mDeformationGradientIncrement::Array{Float64, 2}

    v3Strain::moduleMath.Vector3D
    v3Stress::moduleMath.Vector3D
    v3StrainIncrement::moduleMath.Vector3D
    v3StressIncrement::moduleMath.Vector3D

    function mpmMaterialPoint()
        new(1.0, 1.0, 1.0,
            moduleMath.Vector2D(0.0, 0.0),
            1.0, 0.3,
            moduleMath.Vector2D(0.0, 0.0),
            moduleMath.Vector2D(0.0, 0.0),
            moduleMath.Vector2D(0.0, 0.0),
            moduleMath.Vector2D(0.0, 0.0),
            moduleMath.Vector2D(0.0, 0.0),
            moduleMath.Vector2D(0.0, 0.0),
            zeros(4,2),
            zeros(4,2),
            ones(2),
            ones(2),
            Matrix{Float64}(I, 2, 2),
            Matrix{Float64}(I, 2, 2),
            moduleMath.Vector3D(0.0, 0.0, 0.0),
            moduleMath.Vector3D(0.0, 0.0, 0.0),
            moduleMath.Vector3D(0.0, 0.0, 0.0),
            moduleMath.Vector3D(0.0, 0.0, 0.0))
    end
    function mpmMaterialPoint(fM::Float64, fV::Float64, fEM::Float64, fPR::Float64, v2P::moduleMath.Vector2D, v2V::moduleMath.Vector2D, v2M::moduleMath.Vector2D, v2ExternalForce::moduleMath.Vector2D, v3Strain::moduleMath.Vector3D, v3Stress::moduleMath.Vector3D)
        new(fM, fV, fV,
            moduleMath.Vector2D(0.0, 0.0),
            fEM, fPR,
            moduleMath.Vector2D(v2P),
            moduleMath.Vector2D(0.0, 0.0),
            moduleMath.Vector2D(v2V),
            moduleMath.Vector2D(v2M),
            moduleMath.Vector2D(v2ExternalForce),
            moduleMath.Vector2D(0.0, 0.0),
            zeros(4,2),
            zeros(4,2),
            ones(2),
            ones(2),
            Matrix{Float64}(I, 2, 2),
            Matrix{Float64}(I, 2, 2),
            moduleMath.Vector3D(0.0, 0.0, 0.0),
            moduleMath.Vector3D(0.0, 0.0, 0.0),
            moduleMath.Vector3D(0.0, 0.0, 0.0),
            moduleMath.Vector3D(0.0, 0.0, 0.0))
    end
end

function createMaterialDomain_Rectangle(fCenter_x::Float64, fCenter_y::Float64, fWidth::Float64, fHeight::Float64, fOffset::Float64)
    thisMaterialDomain = Vector{mpmMaterialPoint}(undef, 0)

    fWidth	= floor(fWidth/fOffset) * fOffset	#just in case width is not a multiple of offset
    fHeight	= floor(fHeight/fOffset) * fOffset	#just in case height is not a multiple of offset

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
    thisMaterialDomain = Array{mpmMaterialPoint}(undef, 0)

    fRadius = floor(fRadius/fOffset) * fOffset	#just in case radius is not a multiple of offset

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