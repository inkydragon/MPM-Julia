module moduleMaterialPoint

using LinearAlgebra # I
export mpmMaterialPoint_2D_Classic,
    createMaterialDomain_Circle

# material point container
mutable struct mpmMaterialPoint_2D_Classic
    fMass            :: Float64
    fVolumeInitial   :: Float64
    fVolume          :: Float64

    fElasticModulus  :: Float64
    fPoissonRatio    :: Float64

    v2Centroid       :: Vector{Float64} # position
    v2Velocity       :: Vector{Float64}
    v2Momentum       :: Vector{Float64}
    v2ExternalForce  :: Vector{Float64}
    v2Restraint      :: Vector{Float64} # 0.0=no restraint, 1.0=fully restrained

    v2Corner         :: Array{Float64,2} # corner position

    m22DeformationGradient          :: Array{Float64,2}
    m22DeformationGradientIncrement :: Array{Float64,2}

    v3Strain         :: Vector{Float64} # xx, yy, zz, xy, yz, zx
    v3Stress         :: Vector{Float64}

    function mpmMaterialPoint_2D_Classic()
        new(
            1.0, 1.0, 1.0, # Mass, initial volume, volume
            1.0, 0.3, # elastic modulus, poisson ratio
            zeros(2), # centroid position
            zeros(2), # velocity
            zeros(2), # momentum
            zeros(2), # external force
            zeros(2), # restraint
            zeros(2,4), # array of corner positions, 3d coordinates
            Matrix{Float64}(I, 2, 2), # deformation gradient
            Matrix{Float64}(I, 2, 2), # deformation gradient increment
            zeros(3), # strain
            zeros(3), # stress
        )
    end
end

function createMaterialDomain_Circle(fCenter::Vector{Float64}, fRadius::Float64, fOffset::Float64)
    thisMaterialDomain = Vector{mpmMaterialPoint_2D_Classic}(undef, 0)

    # just in case radius is not a multiple of offset
    fRadius = floor(fRadius/fOffset) * fOffset    

    start = -fRadius+0.5*fOffset
    xy_range = start:fOffset:-start
    for fy in xy_range, fx in xy_range
        if(fx^2 + fy^2 < fRadius^2)
            thisMaterialPoint = mpmMaterialPoint_2D_Classic()
            thisMaterialPoint.v2Centroid = [fCenter[1] + fx; fCenter[2] + fy]
            push!(thisMaterialDomain, thisMaterialPoint)
        end
    end

    thisMaterialDomain
end

end # module moduleMaterialPoint