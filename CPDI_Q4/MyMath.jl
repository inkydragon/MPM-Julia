module moduleMath

mutable struct Vector2D   #node container
    fx::Float64
    fy::Float64

    function Vector2D(x::Float64, y::Float64)
        new(x, y)
    end
end
Vector2D() = Vector2D(0.0, 0.0)

mutable struct Vector3D   #node container
    f1::Float64
    f2::Float64
    f3::Float64

    function Vector3D(x::Float64, y::Float64, z::Float64)
        new(x, y, z)
    end
end
Vector3D() = Vector3D(0.0, 0.0, 0.0)

end # module moduleMath