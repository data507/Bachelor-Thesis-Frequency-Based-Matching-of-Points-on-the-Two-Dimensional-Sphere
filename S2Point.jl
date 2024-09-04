## Datatype S2Point for easier calculations on the sphere

struct S2Point
    #Cartesian
    x::Float64
    y::Float64
    z::Float64

    #Spherical
    θ::Float64
    ϕ::Float64
end

Base.show(io::IO, x::S2Point) = print(io, "S2Point([x,y,z] = [$(x.x), $(x.y), $(x.z)]; [θ, ϕ] = [$(x.θ/pi)π, $(x.ϕ/pi)π] rad)")

## Constructor

#Cartesian to S2Point
function S2Point(x,y,z)

    (x2, y2, z2) = Float64.([x,y,z])

    θ = acos(z2)
    ϕ = sign(y2) * acos(x2 / sqrt(x2.^2 + y2.^2))
    if (x2 == 0 && y2 == 0) ϕ = 0 end
    
    return S2Point(x2,y2,z2,θ,ϕ)
end

#Spherical to S2Point
function S2Point(θ,ϕ)
    
    (θ2,ϕ2) = Float64.([θ,ϕ])
    
    sθ, cθ = sincos(θ2)
    sϕ, cϕ = sincos(ϕ2)
    
    x = sθ * cϕ
    y = sθ * sϕ
    z = cθ
    
    return S2Point(x,y,z,θ2,ϕ2)
end

# Vector to S2Point
function S2Point(v::SVector{3, Float64})
    return S2Point(v.x,v.y,v.z)
end

function S2Point(v::Vector)
    if length(v) == 3
        return S2Point(v[1],v[2],v[3])
    elseif length(v) == 2
        return S2Point(v[1],v[2])
    else 
        throw("Dimension Mismatch: Argument of S2Point must have dimension 2 or 3 but is $(length(v))")
    end
end

# Getter
getx(v::S2Point) = v.x;
gety(v::S2Point) = v.y;
getz(v::S2Point) = v.z;
getθ(v::S2Point) = v.θ;
getϕ(v::S2Point) = v.ϕ;

## Operations

#Rotation Operations
import Base.*, Base.+

function Base.:+(u::S2Point, v::S2Point)
    return [u.x + v.x, u.y + v.y, u.z + v.z]
end

function Base.:+(v::Vector{Float64}, u::S2Point)
    return v + [u.x, u.y, u.z]
end

function Base.:*(R::RotZYZ{Float64}, S2pt::S2Point)
    return S2Point(normalize(R * SVector(S2pt.x, S2pt.y, S2pt.z)))
end

function Base.:*(R::RotZYZ{Float64}, S2ptVec::Vector{S2Point})
    return [R * S2ptVec[i] for i =1:length(S2ptVec)]
end

function Base.:*(R::RotMatrix3{Float64}, S2pt::S2Point)
    return S2Point(normalize(R * SVector(S2pt.x, S2pt.y, S2pt.z)))
end

function Base.:*(R::RotMatrix3{Float64}, S2ptVec::Vector{S2Point})
    return [R * S2ptVec[i] for i =1:length(S2ptVec)]
end


#Other functions

function randS2(n)
    return [S2Point(normalize(randn(3))) for i = 1:n]
end