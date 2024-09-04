import LinearAlgebra, SphericalHarmonics
import SphericalHarmonics.SphericalHarmonicArrays, SphericalHarmonics.SphericalHarmonicModes
import WignerD

include("S2Point.jl")


## useful Rotations

function rot_to_northpole(v::S2Point)
    return RotZYZ(0, -v.θ, -v.ϕ)
end

function rot_u_to_v(x::S2Point,y::S2Point,t::Float64)
    return inv(rot_to_northpole(y)) * RotZ(2*pi*t) * rot_to_northpole(x)
end

function euler_zyz(R::RotMatrix3{Float64})
    α = atan(R[2, 3], R[1, 3])
    β = atan(sqrt(1 - R[3, 3].^2), R[3,3])
    γ = atan(R[3, 2], -R[3, 1])
    return α, β, γ
end

## S2Transformation

function S2TransformS2Point(U::Vector{S2Point}, B::Int64)
    if (B==1) ## Faster evaluation for B=1
        return sum([conj(SphericalHarmonicArrays.SHArray([sqrt(1/(4*pi)), 1/2 * sqrt(3/(2*pi)) * sin(u.θ) * exp(-im * u.ϕ), 1/2 * sqrt(3/pi) * cos(u.θ), -1/2 * sqrt(3/(2*pi)) * sin(u.θ) * exp(im * u.ϕ)], (SphericalHarmonicModes.ML(0:1, -1:1),))) for u in U])
    else 
        return sum([conj(computeYlm(u.θ, u.ϕ, lmax=B)) for u in U])
    end
end


## Evaluation of Convolution

function wigner_eval(V::Vector{S2Point}, U::Vector{S2Point}, B::Int64, Wigner_Matrices) # evaluates the convolution (V conv U) at the given Wigner_Matrices
    u_hat = S2TransformS2Point(U, B)
    v_hat = S2TransformS2Point(V, B)

    coeff = Vector{Any}(undef, B+1)
    for b = 0:B
        u_hat_b = [u_hat[(b,i)] for i in -b:b]
        v_hat_b = [v_hat[(b,i)] for i in -b:b]
        coeff[b+1] = conj(v_hat_b) * transpose(u_hat_b)
    end

    return [eval_SO3Fun(coeff, Wigner_D) for Wigner_D in Wigner_Matrices]
end

function eval_SO3Fun(SO3coeff, Wigner_D) # sums up SO3coefficients with the corresponding Wigner-D Function
    M = broadcast(.*, SO3coeff, Wigner_D)
    return sum(sum(m) for m in M);
end


## Algorithms

function FindRotation_general(U::Vector{S2Point},V::Vector{S2Point},B::Int64, S::Int64) # finds maximum of convolution (V conv U) by sampling S points
    ũ = U[1]
    
    mu_ũ_hat =  S2TransformS2Point(rot_to_northpole(ũ) * U, B)
    nu_Ṽ_hat = [S2TransformS2Point(rot_to_northpole(v) * V, B) for v in V]
    
    c_m_Ṽ = [[sum([conj(nu_ṽ_hat[(j,m)]) * mu_ũ_hat[(j,m)] for j = abs(m):B]) for m = -B:B] for nu_ṽ_hat in nu_Ṽ_hat]
    
    T = range(-1/2,1/2,S)
    p = plan_nfft(T, 2*B+1, reltol=1e-12)
    evaluated_convolution = [real.(p * c_m_ṽ) for c_m_ṽ in c_m_Ṽ]

    I = findmax(hcat(evaluated_convolution...))[2]
    ṽ_max = V[I[2]]

    T_fine = range(T[I[1]-1], T[I[1]+1], S)
    p_fine = plan_nfft(T_fine, 2*B+1, reltol=1e-12)
    evaluated_convolution_fine = real.(p_fine * c_m_Ṽ[I[2]])
    
    I_fine = findmax(evaluated_convolution_fine)
    t_max = T_fine[I_fine[2]]
    
    return rot_u_to_v(ũ, ṽ_max, t_max)
end

function choose_ũ_that_fulfills_polarangle_condition(U::Vector{S2Point})
    N = length(U)
    U_rotated = rot_to_northpole(S2Point(normalize(sum(U)))) * U
    for i in 1:N
        u_rotated = U_rotated[i]
        number_of_identical_polar_angles = sum([u_rotated.θ ≈ u2_rotated.θ for u2_rotated in U_rotated]) - 1
        if (number_of_identical_polar_angles == 0) && !((u_rotated.θ ≈ 0) || (u_rotated.θ ≈ pi))
            return U[i]
        end
    end
    error("U does not fullfill the Polarangel-Condition.")
end

function FindRotation(U::Vector{S2Point}, V::Vector{S2Point}) # finds Rotation that fullfills V = R * U
    ũ = choose_ũ_that_fulfills_polarangle_condition(U)
    
    mu_ũ_hat =  S2TransformS2Point(rot_to_northpole(ũ) * U, 1)
    nu_Ṽ_hat = [S2TransformS2Point(rot_to_northpole(v) * V, 1) for v in V]

    c_0_Ṽ = [conj(nu_ṽ_hat[(0,0)]) * mu_ũ_hat[(0,0)] + conj(nu_ṽ_hat[(1,0)]) * mu_ũ_hat[(1,0)] for nu_ṽ_hat in nu_Ṽ_hat]
    c_1_Ṽ = [conj(nu_ṽ_hat[(1,1)]) * mu_ũ_hat[(1,1)] for nu_ṽ_hat in nu_Ṽ_hat]
    
    Max_on_supp_Ṽ = 2 * abs.(c_1_Ṽ) + real.(c_0_Ṽ)
    index_ṽ_max = findmax(Max_on_supp_Ṽ)[2]
    
    c_1_ṽ_max = c_1_Ṽ[index_ṽ_max]
    t_max = atan(imag(c_1_ṽ_max),real(c_1_ṽ_max)) / (2*pi) 
    ṽ_max = V[index_ṽ_max]

    return rot_u_to_v(ũ, ṽ_max, t_max)
end

function FastFindRotation(U::Vector{S2Point}, V::Vector{S2Point}) # finds Rotation that fullfills V = R * U, but faster
    ũ = choose_ũ_that_fulfills_polarangle_condition(U)
    
    mu_ũ_hat = S2TransformS2Point(rot_to_northpole(ũ) * U, 1)
    nu_hat = S2TransformS2Point(V, 1)

    nu_Ṽ_hat = [SphericalHarmonicArrays.SHArray([nu_hat[1]; wignerD(1, ṽ.ϕ, ṽ.θ, 0)' * nu_hat[2:4]], (SphericalHarmonicModes.ML(0:1, -1:1),)) for ṽ in V]

    c_0_Ṽ = [conj(nu_ṽ_hat[(0,0)]) * mu_ũ_hat[(0,0)] + conj(nu_ṽ_hat[(1,0)]) * mu_ũ_hat[(1,0)] for nu_ṽ_hat in nu_Ṽ_hat]
    c_1_Ṽ = [conj(nu_ṽ_hat[(1,1)]) * mu_ũ_hat[(1,1)] for nu_ṽ_hat in nu_Ṽ_hat]
    
    Max_on_supp_Ṽ = 2 * abs.(c_1_Ṽ) + real.(c_0_Ṽ)
    index_ṽ_max = findmax(Max_on_supp_Ṽ)[2]
    
    c_1_ṽ_max = c_1_Ṽ[index_ṽ_max]
    t_max = atan(imag(c_1_ṽ_max),real(c_1_ṽ_max)) / (2*pi)
    ṽ_max = V[index_ṽ_max] 

    return rot_u_to_v(ũ, ṽ_max, t_max)
end;