## Functions for plotting the support of (V conv U)

function mark_jump(x, y, z, threshold)
    new_x, new_y, new_z = [x[1]], [y[1]], [z[1]]
    for i in 2:length(x)
        distance = sqrt((x[i] - x[i-1])^2 + (y[i] - y[i-1])^2 + (z[i] - z[i-1])^2)
        if distance > threshold
            push!(new_x, NaN)
            push!(new_y, NaN)
            push!(new_z, NaN)
        end
        push!(new_x, x[i])
        push!(new_y, y[i])
        push!(new_z, z[i])
    end    
    return new_x, new_y, new_z
end

function draw_support(U, V, n, ax)

    for (v, i_v) in zip(V,1:length(V))
        for (u,i_u) in zip(U,1:length(U))

            a = u.ϕ;
            b = u.θ;
            x = v.ϕ;
            y = v.θ;
            c1 = cos(b) * cos(x) * sin(y);
            c2 = - sin(b) * cos(x) * cos(y);
            c3 = sin(b) * sin(x);
            c4 = cos(b) * sin(x) * sin(y);
            c5 = - sin(b) * sin(x) * cos(y);
            c6 = - sin(b) * cos(x);
            c7 = sin(b) * sin(y);
            c8 = cos(b) * cos(y);
            c9 = -cos(a) * cos(b) * sin(y);
            c10 = cos(a) * sin(b) * cos(y);
            c11 = - sin(a) * sin(y);
            c12 = -sin(a) * cos(b) * sin(y);
            c13 = sin(a) * sin(b) * cos(y);
            c14 = cos(a) * sin(y);

            
            R_13(t) = c1 + c2 * cos(t) + c3 * sin(t);
            R_23(t) = c4 + c5 * cos(t) + c6 * sin(t);
            R_33(t) = c7 * cos(t) + c8;
            R_31(t) = c9 * cos(t)+ c10 + c11 * sin(t);
            R_32(t) = c12 * cos(t) + c13 + c14 * sin(t);

            α(t) = mod(atan(R_23(t), R_13(t)), 2*pi)
            β(t) = mod(atan(sqrt(1 - (R_33(t))^2), R_33(t)), 2*pi)
            γ(t) = mod(atan(R_32(t), -R_31(t)), 2*pi)

            T = range(0,2*pi, n)
            α, β, γ = mark_jump(α.(T), β.(T), γ.(T), .1)
            ax.plot(α, β, γ, label=latexstring("\$\\mathcal{R}_{u_$(i_u), v_$(i_v)}\$"))
        end
    end
end

function find_support(u, v, n)

    a = u.ϕ;
    b = u.θ;
    x = v.ϕ;
    y = v.θ;
    c1 = cos(b) * cos(x) * sin(y);
    c2 = - sin(b) * cos(x) * cos(y);
    c3 = sin(b) * sin(x);
    c4 = cos(b) * sin(x) * sin(y);
    c5 = - sin(b) * sin(x) * cos(y);
    c6 = - sin(b) * cos(x);
    c7 = sin(b) * sin(y);
    c8 = cos(b) * cos(y);
    c9 = -cos(a) * cos(b) * sin(y);
    c10 = cos(a) * sin(b) * cos(y);
    c11 = - sin(a) * sin(y);
    c12 = -sin(a) * cos(b) * sin(y);
    c13 = sin(a) * sin(b) * cos(y);
    c14 = cos(a) * sin(y);
    
    R_13(t) = c1 + c2 * cos(t) + c3 * sin(t);
    R_23(t) = c4 + c5 * cos(t) + c6 * sin(t);
    R_33(t) = c7 * cos(t) + c8;
    R_31(t) = c9 * cos(t)+ c10 + c11 * sin(t);
    R_32(t) = c12 * cos(t) + c13 + c14 * sin(t);

    α(t) = mod(atan(R_23(t), R_13(t)), 2*pi)
    β(t) = mod(atan(sqrt(1 - (R_33(t))^2), R_33(t)), 2*pi)
    γ(t) = mod(atan(R_32(t), -R_31(t)), 2*pi)

    T = range(0,2*pi, n)

    return [[α(t), β(t), γ(t)] for t in T]
end;