@doc raw"""
    Type: ComplexZonotope{N<:Real, M<:Number}(center::Vector{N}, generators::Matrix{M}, scale_vector::Vector{N})
Complex Zonotope is a set representation for some convex sets containing complex valued vectors.\
It is an extension of the simple zonotope representation to complex numbers.\
Its real projection can however be a non-polytopic set, unlike a simple zonotope which is polytopic.

Let c be a real vector, G be a complex matrix and s be a positive real vector. Then\
ComplexZonotope(c,G,s) = Set{c + G\zeta: \zeta is complex vector, |\zeta| <= s}

Constructor: ComplexZonotope(center,generators,scale_vector)
"""
struct ComplexZonotope{N<:Real, M<:Number}
    center::Vector{N} # center
    generators::Matrix{M} # generator matrix
    scale_vector::Vector{N} # scale vector
end

@doc raw"""
    rand(::ComplexZonotope, dim::Integer, numgen::Integer) 
computes a random complex zonotope of dimension <dim> and number of generators <numgen>.
"""
function Base.rand(::Type{ComplexZonotope}, dim::Integer, numgen::Integer)
    c = rand(dim) .- 0.5 # center
    G = (rand(dim, numgen) .- 0.5) + (rand(dim,numgen)*1im .-0.5im) # generator matrix
    s = abs.(rand(numgen)) # scale vector
    return ComplexZonotope(c,G,s)
end

@doc raw"""
    A*Z: where A::Matrix{L}, Z::ComplexZonotope{N,M}
computes the complex zonotope obtained by pre-multiplying the complex zonotope Z with A.
"""
function Base. *(A::Matrix{L}, Z::ComplexZonotope{N,M}) where {L<:Number, N<:Real, M<:Number}
    return ComplexZonotope(A*Z.center, A*Z.generators, Z.scale_vector)
end

function Base. *(lambda::Number, Z::ComplexZonotope{N,M}) where {N<:Real, M<:Number}
    return ComplexZonotope(lambda*Z.center, lambda*Z.generators, Z.scale_vector)
end

@doc raw"""
    inclusion_scale(Z1::ComplexZonotope{N1, M1}, Z2::ComplexZonotope{N2,M2}) where {N1<:Real, N2<:Real, M1<:Number, M2<:Number}
function computes the minimum amount of expansion required for a complex zonotope to contain another complex zonotope.

Let Z1::ComplexZonotope(G1,c1,s1), Z2::ComplexZonotope(G2,c2,s2)\
Then inclusion\_scale(Z1,c1,s1) = min lambda \in \reals\_{>=0} s.t.\
Z2 \subseteq ComplexZonotope(G1,c1,lambda*s1)
"""
function inclusion_scale(Z1::ComplexZonotope{N1, M1}, Z2::ComplexZonotope{N2,M2}; optimizer = Mosek.Optimizer) where {N1<:Real, N2<:Real, M1<:Number, M2<:Number}
    # get number of generators
    m1 = length(Z1.scale_vector) # number of generators in Z1
    m2 = length(Z2.scale_vector) # number of generators in Z2
    # initialize optimization problem
    problem = Model(optimizer)
    set_silent(problem) # turn off display
    # declare variables in optimization
    @variable(problem, Xre[1:m1, 1:m2])
    @variable(problem, Xim[1:m1, 1:m2])
    @variable(problem, Xabs[1:m1, 1:m2])
    @variable(problem, yre[1:m1])
    @variable(problem, yim[1:m1])
    @variable(problem, yabs[1:m1])
    @variable(problem, lambda)
    # declare constraints in optimization
    G1 = Z1.generators
    G2 = Z2.generators
    @constraint(problem, con1, real(G1)*Xre - imag(G1)*Xim .== real(G2)*diagm(Z2.scale_vector))
    @constraint(problem, con2, real(G1)*Xim + imag(G1)*Xre .== imag(G2)*diagm(Z2.scale_vector))
    @constraint(problem, con3[i=1:m1, j=1:m2], [Xabs[i,j], Xre[i,j], Xim[i,j]] in SecondOrderCone())
    @constraint(problem, con4[i=1:m1], [yabs[i], yre[i], yim[i]] in SecondOrderCone())
    @constraint(problem, con5, Xabs*ones(m2) + yabs .<= lambda*Z1.scale_vector )
    @constraint(problem, con6, lambda >= 0)
    # declare objective
    @objective(problem, Min, lambda)
    # optimize
    optimize!(problem)
    # return expansion value
    if termination_status(problem) == OPTIMAL
        return objective_value(problem)
    else
        println()
        error("problem status is ", termination_status(problem), ". Expansion of the complex zonotope can not satisfy the inclusion.\n Possible near singluar scaled generator matrix.  Try changing complex zonotope number generators.")
        return NaN
    end
end

@doc raw"""
    fast_inclusion_scale(Z1::ComplexZonotope{N1, M1}, Z2::ComplexZonotope{N2,M2}) where {N1<:Real, N2<:Real, M1<:Number, M2<:Number}
function computes the minimum amount of expansion required for a complex zonotope to contain another complex zonotope.

Let Z1::ComplexZonotope(G1,c1,s1), Z2::ComplexZonotope(G2,c2,s2)\
Then inclusion\_scale(Z1,c1,s1) = min lambda \in \reals\_{>=0} s.t.\
Z2 \subseteq ComplexZonotope(G1,c1,lambda*s1)
"""
function fast_inclusion_scale(Z1::ComplexZonotope{N1, M1}, Z2::ComplexZonotope{N2,M2}) where {N1<:Real, N2<:Real, M1<:Number, M2<:Number}
    transfermat = pinv(Z1.generators*diagm(Z1.scale_vector))*Z2.generators*diagm(Z2.scale_vector)
    return opnorm(transfermat, Inf)
end






