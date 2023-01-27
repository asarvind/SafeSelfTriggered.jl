module SafeSelfTriggered
#------------------------------------------------------------------------------------------------------------

export SolverOptions, SelfTriggeredLinearControl,  OfflineQuantities

using LinearAlgebra, LazySets # linear algebra and set based calculations
using JuMP, SCS, MathOptInterface, Mosek, MosekTools # optimization packages
using Random # random number generation
include("ComplexZonotope.jl") # calculations with complex zonotope

@doc raw"""

"""
struct SelfTriggeredLinearControl{N<:Real}
    A::Matrix{N} # state action matrix
    B::Matrix{N} # control action matrix 
    E::Matrix{N} # disturbance action matrix 
    K::Matrix{N} # feedback gain matrix 
    tmin::N # minimum time elapsed before trigger
end

@doc raw"""
Description of system:\
(t\_i)\_{i=0}^{\infty}: Sequence of impulse times\
t\_{i+1}-t\_{i} \in [tau_{\min},tau_{\max}]\
x(t): state of system at time t\
dx(t)/dt = Ax(t) + Ev(t) : \max_i |v(t)|_i \leq 1,  if t\in [t\_i,t\_{i+1})\
x(t\_{i}) = Rt\_{i}^{-} if i>1

Code macros:\
"A" -> A,  "R"->R, "E" -> E, "\mu" -> mu, "tau_{\max}" -> tmax, "tau_{\min}" -> tmin

Type: SelftriggeredLinearImpusive{N<:Number}[A::Matrix{N}, R::Matrix{N}, E::Matrix{N}, mu:Vector{N}, tmax::N, tmin::N]

Constructor:\
SelftriggeredLinearImpulsive(A, T, E, tmin)
"""
struct SelfTriggeredLinearImpulsive{N<:Real}
    A::Matrix{N}
    R::Matrix{N}
    E::Matrix{N}
    tmin::N
end

function convert(::Type{SelfTriggeredLinearImpulsive{<:Real}}, L::SelfTriggeredLinearControl{<:Real})
    n, m = size(L.B)
    l = size(L.E, 2)
    Ai = [L.A L.B; zeros(m, n) zeros(m, m)]
    Ri = [Matrix(1.0I, n, n) zeros(n, m); L.K zeros(m, m)]
    Ei = [L.E; zeros(m, l)]
    return SelfTriggeredLinearImpulsive(Ai, Ri, Ei, L.tmin)
end

@doc raw"""

"""
function bloat_input(L::SelfTriggeredLinearImpulsive{N}, t) where {N<:Real}
    # get dimensions
    n = size(L.A,1)
    m = size(L.E,2)
    # vector radii for computing box hull
    r1 = t*abs.(L.E)*ones(m)
    r2 = 0.5*(t^2)*abs.(L.A*L.E)*ones(m)
    r3 = 0.5*(t^3)*abs.(((L.A)^2)*L.E)*ones(m)
    # radii of box hull
    r = max.(r1, r2, r3)
    return Zonotope(zeros(n), diagm(r))
end

@doc raw"""

"""
mutable struct OfflineComputationOptions
    upper_bound_expansion_rate::Real # initial bound on invariant expansion rate w.t.t state action while doing binary search
    samples_time_gap::Real # same as time step below
    number_samples::Integer # order of complex zonotope
    time_step::Real # time gap between reachability sequences
    number_steps::Integer # lenght of sequence length
    error_zonotope_order::Integer # reduction order of error zonotope during scalar sequence computation
    minimum_expansion_rate_search_gap::Real # tolerance stop criterion while minimizing expansion rate
    maximum_iterations::Integer # maximum number of iterations in invariant computation for minimizing expansion rate
    expansion_index::Integer # while minimizing expansion rate in invariant computation, modify the time step
end

@doc raw"""

"""
function OfflineComputationOptions(L::SelfTriggeredLinearImpulsive{<:Real})
    upper_bound_expansion_rate = 100.0
    samples_time_gap = L.tmin
    number_samples = 2
    time_step = L.tmin
    number_steps = 100
    error_zonotope_order = 10
    minimum_expansion_rate_search_gap = 0.01
    maximum_iterations = 100
    expansion_index = 1
    return OfflineComputationOptions(upper_bound_expansion_rate, samples_time_gap, number_samples, time_step, number_steps, error_zonotope_order, minimum_expansion_rate_search_gap, maximum_iterations, expansion_index)
end

mutable struct SolverOptions
    sequence_length::Real
    order_invariant_zonotope::Integer
    order_error_zonotope::Integer
    optimizer::Function
    function SolverOptions(; sequence_length::Integer=100, order_invariant_zonotope::Integer=3, order_error_zonotope::Integer=10, optimizer::Function=Mosek.Optimizer)
        return new(sequence_length, order_invariant_zonotope, order_error_zonotope, optimizer)
    end
end

function OfflineComputationOptions(L::SelfTriggeredLinearImpulsive, sopt::SolverOptions)
    upper_bound_expansion_rate = 100.0
    samples_time_gap = L.tmin
    number_samples = min(sopt.order_invariant_zonotope, 1)
    time_step = L.tmin
    number_steps = sopt.sequence_length
    error_zonotope_order = sopt.order_error_zonotope
    minimum_expansion_rate_search_gap = 0.01
    maximum_iterations = 100
    expansion_index = 1
    return OfflineComputationOptions(upper_bound_expansion_rate, samples_time_gap, number_samples, time_step, number_steps, error_zonotope_order, minimum_expansion_rate_search_gap, maximum_iterations, expansion_index)   
end

function check_invariant_zero_centered(L::SelfTriggeredLinearImpulsive{<:Real}, Zinit::ComplexZonotope{N,M}, T::Matrix{N}; options::OfflineComputationOptions = OfflineComputationOptions(L), optimizer = Mosek.Optimizer) where {N<:Real, M<:Number}
    # compute state action and input action matrces for invariant computation
    Zbloatinp = bloat_input(L, L.tmin) # input added in each iteration
    W = Zbloatinp.generators # input action matrix
    # compute template of complex zonotope
    n = size(L.A, 1)
    G = Matrix(1.0I, n, n)
    for i = 0:options.number_samples - 1        
        G = [G eigvecs(L.R*exp(L.A*(i*options.samples_time_gap + L.tmin)))]
    end
    # get dimensions in optimization
    m1 = size(G, 2)
    m2 = size(G, 2) + size(W, 2)
    m3 = size(Zinit.generators, 2)
    # create variables
    problem = Model(optimizer) # instantiate problem
    set_silent(problem) # turn off display
    @variable(problem, s[1:m1])
    @variable(problem, Xre[1:m1,1:m2])
    @variable(problem, Xim[1:m1,1:m2])
    @variable(problem, Xabs[1:m1,1:m2])
    @variable(problem, Yre[1:m1,1:m3])
    @variable(problem, Yim[1:m1,1:m3])
    @variable(problem, Yabs[1:m1,1:m3])
    @variable(problem, yre[1:m1])
    @variable(problem, yim[1:m1])
    @variable(problem, yabs[1:m1])
    @variable(problem, Zre[1:m1,1:m1])
    @variable(problem, Zim[1:m1,1:m1])
    @variable(problem, Zabs[1:m1,1:m1])
    @variable(problem, lambda)
    # create constraints
    M1 = L.R*exp(L.A*L.tmin)*G # renaming
    M2 = exp(L.A*L.tmin*options.expansion_index)*G # renaming
    M3 = Zinit.generators*diagm(Zinit.scale_vector)  # renaming
    @constraint(problem, con1re, [real(M1)*diagm(s) W] .== real(G)*Xre .- imag(G)*Xim )
    @constraint(problem, con1im, [imag(M1)*diagm(s) imag(W)] .== real(G)*Xim .+ imag(G)*Xre )
    @constraint(problem, con2[i=1:m1, j=1:m2], [Xabs[i,j], Xre[i,j], Xim[i,j]] in SecondOrderCone())
    @constraint(problem, con3, Xabs*ones(m2) .<= s)
    @constraint(problem, con4re, Zinit.center .== real(G)*yre)
    @constraint(problem, con4im, imag(Zinit.center) .== imag(G)*yim)
    @constraint(problem, con5re, real(M3) .== real(G)*Yre .- imag(G)*Yim)
    @constraint(problem, con5im, imag(M3) .== real(G)*Yim .+ imag(G)*Yre)
    @constraint(problem, con6[i=1:m1, j=1:m3], [Yabs[i,j], Yre[i,j], Yim[i,j]] in SecondOrderCone())
    @constraint(problem, con7[i=1:m1], [yabs[i], yre[i], yim[i]] in SecondOrderCone())
    @constraint(problem, con8, Yabs*ones(m3) .+ yabs .<= s) 
    @constraint(problem, con9, abs.(T*G)*s .+ L.tmin*abs.(T*L.A*G)*s .<= ones(size(T,1))) # safety bounds
    @constraint(problem, con10re, real(M2)*diagm(s) .== real(G)*Zre .- imag(G)*Zim )
    @constraint(problem, con10im, imag(M2)*diagm(s) .== real(G)*Zim .+ imag(G)*Zre )
    @constraint(problem, con11[i=1:m1, j=1:m1], [Zabs[i,j], Zre[i,j], Zim[i,j]] in SecondOrderCone())
    @constraint(problem, con12, Zabs*ones(m1) .<= options.upper_bound_expansion_rate*s)  
    optimize!(problem)
    return termination_status(problem), ComplexZonotope(zeros(n), G, value.(s))
end

@doc raw"""

"""
function large_invariant_zero_centered(L::SelfTriggeredLinearImpulsive{<:Real}, Zinit::ComplexZonotope{N,M}, T::Matrix{N}; options::OfflineComputationOptions = OfflineComputationOptions(L), optimizer = Mosek.Optimizer) where {N<:Real, M<:Number}
    # minimize expansion rate with respect to action with exp(L.A*L.tmin) while guaranteeing invariance with respect to action with action with L.R*exp(L.A*L.tmin)
    iter = 0
    lowrate = 0
    uprate = options.upper_bound_expansion_rate
    currate = uprate
    prec = uprate - lowrate
    newoptions = deepcopy(options)
    while iter <= options.maximum_iterations && prec >= options.minimum_expansion_rate_search_gap
        status, _ = check_invariant_zero_centered(L, Zinit, T; options=newoptions, optimizer = optimizer)
        if status ==  OPTIMAL 
            uprate = currate
            currate = 0.5*(currate + lowrate)
            prec = uprate - lowrate
        else
            lowrate = currate
            currate = 0.5*(currate + uprate)
            prec = uprate - lowrate
        end
        iter += 1
        newoptions.upper_bound_expansion_rate = currate
        println("iteration = ", iter, " status = ", status, " search_gap = ", prec)
    end
    newoptions = deepcopy(options)
    newoptions.upper_bound_expansion_rate = uprate # new options with minimized expansion rate
    # compute invariant with the minimized expansion rate
    status, Zinv = check_invariant_zero_centered(L, Zinit, T; options = newoptions, optimizer = optimizer)
    # scale up the minimized invariant to compute a large invariant
    if status ==  OPTIMAL
        lambda = 1/maximum(abs.(T*Zinv.generators)*Zinv.scale_vector)
        println("invariant is scaled up ", lambda, " times")
        return ComplexZonotope(Zinv.center, Zinv.generators, Zinv.scale_vector*lambda)
    else
        error("could not find invariant complex zonotope")
    end
end

@doc raw"""
    error_scale_sequence(L::SelfTriggeredLinearImpulsive, Zinv::ComplexZonotope{N,M}, t::Real, l::Integer; order = 10) where {N<:Real, M<:Number}
computes sequence of contractions required for inclusion inside an invariant complex zonotope of a sequence of disturbance input sets involved in computing reachable sets in successive time intervals.
"""
function error_scale_sequence(L::SelfTriggeredLinearImpulsive{N}, Zinv::ComplexZonotope{N,M}; options::OfflineComputationOptions = OfflineComputationOptions(L), optimizer = Mosek.Optimizer) where {N<:Real, M<:Number}
    out = AbstractFloat[] # initialize vector of output sequence
    Zadd = bloat_input(L,options.time_step) # amount of input added in each iteration
    Zinp = (0*L.A)*Zadd # initialize current input set to origin
    X = exp(L.A*options.time_step) # operator at each step
    for i = 1:options.number_steps
        Zinp = concretize(X*Zinp + Zadd) # disturbance input between it and (i+1)t
        Zred = reduce_order(Zinp, options.error_zonotope_order) # reduce order of zonotope representing disturbance input
        Zcompinp = ComplexZonotope(Zred.center, Zred.generators, ones(size(Zred.generators,2))) # convert input zonotope to complex zonotope
        beta = inclusion_scale(Zinv, Zcompinp) # minimum contraction of invariant zonotope required for containment of input zonotope
        println("error scaling iteration =", i, " remaining = ", options.number_steps - i)
        if i == 1
            push!(out, beta)
        else
            push!(out, max(out[end], beta)) # the max operator ensures that sequence is increasing
        end
    end
    return out
end

@doc raw"""
    state_scale_sequence(L::SelfTriggeredLinearImpulsive, Zinv::ComplexZonotope{N,M}, t::Real, l::Integer) where {N<:Real, M<:Number}
computes sequence of contractions required for inclusion inside an invariant complex zonotope after transformation by a sequence of sets of transformation operators as {exp(A tau) | (i-1)t <= tau <= itau}_{i=1}^{infinity}. 
"""
function state_scale_sequence(L::SelfTriggeredLinearImpulsive{N}, Zinv::ComplexZonotope{N,M}; options::OfflineComputationOptions = OfflineComputationOptions(L), optimizer = Mosek.Optimizer) where {N<:Real, M<:Number}
    out = AbstractFloat[] # initialize vector of output sequence
    # compute operator vertices and error zonotope of approximate convex hull for exp(L.A*tau)x for 0<=tau<=t, x \in Zinv
    M1 = exp(L.A*options.time_step) # approximate operator 
    alphaerr = norm(transpose(pinv(Zinv.generators*Zinv.scale_vector)*(options.time_step*L.A*Zinv.generators*diagm(Zinv.scale_vector))), 1) # error to be added to approximation inclusion scale
    # iteration
    for i = 1:options.number_steps
        mat = exp(L.A*(i-1)*options.time_step)
        alpha = inclusion_scale(Zinv, mat*M1*Zinv, optimizer = optimizer) + alphaerr
        println("state action scaling iteration =", i, " remaining = ", options.number_steps - i)
        if i == 1
            push!(out, alpha)
        else
            push!(out, max(out[end], alpha)) # the max operator ensures that sequence is increasing
        end
    end
    return out
end

struct OfflineQuantities
    time_step::Real
    invariant::ComplexZonotope{<:Real, <:Number}
    scaling_matrix::Matrix{<:Number}
    state_scales::Vector{<:Real}
    error_scales::Vector{<:Real}
end

function OfflineQuantities(L::SelfTriggeredLinearImpulsive{N}, Zinit::ComplexZonotope{N,M}, T::Matrix{N}; options::OfflineComputationOptions = OfflineComputationOptions(L), optimizer = Mosek.Optimizer) where {N<:Real, M<:Number}
    # compute large invariant
    Zinv = large_invariant_zero_centered(L, Zinit, T; options)
    # compute state action contraction sequence
    statescales = state_scale_sequence(L, Zinv; options = options, optimizer = optimizer)
    # compute input error action contraction sequence
    errorscales = error_scale_sequence(L, Zinv; options = options, optimizer = optimizer)
    # compute matrix used for checking inclusion
    scaling_matrix = pinv(Zinv.generators*diagm(Zinv.scale_vector))
    # construct dictionary of offline quantities
    return OfflineQuantities(options.time_step, Zinv, scaling_matrix, statescales, errorscales)
end

@doc raw"""
    OfflineQuantities(L::SelfTriggeredLinearControl{N}, Zinit::ComplexZonotope{N,M}, T::Matrix{N}, solver_options::SolverOptions=SolverOptions()) where {N<:Real, M<:Number} -> ::OfflineQuantitites
offline pre-computation of quantities required for computing latter the online triggering time upper bound.
"""
function OfflineQuantities(Lc::SelfTriggeredLinearControl{N}, T::Matrix{N}, init_bound::Vector{N}=zeros(size(Lc.A,1)); solver_options::SolverOptions=SolverOptions()) where {N<:Real}
    L = convert(SelfTriggeredLinearImpulsive, Lc)
    Gentop = diagm(init_bound)
    Genbot = Lc.K*Gentop
    generators = vcat(Gentop, Genbot)
    Zinit = ComplexZonotope(zeros(size(generators,1)), generators, ones(size(generators,2)))
    options = OfflineComputationOptions(L, solver_options)
    return OfflineQuantities(L, Zinit, T; options=options, optimizer=solver_options.optimizer)
end


@doc raw"""
    trigger_time_upper_bound(L::SelfTriggeredLinearImpulsive{<:Real}, x::Vector{<:Real}, inctesstateaction::Matrix{<:Number}, t::Real, statescales::Vector{<:Real}, errorscales::Vector{<:Real})
computes the upper bound on the time for next impulse
"""
function trigger_time_upper_bound(L::SelfTriggeredLinearImpulsive{<:Real}, x::Vector{<:Real}, Q::OfflineQuantities)
    # contraction required for containing x
    cont = norm(Q.scaling_matrix*x, Inf)
    # upper bound on trigger time using line search
    upind = length(Q.state_scales)  # initialize upper index for line search
    lowind = 0 # initialize lower index for line search
    curind = 1  # initialize current index for line search
    while lowind < upind - 1 && curind > 0
        if Q.state_scales[curind]*cont + Q.error_scales[curind] <= 1
            lowind = curind
            curind = Int(floor(0.5*(lowind + upind)))
        else
            upind = curind
            curind = lowind
        end
    end
    return lowind*Q.time_step + L.tmin
end

@doc raw"""

"""
function verify(x::Vector{N}, stateaction::Matrix{N}, inputgens::Matrix{N}, num_steps::Integer,
                 T::Matrix{N}) where {N<:Number}
    adddirbounds = zeros(size(T,1)) # bounds to be added in each iteration along the directions of T
    addinpgens = deepcopy(inputgens) # initialize generators of additional input in each Girard style iteration 
    addstate = stateaction*x # initialize state in each Girard style iteration
    ones_vec = ones(size(inputgens,2)) # pre-iteration initialization of ones vector for computational efficiency
    out = true # initialize output 
    for i = 1:num_steps
        # add to previous direction bounds to compute new direction bounds
        adddirbounds += abs.(T*addinpgens)*ones_vec
        dirbounds = T*addstate + adddirbounds
        addstate = stateaction*addstate # new state to add
        addinpgens = stateaction*addinpgens # generators of new input to add
        if norm(dirbounds, Inf) > 1
            out = false
            break
        end
    end
    return out
end

@doc raw"""

"""
function compare_computation_times(number_points::Integer, L::SelfTriggeredLinearImpulsive{<:Real}, T::Matrix{<:Real}, Q::OfflineQuantities; sampling_region_scale::Real = 1.0, seed = 1)
    stateaction = exp(L.A*L.tmin) # state action matrix used in discrete time reachability
    Zinp = bloat_input(L, L.tmin) # input zonotope in discrete time reachability
    inputgens = Zinp.generators # generators of above input zonotope
    # generate random points inside Zinv
    randpoints = []
    for i = 1:number_points
        randpoints = push!(randpoints, sampling_region_scale*real.(Q.invariant.generators*(rand(MersenneTwister(seed), Complex{Float64}, size(Q.invariant.scale_vector)) .*Q.invariant.scale_vector)))
    end
    # compute schedules and number of multiples of minimum impulse time with precomputed reachability scales 
    schedules = []
    number_steps = []
    for x in randpoints
        t = trigger_time_upper_bound(L, x, Q)
        push!(schedules, t)
        push!(number_steps, Int(floor(t/L.tmin)))
    end
     # compute time required for computing schedules for all points by two approaches
    comp_time0 = @elapsed begin
        for x in randpoints
            trigger_time_upper_bound(L, x, Q)
        end
    end
    # compute time required for computing schedules for all points by two approaches
    comp_time1 = @elapsed begin
        for x in randpoints
            trigger_time_upper_bound(L, x, Q)
        end
    end
    comp_time2 = @elapsed begin
        for i in 1:number_points
            verify(randpoints[i], stateaction, inputgens, number_steps[i],
                 T)
        end
    end
    return (schedules, number_steps, min(comp_time0, comp_time1), comp_time2)
end

# end module------------------------------------------------------------------------------------------------
end
