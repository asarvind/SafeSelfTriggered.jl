A = [0 0 1 0;
0.0175 -1.273 -3.559 0;
-0.052 1.273 -2.661 0;
-5 1 0 0 ]

B = Matrix(transpose([0 0.085 21.79 0]))

K = [-0.7214 0.0445 -0.1873 0.2292]
E = 0.02*B
tmin = 0.001

Luw = SampledLinearSystem(A, B, E, K, tmin)

a = 0.5
T = [0.0 0.0 0 a 0; 0.0 0.0 0 -a 0]

# solver options instantiation
sopt = SolverOptions(order_invariant_zonotope=2, sequence_length = 20)
sopt.order_error_zonotope = 10
# compute scheduler and the time required for the offline computation.  # For the first time execution, compilation time is added. 
# To get the correct computation time, repeat the execution of this Jupyter cell.  
@time sdl = Scheduler(Luw, T; solver_options = sopt);

# randmly sampling a point inside the invariant set
sampling_scale = rand();
comb_state = sampling_scale*real.(sdl.invariant.generators*(rand(MersenneTwister(0), Complex{Float64}, size(sdl.invariant.scale_vector)) .*sdl.invariant.scale_vector))
x = comb_state[1:4]
u = [comb_state[5]]
# compute scheduling bound
tmax = sdl(x, u)
println("scheduled upper bound on update time = ", tmax)

verifun = generate_verifier(Luw, T, sdl.time_step)
verifun(x, u, tmax)
schfun_plot = plot_scheduling_function(sdl; title = "Underwater")

(schedules, t1, t2) = compare_computation_times(100, Luw, T, sdl; sampling_region_scale = 1);