# specify system 
A = [0 0 1 0;
     0.0175 -1.273 -3.559 0;
     -0.052 1.273 -2.661 0;
     -5 1 0 0 ]

B = Matrix(transpose([0 0.085 21.79 0]))

K = [-0.7214 0.0445 -0.1873 0.2292]
E = 0.02*B
tmin = 0.001
Lc = SelfTriggeredLinearControl(A, B, E, K, tmin)

# specify safety constraints
a = 0.5
T = [0.0 0.0 0 a 0; 0.0 0.0 0 -a 0]; 

# specify initial set
init_bound = zeros(5)

# compute offline quantities
sopt = SolverOptions(order_invariant_zonotope=2, sequence_length = 400)
sopt.order_error_zonotope = 10
@time offquants = OfflineQuantities(Lc, T; solver_options = sopt);
(schedules, number_steps,t1,t2) = compare_computation_times(100000, Lc, T, offquants; sampling_region_scale = 0.01);
n, m = size(Lc.B)
ngen = length(offquants.invariant.scale_vector)
M = length(offquants.state_scales)
println("no. arithmetic operations = ", (8*(n+m)+3)*ngen + 2*(log(2,M)+1))
println(t1, " ", t2, " ", t2/t1)
println()

