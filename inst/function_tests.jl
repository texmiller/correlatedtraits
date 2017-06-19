# Simple tests for the functions used in the simulation
cd("D:\\GIT\\correlatedtraits\\inst")
@everywhere include("simulation_2.jl")
using Base.Test

# set some parameters for testing
n = 20                 # number of individuals
M = [log(6.1), -0.97]  # phenotype means
V = [22.0, 1.92]       # total phenotypic variances
ρ = [0.9, 0.0]         # additive and environmental covariances
H = [0.8, 0.8]         # trait heritabilities
K = 40.0               # patch carrying capacity
a = 1.0                # slope of logistic function
s = 2.2                # negative binomial scale parameter
p = 3                  # population number
r = 5                  # replicate number
ngens = 20             # number of generations to simulate
bat_time = time()      # the time that the batch process started

# set a random seed
srand(42)

popmatrix = init_inds(n, M, V, ρ, H)


# ------------------------------------------------------------------------------
# UTILITY FUNCTION TESTS -------------------------------------------------------
# ------------------------------------------------------------------------------

@test clocktime(4242.676) == "01:10:42:676"

# rep_each functions -----------------------------------------------------------
A = Vector{Int64}(7)

rep_each!(A, [1, 2, 3, 4], [1, 2, 0, 4]); @test A == [1, 2, 2, 4, 4, 4, 4]
@test rep_each([1.1, 2.2, 3.3], [1, 0, 3])        == [1.1, 3.3, 3.3, 3.3]
@test rep_eachi([1, 2, 3], [1, 0, 3])             == [1, 3, 3, 3]

# vector math functions --------------------------------------------------------
A = rand(100); B = rand(100)
C = Vector{Float64}(100); C[:] = A

@test vector_add(A, B)  == A .+ B
@test vector_mult(A, B) == A .* B
vector_scalar_mult!(A, 0.5); @test A == C .* 0.5

μ1  = rand(50) .+ 5.0
μ2  = rand(50) .+ 3.0
ix  = vcat(collect(1:50), collect(1:50))
D   = zeros(Float64, 50)

sum_prod_deviates!(D, μ1, μ2, A, B, ix);
@test D == sum(reshape((A .- vcat(μ1, μ1)) .* (B .- vcat(μ2, μ2)), 50, 2), 2)[:]

# functions that refer to reference vectors to make calculations ---------------
A = zeros(Int64, 3); B = zeros(Float64, 3)

running_count!(A, [1, 1, 1, 2, 2, 3]);               @test A == [3, 2, 1]
running_sum!(B, [1.1, 2.2, 1.1, 3.3], [1, 2, 1, 3]); @test B == [2.2, 2.2, 3.3]

map_vals!(B, [3.0, 2.0, 1.0], [3, 2, 1]);            @test B == [1.0, 2.0, 3.0]
@test map_vals([3.0, 2.0, 1.0], [3, 2, 1])                   == [1.0, 2.0, 3.0]
map_valsi!(A, [3, 2, 1], [3, 2, 1]);                 @test A == [1, 2, 3]
@test map_valsi([3, 2, 1], [3, 2, 1])                        == [1, 2, 3]

ix = Vector{Int64}(5); dict = Dict(1=>4, -1=>2, 2=>3, 0=>1)
# x only includes entries in dictionary
map_dict!(ix, [-1, -1, 0, 1, 2], dict); @test ix == [2, 2, 1, 4, 3]
# x includes some entries not in dictionary
map_dict!(ix, [-1, -3, 0, 1, 2], dict); @test ix == [2, 0, 1, 4, 3]

# functions that operate on arrays ---------------------------------------------
A = cor2cov([1.0 1.0; 1.0 1.0], [1, 2])
pos_def_check!(A); @test A == [1.0 1.996; 1.996 4.0] # A should changes
pos_def_check!(A); @test A == [1.0 1.996; 1.996 4.0] # A should not change

A = [3 1; 1 2; 2 4; 2 3; 4 5]
sortrowsi!(A); @test A == [1 2; 2 3; 2 4; 3 1; 4 5]

B = Array{Int64}(4, 2)
group_index!(B, A, collect(1:4), 4, 5); @test B == [1 1; 2 3; 4 4; 5 5]
