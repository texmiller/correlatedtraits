# Simple tests for the functions used in the simulation
cd("D:\\GIT\\correlatedtraits\\inst")
@everywhere include("simulation.jl")
using Base.Test

# set some parameters for testing
n = 20                 # number of individuals
M = [log(6.1), -0.97]  # phenotype means
V = [0.2, 1.9]         # total phenotypic variances
ρ = [0.9, 0.0]         # additive and environmental covariances
H = [0.8, 0.8]         # trait heritabilities
K = 40.0               # patch carrying capacity
a = 1.0                # slope of logistic function
s = 2.2                # negative binomial scale parameter
p = 3                  # population number
r = 5                  # replicate number
ngens = 20             # number of generations to simulate

srand(42)
# ------------------------------------------------------------------------------
# from utilities.jl ------------------------------------------------------------
# ------------------------------------------------------------------------------

# rep_each functions ###########################################################
A = Vector{Int64}(7)

rep_each!(A, [1, 2, 3, 4], [1, 2, 0, 4]); @test A == [1, 2, 2, 4, 4, 4, 4]
@test rep_each([1.1, 2.2, 3.3], [1, 0, 3])        == [1.1, 3.3, 3.3, 3.3]
@test rep_eachi([1, 2, 3], [1, 0, 3])             == [1, 3, 3, 3]

# vector math functions ########################################################
A = rand(100); B = rand(100)
C = Vector{Float64}(100); C[:] = A

@test vector_add(A, B)  == A .+ B
vector_add!(A,  B); @test A == C .+ B
vector_add!(A, -B); @test A == C
vector_scalar_add!(A,  0.5); @test A == C .+ 0.5
vector_scalar_add!(A, -0.5); @test A == C
@test vector_mult(A, B) == A .* B
vector_scalar_mult!(A, 0.5); @test A == C .* 0.5

μ1  = rand(50) .+ 5.0
μ2  = rand(50) .+ 3.0
ix  = vcat(collect(1:50), collect(1:50))
D   = zeros(Float64, 50)

sum_prod_deviates!(D, μ1, μ2, A, B, ix);
@test D == sum(reshape((A .- vcat(μ1, μ1)) .* (B .- vcat(μ2, μ2)), 50, 2), 2)[:]

# functions that refer to reference vectors to make calculations ###############
A = zeros(Int64, 3); B = zeros(Float64, 3)

running_count!(A, [1, 1, 1, 2, 2, 3]);               @test A == [3, 2, 1]
running_sum!(B, [1.1, 2.2, 1.1, 3.3], [1, 2, 1, 3]); @test B == [2.2, 2.2, 3.3]

map_vals!(B, [3.0, 2.0, 1.0], [3, 2, 1]);  @test B == [1.0, 2.0, 3.0]
@test map_vals([3.0, 2.0, 1.0], [3, 2, 1])         == [1.0, 2.0, 3.0]
map_valsi!(A, [3, 2, 1], [3, 2, 1]);       @test A == [1, 2, 3]
@test map_valsi([3, 2, 1], [3, 2, 1])              == [1, 2, 3]

ix = Vector{Int64}(5); dict = Dict(1=>4, -1=>2, 2=>3, 0=>1)
# x only includes entries in dictionary
map_dict!(ix, [-1, -1, 0, 1, 2], dict); @test ix == [2, 2, 1, 4, 3]
# x includes some entries not in dictionary
map_dict!(ix, [-1, -3, 0, 1, 2], dict); @test ix == [2, 0, 1, 4, 3]

# functions that operate on arrays #############################################
A = cor2cov([1.0 1.0; 1.0 1.0], [1, 2])
pos_def_check!(A); @test A == [1.0 1.996; 1.996 4.0] # A should change
pos_def_check!(A); @test A == [1.0 1.996; 1.996 4.0] # A should not change

A = [3 1; 1 2; 2 4; 2 3; 4 5]
sortrowsi!(A); @test A == [1 2; 2 3; 2 4; 3 1; 4 5]

B = Array{Int64}(4, 2)
group_index!(B, A, collect(1:4), 4, 5); @test B == [1 1; 2 3; 4 4; 5 5]

# Et cetera ####################################################################
@test clocktime(4242.676) == "01:10:42:676"

@test r1(pi) == 3.1
@test r1(-pi) == -3.1
@test r1(0) == 0

# test random seed generator
seed_test = Vector{Int64}(21 * 21 * 11 * 11 * 10 * 20)
c = 0
for i in vcat(-0.99, -0.9:0.1:0.9, 0.99)
  for j in vcat(-0.99, -0.9:0.1:0.9, 0.99)
    for k in vcat(0.01, 0.1:0.1:0.9, 0.99)
      for l in vcat(0.01, 0.1:0.1:0.9, 0.99)
        for m in 1:10
          for n in 1:20
            c += 1
            seed_test[c] = plant_seed([i j], [k l], m, n)
          end
        end
      end
    end
  end
end
@test length(seed_test) == length(unique(seed_test))

# ------------------------------------------------------------------------------
# from simulation.jl -----------------------------------------------------------
# ------------------------------------------------------------------------------

# helper functions #############################################################
@test outname(p, r) == "correlated_traits_3_5.csv"

# growth functions #############################################################
@test growth(0, K, M[2], a) == 0
@test growth(3, K, M[2], a)  ≈ 8.24977881
@test growth(10, K, M[2], a) ≈ 3.83924821
@test growth(50, K, M[2], a) ≈ 0.79999992

# genetic variance functions ###################################################

# test the popmatrix function --------------------------------------------------
pmx = init_inds(10^6, M, V, ρ, H)

# check basic qualities of the data frame
@test nrow(pmx) == 10^6
@test ncol(pmx) == 5
@test unique(pmx[:X]) == [0]

# check that mean phenotypes are correct
@test r1(mean(pmx[:PD])) ≈ r1(M[1])
@test r1(mean(pmx[:Pr])) ≈ r1(M[2])

# check that variances are correct
@test r1(var(pmx[:PD])) ≈ V[1]                            # total
@test r1(var(pmx[:Pr])) ≈ V[2]
@test r1(var(pmx[:D])) ≈ r1(V[1] * H[1])                  # additive
@test r1(var(pmx[:r])) ≈ r1(V[2] * H[2])
@test r1(var(pmx[:PD] - pmx[:D])) ≈ r1(V[1] * (1 - H[1])) # environmental
@test r1(var(pmx[:Pr] - pmx[:r])) ≈ r1(V[2] * (1 - H[2]))

# check that correlations are correct
@test r1(cor(pmx[:D], pmx[:r])) ≈ ρ[1]
@test r1(cor(pmx[:PD] - pmx[:D], pmx[:Pr] - pmx[:r])) ≈ ρ[2]

# test the covariance matrix functions -----------------------------------------
@test    VE(V, 0, [1, 1]) == zeros(Float64, 2, 2)
@test    VE(V, 1, [1, 1]) == zeros(Float64, 2, 2)
@test    VE(V, 0, [0, 0])  ≈ [V[1]   0; 0   V[2]]
@test r1(VE(V, 1, [0, 0])) ≈ [V[1] 0.6; 0.6 V[2]]

@test    VA(V, 0, [0, 0]) == zeros(Float64, 2, 2)
@test    VA(V, 1, [0, 0]) == zeros(Float64, 2, 2)
@test    VA(V, 0, [1, 1])  ≈  [V[1]   0; 0   V[2]]
@test r1(VA(V, 1, [1, 1])) == [V[1] 0.6; 0.6 V[2]]

# accordian_bins ###############################################################
@test length(accordian_bins(collect(-2000:2000))) == 2301
@test accordian_bins(collect(-2000:2000))[1:100]     == collect(-2000:-1901)
@test accordian_bins(collect(-2000:2000))[2202:2301] == collect(1901:2000)
@test accordian_bins(collect(-2000:2000))[1101:1201] == collect(-50:50)
@test accordian_bins(collect(-1:1)) == collect(-1150:1150)

# metrics ######################################################################

# occupied patches, genotype and phenotype metrics -----------------------------
pmx = init_inds(10^6, M, V, ρ, H)
A = metrics(pmx, M, false, "occupied")
@test nrow(A)   == 1
@test A[1, :X]  == 0
@test A[1, :Nx] == 10^6
@test r1(A[1, :Mean_D]) == r1(M[1])
@test r1(A[1, :Mean_r]) == r1(M[2])
@test r1(A[1, :V_D])    == r1(V[1] * H[1])
@test r1(A[1, :V_r])    == r1(V[2] * H[2])
@test r1(A[1, :C_Dr])   == r1(cov(pmx[:D], pmx[:r]))

A = metrics(pmx, M, true, "occupied")
@test nrow(A)   == 1
@test A[1, :X]  == 0
@test A[1, :Nx] == 10^6
@test r1(A[1, :Mean_D]) == r1(M[1])
@test r1(A[1, :Mean_r]) == r1(M[2])
@test r1(A[1, :V_D])    == r1(V[1] * (1 - H[1]))
@test r1(A[1, :V_r])    == r1(V[2] * (1 - H[2]))
@test r1(A[1, :C_Dr])   == r1(cov(pmx[:PD] - pmx[:D], pmx[:Pr] - pmx[:r]))

# accordian patches, genotype and phenotype metrics ----------------------------
A = metrics(pmx, M, false, "accordian")
@test nrow(A)   == 2301
@test A[:X]  == collect(-1150:1150)
@test A[1151, :Nx] == 10^6
@test A[vcat(collect(1:1150), collect(1152:2301)), :Nx] == zeros(Int64, 2300)
@test r1(A[1151, :Mean_D]) == r1(M[1])
@test r1(A[1151, :Mean_r]) == r1(M[2])
@test r1(A[1151, :V_D])    == r1(V[1] * H[1])
@test r1(A[1151, :V_r])    == r1(V[2] * H[2])
@test r1(A[1151, :C_Dr])   == r1(cov(pmx[:D], pmx[:r]))


A = metrics(pmx, M, true, "accordian")
@test nrow(A)   == 2301
@test A[:X]  == collect(-1150:1150)
@test A[1151, :Nx] == 10^6
@test A[vcat(collect(1:1150), collect(1152:2301)), :Nx] == zeros(Int64, 2300)
@test r1(A[1151, :Mean_D]) == r1(M[1])
@test r1(A[1151, :Mean_r]) == r1(M[2])
@test r1(A[1151, :V_D])    == r1(V[1] * (1 - H[1]))
@test r1(A[1151, :V_r])    == r1(V[2] * (1 - H[2]))
@test r1(A[1151, :C_Dr])   == r1(cov(pmx[:PD] - pmx[:D], pmx[:Pr] - pmx[:r]))

# edge patches, genotype and phenotype metrics ---------------------------------
A = metrics(pmx, M, false, "edges")
@test nrow(A)   == 3
@test A[:X]  == collect(-1:1)
@test A[2, :Nx] == 10^6
@test A[[1, 3], :Nx] == zeros(Int64, 2)
@test r1(A[2, :Mean_D]) == r1(M[1])
@test r1(A[2, :Mean_r]) == r1(M[2])
@test r1(A[2, :V_D])    == r1(V[1] * H[1])
@test r1(A[2, :V_r])    == r1(V[2] * H[2])
@test r1(A[2, :C_Dr])   == r1(cov(pmx[:D], pmx[:r]))

A = metrics(pmx, M, true, "edges")
@test nrow(A)   == 3
@test A[:X]  == collect(-1:1)
@test A[2, :Nx] == 10^6
@test A[[1, 3], :Nx] == zeros(Int64, 2)
@test r1(A[2, :Mean_D]) == r1(M[1])
@test r1(A[2, :Mean_r]) == r1(M[2])
@test r1(A[2, :V_D])    == r1(V[1] * (1 - H[1]))
@test r1(A[2, :V_r])    == r1(V[2] * (1 - H[2]))
@test r1(A[2, :C_Dr])   == r1(cov(pmx[:PD] - pmx[:D], pmx[:Pr] - pmx[:r]))

### LIFE HISTORY FUNCTIONS #####################################################

# mate pairing -----------------------------------------------------------------
pai = [1 1; 2 3; 4 4; 5 5]
sx  = [1 2; 2 3; 2 4; 3 1; 4 5]
m = Vector{Int64}(5)
find_mates!(m, pai, sx, 4); @test m == [-99, -99, 4, 3, -99]

pmx = init_inds(5, M, V, ρ, H)
pmx[:X] = [3, 1, 2, 2, 4]
@test mate(pmx) == [-99, -99, 4, 3, -99]

# reproduction -----------------------------------------------------------------
Nx = ones(Int64, 5)
r  = zeros(Float64, 5)
@test reproduce(Nx, r, m, K, a)[[1, 2, 5]] == [0, 0, 0]

Nx = ones(Int64, 10^6)
r = zeros(Float64, 10^6)
m = ones(Int64, 10^6)
@test r1(mean(reproduce(Nx, r, m, K, a))) ≈ r1(K ./ 2)
@test r1( var(reproduce(Nx, r, m, K, a))) ≈ r1(K ./ 2)

Nx = rep_eachi([2], [10^6])
r  = rep_each([M[2]], [10^6])
m  = ones(Int64, 10^6)
@test r1(mean(reproduce(Nx, r, m, K, a))) ≈ 9.5
@test r1( var(reproduce(Nx, r, m, K, a))) ≈ 9.5

# dispersal --------------------------------------------------------------------
# this requires fairly large sample sizes for the variance to work out correctly
x = zeros(Int64, 10^7)
μ = rep_each([exp.(M[1])], [10^7])
@test r1(mean(    disperse(x, μ, s)))  == 0
@test r1(mean(abs(disperse(x, μ, s)))) == r1(exp.(M[1]))
@test r1( var(abs(disperse(x, μ, s)))) == r1(exp.(M[1]) + exp.(M[1])^2 / s)

@test abs(minimum(disperse(x, μ, s))) < 75
@test abs(maximum(disperse(x, μ, s))) < 75

x = floor(Int64, vcat(1:(10^6)/2, 1:(10^6)/2))
μ = rep_each([exp.(M[1])], [10^6])
@test r1(mean(abs(disperse(x, μ, s)  .- x))) == r1(exp.(M[1]))

# combined reproduction and dispersal ------------------------------------------

# reduce sources of variance for testing
V1 = [0.0001, 0.0001]
H1 = [0.01,   0.01]
ρ1 = [0.0,    0.0]
pmx = init_inds(10^6, M, V1, ρ1, H1)
# ensure 2 individuals/patch so matings are successful
pmx[:X] = vcat(1:(nrow(pmx)/2), 1:(nrow(pmx)/2))
pmx1 = repro_disp(pmx, M, V1, ρ1, H1, K, a, s, true)

# test traits
@test r1(mean(pmx1[:PD])) == r1(mean(pmx[:PD]))
@test r1(mean(pmx1[:Pr])) == r1(mean(pmx[:Pr]))
@test r1(mean(pmx1[:D]))  == r1(mean(pmx[:D]))
@test r1(mean(pmx1[:r]))  == r1(mean(pmx[:r]))

# test reproduction
@test r1(nrow(pmx1) / nrow(pmx)) ≈ 9.5

# test dispersal
@test r1(mean(abs(pmx1[:X] - pmx1[:dX]))) ≈ r1(exp.(M[1]))
@test r1( var(abs(pmx1[:X] - pmx1[:dX]))) ≈ r1(exp.(M[1]) + exp.(M[1])^2 / s)
@test abs(maximum(pmx1[:X]) -  maximum(pmx[:X])) < 75
@test abs(minimum(pmx1[:X]) -  minimum(pmx[:X])) < 75

# more detailed test of traits using the usual parameters
pmx = init_inds(40*10000, M, V, ρ, H)
pmx[:X] = rep_eachi(vcat(1:10000), rep_eachi([40], [10000]))
pmx1 = repro_disp(pmx, M, V, ρ, H, K, a, s)

# test (co)variances
CA = cov(hcat(convert(Vector{Float64}, pmx1[:D]),
              convert(Vector{Float64}, pmx1[:r])))
ΣA = VA(V, ρ[1], H)
@test r1(CA[1]) == r1(ΣA[1])
@test r1(CA[2]) == r1(ΣA[2])
@test r1(CA[4]) == r1(ΣA[4])

CE = cov(hcat(convert(Vector{Float64}, pmx1[:PD] - pmx1[:D]),
              convert(Vector{Float64}, pmx1[:Pr] - pmx1[:r])))
ΣE = VE(V, ρ[2], H)
@test r1(CE[1]) == r1(ΣE[1])
@test r1(CE[2]) == r1(ΣE[2])
@test r1(CE[4]) == r1(ΣE[4])

@test r1(ΣA  ./ (ΣA .+ ΣE)) == r1(CA  ./ (CA .+ CE))

### SIMULATION FUNCTIONS #######################################################

# generate output --------------------------------------------------------------
i = p = r = 2
pmx = init_inds(20, M, V, ρ, H)
pars = [p, r, M[1], M[2], H[1], H[2], ρ[1], ρ[2], i]

A  = generate_output(pmx, p, r, M, H, ρ, i, "full")
@test size(A)    == (12, 2315)
@test A[1, 2:11] == vcat(pars, nrow(pmx))
@test A[2, 1165] == nrow(pmx)
@test A[2, 13]   == 20

A  = generate_output(pmx, p, r, M, H, ρ, i, "zero")
@test size(A)       == (12, 2315)
@test A[1, 2:11]    == vcat(pars, 0)
@test A[2, 12:2315] == rep_each([0.0], [2304])
@test A[2, 13]      == 0

A  = generate_output(pmx, p, r, M, H, ρ, i, "timeout")
@test size(A)       == (12, 2315)
@test A[1, 2:11]    == vcat(pars, -i)
@test A[2, 15:2315] == rep_each([0.0], [2301])
@test A[2, 13]      == -999

#runsim ------------------------------------------------------------------------
p = r = 2
pmx = init_inds(20, M, V, ρ, H)
pars = [p, r, M[1], M[2], H[1], H[2], ρ[1], ρ[2]]

A = runsim(pmx, 5, M, V, ρ, H, K, a, s, p, r, time())
@test size(A)          == (60, 2315)
@test A[1:12:60, 2:9]  == repmat(pars', 5, 1)
@test A[1:12:60, 10]   == vcat(1.0:5.0)

A = runsim(pmx, 5, M, V, ρ, H, K, a, s, p, r, time() - 3600*40)
@test size(A)          == (60, 2315)
@test A[1:12:60, 2:9]  == repmat(pars', 5, 1)
@test A[1:12:60, 10]   == vcat(1.0:5.0)
@test A[1:12:60, 11]   == -vcat(1.0:5.0)
@test A[2, 15:2315]    == rep_each([0.0], [2301])
@test A[2, 12:14]      == [-999, -999, -999]

pmx = init_inds(0, M, V, ρ, H)
A = runsim(pmx, 5, M, V, ρ, H, K, a, s, p, r, time())
@test size(A)          == (60, 2315)
@test A[1:12:60, 2:9]  == repmat(pars', 5, 1)
@test A[1:12:60, 10]   == vcat(1.0:5.0)
@test A[1:12:60, 11]   == rep_each([0.0], [5])
@test A[2, 15:2315]    == rep_each([0.0], [2301])
@test A[2, 12:14]      == [0, 0, 0]
