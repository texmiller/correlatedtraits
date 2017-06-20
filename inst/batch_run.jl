import Distributions
import DataFrames

println("# of cores: ", Sys.CPU_CORES, "; # of workers: ", nworkers())

# load simulation functions
@everywhere include("simulation.jl")

## PARAMETER VALUES ############################################################
tic()

# parse input arguments
p = Int(eval(parse(ARGS[1])))  # population number
r = Int(eval(parse(ARGS[2])))  # replicate number

# varying parameters -----------------------------------------------------------
H = [0.1, 0.5, 0.9]                        # heritabilities
ρ = [-0.9, -0.5, -0.1, 0.0, 0.1, 0.5, 0.9] # covariances

# fixed parameters -------------------------------------------------------------
n = 20                 # number of individuals
M = [log(6.1), -0.35]  # phenotype means
V = [0.2, 1.9]         # total phenotypic variances
K = 25.4               # patch carrying capacity
a = 1.0                # slope of logistic function
s = 2.2                # negative binomial scale parameter
ngens = 20             # number of generations to simulate
bat_time = time()      # batch processing time
c = 0                  # loop counter

# Vectors to store parameter combinations
vl = length(H)^2 * length(ρ)^2

P0 = Array{Any}(vl)  # popmatrix
P1 = Array{Any}(vl)  # ngens
P2 = Array{Any}(vl)  # phenotype means
P3 = Array{Any}(vl)  # variances
P4 = Array{Any}(vl)  # correlations
P5 = Array{Any}(vl)  # heritabilities
P6 = Array{Any}(vl)  # carrying capacity
P7 = Array{Any}(vl)  # slope of logistic
P8 = Array{Any}(vl)  # negative binomial scale parameter
P9 = Array{Any}(vl)  # population
PX = Array{Any}(vl)  # rep
PT = Array{Any}(vl)  # bat_time

for h1 in H
	for h2 in H
		for ρ1 in ρ
			for ρ2 in ρ

				c += 1

				# set seed so that founding population can be replicated later
				srand(plant_seed([ρ1, ρ2], [h1, h2], r, p))

				P0[c] = init_inds(n, M, V, [ρ1, ρ2], [h1, h2]) # founding population
				P1[c] = ngens        # ngens
				P2[c] = M            # phenotype means
				P3[c] = V            # variances
				P4[c] = [ρ1, ρ2]     # correlations
				P5[c] = [h1, h2]     # heritabilities
				P6[c] = K            # carrying capacity
				P7[c] = a            # slope of logistic
				P8[c] = s            # negative binomial scale parameter
				P9[c] = p            # population
				PX[c] = r            # rep
				PT[c] = bat_time     # bat_time

			end
		end
	end
end

## SIMULATION ##################################################################

# run all parameter combinations in parallel
R = pmap((p0,p1,p2,p3,p4,p5,p6,p7,p8,p9,pX,PT) ->
   runsim(p0,p1,p2,p3,p4,p5,p6,p7,p8,p9,pX,PT),
	        P0,P1,P2,P3,P4,P5,P6,P7,P8,P9,PX,PT);

## RECORD OUTPUT ###############################################################

# convert R to a single array (could be a better way to do this)
out = Array{Any}(size(R, 1) * ngens * 12, size(R[1], 2))

# fill 'out' with data from R
for i in 1:size(R, 1)
	out[((i -1) * ngens * 12 + 1) : i * ngens * 12, :] = R[i]
end

# write output
writedlm(outname(p, r), out, '\,')
toc()
