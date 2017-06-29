import Distributions
import DataFrames

println("# of cores: ", Sys.CPU_CORES, "; # of workers: ", nworkers())

# load simulation functions
@everywhere include("simulation.jl")

## PARAMETER VALUES #####################################################
tic()

# parse input arguments
p = Int(eval(parse(ARGS[1])))  # population number
r = Int(eval(parse(ARGS[2])))  # replicate number

# varying parameters ----------------------------------------------------
H = [0.1, 0.5, 0.9]                        # heritabilities                     TODO: run beetle values
ρ = [-0.9, -0.5, -0.1, 0.0, 0.1, 0.5, 0.9] # covariances                        TODO: run beetle values

# fixed parameters ------------------------------------------------------
n = 20               # number of individuals
M = [1.63, 2.74]     # phenotype means
V = [0.40, 0.35]     # total phenotypic variances                               TODO: add additional fixed variances for appendix
b = 4.05/10.0        # Beverton-Holt parameter (re-scaled for simulation)       TODO: add additional values for appendix
ngens = 20           # number of generations to simulate
bat_time = time()    # batch processing time
c = 0                # loop counter

# Vectors to store parameter combinations
vl = length(H)^2 * length(ρ)^2

P0 = Array{Any}(vl)  # popmatrix
P1 = Array{Any}(vl)  # ngens
P2 = Array{Any}(vl)  # phenotype means
P3 = Array{Any}(vl)  # variances
P4 = Array{Any}(vl)  # correlations
P5 = Array{Any}(vl)  # heritabilities
P6 = Array{Any}(vl)  # Beverton-Holt parameter
P7 = Array{Any}(vl)  # population
P8 = Array{Any}(vl)  # rep
PT = Array{Any}(vl)  # bat_time

for h1 in H
	for h2 in H
		for ρ1 in ρ
			for ρ2 in ρ

				c += 1

				# set seed so that founding population can be replicated later
				srand(plant_seed([ρ1, ρ2], [h1, h2], r, p))

				P0[c] = init_inds(n, M, V, [ρ1, ρ2], [h1, h2]) # founding pop'n
				P1[c] = ngens        # ngens
				P2[c] = M            # phenotype means
				P3[c] = V            # variances
				P4[c] = [ρ1, ρ2]     # correlations
				P5[c] = [h1, h2]     # heritabilities
				P6[c] = b            # Beverton-Holt parameter
				P7[c] = p            # population
				P8[c] = r            # rep
				PT[c] = bat_time     # bat_time
			end
		end
	end
end

## SIMULATION ##################################################################

# run all parameter combinations in parallel
R = pmap((p0,p1,p2,p3,p4,p5,p6,p7,p8,PT) ->
   runsim(p0,p1,p2,p3,p4,p5,p6,p7,p8,PT),
	        P0,P1,P2,P3,P4,P5,P6,P7,P8,PT);

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
