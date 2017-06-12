import Distributions
import DataFrames

println("# of cores: ", Base.CPU_CORES, "; # of workers: ", nworkers())

# load simulation functions
@everywhere include("simulation.jl")

TIME = time()
tic()

## PARAMETER VALUES ############################################################

# parse input arguments
h2D = eval(parse(ARGS[1]))  # dispersal heritability parameter
pop = eval(parse(ARGS[2]))  # population number
rep = eval(parse(ARGS[3]))  # replicate number

# fixed parameters
μD   = log(2) #log(5)# mean dispersal trait value
n    = 20     # starting N of individs
spX  = 0.01   # starting +/- extent
VPD  = 0.5    # total phenotypic variance
ngen = 20     # number of generations
K    = 40     # carrying capacity

# Make vectors to store parameter combinations
vl = 70#105   # vector length

P0 = Array{Any}(vl)
P1 = Array{Any}(vl)
P2 = Array{Any}(vl)
P3 = Array{Any}(vl)
P4 = Array{Any}(vl)
P5 = Array{Any}(vl)
P6 = Array{Any}(vl)
P7 = Array{Any}(vl)
P8 = Array{Any}(vl)
P9 = Array{Any}(vl)
PX = Array{Any}(vl)
PY = Array{Any}(vl)
PT = Array{Any}(vl)

# Initialize a counter
c = 0

# set a unique seed based on h2D and pop, so simulations can be repeated with
# the same popmatrix
srand(Int(h2D*10000 + 42*pop))

# generate a single founding population
popmatrix = init_inds(n, spX, μD, h2D, VPD)

# loop over Allee threshold
for Nt in [5] #[2.5, 2, 3, 5]

	# set low-density growth rates
	ldgr = [ K/Nt*1.25, K/Nt, K/Nt*0.5, 1, 2/3, 1/3, 0.00]

	# loop over kurtosis
	for k in [10, 5, 0, -0.6, -1.2]
		# choose dispersal kernel
		if k >= 0
			dk = NSt(0,1,k) # normal & non-std. t-dist
		elseif k < 0
			dk = ssB(0,1,k) # beta distribution
		end

		# loop over low density growth rates
		for l in ldgr
			# set population growth rate
			gp = growth_paramsHS(l, Nt, K)

			# loop over bin width
			for bw in [0.2, 1]#[0.1, 0.2, 1]

				# update counter
				c += 1

				# make vectors of parameter combinations to simulate
				P0[c] = popmatrix
				P1[c] = n
				P2[c] = spX
				P3[c] = gp
				P4[c] = dk
				P5[c] = ngen
				P6[c] = μD
				P7[c] = h2D
				P8[c] = VPD
				P9[c] = bw
				PX[c] = pop
				PY[c] = rep
				PT[c] = TIME
			end
		end
	end
end

## SIMULATION ##################################################################

# run all parameter combinations in parallel
R = pmap((p0,p1,p2,p3,p4,p5,p6,p7,p8,p9,pX,pY,PT) ->
   runsim(p0,p1,p2,p3,p4,p5,p6,p7,p8,p9,pX,pY,PT),
	        P0,P1,P2,P3,P4,P5,P6,P7,P8,P9,PX,PY,PT);

## RECORD OUTPUT ###############################################################

# convert R to a single array (could be a better way to do this)
out = Array{Any}(size(R,1) * ngen * 5, size(R[1],2))

# fill 'out' with data from R
for i in 1:size(R,1)
	out[((i-1)*ngen*5+1) : i*ngen*5, :] = R[i]
end

# write output
writedlm(outname(h2D, pop, rep), out, '\,')
toc()
