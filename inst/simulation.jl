# ------------------------------------------------------------------------------
# --- DESCRIPTION OF MODEL -----------------------------------------------------
# ------------------------------------------------------------------------------
# Discrete space, discrete time model for examining stochastic evolutionary
# processes during range expansion.
#
# Original R Code by Ben Phillips (2015)
# Translated to Julia by Brad Ochocki (2016)
# Adapted for current project by Brad Ochocki (2017)

# ------------------------------------------------------------------------------
# --- PACKAGE MANAGEMENT -------------------------------------------------------
# ------------------------------------------------------------------------------
installPackages = false
updatePackages  = false
loadPackages    = true

# Install required packages ----------------------------------------------------
if installPackages
	Pkg.add("Distributions")
	Pkg.add("DataFrames")
end

# Update packages --------------------------------------------------------------
if updatePackages
	Pkg.update()
end

# Set packages to use ----------------------------------------------------------
if loadPackages
	using Distributions
	using DataFrames
end

# Include utility functions ----------------------------------------------------
include("utilities.jl")

# ------------------------------------------------------------------------------
# --- FUNCTION DEFINITIONS -----------------------------------------------------
# ------------------------------------------------------------------------------

### HELPER FUNCTIONS ###########################################################

# function for naming output files
function outname(p, r)
	return string("correlated_traits_", p, "_" , r, ".csv")
end

### GROWTH FUNCTIONS ###########################################################

# Beverton-Holt per-capits population growth -----------------------------------
function growth(N::Int64, r::Float64, b::Float64)
	# Given [N, r, b] returns expected [per-capita growth rate]
  # N : number of individuals
	# r : low-density per-capita growth rate phenotype
	# b : Beverton-Holt parameter

	# use a logistic link function to map growth phenotype to a growth rate that
	# is bounded between 0 and K.
  if N == 0
    return 0
  else
    return r / (1 + b .* N)
  end
end


### GENETIC VARIANCE FUNCTIONS #################################################

# Initialize population --------------------------------------------------------
function init_inds(n, M, V, ρ, H)
	# Given [n, M, V, C, H] returns [X, D, PD, K, PK, r, Pr], which is
	# referred to in other functions as 'popmatrix'.
  # n   : number of starting individuals
	# M   : vector of mean phenotypes [μD μr]
	# V   : vector of total phenotypic variance [VD Vr]
	# ρ   : vector of correlations [ρA ρE]
	# H   : vector of heritabilities [h²D h²r]

	# X   : spatial location of each individual
	# D   : dispersal genotype
	# PD  : dispersal phenotype (accounts for environmental variance)
	# r   : low-density growth rate (r) genotype
	# Pr  : r phenotype (accounts for environmental variance)
	# ΣA  : additive genetic covariance matrix
	# ΣE  : environmental covariance matrix

	# set all starting locations to zero
	X = zeros(Int64, n)

	ΣA = VA(V, ρ[1], H)
	ΣE = VE(V, ρ[2], H)
  G = rand(MvNormal(M, ΣA), n)'
	P = rand(MvNormal([0.0, 0.0], ΣE), n)' + G

	return DataFrame(X=X, D=G[:, 1], PD=P[:, 1], r=G[:, 2], Pr=P[:, 2])
end

# Calculate non-additive (environmental) phenotypic covariance matrix ----------
function VE(V, ρ, H)
	# Given [V, ρ, H], returns ΣE
	# V  : vector of total phenotypic variance [VD Vr]
  # ρ  : covariance between D and r
	# H  : vector of heritabilities [h²D h²r]

	# ΣE : environmental variance matrix

	return sqrt.((1.0 - H) .* V) * sqrt.((1.0 - H) .* V)' .* [1.0 ρ; ρ 1.0]
end

# Calculate additive genetic covariance matrix ---------------------------------
function VA(V, ρ, H)
	# Given [V, ρ, H], returns ΣA
	# V  : vector of total phenotypic variance [VD Vr]
  # ρ  : covariance between D and r
	# H  : vector of heritabilities [h²D h²r]

	# ΣA : environmental variance matrix

  return sqrt.(H .* V) * sqrt.(H .* V)' .* [1.0 ρ; ρ 1.0]
end

### FUNCTIONS FOR PATCH-RELATED CALCULATIONS ###################################

# Determine which spatial locations to sample from -----------------------------
function accordian_bins(A::Vector{Int64})

	ls = minimum(A)
	rs = maximum(A)

	# center
	C = collect(-50:1:50)

	# left side
	if ls < -1150
		L  = collect(ls:1:(ls + 99))
		LC = floor(Int64, linspace(L[100], -50, 1002)[2:1001])
	else
		L  = collect(-1150:1:-1051)
		LC = collect(-1050:1:-51)
	end

	# right side
	if rs > 1150
		R  = collect((rs - 99):1:rs)
		CR = floor(Int64, linspace(50, R[1], 1002)[2:1001])
	else
		R  = collect(1051:1:1150)
		CR = collect(51:1:1050)
	end

  pa    = Vector{Int64}(2301)
	pa[:] = vcat(L, LC, C, CR, R)

	return pa
end

### METRICS FUNCTIONS ##########################################################

# Calculate breeding value metrics wrt individuals -----------------------------
function metrics(popmatrix, M, phenotypes = false, patches = "occupied")
	# Given   [popmatrix, phenotypes, occupied]
	# returns [Nx MeanD Meanr V_D V_r C_Dr pID]
	# popmatrix : see init_inds
	# phenotypes: should metrics calculate phenotypes instead of genotypes?
	# patches   : one of the following:
	#             "occupied" : all patches occupied by at least one individual.
	#             "accordian": left edge, center, right edge, & some evenly-
	#                          distributed intermediate patches.
	#             "edges"    : left-most, center-most, & right-most patches only.

	# Nx        : population density in each patch
	# MeanD			: mean dispersal genotype/phenotype in each individual's patch
	# Meanr			: mean low-density growth rate genotype/phenotype
	# V_D       : variance in D
	# V_r       : variance in r
	# C_Dr      : covariance between D and r

  # If the population is extinct, return NULL
  ni = nrow(popmatrix)
	if ni < 1
   	return []
  end

	# get individual-level data
	x  = convert(Vector{Int64}, popmatrix[:X])

  # center traits at zero
	if phenotypes
    D = vector_add(convert( Vector{Float64}, popmatrix[:PD]),
		               -convert(Vector{Float64}, popmatrix[ :D]))
		r = vector_add(convert( Vector{Float64}, popmatrix[:Pr]),
		               -convert(Vector{Float64}, popmatrix[ :r]))
	else
		D = vector_scalar_add(convert(Vector{Float64}, popmatrix[:D]), -M[1])
		r = vector_scalar_add(convert(Vector{Float64}, popmatrix[:r]), -M[2])
	end

	# which patches should be analyzed?
	if patches == "occupied"
    pa = Vector{Int64}()
    sizehint!(pa, ni)
		pa = unique(x)
	elseif patches == "accordian"
		pa = accordian_bins(x)
	elseif patches == "edges"
		pa  = [minimum([minimum(x), -1.0]), 0, maximum([maximum(x), 1.0])]
	end

  # total number of patches
  np = length(pa)

	# make a patch ID dictionary to denote array indices for each patch location
  pID = Dict{Int64, Int64}()
  sizehint!(pID, np)
  for i = 1:np
    pID[pa[i]] = i
  end

  # map each individual's patch to an appropriate array index using pID
	ix = Vector{Int64}(ni)
  map_dict!(ix, x, pID)

  # arrays to store patch-specific metrics
  Nx     = zeros(Int64,   np)
  sumD   = zeros(Float64, np)
  sumr   = zeros(Float64, np)
  spd_D  = zeros(Float64, np)
  spd_r  = zeros(Float64, np)
  spd_Dr = zeros(Float64, np)

  # get patch-specific totals
  running_sum!(sumD, D, ix)
  running_sum!(sumr, r, ix)
  running_count!(Nx, ix)

  # calculate patch-specific means
  Mean_D = sumD ./ Nx
  Mean_r = sumr ./ Nx

  # patch-specific sums of the products of deviations from means
  # (used for calculating variances and covariances)
  sum_prod_deviates!(spd_D,  Mean_D, Mean_D, D, D, ix)
  sum_prod_deviates!(spd_r,  Mean_r, Mean_r, r, r, ix)
  sum_prod_deviates!(spd_Dr, Mean_D, Mean_r, D, r, ix)

  # calculate the uncorrected (co)variances.
  # using uncorrected formula since population mean is (exactly) known.
  V_D  = spd_D  ./ Nx
  V_r  = spd_r  ./ Nx
  C_Dr = spd_Dr ./ Nx

	# center traits on (expected) population mean
	if phenotypes
		genotypes = metrics(popmatrix, M, false, patches)
		vector_add!(Mean_D, convert(Vector{Float64}, genotypes[:Mean_D]))
		vector_add!(Mean_r, convert(Vector{Float64}, genotypes[:Mean_r]))
	else
		vector_scalar_add!(Mean_D, M[1])
		vector_scalar_add!(Mean_r, M[2])
	end

	return DataFrame(X=pa, Nx=Nx, Mean_D=Mean_D, Mean_r=Mean_r, V_D=V_D, V_r=V_r,
	                 C_Dr=C_Dr, pID=collect(pID))
end

### LIFE HISTORY FUNCTIONS #####################################################

# Random mating among individuals that share a patch ---------------------------
function find_mates!(m::Vector{Int64}, pai::Array{Int64, 2},
	                   sx::Array{Int64, 2}, np::Int64)
	# given [m, pai, sx, np] returns [m]
	# m   : a vector of mate matches for each individual
	# pai : an array to track starting and ending indices for individuals that
	#       share a patch. each row corresponds to a row in pa.
	# sx  : and array where the first column is a sorted vector of x and the
	#       second column provides the entry's original index in x.
	# np  : the number of unique occupied patches

	@inbounds for i = 1:np

		# potential mates in patch i
		pm  = sx[pai[i, 1] : pai[i, 2], 2]
		lpm = length(pm)

		if lpm == 1
			# lone individuals do not find mates
			m[sx[pai[i, 1], 2]] = -99
		else
			# individuals sample mates with replacement
			for j = 1:lpm
				m[sx[pai[i, 1] + (j - 1), 2]] = rand(pm[vcat(1:(j - 1), (j + 1) : lpm)])
			end
		end
	end
	nothing
end

# Sexual mating ----------------------------------------------------------------
function mate(popmatrix)
	# Given [popmatrix], returns [m].
	# popmatrix : see init_inds
	# m         : a vector of mate matches (sampled with replacement)

	# If the population is extinct, return NULL
	ni = nrow(popmatrix)
  if ni < 1
  	return []
  end

  # get individual locations
  x  = convert(Vector{Int64}, popmatrix[:X])

	# get a sorted vector of unique occupied patches
	pa = Vector{Int64}()
	sizehint!(pa, ni)
	pa = sort(unique(x))

  # total number of patches
  np = length(pa)

	# make a sorted vector of x and track the original indices
	sx = [x collect(1:ni)]
	sortrowsi!(sx)

	# an array to track starting and ending indices for individuals that share a
	# patch. each row corresponds to a row in pa.
	pai = Array{Int}(np, 2)
	group_index!(pai, sx, pa, np, ni)

  # mate sampling
	m  = Array{Int}(ni)
	find_mates!(m, pai, sx, np)

  return m
end

# Reproduction and dispersal ---------------------------------------------------
function reproduce(Nx::Vector{Int64}, r::Vector{Float64}, m::Vector{Int64},
	                 b::Float64)
  # Given [Nx, r, m, K, a], return [Noff]
	# Nx   : patch density for individual i
	# r    : growth rate for individual i
	# m    : vector of mates
	# b    : Beverton-Holt parameter

	# Noff : the number of offspring produced by individual i

	Noff = Vector{Int64}(length(Nx))
  for i = 1:length(Nx)
		if m[i] == -99
			Noff[i] = 0      # no mate, no offspring
		else
			Noff[i] = rand(Poisson(growth(Nx[i], exp.(r[i]), b)))
		end
  end
	return Noff
end

function disperse(x::Vector{Int64}, λ::Vector{Float64})
	# Given [x, λ], returns [y].
	# x : an integer vector of each individual's current spatial location
	# λ : the Poisson parameter

	# y : an integer vector of new spatial locations for each individual

	lx = length(x)
	y = Vector{Int64}(lx)
	n = Vector{Int64}(lx)

	# since the Poisson only gives positive numbers, flip a coin to determine
	# which direction each individual moves
	d = rand(Bernoulli(), lx)

	for i = 1:lx
		n[i] = rand(Poisson(λ[i]))
		y[i] = x[i] + (d[i] * 2 - 1) * n[i]
	end
	return y
end

function repro_disp(popmatrix, M, V, ρ, H, b, return_dx = false)
	# Given [popmatrix(t), H, V, bw], returns [popmatrix(t+1)]
	# popmatrix : see init_inds
  # V : vector of total phenotypic variance [VD Vr]
  # ρ : vector of additive and non-additive genetic correlations between D and r
	# H : vector of heritabilities [h²D h²r]
	# b : Beverton-Holt parameter

	# If the population is extinct, or if there is only 1 individual, return NULL
	ni = nrow(popmatrix)
	if ni < 2
		return []
  end

	# convert some popmatrix columns to standalone vectors
	x  = convert(Vector{Int64},   popmatrix[:X])
	D  = convert(Vector{Float64}, popmatrix[:D])
	r  = convert(Vector{Float64}, popmatrix[:r])
	Pr = convert(Vector{Float64}, popmatrix[:Pr])

	# Calculate traits through space
  mets = metrics(popmatrix, M)

	# map patch-specific metrics to individuals
  ix  = Vector{Int64}(ni)
	pID = Dict(convert(Vector, mets[:pID]))
	map_dict!(ix, x, pID)

	Nx   = map_valsi(convert(Vector{Int64},  mets[:Nx]),   ix)
	V_D  = map_vals(convert(Vector{Float64}, mets[:V_D]),  ix)
	V_r  = map_vals(convert(Vector{Float64}, mets[:V_r]),  ix)
	C_Dr = map_vals(convert(Vector{Float64}, mets[:C_Dr]), ix)

	# (1) Density-dependent reproduction -----------------------------------------

	m    = mate(popmatrix)         # find a mate (or not) for each individual
	Noff = reproduce(Nx, Pr, m, b) # number of offspring for each parent
	no   = sum(Noff)               # total number of offspring

  # If there are no offspring, return NULL
  if no == 0
		return []
  end

	# (2) Inheritence of genotypes from parents to offspring ---------------------

	# make vectors of values to describe the parents of each offspring
	did   = rep_eachi(collect(1:ni), Noff)  # 'dam' ids
	sid   = rep_eachi(m, Noff)              # 'sire' ids
	dx    = map_valsi(x, did)               # dam location
	dD    = map_vals(D,  did)               # dam dispersal genotype
	sD    = map_vals(D,  sid)               # sire dispersal genotypes
	dr    = map_vals(r,  did)               # dam growth genotype
	sr    = map_vals(r,  sid)               # sire growth genotype

	# get patch-level additive genetic covariance matrix elements for offspring
	oV_D  = rep_each(V_D,  Noff)            # patch-level variance in D
	oV_r  = rep_each(V_r,  Noff)            # patch-level variance in r
	oC_Dr = rep_each(C_Dr, Noff)            # patch-level covariance in D and r

	# multiply elements of additive genetic covariance matrix by 0.5; this is the
	# variance among full-siblings when sampling traits (Roughgarden 1979)
	vector_scalar_mult!(oV_D,  0.5)
	vector_scalar_mult!(oV_r,  0.5)
	vector_scalar_mult!(oC_Dr, 0.5)

  # calculate midparent values for each offspring
	mpvD = vector_add(dD, sD)
	mpvr = vector_add(dr, sr)
	vector_scalar_mult!(mpvD, 0.5)
	vector_scalar_mult!(mpvr, 0.5)

	# initiate some variables
	ΣE   = VE(V, ρ[2], H)        # non-additive (environmental) covariance
	ΣA   = Array{Float64}(2,  2) # additive genetic covariance (placeholder)
	oDr  = Array{Float64}(no, 2) # genotype deviates from mean
	oPDr = Array{Float64}(no, 2) # phenotype deviates from mean

	# loop over all offspring
	for i = 1:no

		ΣA[1, 1] = oV_D[i]
		ΣA[1, 2] = oC_Dr[i]
		ΣA[2, 1] = oC_Dr[i]
		ΣA[2, 2] = oV_r[i]

		# make sure ΣA is positive definite. If not, it's most likely because
		# the matrix is singular, so nudge the corresponding correlation matrix
		# slightly to make the matrix positive definite.
		pos_def_check!(ΣA)

		# sample additive and non-additive deviates
		oDr[i,  1:2] = rand(MvNormal(ΣA))
		oPDr[i, 1:2] = rand(MvNormal(ΣE))
	end

	# sum means and deviates to get trait values
	oD  = vector_add( oDr[:, 1], mpvD)
	or  = vector_add( oDr[:, 2], mpvr)
	oPD = vector_add(oPDr[:, 1], oD)
	oPr = vector_add(oPDr[:, 2], or)

	# (3) Offspring dispersal ----------------------------------------------------
  ox = disperse(dx, exp.(oPD))

	if return_dx
		return DataFrame(X=ox, dX=dx, D=oD, PD=oPD, r=or, Pr=oPr)
	else
		return DataFrame(X=ox, D=oD, PD=oPD, r=or, Pr=oPr)
	end
end

### SIMULATION FUNCTIONS #######################################################

# Save output ------------------------------------------------------------------
function generate_output(popmatrix, p, r, M, H, ρ, i, flag = "full")
	# Given [popmatrix, p, r, H, C, ρ, i, bw], generate summary stats for
	# generation i
	# popmatrix : see init_inds
	# p         : population number for this simulation
	# r         : replicate number for this simulation
	# H         : heritability vector
	# ρ         : correlation vector
	# i         : current generation

	# make blank dataframes for 'zero' and 'time' flags
	zero_DFa = DataFrame(X=-1150:1150,
	                     Nx=0,    Mean_D=0, Mean_r=0, V_D=0, V_r=0, C_Dr=0)
	zero_DFe = DataFrame(X=-1:1,
	                     Nx=0,    Mean_D=0, Mean_r=0, V_D=0, V_r=0, C_Dr=0)
	time_DFe = DataFrame(X=-1:1,
	                     Nx=-999, Mean_D=0, Mean_r=0, V_D=0, V_r=0, C_Dr=0)

	if flag == "full"
    # simulation is running as planned, get metrics
		ni = nrow(popmatrix)                           # population size
	  ag = metrics(popmatrix, M, false, "accordian") # genotype summary metrics
		ap = metrics(popmatrix, M, true,  "accordian") # phenotype summary metrics
		eg = metrics(popmatrix, M, false, "edges")     # edge patches gen. summary
		ep = metrics(popmatrix, M, true,  "edges")     # edge patches phen. summary
	elseif flag == "zero"
		# population has gone extinct, write a bunch of zeros
		ni = 0
	  ag = zero_DFa
		ap = zero_DFa
		eg = zero_DFe
		ep = zero_DFe
	elseif flag == "timeout"
		# simulation has timed out, write some nonsense-values
		ni = -i
	  ag = zero_DFa
		ap = zero_DFa
		eg = time_DFe
		ep = time_DFe
	end

	# parameter info
	parameters = repmat([p, r, M[1], M[2], H[1], H[2], ρ[1], ρ[2], i, ni]', 12, 1)

	# reshaped summary statistics from sum_metrics and LE_sum_metrics
  row_id    = ["patch", "N",
	             "Mean_D",  "Mean_r",  "V_D",  "V_r",  "C_Dr",
							 "Mean_PD", "Mean_Pr", "V_PD", "V_Pr", "C_PDr"]
	edges     = vcat(eg[:X]', eg[:Nx]',
	                 eg[:Mean_D]', eg[:Mean_r]', eg[:V_D]', eg[:V_r]', eg[:C_Dr]',
								   ep[:Mean_D]', ep[:Mean_r]', ep[:V_D]', ep[:V_r]', ep[:C_Dr]')
	accordian = vcat(ag[:X]', ag[:Nx]',
	                 ag[:Mean_D]', ag[:Mean_r]', ag[:V_D]', ag[:V_r]', ag[:C_Dr]',
								   ap[:Mean_D]', ap[:Mean_r]', ap[:V_D]', ap[:V_r]', ap[:C_Dr]')

	return hcat(row_id, parameters, edges, accordian)
end

# Run the simulation -----------------------------------------------------------
function runsim(popmatrix, ngens, M, V, ρ, H, b, p, r, bat_time)
	# Given [popmatrix, ngens, M, V, ρ, H, b, p, r, bat_time]
	# returns changes in population size, extent, and genetic variance over time.
  # popmatrix : see init_inds
	# n  	      : initial number of starting individuals
	# spX       : initial spatial locations
	# ngens     : number of generations to simulate
	# M         : vector of mean phenotypes [μD μK μr]
	# V         : vector of total phenotypic variance [VD VK Vr]
	# ρ         : vector of correlations [ρA ρE]
	# H         : vector of heritabilities [h²D h²r]
	# b         : Beverton-Holt parameter
	# p         : population number for this generation
	# r         : replicate number for this simulation
	# bat_time  : time that the current batch process started

	# start timers
	tic()
	sim_time = time()
	sim_flag = ""
	met_flag = "full"

	# set a unique seed based on H, p, and r, so simulations can be repeated.
	srand(plant_seed(ρ, H, r, p))

	# initialize output
	batchoutput = Array{Any}(ngens * 12, 2315)

	# loop over generations
  for i = 1:ngens

		# an array index for the output
		ind = (1:12) + (12) * (i - 1)

		if met_flag == "full"
			# do reproduction and dispersal
			popmatrix = repro_disp(popmatrix, M, V, ρ, H, b)
		end

		# check if future generations should be simulated
		if popmatrix == []
			# stop simulating if population went extinct
			met_flag = "zero"
		elseif (time() - bat_time) > 16.0 * 3600.0
			# stop simulating if the batch has run longer than 16 hours
			met_flag = "timeout"
			sim_flag = "!BATCH_RUN_TIMEOUT!"
		elseif (time() - sim_time) > 10.0 * 3600.0
			# stop simulating if the simulation has run longer than 10 hours
			met_flag = "timeout"
			sim_flag = "!SIM_RUN_TIMEOUT!"
		else
		end

		# record data
		batchoutput[ind, :] = generate_output(popmatrix, p, r, M, H, ρ, i, met_flag)
	end

	# track simulation info and run times
	print(  "rep ",       (@sprintf "%2.0f" r),
				", pop ",       (@sprintf "%2.0f" p),
				", h²D ",       (@sprintf "%0.2f" H[1]),
				", h²r ",       (@sprintf "%0.2f" H[2]),
				", ρA ",        (@sprintf "%5.2f" ρ[1]),
				", ρE ",        (@sprintf "%5.2f" ρ[2]),
				", time ",      (clocktime(toq())),
				", final pop ", (@sprintf "%-9.0f" batchoutput[end, 11]),
				" | total time ", (clocktime(time() - bat_time)),
				" | ", sim_flag, "\n")

	return batchoutput
end
