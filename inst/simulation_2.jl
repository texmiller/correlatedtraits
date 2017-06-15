# ------------------------------------------------------------------------------
# --- DESCRIPTION OF MODEL -----------------------------------------------------
# ------------------------------------------------------------------------------
# Continuous space, discrete time model for examining stochastic evolutionary
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
	using Base.LinAlg.BLAS
end

# ------------------------------------------------------------------------------
# --- FUNCTION DEFINITIONS -----------------------------------------------------
# ------------------------------------------------------------------------------

### HELPER FUNCTIONS ###########################################################

# function for naming output files
function outname(H, ρ, p, r)
	return string(H[1], "_", H[2], "_", ρ[1], "_", Int(p), "_" , Int(r), "_R.csv")
end

# function for turning seconds into a HH:MM:SS:mmm format
function clocktime(elapsed)
	hh   = trunc(elapsed / 3600)                  # hours
	mm   = trunc((elapsed - hh * 3600) / 60)      # minutes
	ss   = trunc((elapsed - hh * 3600 - mm * 60)) # seconds
	mill = (elapsed - trunc(elapsed)) * 1000      # milliseconds
	out  = string((@sprintf "%02.0f" hh),  ":",
	              (@sprintf "%02.0f" mm),  ":",
							  (@sprintf "%02.0f" ss),  ":",
							  (@sprintf "%03.0f" mill));
	return out
end

### GROWTH FUNCTIONS ###########################################################

# Modified 'Skellam' (exponential) population growth ---------------------------
function growth(N::Float64, K::Float64, r::Float64)
	# Given [N, K, r] returns expected [per-capita growth rate]
  # N : number of individuals
	# K : carrying capacity
	# r : low-density per-capita growth rate

	return (K - K * (1.0 - r / K).^N) ./ N
end

### GENETIC VARIANCE FUNCTIONS #################################################

# Initialize population --------------------------------------------------------
function init_inds(n, spX, M, V, C, H)
	# Given [n, spX, M, V, C, H] returns [X, D, PD, K, PK, r, Pr], which is
	# referred to in other functions as 'popmatrix'.
  # n   : number of starting individuals
	# spX : absolute value of starting population's spatial boundaries
	# M   : vector of mean phenotypes [μD μr]
	# V   : vector of total phenotypic variance [VD Vr]
	# C   : vector of covariances [C_Dr]
	# H   : vector of heritabilities [h²D h²r]

	# X   : spatial location of each individual
	# D   : dispersal genotype
	# PD  : dispersal phenotype (accounts for environmental variance)
	# r   : low-density growth rate (r) genotype
	# Pr  : r phenotype (accounts for environmental variance)
	# ΣA  : additive genetic covariance matrix
	# ΣE  : environmental covariance matrix

	X  = rand(Uniform(-spX, spX), n)

  # additive genetic covariances
	CV = [0 C;
	      C 0]
	ΣA = eye(2) .* (V .* H) + CV           # additive genetic covariance matrix
	ΣE = eye(2) .* (V .* (1 - H))          # environmental covariance matrix
  G = rand(MvNormal(M, ΣA), n)'          # genotypes
	P = rand(MvNormal([0, 0], ΣE), n)' + G # phenotypes

	return(DataFrame(X=X, D=G[:, 1], PD=P[:, 1],
												r=G[:, 2], Pr=P[:, 2]))
end

# Calculate non-additive (environmental) phenotypic variance -------------------
function VE(H, V)
	# Given [H and V], returns [ΣE]
	# V   : vector of total phenotypic variance [VD VK Vr]
	# H   : vector of heritabilities [h²D h²K h²r]

	# ΣE  : environmental variance matrix

	return eye(2) .* (V .* (1 - H))
end

### FUNCTIONS FOR NEIGHBORHOOD-RELATED CALCULATIONS ############################

# Determine which spatial locations to sample from -----------------------------
function accordian_bins(A::Array{Float64,1}, bw::Float64)
	x = vcat(minimum(A), maximum(A))

	# make a scale so that accordian_bins changes appropriately with bw; when
	# bw = 0.25, "fixed" bins should be spaced 1 unit apart.
	scl = bw * 4

	# center
	C = collect(linspace(-50.0, 50.0, 101) * scl)

	# left side
	if x[1] < -1150 * scl
	  L  = collect(linspace(  x[1], x[1] + 99 * scl,  100))
	  LC = collect(linspace(L[100],            C[1], 1002))[2:1001]
	else
	  L  = collect(linspace(-1150 * scl, -1051 * scl,  100))
	  LC = collect(linspace(-1050 * scl,   -51 * scl, 1000))
	end

	# right side
	if x[2] > 1150 * scl
	  R  = collect(linspace(x[2] - 99 * scl, x[2],  100))
	  CR = collect(linspace(         C[101], R[1], 1002))[2:1001]
	else
	  R  = collect(linspace(1051 * scl, 1150 * scl,  100))
	  CR = collect(linspace(  51 * scl, 1050 * scl, 1000))
	end

	return vcat(L, LC, C, CR, R)
end

# Calculate distance between one location (b) and all other locations (A) ------
function x_calc!(tmpX::Array{Float64,1}, A::Array{Float64,1}, b::Float64)
	# tmpX : a pre-allocated array to store the result
	# A    : an array of all locations
	# b    : the focal location

	@fastmath @inbounds @simd for j = 1:length(tmpX)
		tmpX[j] = A[j] - b
	end
	nothing
end

# Calculate the weighted density of individuals wrt distance -------------------
function w_calc!(w::Array{Float64,1}, A::Array{Float64,1}, sd::Float64,
	               scl::Float64)
	# w   : the neighborhood-weighted density of each location in A
	# A   : distance between one location and all other locations (from x_calc!)
	# sd  : the bin-width of the Gaussian neighborhood
	# scl : a scaling factor for the neighborhood, given by pdf(Normal(0,sd),0)

	# The original code, which was simpler to write but much slower:
	#  Distributions.pdf!(w, Normal(0,sd), A)
  #  scale!(w,scl)

	# Manually calculating the normal pdf for A:
	C1 = -1 / (2 * sd * sd)            # A constant for the pdf's exponent
	C2 = scl / sqrt(2 * sd * sd * pi)  # A constant for the first term in the pdf
	@fastmath @inbounds @simd for j = 1:length(w)
		w[j] = A[j] * A[j] * C1          # calculate (A^2) * C1
		if w[j] < -15.5                  # 0 when density will be < than exp(-15.5)
			w[j] = 0.0
		else
			w[j] = exp(w[j]) * C2          # calculate C2 * exp((A^2) * C1)
		end
	end
	nothing
end

# Calculate the total density at a location ------------------------------------
function Nx_calc!(Nx::Array{Float64,1}, w::Array{Float64,1}, i::Int64)
	# Nx : the density at the location of individual i
	# w  : the weighted density of all individuals in the population wrt ind. i
	# i  : the index of the individual

	Nx[i] = sum(w)
	nothing
end

# Calculate the density at a location at a location, minus the individual itself
function D_calc!(D::Array{Float64,1}, w::Array{Float64,1}, i::Int64)
	# D : the density at the location of individual i minus the individual itself
	# w : the weighted density of all individuals in the population wrt ind. i
	# i : the index of the individual

	D[i] = sum(w) - w[i]
	nothing
end

# Calculate the mean genotype for each location --------------------------------
function Mean_calc!(Mean::Array{Float64,1}, R::Array{Float64,1},
	                   w::Array{Float64,1}, Nx::Float64, i::Int64)
	# Mean  : the average genotype in a location given by index i
	# R     : the value of genotypes for each index
	# w     : the weighted density of all individuals in the population wrt ind. i
	# Nx    : the density at the location of ind. i
	# i     : the index of the invididual

	Mean[i]  = (*(R', w) / Nx)[1] # sum of weighted RD, scaled by patch density
	nothing
end

# Calculate the mean phenotype for each location -------------------------------
function MeanP_calc!(MeanP::Array{Float64,1}, P::Array{Float64,1},
	                    w::Array{Float64,1}, Nx::Float64, i::Int64)
	# MeanP  : the average phenotype in a location given by index i
	# P      : the value of phenotypes for each index
	# w      : the weighted dens. of all individuals in the population wrt ind. i
	# Nx     : the density at the location of ind. i
	# i      : the index of the invididual

	MeanP[i]  = (*(P', w) / Nx)[1] # sum of weighted P, scaled by patch density
	nothing
end

# Calculate the covariance among breeding values for a location ----------------
function COV_calc!(COV::Array{Float64,1}, w::Array{Float64,1},
	                 R1::Array{Float64,1}, R2::Array{Float64,1},
									 Mean1::Float64, Mean2::Float64,
									 Nx::Float64, i::Int64, tmp::Array{Float64,1})
	# COV   : covariance among breeding values at location i
	# w     : the weighted density of all individuals in the population wrt ind. i
	# R1    : the genotype of trait 1 at index i
	# R2    : the genotype of trait 2 at index i
	# Mean1 : the average genotype of trait 1 in a location given by index i
	# Mean2 : the average genotype of trait 1 in a location given by index i
	# Nx    : the density at the location of ind. i
	# i     : the index of the invididual
	# tmp   : vector to store the output

	@fastmath @inbounds @simd for j = 1:length(tmp)
		tmp[j] = (R1[j] - Mean1) * (R2[j] - Mean2)
		tmp[j] = tmp[j] * w[j]
	end

	COV[i] = sum(tmp) / Nx
	nothing
end

### METRICS FUNCTIONS ##########################################################

# Calculate breeding value metrics wrt individuals -----------------------------
function metrics(popmatrix, bw)
	# Given   [popmatrix, bw]
	# returns [Nx MeanD Meanr V_D V_r C_Dr]
	# popmatrix : see init_inds
	# bw        : neighborhood (bin) width, Gaussian σ

	# Nx        : population density in continuous, Gaussian neighborhood
	# MeanD			: mean dispersal genotype in neighborhood
	# Meanr			: mean low-density growth rate genotype in neighborhood
	# V_D       : additive genetic variance in D
	# V_r       : additive genetic variance in r
	# C_Dr      : additive genetic covariance between D and r

  # If the population is extinct, return NULL
  if size(popmatrix, 1) < 1
  	return []
  end

  # Define within-function variables
  n   = size(popmatrix, 1)              # number of individuals in population
  x   = convert(Array, popmatrix[:X])   # individual locations
  RD  = convert(Array, popmatrix[:D])
	Rr  = convert(Array, popmatrix[:r])
	scl = 1 / pdf(Normal(0, bw), 0)       # a scaling factor for n'hood size

	# NOTE: For Gaussian distribution, ~95% of values are within 2σ of the mean.
	# When bw = 0.25, an individual at spatial location 0.5 will have a
	# neighborhood where ~95% of its neighborhood is within the bounds 0:1, which
	# roughly translates to a discrete-patch system.

  # Initialize pre-allocated arrays
	w     = Array{Float64}(n)             # density wrt individual
	Nx    = Array{Float64}(n)
	MeanD = Array{Float64}(n)
	Meanr = Array{Float64}(n)
	V_D    = Array{Float64}(n)
	V_r    = Array{Float64}(n)
	C_Dr   = Array{Float64}(n)
	tmpX  = Array{Float64}(length(x))     # temporary array for calculations

	# Calculations, looped over individuals
	# NOTE: densities (w) are calculated st solitary individuals experience w = 1
	for i = 1:n
		x_calc!(tmpX, x, x[i])
		w_calc!(w, tmpX, bw, scl)
		Nx_calc!(Nx, w, i)

		Mean_calc!(MeanD, RD, w, Nx[i], i)
	  Mean_calc!(Meanr, Rr, w, Nx[i], i)

	  COV_calc!(V_D,  w, RD, RD, MeanD[i], MeanD[i], Nx[i], i, tmpX)
	  COV_calc!(V_r,  w, Rr, Rr, Meanr[i], Meanr[i], Nx[i], i, tmpX)
	  COV_calc!(C_Dr, w, RD, Rr, MeanD[i], Meanr[i], Nx[i], i, tmpX)
	end

	out = convert(DataFrame, [Nx MeanD Meanr V_D V_r C_Dr])
  names!(out, [:Nx, :MeanD, :Meanr, :V_D, :V_r, :C_Dr])
	return out
end

# Calculate breeding value metrics wrt spatial location ------------------------
function sum_metrics(popmatrix, bw)
	# Given   [popmatrix, bw]
	# returns [b Nx MeanD Meanr MeanPD MeanPr V_D V_r C_Dr]
	# popmatrix : see init_inds
	# bw        : neighborhood (bin) width, Gaussian σ

	# Nx        : population density in continuous, Gaussian neighborhood
	# MeanD			: mean dispersal genotype in neighborhood
	# Meanr			: mean low-density growth rate genotype in neighborhood
	# MeanPD		: mean dispersal phenotype in neighborhood
	# MeanPr		: mean low-density growth rate phenotype in neighborhood
	# V_D       : additive genetic variance in D
	# V_r       : additive genetic variance in r
	# C_Dr      : additive genetic covariance between D and r

  # If the population is extinct, return NULL
  if(size(popmatrix, 1) < 1)
  	return []
  end

  # Define within-function variables
  n   = size(popmatrix,1)               # number of individuals in population
  x   = convert(Array, popmatrix[:X])   # individual locations
	RD  = convert(Array, popmatrix[:D])
	Rr  = convert(Array, popmatrix[:r])
	PD  = convert(Array, popmatrix[:PD])
	Pr  = convert(Array, popmatrix[:Pr])
	b   = accordian_bins(x, bw)           # spatial extent
	nb  = length(b)                       # number of bins
	scl = 1 / pdf(Normal(0, bw), 0)       # a scaling factor for n'hood size

  # Initialize pre-allocated arrays
	w      = Array{Float64}(length(x))    # density wrt bin
	Nx     = Array{Float64}(nb)
	MeanD  = Array{Float64}(nb)
	Meanr  = Array{Float64}(nb)
	MeanPD = Array{Float64}(nb)
	MeanPr = Array{Float64}(nb)
	V_D    = Array{Float64}(nb)
	V_r    = Array{Float64}(nb)
	C_Dr   = Array{Float64}(nb)
	tmpX   = Array(Float64, length(x))    # temporary array for calculations

	# Calculations, looped over bin
	for i = 1:nb
		x_calc!(tmpX, x, b[i])
		w_calc!(w, tmpX, bw, scl)
		Nx_calc!(Nx, w, i)

		Mean_calc!(MeanD, RD, w, Nx[i], i)
		Mean_calc!(Meanr, Rr, w, Nx[i], i)

		MeanP_calc!(MeanPD, PD, w, Nx[i], i)
		MeanP_calc!(MeanPr, Pr, w, Nx[i], i)

		COV_calc!(V_D,  w, RD, RD, MeanD[i], MeanD[i], Nx[i], i, tmpX)
		COV_calc!(V_r,  w, Rr, Rr, Meanr[i], Meanr[i], Nx[i], i, tmpX)
		COV_calc!(C_Dr, w, RD, Rr, MeanD[i], Meanr[i], Nx[i], i, tmpX)
	end

	out = convert(DataFrame, [b Nx MeanD Meanr MeanPD MeanPr V_D V_r C_Dr])
  names!(out, [:b, :Nx, :MeanD, :Meanr, :MeanPD, :MeanPr, :V_D, :V_r, :C_Dr])

  return out
end

# Calculate breeding value metrics wrt leading-edge (and zero) locations -------
function LE_sum_metrics(popmatrix, bw)
	# Given   [popmatrix, bw]
	# returns [b Nx MeanD Meanr MeanPD MeanPr V_D V_r C_Dr]
	# popmatrix : see init_inds
	# bw        : neighborhood (bin) width, Gaussian σ

	# Nx        : population density in continuous, Gaussian neighborhood
	# MeanD			: mean dispersal genotype in neighborhood
	# Meanr			: mean low-density growth rate genotype in neighborhood
	# MeanPD		: mean dispersal phenotype in neighborhood
	# MeanPr		: mean low-density growth rate phenotype in neighborhood
	# V_D       : additive genetic variance in D
	# V_r       : additive genetic variance in r
	# C_Dr      : additive genetic covariance between D and r

  # If the population is extinct, return NULL
  if(size(popmatrix, 1) < 1)
  	return []
  end

  # Define within-function variables
  n   = size(popmatrix, 1)               # number of individuals in population
  x   = convert(Array, popmatrix[:X])    # individual locations
	RD  = convert(Array, popmatrix[:D])
	Rr  = convert(Array, popmatrix[:r])
  PD = convert(Array, popmatrix[:PD])
	Pr = convert(Array, popmatrix[:Pr])
	b  = [findmin(popmatrix[:X])[1],       # leftmost density > 0
				findmin(abs(popmatrix[:X]))[1],  # density closest to X=0 (center)
				findmax(popmatrix[:X])[1]]       # rightmost density > 0
	nb  = 3                                # number of bins
	scl = 1 / pdf(Normal(0, bw), 0)        # a scaling factor for n'hood size

  # Initialize pre-allocated arrays
	w      = Array{Float64}(length(x))    # density wrt bin
	Nx     = Array{Float64}(nb)
	MeanD  = Array{Float64}(nb)
	Meanr  = Array{Float64}(nb)
	MeanPD = Array{Float64}(nb)
	MeanPr = Array{Float64}(nb)
	V_D    = Array{Float64}(nb)
	V_r    = Array{Float64}(nb)
	C_Dr   = Array{Float64}(nb)
	tmpX   = Array(Float64,length(x))     # temporary array for calculations

	# Calculations, looped over bin
	for i = 1:nb
		x_calc!(tmpX, x, b[i])
		w_calc!(w, tmpX, bw, scl)
		Nx_calc!(Nx, w, i)

		Mean_calc!(MeanD, RD, w, Nx[i], i)
		Mean_calc!(Meanr, Rr, w, Nx[i], i)

		MeanP_calc!(MeanPD, PD, w, Nx[i], i)
		MeanP_calc!(MeanPr, Pr, w, Nx[i], i)

		COV_calc!(V_D,  w, RD, RD, MeanD[i], MeanD[i], Nx[i], i, tmpX)
		COV_calc!(V_r,  w, Rr, Rr, Meanr[i], Meanr[i], Nx[i], i, tmpX)
		COV_calc!(C_Dr, w, RD, Rr, MeanD[i], Meanr[i], Nx[i], i, tmpX)
	end

	out = convert(DataFrame, [b Nx MeanD Meanr MeanPD MeanPr V_D V_r C_Dr])
	names!(out, [:b, :Nx, :MeanD, :Meanr, :MeanPD, :MeanPr, :V_D, :V_r, :C_Dr])
  return out
end

### LIFE HISTORY FUNCTIONS #####################################################

# Sexual mating ----------------------------------------------------------------
# NOTE: This mating assumes sampling with replacement
function mate(popmatrix, bw)
	# Given [popmatrix, bw], returns [M].
	# popmatrix : see init_inds
	# bw        : neighborhood (bin) width, Gaussian σ

	# M         : a mate for each individual

	# If the population is extinct, return NULL
  if size(popmatrix, 1) < 1
  	return []
  end

  # Define within-function variables
  n   = size(popmatrix, 1)              # number of individuals in pop
  x   = convert(Array,popmatrix[:X])    # individual locations
  scl = 1 / pdf(Normal(0, bw), 0)       # a scaling factor for n'hood size

	# Initialize pre-allocated arrays
	w      = Array{Float64, 1}(n)         # weighted dens. to all other inds.
	D      = Array{Float64}(n)            # n'hood density (excluding yourself)
	M      = Array{Int}(n)                # mate matches
	tmpX   = Array{Float64}(length(x))    # temporary array for calculations
	tmpvec = Vector{Float64}(n)           # temp. vector to store cumulative sums

	# Draw a random number for each individual
  rnum = rand(Uniform(0, 1), n)

	# Loop over all invividuals
	for i = 1:n
		x_calc!(tmpX, x, x[i])
		w_calc!(w, tmpX, bw, scl)
		D_calc!(D, w, i)

		# Matchmaker, matchmaker, make me a match
		if D[i] == 0   # if there's nobody near you...
			M[i] = -99   # ...no mate (results in mate-finding Allee effect)
		else
			scale!(w, 1 / D[i])                       # scale w by density for ind. i
			cumsum!(tmpvec, w, 1)                     # calculate cumsum of densities
			M[i] = searchsortedfirst(tmpvec, rnum[i]) # choose mate according to rnum
		end
	end
  return M
end

# Reproduction and dispersal ---------------------------------------------------
function repro_disp(popmatrix, H, V, K, bw)
	# Given [popmatrix(t), H, V, bw], returns [popmatrix(t+1)]
	# popmatrix : see init_inds
	# H   : vector of heritabilities [h²D h²r]
	# V   : vector of total phenotypic variance [VD Vr]
	# bw  : neighborhood (bin) width, Gaussian σ

	# repro_disp simulates:
	# (1) Density-dependent reproduction
	#	(2) Inheritence of dispersal genotype from parents to offspring
	#	(3) Offsrping dispersal

	# If the population is extinct, or if there is only 1 individual, return NULL
  if size(popmatrix,1) < 2
		return []
  end

	# (1) Density-dependent reproduction -----------------------------------------
  # a place to store the number of offspring produced by each parent
  Noff = Array{Int}(size(popmatrix, 1))

  # Calculate non-additive (environmental) phenotypic variances
  ΣE = VE(H, V)

  # Calculate traits through space
  mets = metrics(popmatrix, bw)

  # Density-dependent population growth (with demographic stochasticity)
  for i = 1:size(popmatrix, 1)
		println([mets[i,:Nx], K,	abs(popmatrix[i,:Pr])])
		Noff[i] = rand(Poisson(growth(mets[i,:Nx], K,	exp(popmatrix[i,:Pr]))))
  end

	println(Noff)
	println(length(Noff))
	println(sum(Noff))

  # If there are no offspring, return NULL
  if sum(Noff) == 0
		return []
  end

	# (2) Inheritence of dispersal phenotype from parents to offspring -----------

	# Find a mate for each individual
  m = mate(popmatrix, bw)

  # Create a new popmatrix (npm), accounting for the number of new offspring
  npm = DataFrame(Float64, sum(Noff), 6)
  names!(npm, [:pX, :D, :PD, :r, :Pr, :X])

	# initiate some variables
	ΣA  = Array{Float64}(2, 2)  # an additive genetic covariance matrix
  mpv = Vector{Float64}(2)    # a vector of midparent values
	c   = 0                     # a counter for indexing

	# loop over all of the parents
	for i = 1:size(popmatrix, 1)
		if m[i] == -99
			continue   # if thre's no mate, there are no offspring
    else
 			# loop over all of the offspring from each parent:
			for j = 1:Noff[i]

				# increment counter
				c += 1

				# offspring inherit location from their parents
        npm[c, :pX] = popmatrix[i, :X]

				# get midparent values (MPV) for each (mating) parent
        mpv[:] = ([popmatrix[  i,  :D] popmatrix[  i,  :r]]  +
				          [popmatrix[m[i], :D] popmatrix[m[i], :r]]) * 0.5

				# get the additive genetic covariance matrix
				ΣA[:, :] = [mets[i, :V_D]  mets[i, :C_Dr]
										mets[i, :C_Dr] mets[i, :V_r]]

				# calculate genotypes
				(npm[c, :D], npm[c,:r]) = rand(MvNormal(mpv, 0.5 * ΣA), 1)

				# calculate phenotypes
				(npm[c, :PD], npm[c,:Pr]) = rand(MvNormal([0, 0], ΣE), 1) +
				                                        [npm[c, :D], npm[c, :r]]

				# (3) Offspring dispersal ----------------------------------------------
				npm[c, :X] = rand(Normal(npm[c, :pX], exp(npm[c, :PD])), 1)[]
      end
    end
  end

  # Delete NA rows that result from building npm without knowing which matings
  # were unsuccessful (i.e., m[i] == -99).
  deleterows!(npm, find(isna(npm[:, Symbol("X")])))
  return npm[:, [:X, :D, :PD, :r, :Pr]]
end

### SIMULATION FUNCTIONS #######################################################

# Save output ------------------------------------------------------------------
function calcoutput(popmatrix, p, r, H, C, ρ, i, bw)
	# Given [popmatrix, p, r, H, C, ρ, i, bw], get summary stats for gen i
	# popmatrix : see init_inds
	# p         : population number for this simulation
	# r         : replicate number for this simulation
	# H         : heritability vector
	# C         : covariance vector
	# ρ         : correlation vector
	# i         : current generation
	# bw        : neighborhood (bin) width, Gaussian σ

	popsize = nrow(popmatrix)                 # population size
  sm      = sum_metrics(popmatrix, bw)      # summary metrics
	le      = LE_sum_metrics(popmatrix, bw)

	# summary data sepecific to simulation, generation, etc.
	tmp = repmat([p, r, bw, H[1], H[2], C, ρ, i, popsize]', 9, 1)

	# reshaped summary statistics from sum_metrics and LE_sum_metrics
	metrix = vcat(sm[:, :b]', sm[:, :Nx]', sm[:,:MeanD]', sm[:,:Meanr]', sm[:,:MeanPD]', sm[:,:MeanPr]', sm[:,:V_D]', sm[:,:V_r]', sm[:,:C_Dr]')
	le_met = vcat(le[:, :b]', le[:, :Nx]', le[:,:MeanD]', le[:,:Meanr]', le[:,:MeanPD]', le[:,:MeanPr]', le[:,:V_D]', le[:,:V_r]', le[:,:C_Dr]')

	return hcat(tmp, le_met, metrix)
end

function blankoutput(popmatrix, p, r, H, C, ρ, i, bw)
	# Given [popmatrix, p, r, H, C, ρ, i, bw], get summary stats for gen i
	# popmatrix : see init_inds
	# p         : population number for this simulation
	# r         : replicate number for this simulation
	# H         : heritability vector
	# C         : covariance vector
	# ρ         : correlation vector
	# i         : current generation
	# bw        : neighborhood (bin) width, Gaussian σ

	popsize = 0                    # population size
  mL      = 0                    # leftmost density > 0
  mR      = 0                    # rightmost density > 0
  m0      = 0                    # density closest to X=0 (center)
  sm      = repmat([0], 2301, 9) # summary metrics

	# summary data specific to simulation, generation, etc.
	tmp = repmat([p, r, bw, H[1], H[2], C, ρ, i, popsize]', 9, 1)

	# reshaped summary statistics
	metrix = vcat(sm')

  return hcat(tmp, metrix)
end

function timeoutoutput(popmatrix, p, r, H, C, ρ, i, bw)
	# Given [popmatrix, p, r, H, C, ρ, i, bw], get summary stats for gen i
	# popmatrix : see init_inds
	# p         : population number for this simulation
	# r         : replicate number for this simulation
	# H         : heritability vector
	# C         : covariance vector
	# ρ         : correlation vector
	# i         : current generation
	# bw        : neighborhood (bin) width, Gaussian σ

	popsize = -i                   # population size
  mL      = -999                 # leftmost density > 0
  mR      = -999                 # rightmost density > 0
  m0      = -999                 # density closest to X=0 (center)
  sm      = repmat([0], 2301, 9) # summary metrics

	# summary data specific to simulation, generation, etc.
	tmp = repmat([p, r, bw, H[1], H[2], C, ρ, i, popsize]', 9, 1)

	# reshaped summary statistics
	metrix = vcat(sm')

  return hcat(tmp, metrix)
end

# Run the simulation -----------------------------------------------------------
function runsim(popmatrix, n, spX, ngens, M, V, C, H, K, ρ, bw, p, r, TIME, outname = "_")
	# Given [popmatrix, n, spX, ngens, M, V, C, H, bw, p, r, TIME, outname]
	# returns changes in population size, extent, and genetic variance over time.
  # popmatrix : see init_inds
	# n  	      : initial number of starting individuals
	# spX       : initial spatial locations
	# ngens     : number of generations to simulate
	# M         : vector of mean phenotypes [μD μK μr]
	# V         : vector of total phenotypic variance [VD VK Vr]
	# C         : vector of covariances [C_DK C_Dr C_Kr]
	# H         : vector of heritabilities [h²D h²K h²r]
	# bw        : bin width
	# p         : population number for this generation
	# r         : replicate number for this simulation
	# TIME      : processing time
	# outname   : the output name for this simulation

	# start timers
	tic()
	simtime = time()
	flag = ""

	# set a unique seed based on H, p, and r, so simulations can be repeated.
	srand(Int(([1 10] .* H + (ρ + 1) * 100)[] * 10000 + 42 * p + r))

	# initialize output
	batchoutput = Array{Any}(ngens * 9, 2313)

	# loop over generations
  for i = 1:ngens

		# do reproduction and dispersal
		popmatrix = repro_disp(popmatrix, H, V, K, bw)

		# if population goes extinct, make output blank for remaining generations
    if size(popmatrix, 1) < 1
      for j = i:ngens
        batchoutput[(1:9) + (9) * (j - 1), :] = blankoutput(popmatrix, p, r, H, C, ρ, i, bw)
      end
		break
		# if rerunwrap.jl has been running for more than 16 hours, get out
	  elseif (time() - TIME) > 16 * 3600
			for j = i:ngens
				batchoutput[(1:9) + (9) * (j - 1), :] = timeoutoutput(popmatrix, p, r, H, C, ρ, i, bw)
				flag = "!BATCHRUNTIMEOUT!"
			end
		break
		# if this simulation has been running for more than 8 hours, stop it
	  elseif (time() - simtime) > 10 * 3600
			for j = i:ngens
				batchoutput[(1:9) + (9) * (j - 1), :] = timeoutoutput(popmatrix, p, r, H, C, ρ, i, bw)
				flag = "!SIMTIMEOUT!"
			end
		break
		# otherwise, calculate relevant output metrics for this generation
		else
      batchoutput[(1:9) + (9) * (i - 1), :] = calcoutput(popmatrix, p, r, H, C, ρ, i, bw)
    end
  end

	# track simulation info and runtimes
	print(  "rep ",       (@sprintf "%2.0f" r),
				", pop ",       (@sprintf "%2.0f" p),
				", h²D ",       (@sprintf "%0.3f" H[1]),
				", h²r ",       (@sprintf "%0.3f" H[2]),
				", ρ ",         (@sprintf "%0.3f" ρ),
				", bw ",        (@sprintf "%4.2f" bw),
				", time ",      (clocktime(toq())),
				", final pop ", (@sprintf "%-6.0f" batchoutput[end,9]),
				" | total time ", (clocktime(time()-TIME)),
				" | ", flag, "\n")

	return batchoutput
end
