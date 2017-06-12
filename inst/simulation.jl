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
	return string(H[1], "_", H[2], "_", H[3], "_",
	              ρ[1], "_", ρ[2], "_", ρ[3], "_",
								Int(p), "_" , Int(r), "_R.csv")
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

# Hockey-stick population growth -----------------------------------------------
# Hockey-stick growth parameters
immutable growth_paramsHS
	lambda::Float64  # low density growth rate
	Nt::Float64      # threshold dens. for low-density -> high density growth rate
	K::Float64       # carrying capacity
end

function hockey_stick(N, gp)
	# Given [N, lambda, Nt, K] returns [per-capita growth rate]
  # N      : number of individuals
	# lambda : low-density per-capita growth rate
  # Nt     : threshold density for low-density -> high density growth rate
  # K      : carrying capacity
	# NOTE: for no Allee effect, lambda = K/Nt

	if N <=0
		return 0
	elseif N <= gp.Nt
    return gp.lambda
  else
    return gp.K/N
  end
end

# Modified 'Skellam' (exponential) population growth ---------------------------
function hockey_stick(N, K, r)
	# Given [N, K, r] returns expected [per-capita growth rate]
  # N : number of individuals
	# K : carrying capacity
	# r : low-density per-capita growth rate

	return K - K * (1 - r / K)^N
end


### GENETIC VARIANCE FUNCTIONS #################################################

# Initialize population --------------------------------------------------------
function init_inds(n, spX, M, V, C, H)
	# Given [n, spX, μD, h2D, VPD] returns [X, D, K, r, PD, PK, Pr], which is
	# referred to in other functions as 'popmatrix'.
  # n   : number of starting individuals
	# spX : absolute value of starting population's spatial boundaries
	# M   : vector of mean phenotypes [μD μK μr]
	# V   : vector of total phenotypic variance [VD VK Vr]
	# C   : vector of covariances [C_DK C_Dr C_Kr]
	# H   : vector of heritabilities [h²D h²K h²r]

	# X   : spatial location of each individual
	# D   : dispersal genotype
	# PD  : dispersal phenotype (accounts for environmental variance)
	# K   : carrying capacity genotype
	# PK  : carrying capacity phenotype (accounts for environmental variance)
	# r   : low-density growth rate (ldgr) genotype
	# Pr  : ldgr phenotype (accounts for environmental variance)
	# ΣA  : additive genetic covariance matrix
	# ΣE  : environmental covariance matrix

	X  = rand(Uniform(-spX,spX), n)

	CV = [  0  C[1] C[2];
	      C[1]   0  C[3];
				C[2] C[3]   0]
	ΣA = eye(3) .* (V .* H) + CV
	ΣE = eye(3) .* (V .* (1 - H))
  G = rand(MvNormal(M, ΣA), n)'
	P = rand(MvNormal([0, 0, 0], ΣE), n)' + G

	return(DataFrame(X=X, D=G[:, 1], PD=P[:, 1],
	                      K=G[:, 2], PK=P[:, 2],
												r=G[:, 3], Pr=P[:, 3]))
end

# function init_inds(n, spX, μD, h2D, VPD)
# 	# Given [n, spX, μD, h2D, VPD] returns [X, D, PD], which is referred to in
# 	# future functions as 'popmatrix'.
#   # n   : number of starting individuals
# 	# spX : absolute value of starting population's spatial boundaries
# 	# μD  : mean dispersal genotype
# 	# h2D : dispersal heritability
# 	# VPD : phenotypic variance in dispersal
# 	# X   : spatial location of each individual
# 	# D   : dispersal genotype
# 	# PD  : dispersal phenotype (accounts for environmental variance)
#
# 	X  = rand(Uniform(-spX,spX), n)
#   D  = rand(Normal(μD,sqrt(h2D*VPD)), n)
#   PD = rand(Normal(0,sqrt((1-h2D)*VPD)), n) + D
# 	return(DataFrame(X=X, D=D, PD=PD))
# end

# Calculate non-additive (environmental) phenotypic variance -------------------
function VE(H,V)
	# Given [H and V], returns [ΣE]
	# V   : vector of total phenotypic variance [VD VK Vr]
	# H   : vector of heritabilities [h²D h²K h²r]
	return eye(3) .* (V .* (1 - H))
end

### FUNCTIONS FOR NEIGHBORHOOD-RELATED CALCULATIONS ############################

# Determine which spatial locations to sample from -----------------------------
function accordian_bins(A::Array{Float64,1}, bw::Float64)
	x = vcat(minimum(A), maximum(A))

	scl = bw * 4
	# make a scale so that accordian_bins changes appropriately with bw; when
	# bw = 0.25, "fixed" bins should be spaced 1 unit apart.

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
function w_calc!(w::Array{Float64,1}, A::Array{Float64,1},
                 sd::Float64, scl::Float64)
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

	Mean[i]  = (*(R',w)/Nx)[1] # sum of weighted RD, scaled by patch density
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
	# Given [popmatrix], returns [Nx, MeanD, and SDD].
	# popmatrix : see init_inds
	# bw        : neighborhood (bin) width, Gaussian σ
	# Nx        : population density in continuous, Gaussian neighborhood
	# MeanD			: mean dispersal genotype in neighborhood
	# SDD       : σ(D), standard deviation in dispersal breeding values

  # If the population is extinct, return NULL
  if size(popmatrix, 1) < 1
  	return []
  end

  # Define within-function variables
  n   = size(popmatrix,1)              # number of individuals in population
  x   = convert(Array,popmatrix[:X])   # individual locations
  RD  = convert(Array,popmatrix[:D])   # individual dispersal genotypes
	RK  = convert(Array,popmatrix[:K])   # individual carrying capacity genotypes
	Rr  = convert(Array,popmatrix[:r])   # individual ldgr genotypes
	scl = 1 / pdf(Normal(0, bw), 0)      # a scaling factor for n'hood size

	# NOTE: For Gaussian distribution, ~95% of values are within 2σ of the mean.
	# When bw=0.25, an individual at spatial location 0.5 will have a neighborhood
	# where ~95% of its neighborhood is within the bounds 0:1, which roughly
	# translates to a discrete-patch system.

  # Initialize pre-allocated arrays
	w     = Array{Float64}(n)            # density wrt individual
	Nx    = Array{Float64}(n)            # total density each ind. experiences
	MeanD = Array{Float64}(n)            # mean dispersal genotype at each ind.'s n'hood
	MeanK = Array{Float64}(n)            # mean carrying capacity genotype
	Meanr = Array{Float64}(n)            # mean ldgr genotype

	tmpX  = Array{Float64}(length(x))    # temporary array for calculations

	# Calculations, looped over individuals
	# NOTE: densities (w) are calculated s.t. solitary individuals experience w=1
	for i = 1:n
		x_calc!(tmpX, x, x[i])
		w_calc!(w, tmpX, bw, scl)
		Nx_calc!(Nx, w, i)

		Mean_calc!(MeanD, RD, w, Nx[i], i)
	  Mean_calc!(MeanK, RK, w, Nx[i], i)
	  Mean_calc!(Meanr, Rr, w, Nx[i], i)

	  COV_calc!(V_D,  w, RD, RD, MeanD[i], MeanD[i], Nx[i], i, tmpX)
	  COV_calc!(V_K,  w, RK, RK, MeanK[i], MeanK[i], Nx[i], i, tmpX)
	  COV_calc!(V_r,  w, Rr, Rr, Meanr[i], Meanr[i], Nx[i], i, tmpX)
	  COV_calc!(C_DK, w, RD, RK, MeanD[i], MeanK[i], Nx[i], i, tmpX)
	  COV_calc!(C_Dr, w, RD, Rr, MeanD[i], Meanr[i], Nx[i], i, tmpX)
	  COV_calc!(C_Kr, w, RK, Rr, MeanK[i], Meanr[i], Nx[i], i, tmpX)
	end

	out = convert(DataFrame,[Nx MeanD MeanK Meanr V_D V_K V_r C_DK C_Dr C_Kr])
  names!(out, [:Nx :MeanD :MeanK :Meanr :V_D :V_K :V_r :C_DK :C_Dr :C_Kr])
	return out
end

# Calculate breeding value metrics wrt spatial location ------------------------
function sum_metrics(popmatrix, bw)
	# Given [popmatrix], returns [b, Nx, MeanD, MeanPD, and SDD].
	# popmatrix : see init_inds
	# bw        : neighborhood (bin) width, Gaussian σ
  # b         : spatial bin
	# Nx        : population density in continuous, Gaussian neighborhood
	# MeanD			: mean dispersal genotype in neighborhood
	# MeanPD    : mean dispersal phenotype in neighborhood
	# SDD       : σ(D), standard deviation in dispersal breeding values

  # If the population is extinct, return NULL
  if(size(popmatrix, 1) < 1)
  	return []
  end

  # Define within-function variables
  n   = size(popmatrix,1)               # number of individuals in population
  x   = convert(Array,popmatrix[:X])    # individual locations
	RD  = convert(Array,popmatrix[:D])    # individual dispersal genotypes
	RK  = convert(Array,popmatrix[:K])    # individual carrying capacity genotypes
	Rr  = convert(Array,popmatrix[:r])    # individual ldgr genotypes
	PD = convert(Array,popmatrix[:PD])    # individual D phenotypes
	PK = convert(Array,popmatrix[:PK])    # individual K phenotypes
	PR = convert(Array,popmatrix[:Pr])    # individual r phenotypes
	b   = accordian_bins(x,bw)            # spatial extent
	nb  = length(b)                       # number of bins
	scl = 1 / pdf(Normal(0, bw), 0)       # a scaling factor for n'hood size

  # Initialize pre-allocated arrays
	w      = Array{Float64}(length(x))    # density wrt bin
	Nx     = Array{Float64}(nb)           # total density each bin experiences
	MeanD  = Array{Float64}(nb)           # mean genotype at each bin
	MeanK  = Array{Float64}(nb)           # mean genotype at each bin
	Meanr  = Array{Float64}(nb)           # mean genotype at each bin
	MeanPD = Array{Float64}(nb)           # mean genotype at each bin
	MeanPK = Array{Float64}(nb)           # mean genotype at each bin
	MeanPr = Array{Float64}(nb)           # mean genotype at each bin
	V_D   = Array{Float64}(n)             # variance in dispersal breeding values
	V_K   = Array{Float64}(n)             # var in carrying capacity b.v.
	V_r   = Array{Float64}(n)             # var in ldgr b.v.
	C_DK  = Array{Float64}(n)             # cov in D and K breeding values
	C_Dr  = Array{Float64}(n)             # cov in D and r b.v.
	C_Kr  = Array{Float64}(n)             # cov in K and r b.v.
	tmpX   = Array(Float64,length(x))     # temporary array for calculations

	# Calculations, looped over bin
	for i = 1:nb
		x_calc!(tmpX, x, b[i])
		w_calc!(w, tmpX, bw, scl)
		Nx_calc!(Nx, w, i)

		Mean_calc!(MeanD, RD, w, Nx[i], i)
		Mean_calc!(MeanK, RK, w, Nx[i], i)
		Mean_calc!(Meanr, Rr, w, Nx[i], i)

		MeanP_calc!(MeanPD, PD, w, Nx[i], i)
		MeanP_calc!(MeanPK, PK, w, Nx[i], i)
		MeanP_calc!(MeanPr, Pr, w, Nx[i], i)

		COV_calc!(V_D,  w, RD, RD, MeanD[i], MeanD[i], Nx[i], i, tmpX)
		COV_calc!(V_K,  w, RK, RK, MeanK[i], MeanK[i], Nx[i], i, tmpX)
		COV_calc!(V_r,  w, Rr, Rr, Meanr[i], Meanr[i], Nx[i], i, tmpX)
		COV_calc!(C_DK, w, RD, RK, MeanD[i], MeanK[i], Nx[i], i, tmpX)
		COV_calc!(C_Dr, w, RD, Rr, MeanD[i], Meanr[i], Nx[i], i, tmpX)
		COV_calc!(C_Kr, w, RK, Rr, MeanK[i], Meanr[i], Nx[i], i, tmpX)
	end

	out = convert(DataFrame,[b Nx MeanD MeanK Meanr MeanPD MeanPK MeanPr V_D V_K V_r C_DK C_Dr C_Kr])
  names!(out, [:b :Nx :MeanD :MeanK :Meanr :MeanPD :MeanPK :MeanPr :V_D :V_K :V_r :C_DK :C_Dr :C_Kr])

  return out
end

# Calculate breeding value metrics wrt leading-edge (and zero) locations -------
function LE_sum_metrics(popmatrix, bw)
	# Given [popmatrix], returns [b, Nx, MeanD, MeanPD, and SDD].
	# popmatrix : see init_inds
	# bw        : neighborhood (bin) width, Gaussian σ
  # b         : spatial bin
	# Nx        : population density in continuous, Gaussian neighborhood
	# MeanD			: mean dispersal genotype in neighborhood
	# MeanPD    : mean dispersal phenotype in neighborhood
	# SDD       : σ(D), standard deviation in dispersal breeding values

  # If the population is extinct, return NULL
  if(size(popmatrix, 1) < 1)
  	return []
  end

  # Define within-function variables
  n   = size(popmatrix,1)               # number of individuals in population
  x   = convert(Array,popmatrix[:X])    # individual locations
	RD  = convert(Array,popmatrix[:D])    # individual dispersal genotypes
	RK  = convert(Array,popmatrix[:K])    # individual carrying capacity genotypes
	Rr  = convert(Array,popmatrix[:r])    # individual ldgr genotypes
  PD = convert(Array,popmatrix[:PD])    # individual D phenotypes
	PK = convert(Array,popmatrix[:PK])    # individual K phenotypes
	PR = convert(Array,popmatrix[:Pr])    # individual r phenotypes
	b   = [findmin(popmatrix[:X])[1],     # leftmost density > 0
				findmin(abs(popmatrix[:X]))[1], # density closest to X=0 (center)
				findmax(popmatrix[:X])[1]]      # rightmost density > 0
	nb  = 3                               # number of bins
	scl = 1 / pdf(Normal(0, bw), 0)       # a scaling factor for n'hood size

  # Initialize pre-allocated arrays
	w      = Array{Float64}(length(x))    # density wrt bin
	Nx     = Array{Float64}(nb)           # total density each bin experiences
	MeanD  = Array{Float64}(nb)           # mean genotype at each bin
	MeanK  = Array{Float64}(nb)           # mean genotype at each bin
	Meanr  = Array{Float64}(nb)           # mean genotype at each bin
	MeanPD = Array{Float64}(nb)           # mean genotype at each bin
	MeanPK = Array{Float64}(nb)           # mean genotype at each bin
	MeanPr = Array{Float64}(nb)           # mean genotype at each bin
	V_D   = Array{Float64}(n)             # variance in dispersal breeding values
	V_K   = Array{Float64}(n)             # var in carrying capacity b.v.
	V_r   = Array{Float64}(n)             # var in ldgr b.v.
	C_DK  = Array{Float64}(n)             # cov in D and K breeding values
	C_Dr  = Array{Float64}(n)             # cov in D and r b.v.
	C_Kr  = Array{Float64}(n)             # cov in K and r b.v.
	tmpX   = Array(Float64,length(x))     # temporary array for calculations

	# Calculations, looped over bin
	for i = 1:nb
		x_calc!(tmpX, x, b[i])
		w_calc!(w, tmpX, bw, scl)
		Nx_calc!(Nx, w, i)

		Mean_calc!(MeanD, RD, w, Nx[i], i)
		Mean_calc!(MeanK, RK, w, Nx[i], i)
		Mean_calc!(Meanr, Rr, w, Nx[i], i)

		MeanP_calc!(MeanPD, PD, w, Nx[i], i)
		MeanP_calc!(MeanPK, PK, w, Nx[i], i)
		MeanP_calc!(MeanPr, Pr, w, Nx[i], i)

		COV_calc!(V_D,  w, RD, RD, MeanD[i], MeanD[i], Nx[i], i, tmpX)
		COV_calc!(V_K,  w, RK, RK, MeanK[i], MeanK[i], Nx[i], i, tmpX)
		COV_calc!(V_r,  w, Rr, Rr, Meanr[i], Meanr[i], Nx[i], i, tmpX)
		COV_calc!(C_DK, w, RD, RK, MeanD[i], MeanK[i], Nx[i], i, tmpX)
		COV_calc!(C_Dr, w, RD, Rr, MeanD[i], Meanr[i], Nx[i], i, tmpX)
		COV_calc!(C_Kr, w, RK, Rr, MeanK[i], Meanr[i], Nx[i], i, tmpX)
	end

	out = convert(DataFrame,[b Nx MeanD MeanK Meanr MeanPD MeanPK MeanPr V_D V_K V_r C_DK C_Dr C_Kr])
	names!(out, [:b :Nx :MeanD :MeanK :Meanr :MeanPD :MeanPK :MeanPr :V_D :V_K :V_r :C_DK :C_Dr :C_Kr])
  return out
end

### LIFE HISTORY FUNCTIONS #####################################################

# Sexual mating ----------------------------------------------------------------
function mate(popmatrix, bw)
	# Given [popmatrix], returns [M].
	# popmatrix : see init_inds
	# bw        : neighborhood (bin) width, Gaussian σ
	# M         : a mate for each individual
	# NOTE: This mating assumes sampling with replacement

	# If the population is extinct, return NULL
  if size(popmatrix,1) < 1
  	return []
  end

  # Define within-function variables
  n   = size(popmatrix,1)               # number of individuals in pop
  x   = convert(Array,popmatrix[:X])    # individual locations
  scl = 1/pdf(Normal(0,bw),0)           # a scaling factor for n'hood size

	# Initialize pre-allocated arrays
	w      = Array{Float64,1}(n)          # weighted dens. to all other inds.
	D      = Array{Float64}(n)            # n'hood density (excluding yourself)
	M      = Array{Int}(n)                # mate matches
	tmpX   = Array{Float64}(length(x))    # temporary array for calculations
	tmpvec = Vector{Float64}(n)           # temp. vector to store cumulative sums

	# Draw a random number for each individual
  rnum = rand(Uniform(0,1), n)

	# Loop over all invividuals
	for i = 1:n
		x_calc!(tmpX, x, x[i])
		w_calc!(w, tmpX, bw, scl)
		D_calc!(D, w, i)

		# Matchmaker, matchmaker, make me a match
		if D[i] == 0   # if there's nobody near you...
			M[i] = -99   # no mate (results in mate-finding Allee effect)
		else
			scale!(w, 1 / D[i])                       # scale w by density for ind. i
			cumsum!(tmpvec, w, 1)                     # calculate cumsum of densities
			M[i] = searchsortedfirst(tmpvec, rnum[i]) # choose mate according to rnum
		end
	end
  return M
end

# Asexual "mating" -------------------------------------------------------------
function clone_mate(popmatrix)
	# Given [popmatrix], returns [M].
	# popmatrix : see init_inds
	# M         : a mate (themselves) for each individual
	# NOTE: This mating assumes each individual mates with itself

  # If the population is extinct, return NULL
  if size(popmatrix,1) < 1
  	return []
  end

  # Define within-function variables
  n = size(popmatrix,1)
  M = Array{Int}(n)

	# Matchmaker, matchmaker, make me my match
  for i = 1:n
  	M[i] = i
  end
  return M
end

# Reproduction and dispersal ---------------------------------------------------
function repro_disp(popmatrix, dk, H, V, bw)
	# Given [popmatrix(t), gp, dk, h2D, VPD], returns [popmatrix(t+1)]
	# popmatrix : see init_inds
	# gp        : see growth_params
	# dk        : dispersal kernel (see NSt)
	# h2D       : dispersal trait heritability
	# VPD       : phenotypic variance in dispersal
	# bw        : neighborhood (bin) width, Gaussian σ

	# repro_disp simulates:
	# (1) Density-dependent reproduction
	#	(2) Inheritence of dispersal genotype from parents to offspring
	#	(3) Offsrping dispersal

	# If the population is extinct, or if there is only 1 individual, return NULL
  if size(popmatrix,1) < 2
		return []
		#return DataFrame(X=[], D=[], H=[], PD=[], PH=[])
  end

	# (1) Density-dependent reproduction -----------------------------------------
  # a place to store the number of offspring produced by each parent
  Noff = Array{Int}(size(popmatrix,1))

  # Calculate non-additive (environmental) phenotypic variances
  ΣE = VE(H, V)

  # Calculate traits through space
  mets = metrics(popmatrix, bw)

  # Density-dependent population growth
  for i = 1:size(popmatrix,1)
		# with demographic stochasticity:
		Noff[i] = rand(Poisson(growth(mets[i,:Nx], popmatrix[i,:PK], popmatrix[i,:Pr])))
  end

  # If there are no offspring, return NULL
  if sum(Noff) == 0
		return []
		#return DataFrame(X=[], D=[], H=[], PD=[], PH=[])
  end

	# (2) Inheritence of dispersal phenotype from parents to offspring -----------

	# Find a mate for each individual
	# m = clone_mate(popmatrix) # asexual
  m = mate(popmatrix, bw)     # sexual

  # Create a new popmatrix (npm), accounting for the number of new offspring
  npm = DataFrame(Float64,sum(Noff),8)
  names!(npm,[:pX, :D, :PD, :K, :PK, :r, :Pr, :X])

	# initiate some variables
	ΣA  = Array{Float64}(3, 3)  # an additive genetic covariance matrix
  mpv = Vector{Float64}(3)    # a vector of midparent values
	c   = 0                     # a counter for indexing

	# loop over all of the parents
	for i = 1:size(popmatrix,1)
		if m[i] == -99
			continue   # if thre's no mate, there are no offspring
    else
 			# loop over all of the offspring from each parent:
			for j = 1:Noff[i]
				c += 1   # increment counter

				# offspring inherit location from their parents
        npm[c,:pX] = popmatrix[i,:X]

				# get midparent values (MPV) for each (mating) parent
        mpv[:] = ([popmatrix[  i, :D] popmatrix[  i, :K] popmatrix[  i, :r]]  +
				          [popmatrix[m[i],:D] popmatrix[m[i],:K] popmatrix[m[i],:r]]) *
									0.5

				# get the additive genetic covariance matrix
				ΣA[:, :] = [mets[i, :V_D]  mets[i, :C_DK] mets[i, :C_Dr]
				            mets[i, :C_DK] mets[i, :V_K]  mets[i, :C_Kr]
										mets[i, :C_Dr] mets[i, :C_Kr] mets[i, :V_r]]

				# calculate genotypes
				(npm[c,:D], npm[c,:K], npm[c,:R]) = rand(MvNormal(mpv, 0.5 * ΣA)), 1)

				# calculate phenotypes
				(npm[c,:PD], npm[c,:PK], npm[c,:PR]) = rand(Normal([0 0 0], ΣE), 1) +
				                                       [npm[c,:D], npm[c,:K], npm[c,:R]]

				# (3) Offspring dispersal ----------------------------------------------
				npm[c,:X] = rand(Normal(npm[c,:pX], exp(npm[c,:PD])),1)[]
      end
    end
  end

  # Delete NA rows that result from building npm without knowing which matings
  # were unsuccessful (i.e., m[i] == -99).
  deleterows!(npm,find(isna(npm[:,Symbol("X")])))
  return npm[:, [:X, :D, :PD, :K, :PK, :r, :Pr]]
end

### SIMULATION FUNCTIONS #######################################################

# Run the simulation -----------------------------------------------------------
function runsim(popmatrix,n,spX,gp,dk,ngens,μD,h2D,VPD,bw,p,r,TIME,outname="_")
	# Given [popmatrix, n, spX, gp, dk, ngens, μD, h2D, VPD, p, r], returns
	# changes in population size, extent, and genetic variance over time.
  # popmatrix : see init_inds
	# n  	      : initial number of starting individuals
	# spX       : initial spatial locations
	# gp        : see growth_params
	# dk        : dispersal kernel (see NSt)
	# ngens     : number of generations to simulate
	# μD	      : initial mean dispersal genotype of population
	# h2D       : dispersal heritability
	# VPD	      : initial total phenotypic variance in dispersal
	# p         : population number for this generation
	# r         : replicate number for this simulation
	# outname   : the output name for this simulation

	# start timers
	tic()
	simtime = time()
	flag = ""

	# set a unique seed based on h2D, p, and r, so simulations can be repeated.
	srand(Int(([1 10 100] * H)[] * 10000 + 42 * p + r))

	# initialize output
	batchoutput = Array{Any}(ngens * 9, 2314)

	# loop over generations
  for i = 1:ngens

		# do reproduction and dispersal
		popmatrix = repro_disp(popmatrix, gp, dk, h2D, VPD, bw)

		# if population goes extinct, make output blank for remaining generations
    if size(popmatrix,1) < 1
      for j = i:ngens
        batchoutput[(1:9)+(9)*(j-1),:] = blankoutput(p,r,gp,dk,h2D,j,bw)
      end
			break
		# if rerunwrap.jl has been running for more than 16 hours, get out
	elseif (time()-TIME) > 16*3600
			for j = i:ngens
				batchoutput[(1:9)+(9)*(j-1),:] = timeoutoutput(p,r,gp,dk,h2D,i,bw)
				flag = "!BATCHRUNTIMEOUT!"
			end
			break
		# if this simulation has been running for more than 8 hours, stop it
	elseif (time()-simtime) > 10*3600
			for j = i:ngens
				batchoutput[(1:9)+(9)*(j-1),:] = timeoutoutput(p,r,gp,dk,h2D,i,bw)
				flag = "!SIMTIMEOUT!"
			end
			break
		# otherwise, calculate relevant output metrics for this generation
		else
      batchoutput[(1:9)+(9)*(i-1),:] = calcoutput(popmatrix,p,r,gp,dk,h2D,i,bw)
    end
		# print("gen:",i, " pop. size:", size(popmatrix)[1],"\n") # track progress
  end

	# track simulation info and runtimes
	print(  "rep ",       (@sprintf "%2.0f" r),
				", pop ",       (@sprintf "%2.0f" p),
				", h2D ",       (@sprintf "%0.3f" h2D),
				", Nt ",        (@sprintf "%2.0f" gp.Nt),
				", bw ",        (@sprintf "%4.2f" bw),
				", kur ",       (@sprintf "%4.1f" dk.kur),
				", lam ",       (@sprintf "%5.2f" gp.lambda),
				", time ",      (clocktime(toq())),
				", final pop ", (@sprintf "%-6.0f" batchoutput[end,10]),
				" | total time ", (clocktime(time()-TIME)),
				" | ", flag, "\n")

	return batchoutput
end

function calcoutput(popmatrix, p, r, gp, dk, h2D, i, bw)
	# Given [popmatrix, p, r, gp, h2D, i], get summary stats for gen i
	# popmatrix : see init_inds
	# p         : population number for this simulation
	# r         : replicate number for this simulation
	# gp        : see growth_params
	# dk        : dispersal kernel for this simulation
	# h2D       : dispersal heritability
	# i         : current generation
	# bw        : neighborhood (bin) width, Gaussian σ

	popsize = nrow(popmatrix)                 # population size
  sm      = sum_metrics(popmatrix, bw)      # summary metrics
	le      = LE_sum_metrics(popmatrix,bw)

	# summary data sepecific to simulation, generation, etc.
	tmp = repmat([p,r,bw,gp.lambda,gp.Nt,gp.K,dk.kur,h2D,i,popsize]',5,1)

	# reshaped summary statistics from sum_metrics and LE_sum_metrics
	metrix = vcat(sm[:,:b]',sm[:,:Nx]',sm[:,:MeanD]',sm[:,:MeanPD]',sm[:,:SDD]')
	le_met = vcat(le[:,:b]',le[:,:Nx]',le[:,:MeanD]',le[:,:MeanPD]',le[:,:SDD]')

	return hcat(tmp,le_met,metrix)
end

function blankoutput(p, r, gp, dk, h2D, i, bw)
	# Given [p, r, gp, dk, h2D, i], get '0' summary stats for for generation i
	# p         : population number for this simulation
	# r         : replicate number for this simulation
	# gp        : see growth_params
	# dk        : see NSt
	# h2D       : dispersal heritability
	# i         : current generation
	# bw        : neighborhood (bin) width, Gaussian σ

	popsize = 0                   # population size
  mL      = 0                   # leftmost density > 0
  mR      = 0                   # rightmost density > 0
  m0      = 0                   # density closest to X=0 (center)
  sm      = repmat([0],2301,5)  # summary metrics

	# summary data specific to simulation, generation, etc.
  tmp = repmat([p,r,bw,gp.lambda,gp.Nt,gp.K,dk.kur,h2D,i,popsize,mL,m0,mR]',5,1)

	# reshaped summary statistics
	metrix = vcat(sm')

  return hcat(tmp,metrix)
end

function timeoutoutput(p, r, gp, dk, h2D, i, bw)
	# Given [p, r, gp, dk, h2D, i], get '0' summary stats for for generation i
	# p         : population number for this simulation
	# r         : replicate number for this simulation
	# gp        : see growth_params
	# dk        : see NSt
	# h2D       : dispersal heritability
	# i         : current generation
	# bw        : neighborhood (bin) width, Gaussian σ

	popsize = -i                  # population size
  mL      = -999                # leftmost density > 0
  mR      = -999                # rightmost density > 0
  m0      = -999                # density closest to X=0 (center)
  sm      = repmat([0],2301,5)  # summary metrics

	# summary data specific to simulation, generation, etc.
  tmp = repmat([p,r,bw,gp.lambda,gp.Nt,gp.K,dk.kur,h2D,i,popsize,mL,m0,mR]',5,1)

	# reshaped summary statistics
	metrix = vcat(sm')

  return hcat(tmp,metrix)
end
