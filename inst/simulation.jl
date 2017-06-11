# ------------------------------------------------------------------------------
# --- DESCRIPTION OF MODEL -----------------------------------------------------
# ------------------------------------------------------------------------------
# Continuous space, discrete time model for examining stochastic evolutionary
# processes during range expansion.
#
# Original R Code by Ben Phillips (2015)
# Translated to Julia by Brad Ochocki (2016)

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
function outname(h,p,r)
	return string(h,"_",Int(p),"_",Int(r),"_R.csv")
end

# function for turning seconds into a HH:MM:SS:mmm format
function clocktime(elapsed)
	hh   = trunc(elapsed/3600)                 # hours
	mm   = trunc((elapsed - hh*3600)/60)       # minutes
	ss   = trunc((elapsed - hh*3600 - mm*60))  # seconds
	mill = (elapsed - trunc(elapsed))*1000     # milliseconds
	out  = string((@sprintf "%02.0f" hh),  ":",
	              (@sprintf "%02.0f" mm),  ":",
							  (@sprintf "%02.0f" ss),  ":",
							  (@sprintf "%03.0f" mill));
	return out
end

### DISPERSAL KERNEL FUNCTIONS #################################################

# Define non-standardized t-distribution ---------------------------------------
# We need to define a few properties of the distribution so that it plays well
# with the Distributions.jl package.

# Define a new univariate Distribution type.
immutable NSt <: ContinuousUnivariateDistribution
	mu::Float64   # mean of the distribution
	sd::Float64   # standard deviation of the distribution
	kur::Float64  # *EXCESS* kurtosis
end

# Minimum and maximum support values
Distributions.minimum(d::NSt)  = -Inf
Distributions.maximum(d::NSt)  = Inf

# Statistics
Distributions.mean(d::NSt)     = d.mu
Distributions.median(d::NSt)   = d.mu
Distributions.mode(d::NSt)     = d.mu
Distributions.var(d::NSt)      = d.sd^2
Distributions.kurtosis(d::NSt) = d.kur + 3

# Evaluation and sampling
# pdf function
function Distributions.pdf(d::NSt, x::Real)
	# note: formulated in terms of excess kurtosis
	nu  = 6/d.kur + 4                                             # d.o.f.
	s   = sqrt((nu-2)*(d.sd^2)/nu)                                # scale
	num = gamma((nu+1)/2) *(1+(1/nu)*((x-d.mu)/s)^2)^(-(nu+1)/2)  # numerator
	den = (gamma(nu/2)*sqrt(pi*nu)*s)                             # denominator
	num/den                                                       # density(s)
end

# quantile function
function Base.quantile(d::NSt, q::Real)
	# note: formulated in terms of excess kurtosis
	nu  = 6/d.kur + 4                # d.o.f.
	s   = sqrt((nu-2)*(d.sd^2)/nu)   # scale
	d.mu + s.*quantile(TDist(nu),q)  # quantile(s)
end

# sampling function
function Distributions.rand(d::NSt, n::Int64)
	# note: formulated in terms of excess kurtosis
	nu  = 6/d.kur + 4                # d.o.f.
	s   = sqrt((nu-2)*(d.sd^2)/nu)   # scale
	r   = rand(Uniform(0,1),n)       # uniform random number(s)
	d.mu + s.*quantile(TDist(nu),r)  # random sample(s)
end

# Test the distribution
#test = NSt(10,3,0)
#rand(test,10)
#quantile(test,[0.42,0.50])
#maximum(rand(test,10000))
#minimum(rand(test,10000))
#mean(rand(test,10000))
#var(rand(test,10000))
#pdf(test,[-1,0,1])

# Define scaled-shifted beta distribution --------------------------------------
# We need to define a few properties of the distribution so that it plays well
# with the Distributions.jl package.

# Define a new univariate Distribution type.
immutable ssB <: ContinuousUnivariateDistribution
	mu::Float64
	sd::Float64   # standard deviation of the distribution
	kur::Float64  # *EXCESS* kurtosis
end

# Statistics
Distributions.mean(d::ssB)     = mu
Distributions.median(d::ssB)   = mu
Distributions.mode(d::ssB)     = mu
Distributions.var(d::ssB)      = d.sd^2
Distributions.kurtosis(d::ssB) = d.kur + 3

# Evaluation and sampling
# quantile function
function Base.quantile(d::ssB, q::Real)
	# note: formulated in terms of excess kurtosis
	alpha = -3/d.kur - 3/2                         # shape parameter
	scale = 2 * d.sd * sqrt(2 * alpha + 1)         # scale
	shift = -0.5                                   # zero-centered
  d.mu + shift*scale  + scale.*quantile(Beta(alpha,alpha),q)  # quantile(s)
end

# sampling function
function Distributions.rand(d::ssB, n::Int64)
	# note: formulated in terms of excess kurtosis
	alpha = -3/d.kur - 3/2                         # shape parameter
	scale = 2 * d.sd * sqrt(2 * alpha + 1)         # scale
	shift = -0.5                                    # zero-centered
	r     = rand(Uniform(0,1),n)                   # uniform random number(s)
	d.mu + shift*scale + scale.*quantile(Beta(alpha,alpha),r)   # random sample(s)
end

# Test the distribution
#test = ssB(5.196,3,-1.2)
#rand(test,10)
#quantile(test,[0.42,0.50])
#maximum(rand(test,1000000))
#minimum(rand(test,1000000))
#mean(rand(test,100000))
#var(rand(test,100000))

### GROWTH FUNCTIONS ###########################################################

# Piecewise population growth --------------------------------------------------
# piecewise growth parameters
immutable growth_paramsPW
	N_A::Float64    # critical Allee threshold
	N_max::Float64  # density for maximum growth rate (L_max)
	K::Float64      # carrying capacity
	L_max::Float64  # maximum growth rate (@ N=N_max)
	L_low::Float64  # low density (@ N=0) growth rate
end

function pwLambda(N, gp)
  # Given [N, gp] returns [lambda]
	# N         : number of individuals
	# gp        : see growth_params
	# gp.N_A    : critical Allee threshold
	# gp.N_max  : density for maximum growth rate (L_max)
	# gp.K      : carrying capacity
	# gp.L_max  : maximum growth rate (@ N=N_max)
	# gp.L_low  : low density (@ N=0) growth rate
  # gp.lambda : per-capita growth rate

	if  N < gp.N_A
		return 0.0
	elseif gp.N_A <= N < gp.N_max
		return (gp.L_max - gp.L_low)/(gp.N_max - gp.N_A)*(N-gp.N_A) + gp.L_low
	elseif gp.N_max <= N < gp.K
		return (1-gp.L_max)/(gp.K-gp.N_max)*(N-gp.N_max) + gp.L_max
	elseif gp.K <= N
		return gp.K/N
	end
end

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

### GENETIC VARIANCE FUNCTIONS #################################################

# Initialize population --------------------------------------------------------
function init_inds(n, spX, μD, h2D, VPD)
	# Given [n, spX, μD, h2D, VPD] returns [X, D, PD], which is referred to in
	# future functions as 'popmatrix'.
  # n   : number of starting individuals
	# spX : absolute value of starting population's spatial boundaries
	# μD  : mean dispersal genotype
	# h2D : dispersal heritability
	# VPD : phenotypic variance in dispersal
	# X   : spatial location of each individual
	# D   : dispersal genotype
	# PD  : dispersal phenotype (accounts for environmental variance)

	X  = rand(Uniform(-spX,spX), n)
  D  = rand(Normal(μD,sqrt(h2D*VPD)), n)
  PD = rand(Normal(0,sqrt((1-h2D)*VPD)), n) + D
	return(DataFrame(X=X, D=D, PD=PD))
end

# Calculate non-additive (environmental) phenotypic variance -------------------
function VE(h2,VP)
	# Given [h2 and VP], returns [VE]
	# h2 : heritability of a trait
	# VP : phenotypic variance of a trait
	# VE : environmental variance of a trait
	return (1-h2)*VP
end

### FUNCTIONS FOR NEIGHBORHOOD-RELATED CALCULATIONS ############################

# Determine which spatial locations to sample from -----------------------------
function accordian_bins(A::Array{Float64,1},bw::Float64)
	x = vcat(minimum(A),maximum(A))

	scl = bw*4
	# make a scale so that accordian_bins changes appropriately with bw; when
	# bw = 0.25, "fixed" bins should be spaced 1 unit apart.

	# center
	C = collect(linspace(-50.0,50.0,101)*scl)

	# left side
	if x[1] < -1150*scl
	  L  = collect(linspace(  x[1], x[1]+99*scl,  100))
	  LC = collect(linspace(L[100],        C[1], 1002))[2:1001]
	else
	  L  = collect(linspace(-1150*scl, -1051*scl,  100))
	  LC = collect(linspace(-1050*scl,   -51*scl, 1000))
	end

	# right side
	if x[2] > 1150*scl
	  R  = collect(linspace(x[2]-99*scl, x[2],  100))
	  CR = collect(linspace(     C[101], R[1], 1002))[2:1001]
	else
	  R  = collect(linspace(1051*scl, 1150*scl,  100))
	  CR = collect(linspace(  51*scl, 1050*scl, 1000))
	end

	return vcat(L,LC,C,CR,R)
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
	C1 = -1/(2*sd*sd)          # A constant for the pdf's exponent
	C2 = scl/sqrt(2*sd*sd*pi)  # A constant for the first term in the pdf
	@fastmath @inbounds @simd for j = 1:length(w)
		w[j] = A[j]*A[j] * C1    # calculate (A^2) * C1
		if w[j] < -15.5            # ignore cases when density will be < than exp(-20)
			w[j] = 0.0
		else
			w[j] = exp(w[j]) * C2  # calculate C2 * exp((A^2) * C1)
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

	D[i] = sum(w)-w[i]
	nothing
end

# Calculate the mean genotype for each location --------------------------------
function MeanD_calc!(MeanD::Array{Float64,1}, RD::Array{Float64,1},
	                   w::Array{Float64,1}, Nx::Float64, i::Int64)
	# MeanD : the average genotype in a location given by index i
	# RD    : the value of genotypes for each index
	# w     : the weighted density of all individuals in the population wrt ind. i
	# Nx    : the density at the location of ind. i
	# i     : the index of the invididual

	MeanD[i]  = (*(RD',w)/Nx)[1] # sum of weighted RD, scaled by patch density
	nothing
end

# Calculate the mean phenotype for each location -------------------------------
function MeanPD_calc!(MeanPD::Array{Float64,1}, RPD::Array{Float64,1},
	                    w::Array{Float64,1}, Nx::Float64, i::Int64)
	# MeanPD : the average phenotype in a location given by index i
	# RPD    : the value of phenotypes for each index
	# w      : the weighted dens. of all individuals in the population wrt ind. i
	# Nx     : the density at the location of ind. i
	# i      : the index of the invididual

	MeanPD[i]  = (*(RPD',w)/Nx)[1] # sum of weighted RPD, scaled by patch density
	nothing
end

# Calculate the standard deviation in breeding values for each location --------
function SDD_calc!(SDD::Array{Float64,1}, w::Array{Float64,1},
	                 RD::Array{Float64,1}, MeanD::Float64, Nx::Float64, i::Int64,
									 tmp::Array{Float64,1})
	# SDD   : standard deviation in breeding values at location i
	# w     : the weighted density of all individuals in the population wrt ind. i
	# RD    : the genotype at index i
	# MeanD : the average genotype in a location given by index i
	# Nx    : the density at the location of ind. i
	# i     : the index of the invididual
	# tmp   : vector to store the output

	# The original code, which was simpler to write but much slower:
	# SDD[i]   = sqrt(sum(w.*(MeanD - RD).^2)/Nx)

	@fastmath @inbounds @simd for j = 1:length(tmp)
		tmp[j] = MeanD - RD[j]             # subtract each element of RD from MeanD
		tmp[j] = tmp[j] * tmp[j] * w[j]    # tmp^2 * w (elementwise)
	end
	SDD[i] = sqrt(sum(tmp)/Nx)           # sqrt(sum(tmp)/(total density))
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
  if size(popmatrix,1) < 1
  	return []
  end

  # Define within-function variables
  n   = size(popmatrix,1)              # number of individuals in population
  x   = convert(Array,popmatrix[:X])   # individual locations
  RD  = convert(Array,popmatrix[:D])   # individual dispersal genotypes
  scl = 1 / pdf(Normal(0, bw), 0)      # a scaling factor for n'hood size

	# NOTE: For Gaussian distribution, ~95% of values are within 2σ of the mean.
	# When bw=0.25, an individual at spatial location 0.5 will have a neighborhood
	# where ~95% of its neighborhood is within the bounds 0:1, which roughly
	# translates to a discrete-patch system.

  # Initialize pre-allocated arrays
	w     = Array{Float64}(n)            # density wrt individual
	Nx    = Array{Float64}(n)            # total density each ind. experiences
	MeanD = Array{Float64}(n)            # mean genotype at each ind.'s n'hood
	SDD   = Array{Float64}(n)            # sd in dispersal breeding values
	tmpX  = Array{Float64}(length(x))    # temporary array for calculations

	# Calculations, looped over individuals
	# NOTE: densities (w) are calculated s.t. solitary individuals experience w=1
	for i = 1:n
		x_calc!(tmpX, x, x[i])
		w_calc!(w, tmpX, bw, scl)
		Nx_calc!(Nx, w, i)
		MeanD_calc!(MeanD, RD, w, Nx[i], i)
		SDD_calc!(SDD, w, RD, MeanD[i], Nx[i], i, tmpX)
	end

	out = convert(DataFrame,[Nx MeanD SDD])  # output as dataframe
  names!(out,[:Nx,:MeanD,:SDD])            # update names
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
  if(size(popmatrix,1)<1)
  	return []
  end

  # Define within-function variables
  n   = size(popmatrix,1)               # number of individuals in population
  x   = convert(Array,popmatrix[:X])    # individual locations
  RD  = convert(Array,popmatrix[:D])    # individual dispersal genotypes
  RPD = convert(Array,popmatrix[:PD])   # individual dispersal phenotypes
	b   = accordian_bins(x,bw)            # spatial extent
	nb  = length(b)                       # number of bins
	scl = 1 / pdf(Normal(0, bw), 0)       # a scaling factor for n'hood size

  # Initialize pre-allocated arrays
	w      = Array{Float64}(length(x))    # density wrt bin
	Nx     = Array{Float64}(nb)           # total density each bin experiences
	MeanD  = Array{Float64}(nb)           # mean genotype at each bin
	MeanPD = Array{Float64}(nb)           # mean genotype at each bin
	SDD    = Array{Float64}(nb)           # sd in dispersal breeding values
	tmpX   = Array(Float64,length(x))     # temporary array for calculations

	# Calculations, looped over bin
	for i = 1:nb
		x_calc!(tmpX, x, b[i])
		w_calc!(w, tmpX, bw, scl)
		Nx_calc!(Nx, w, i)
		MeanD_calc!(MeanD, RD, w, Nx[i], i)
		MeanPD_calc!(MeanPD, RPD, w, Nx[i], i)
		SDD_calc!(SDD, w, RD, MeanD[i], Nx[i], i, tmpX)
	end

	out = convert(DataFrame,[b Nx MeanD MeanPD SDD]) # output as dataframe
	names!(out,[:b, :Nx, :MeanD, :MeanPD, :SDD])     # update names
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
  if(size(popmatrix,1)<1)
  	return []
  end

  # Define within-function variables
  n   = size(popmatrix,1)               # number of individuals in population
  x   = convert(Array,popmatrix[:X])    # individual locations
  RD  = convert(Array,popmatrix[:D])    # individual dispersal genotypes
  RPD = convert(Array,popmatrix[:PD])   # individual dispersal phenotypes
	b   = [findmin(popmatrix[:X])[1],     # leftmost density > 0
				findmin(abs(popmatrix[:X]))[1], # density closest to X=0 (center)
				findmax(popmatrix[:X])[1]]      # rightmost density > 0
	nb  = 3                               # number of bins
	scl = 1 / pdf(Normal(0, bw), 0)       # a scaling factor for n'hood size

  # Initialize pre-allocated arrays
	w      = Array{Float64}(length(x))    # density wrt bin
	Nx     = Array{Float64}(nb)           # total density each bin experiences
	MeanD  = Array{Float64}(nb)           # mean genotype at each bin
	MeanPD = Array{Float64}(nb)           # mean genotype at each bin
	SDD    = Array{Float64}(nb)           # sd in dispersal breeding values
	tmpX   = Array(Float64,length(x))     # temporary array for calculations

	# Calculations, looped over bin
	for i = 1:nb
		x_calc!(tmpX, x, b[i])
		w_calc!(w, tmpX, bw, scl)
		Nx_calc!(Nx, w, i)
		MeanD_calc!(MeanD, RD, w, Nx[i], i)
		MeanPD_calc!(MeanPD, RPD, w, Nx[i], i)
		SDD_calc!(SDD, w, RD, MeanD[i], Nx[i], i, tmpX)
	end

	out = convert(DataFrame,[b Nx MeanD MeanPD SDD]) # output as dataframe
	names!(out,[:b, :Nx, :MeanD, :MeanPD, :SDD])     # update names
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
function repro_disp(popmatrix, gp, dk, h2D, VPD, bw)
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
  	return DataFrame(X=[], D=[], H=[], PD=[], PH=[])
  end

	# (1) Density-dependent reproduction -----------------------------------------
  # a place to store the number of offspring produced by each parent
  Noff = Array{Int}(size(popmatrix,1))

  # Calculate non-additive (environmental) phenotypic variances
  sdVED = sqrt(VE(h2D,VPD))

  # Calculate traits through space
  mets = metrics(popmatrix, bw)

  # Density-dependent population growth
  for i = 1:size(popmatrix,1)
		# with demographic stochasticity:
  	# Noff[i] = rand(Poisson(pwLambda(mets[i,:Nx],gp)),1)
		# Noff[i] = ceil(Allee_oop(mets[i,:Nx],0,10))
		Noff[i] = rand(Poisson(hockey_stick(mets[i,:Nx],gp)))

		# without demographic stochasticity (likely messes up growth function):
		#Noff[i] = round((pwLambda(mets[i,:Nx],gp)))
  end

  # If there are no offspring, return NULL
  if sum(Noff) == 0
  	return DataFrame(X=[], D=[], H=[], PD=[], PH=[])
  end

	# (2) Inheritence of dispersal phenotype from parents to offspring -----------
	# get standard deviation in dispersal breeding values
	SDD = mets[:SDD]

	# Find a mate for each individual
	# m = clone_mate(popmatrix) # asexual
  m = mate(popmatrix, bw)     # sexual

  # Create a new popmatrix (npm), accounting for the number of new offspring
  npm = DataFrame(Float64,sum(Noff),5)
  names!(npm,[:mD,:pX,:D,:PD,:X])

	# start a counter for indexing
	c = 0

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
        npm[c,:mD] = 0.5 * (popmatrix[i,:D] + popmatrix[m[i],:D])

				# calculate dispersal genotype
				if SDD[i] == 0
          npm[c,:D] = npm[c,:mD][] # if no genetic variance, genotype = MPV
        else
          npm[c,:D] = rand(Normal(npm[c,:mD], sqrt(0.5 * SDD[i]^2)), 1)[]
        end

				# calculate dispersal phenotype
				npm[c,:PD] = npm[c,:D] + rand(Normal(0,sdVED),1)[]

				# (3) Offspring dispersal ----------------------------------------------
				if dk.kur == 0
					# if excess kurtosis is 0, use Gaussian dispersal
					npm[c,:X] = rand(Normal(npm[c,:pX], exp(npm[c,:PD])),1)[]
				elseif dk.kur > 0
					# if excess kurtosis is >0, use NSt dispersal
					npm[c,:X] = rand(NSt(npm[c,:pX], exp(npm[c,:PD]), dk.kur),1)[]
				elseif dk.kur < 0
					# if excess kurtosis is <0, use ssB dispersal
					npm[c,:X] = rand(ssB(npm[c,:pX], exp(npm[c,:PD]), dk.kur),1)[]
				end
      end
    end
  end

  # Delete NA rows that result from building npm without knowing which matings
  # were unsuccessful (i.e., m[i] == -99).
  deleterows!(npm,find(isna(npm[:,Symbol("X")])))
  return npm[:,[:X,:D,:PD]]
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
	srand(Int(h2D*10000 + 42*p + r))

	# initialize output
	batchoutput = Array{Any}(ngens*5, 2314)

	# loop over generations
  for i = 1:ngens

		# do reproduction and dispersal
		popmatrix = repro_disp(popmatrix, gp, dk, h2D, VPD, bw)

		# if population goes extinct, make output blank for remaining generations
    if size(popmatrix,1) < 1
      for j = i:ngens
        batchoutput[(1:5)+(5)*(j-1),:] = blankoutput(p,r,gp,dk,h2D,j,bw)
      end
			break
		# if rerunwrap.jl has been running for more than 16 hours, get out
	elseif (time()-TIME) > 16*3600
			for j = i:ngens
				batchoutput[(1:5)+(5)*(j-1),:] = timeoutoutput(p,r,gp,dk,h2D,i,bw)
				flag = "!BATCHRUNTIMEOUT!"
			end
			break
		# if this simulation has been running for more than 8 hours, stop it
	elseif (time()-simtime) > 10*3600
			for j = i:ngens
				batchoutput[(1:5)+(5)*(j-1),:] = timeoutoutput(p,r,gp,dk,h2D,i,bw)
				flag = "!SIMTIMEOUT!"
			end
			break
		# otherwise, calculate relevant output metrics for this generation
		else
      batchoutput[(1:5)+(5)*(i-1),:] = calcoutput(popmatrix,p,r,gp,dk,h2D,i,bw)
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
