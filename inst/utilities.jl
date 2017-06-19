

# Turn seconds into a HH:MM:SS:mmm format --------------------------------------
function clocktime(elapsed)
	hh   = trunc(elapsed / 3600.0)                  # hours
	mm   = trunc((elapsed - hh * 3600.0) / 60.0)      # minutes
	ss   = trunc((elapsed - hh * 3600.0 - mm * 60.0)) # seconds
	mill = (elapsed - trunc(elapsed)) * 1000.0      # milliseconds
	out  = string((@sprintf "%02.0f" hh),  ":",
	              (@sprintf "%02.0f" mm),  ":",
							  (@sprintf "%02.0f" ss),  ":",
							  (@sprintf "%03.0f" mill));
	return out
end

# Replicate each element in an array a number of times (in place) --------------
function rep_each!(A::Vector{Int64}, B::Vector{Int64}, C::Vector{Int64})
	# Given [A, B, C], returns A
	# A A vector where each element B_i is repeated C_i times.
	# B vector whose elements should be repeated
	# C the number of times each element in A should be repeated

	c = 1
	for i = 1:length(B)
		A[c:(c + C[i] - 1)] = B[i]
		c += C[i]
	end
	nothing
end

# Replicate each element in an array a number of times -------------------------
function rep_each(B::Vector{Float64}, C::Vector{Int64})
	# Given [B, C], returns A
	# A A vector where each element B_i is repeated C_i times.
	# B vector whose elements should be repeated
	# C the number of times each element in A should be repeated

	A = Vector{Float64}(sum(C))

	c = 1
	for i = 1:length(B)
		A[c:(c + C[i] - 1)] = B[i]
		c += C[i]
	end
	return A
end

# Replicate each element in an integer array a number of times -----------------
function rep_eachi(B::Vector{Int64}, C::Vector{Int64})
	# Given [B, C], returns A
	# A A vector where each element B_i is repeated C_i times.
	# B vector whose elements should be repeated
	# C the number of times each element in A should be repeated

	A = Vector{Int64}(sum(C))

	c = 1
	for i = 1:length(B)
		A[c:(c + C[i] - 1)] = B[i]
		c += C[i]
	end
	return A
end

# Calculate sum of the product of the deviates from mean trait values ----------
function sum_prod_deviates!(spd::Vector{Float64},
                            μ1::Vector{Float64}, μ2::Vector{Float64},
                            x1::Vector{Float64}, x2::Vector{Float64},
                            ix::Vector{Int64})
	# Given [spd, μ1, μ2, x1, x2, ix], returns spd.
	# spd : vector of the sum of the product of the deviates for each patch
	# μ1  : a vector of mean trait 1 values for each patch
	# μ2  : a vector of mean trait 2 values for each patch
	# x1  : a vector of each individual's trait 1 values
	# x2  : a vector of each individual's trait 1 values
	# ix  : a vector to map each individual's location to an index in μ1 & μ2

  for i = 1:length(ix)
    if ix[i] > 0
      d1 = x1[i] - μ1[ix[i]]
      d2 = x2[i] - μ2[ix[i]]
      spd[ix[i]] += d1 * d2
    end
  end
  nothing
end

# Sum each occurance of A_i in vector B ----------------------------------------
function running_sum!(A::Vector{Float64}, B::Vector{Float64},
                      ix::Vector{Int64})
	# given [A, B, ix], returns A.
	# A  : a vector that contains the sum of elements of B as they relate to ix
	# B  : a vector of values
	# ix : a vector to map each element of B to an index in A
  for i = 1:length(ix)
    if ix[i] > 0
      A[ix[i]] += B[i]
    end
  end
  nothing
end

# Count number of occurances ---------------------------------------------------
function running_count!(A::Vector{Int64}, ix::Vector{Int64})
	# given [A, ix], returns A.
	# A  : a vector to count the number of occurances of each element of ix
	# ix : a vector integers
  for i = 1:length(ix)
    if ix[i] > 0
      A[ix[i]] += 1
    end
  end
  nothing
end

# Sort rows in place (for an array of integers) --------------------------------
function sortrowsi!(A::Array{Int64, 2})
	# A : an integer array to be sorted.
	r = size(A)[1]
	c = size(A)[2]
	A[1:r, 1:c] = sortrows(A)
end

# Map each individual's patch id to an appropriate array index -----------------
function map_dict!(ix::Vector{Int64}, x::Vector{Int64}, dict::Dict{Int64,Int64})
  for i = 1:length(ix)
    #ix[i] = try dict[x[i]] catch; 0 end
    ix[i] = get(dict, x[i], 0)
  end
  nothing
end

# Map values of a Float64 vector to another Float64 vector ---------------------
function map_vals!(A::Vector{Float64}, B::Vector{Float64}, ix::Vector{Int64})
	# Given [A, B, ix], returns [A]
	# A  : a vector of values mapped from B according to ix
	# B  : a vector of values
	# ix : a vector of indices in B
  for i = 1:length(ix)
    A[i] = B[ix[i]]
  end
  nothing
end

# Map values of a Float64 vector to another Float64 vector ---------------------
function map_vals(A::Vector{Float64}, ix::Vector{Int64})
	# Given [A, ix], returns [B]
	# A  : a vector of values
	# ix : a vector of indices in A

	# B  : a vector of values mapped from A according to ix

	B = Vector{Float64}(length(ix))
  for i = 1:length(ix)
    B[i] = A[ix[i]]
  end
  return B
end

# Map values of an Int64 vector to another Int64 vector ------------------------
function map_valsi!(A::Vector{Int64}, B::Vector{Int64}, ix::Vector{Int64})
	# Given [A, B, ix], returns [A]
	# A  : a vector of values mapped from B according to ix
	# B  : a vector of values
	# ix : a vector of indices in B
	for i = 1:length(ix)
    A[i] = B[ix[i]]
  end
  nothing
end

# Map values of a Int64 vector to another Float64 vector ---------------------
function map_valsi(A::Vector{Int64}, ix::Vector{Int64})
	# Given [A, ix], returns [B]
	# A  : a vector of values
	# ix : a vector of indices in A

	# B  : a vector of values mapped from A according to ix

	B = Vector{Int64}(length(ix))
  for i = 1:length(ix)
    B[i] = A[ix[i]]
  end
  return B
end

# Find indices of individuals in the same patch --------------------------------
function group_index!(pai::Array{Int64, 2}, sx::Array{Int64, 2},
	                    pa::Vector{Int64}, np::Int64, ni::Int64)
	# pai : an array to track starting and ending indices for individuals that
  #       share a patch. each row corresponds to a row in pa.
	# sx  : and array where the first column is a sorted vector of x and the
	#       second column provides the entry's original index in x.
	# pa  : a sorted vector of unique occupied patches
	# np  : the number of unique occupied patches
	# ni  : the number of individuals
	c = 1
	for i = 1:np
		pai[i, 1] = c
		while pa[i] == sx[c, 1]
			c += 1
			if c == ni + 1
				break
			end
		end
		pai[i, 2] = c - 1
	end
	nothing
end

# Check if a function is positive-definite. If not, tweak it a bit -------------
function pos_def_check!(A::Array{Float64, 2})
	# given [A] returns [A]
	# A : a matrix that should be positive definite.
	#     If A is already positive definite, A is returned unchanged.
	#     If A is not positive definite, the strength of the correlation between
	#     traits is reduced by 0.2% to make A positive definite.

	if !isposdef(A)
		# get standard deviations
		σ = sqrt(diag(A))

		# get correlation
		C = cov2cor(A, σ)

		C[2] *= 0.998
		C[3] *= 0.998

		# slightly weaken correlation, make new covariance matrix
		A[1:2, 1:2] = cor2cov(C, σ)

	end
	nothing
end

# Add two vectors element-wise (without using broadcast) -----------------------
function vector_add(A::Vector{Float64}, B::Vector{Float64})
	# Given [A, B], returns [C]
	# A : a vector of the same length as B
	# B : a vector of the same length as A

	C = Vector{Float64}(length(A))
	for i = 1:length(A)
		C[i] = A[i] + B[i]
	end
	return C
end

# Multiply two vectors element-wise (without using broadcast) ------------------
function vector_mult(A::Vector{Float64}, B::Vector{Float64})
	# Given [A, B], returns [C]
	# A : a vector of the same length as B
	# B : a vector of the same length as A

	C = Vector{Float64}(length(A))
	for i = 1:length(A)
		C[i] = A[i] * B[i]
	end
	return C
end

# Multiply a vector by a scalar, element-wise (without using broadcast) --------
function vector_scalar_mult!(A::Vector{Float64}, B::Float64)
	# Given [A, B], returns [A]
	# A : a vector of the same length as B
	# B : a scalar

	for i = 1:length(A)
		A[i] *= B
	end
	nothing
end
