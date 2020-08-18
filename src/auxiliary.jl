# Helper functions for PCHIP interpolation

"""
    checkinput(x::Vector{<:Real}, y::AbstractArray{T,N} where T<:Real where N)
      -> x, y, (x = xsize, y = ysize), T

Check and/or correct the input data: `x` must be a vector with monotonic ascending data,
`y` a vector or array with the first dimension of `length(x)`; `y` data of higher dimensions
than `2` are reshaped to a `Matrix`.

Return the corrected `x` and `y`-data, a tuple of the dimensions of `x` and `y`,
and the inferred type `T` of the data, which must be of type float.
"""
function checkinput(x::Vector{<:Real}, y::AbstractArray{T,N} where T<:Real where N)
  # Infer output type (cannot be integer)
  T = promote_type(eltype(x), eltype(y))
  T<:Integer && (T = float(T))
  x, y = T.(x), T.(y)

  # Get array sizes
  xsize, ysize = length(x), size(y)
  # Check for minimum size and correct corresponding lengths of x and y data
  if xsize < 2
    throw(DataError("PCHIP requires at least 2 points for interpolation", x))
  elseif xsize ≠ ysize[1]
    ydim = findfirst(ysize .== xsize)
    if ydim === nothing
      throw(DataError(string("`x` and first dimension of `y` data must be of same length; ",
        "got $xsize, $(size(y, 1))")))
    else
      ydims = [ydim; setdiff(1:length(ysize), ydim)...]
      @warn("first dimension must be of length `x`; dimension reordered", ydims)
      y = permutedims(y, ydims)
      ysize = size(y)
    end
  end
  # Ensure maximum of 2 dimensions for interpolation
  length(ysize) ≥ 2 && (y = reshape(y, ysize[1], prod(ysize[2:end])))
  # Check for monotonic data and make sure data is ascending
  if x ≠ sort(x)
    reverse!(x)
    if x == sort(x)
      y = reverse(y, dims = 1)
    else
      reverse!(x)
      throw(DataError("unsorted x data; monotonic data needed", x))
    end
  end

  return x, y, (x = xsize, y = length(ysize) > 1 ? ysize[2:end] : 1), T
end #function checkinput


"""
    lindex(x::Array{Float64,1}, v::Float64)

Find the index of next lower `x` value in the original data compared to value `v`.
If `v` exists as `x` value, return the index of `x`.
"""
function lindex(x::Vector{<:Real}, v::Real)
    # Binary search
    len = size(x,1)
    li = 1
    ui = len
    mi = div(li+ui,2)
    done = false
    while !done
        if v<x[mi]
            ui = mi
            mi = div(li+ui,2)
        elseif v>x[mi+1]
            li = mi
            mi = div(li+ui,2)
        else
            done = true
        end
        if mi == li
            done = true
        end
    end
    return mi
end #function lindex


"""
    store(
      x::Vector{T},
      y::AbstractArray{T,N},
      s::AbstractArray{T,N},
      dx::Vector{T},
      divdif::AbstractArray{T,N},
      dim::Union{Int,Tuple{Vararg{Int}}}
    ) where T where N -> Polynomial{T}

Construct a `Polynomial` with all necessary data needed for `PCHIP` interpolation
from the `x` and `y` input data, the derivates (or slopes) at each point `s`,
the spaces in x-directon between points `dx`, `divdif` := dy / dx (where dy = diff(y, dims=1)),
and the dimension of the y-data in the first dimension `dim := size(y, 1)`.
"""
function store(
  x::Vector{T},
  y::AbstractArray{T,N} where N,
  s::AbstractArray{T,N} where N,
  dx::Vector{T},
  divdif::AbstractArray{T,N} where N,
  dim::Union{Int,Tuple{Vararg{Int}}}
) where T

  # Get dimensions
  n = length(x)
  d = size(y,2)
  # Duplicate dx to fill all dimensions
  dxd = repeat(dx, outer=[1,d]);

  # Calculate coefficients and instantiate the Polynomial
  dzzdx, dzdxdx = (divdif .- s[1:n-1,:])./dxd, (s[2:n,:] .- divdif) ./ dxd
  dnm1 = (n-1)d
  polynomial(x, [reshape(((dzdxdx .- dzzdx) ./ dxd)', dnm1, 1);
    reshape((2dzzdx .- dzdxdx)', dnm1, 1); reshape(s[1:n-1,:]', dnm1, 1);
    reshape(y[1:n-1,:]', dnm1, 1)], dim)
end #function store


"""
    pchipslopes(x,y,del)

Calculate the first derivatives `d(i) = P'(x(i))` at each point (`x`, `y`) for
the shape-preserving Piecewise Cubic Hermite Interpolation.

`d(i)` is the weighted avaerage of `del(i-1`) and `del(i)` when they have the same sign.
`d(i) = 0` when `del(i-1`) and `del(i)` have opposites signs or either is zero.
"""
function pchipslopes(x,y,del)
  # Special case n=2, use linear interpolation.
  n = length(x)
  if n==2
    return repeat(del[1:1], outer=[size(y)...])
  end

  # Preallocate d with zeros
  d = zeros(size(y))

  # Slopes at interior points.
  # d(k) = weighted average of del(k-1) and del(k) when they have the same sign.
  # d(k) = 0 when del(k-1) and del(k) have opposites signs or either is zero.

  # Interior indices
  k = findall(sign.(del[1:n-2]).*sign.(del[2:n-1]) .> 0)
  # Compute weights, min/max
  h = diff(x)
  hs = h[k]+h[k.+1]
  w1 = (h[k] .+ hs) ./ 3hs
  w2 = (hs .+ h[k.+1]) ./ 3hs
  dmax = max(abs.(del[k]), abs.(del[k.+1]))
  dmin = min(abs.(del[k]), abs.(del[k.+1]))
  # Compute slopes at interior points
  d[k.+1] = dmin./conj.(w1.*(del[k]./dmax) .+ w2.*(del[k.+1]./dmax));

  # Compute slopes at end points
  # Set d(1) and d(n) via non-centered, shape-preserving three-point formulae
  d[1] = ((2h[1] + h[2]) * del[1] - h[1] * del[2]) / (h[1] + h[2])
  if sign(d[1]) ≠ sign(del[1])
    d[1] = 0
  elseif sign(del[1]) ≠ sign(del[2]) && abs(d[1]) > abs(3del[1])
    d[1] = 3del[1]
  end
  d[n] = ((2h[n-1] + h[n-2]) * del[n-1] - h[n-1] * del[n-2]) / (h[n-1] + h[n-2])
  if sign(d[n]) ≠ sign(del[n-1])
    d[n] = 0
  elseif sign(del[n-1]) ≠ sign(del[n-2]) && abs(d[n]) > abs(3del[n-1])
    d[n] = 3del[n-1]
  end

  return d
end #function pchipslopes


"""
    polynomial(
      breaks::Vector{<:Real},
      coeffs::AbstractArray{T,N} where T<:Real where N,
      d::Int=1
    )

Construct a `Polynomial{T}` from the `breaks`, `coeffs`, and y-dimensions `d`.
"""
function polynomial(
  breaks::Vector{T},
  coeffs::AbstractArray{T,N} where N,
  d::Union{Int,Tuple{Vararg{Int}}}=1
) where T
  # Calculate order of Polynomial, number of breaks and dimensions to reshape coeff matrix
  dlk=length(coeffs); l=length(breaks)-1; dl=prod(d)*l
  k = dlk/dl+100eps(); k = k < 0 ? ceil(Int, k) : floor(Int, k)
  k ≤ 0 || dl*k ≠ dlk &&  @error "Number mismatch in coefficents" l d dlk

  # Instantiate and return the Polynomial
  Polynomial{T}(breaks, reshape(coeffs, dl ,k), l, k, d)
end #function polynomial
