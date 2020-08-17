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
  @show ysize
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
  length(ysize) > 2 && (y = reshape(y, ysize[1], prod(ysize[2:end])))
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
    pchipslopes(x,y,del)

Calculate the first derivatives `d(k) = P'(x(k))` at each point (`x`, `y`) for
the shape-preserving Piecewise Cubic Hermite Interpolation.

`d(k)` is the weighted avaerage of `del(k-1`) and `del(k)` when they have the same sign.
`d(k) = 0` when `del(k-1`) and `del(k)` have opposites signs or either is zero.
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
