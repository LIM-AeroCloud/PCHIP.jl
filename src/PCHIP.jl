"""
# Module PCHIP

Piecewise Cubic Hermite Interpolating Polynomial (PCHIP) for interpolation of data
with monotonic x data.

The module is based on MATLAB's PCHIP routine and the works of:

F. N. Fritsch and R. E. Carlson, "Monotone Piecewise Cubic Interpolation",
SIAM J. Numerical Analysis 17, 1980, 238-246.

and

David Kahaner, Cleve Moler and Stephen Nash, "Numerical Methods and Software",
Prentice Hall, 1988.
"""
module PCHIP

## Export functions
export Polynomial, pchip, interpolate


## Define Structs
"""
    Polynomial{T<:Real}

`Polynomial` that stores all necessary data needed for peacewise cubic Hermite
interpolation. The structs holds the following fields:

- `breaks::Vector{T}`: Vector with breaks between the start and the end point of
  the data range, where the second derivative is allowed to be discontinuous
  corresponding to the `x` values in true or measured data. Interpolation occurs
  in the intervals between the `breaks`.
- `coeffs::AbstractArray{T,N} where N`: matrix with coefficients of the polynomial
  at any given point in the `x` range; matrix is of size `pieces⋅dim` × `order`
- `pieces::Int`: number of intervals between `breaks`
- `order::Int`: order of the `Polynomial`
- `dim::Int`: dimension of the y data input in x-direction
"""
struct Polynomial{T<:Real}
  breaks::Vector{T}
  coeffs::AbstractArray{T,N} where N
  pieces::Int
  order::Int
  dim::Union{Int,Tuple{Vararg{Int}}}
end

# Ensure Polynomial is seen as scaler during broadcasting
Broadcast.broadcastable(p::Polynomial) = Ref(p)


## Exception handling

"""
    RangeError(val::Real, pc::PCHIPdata{T}) where {T}

Warn when `val` is out of x range in `pc`.
"""
struct RangeError <: Exception
  val::Real
  range::NamedTuple{(:min,:max),Tuple{Real,Real}}

  function RangeError(val::Real, pc::Polynomial{T}) where {T}
    range = (min=pc.breaks[1], max=pc.breaks[end])
    new(val, range)
  end
end
# Format Error message
Base.showerror(io::IO, e::RangeError) = print(io, typeof(e), ": $(e.val) not within data range $(e.range.min)..$(e.range.max)")


"""
    DataError(msg, data)

Warn of any misfits or errors in `data` with a message (`msg`).
Warnings can include unsorted data or data that is not monotonic.
"""
struct DataError <: Exception
  msg::String
end
# Format Error message
Base.showerror(io::IO, e::DataError) = print(io, typeof(e), ": ", e.msg)
# Define default alarm and info message

## Include files
include("auxiliary.jl") # Helper functions for PCHIP interpolation

## Public functions

"""
    pchip(x::Vector{<:Real}, y::AbstractArray{T,N} where T<:Real where N) -> PCHIPdata{T}

Create the `Polynomial{T}` needed for piecewise cubic Hermite interpolation

# Arguments
- `x`: a vector of monotonic x values at which the function is known
- `y`: an vector or n-dimensional matrix of y values corresonding to these x values
"""
function pchip(x::Vector{<:Real}, y::AbstractArray{T,N} where T<:Real where N)
  # Check input data and attempt to correct faulty data
  x, y, dim, T = checkinput(x, y)
  # Get slopes at each point
  h, m = diff(x), prod(dim.y)
  del = diff(y,dims=1) ./ repeat(h, outer=[1, m])

  # Compute first derivatives at each break point
  slopes = zeros(T, size(y))
  for r = 1:m
    if isreal(del)
      slopes[:,r] = pchipslopes(x, y[:,r], del[:,r])
    else
      realslopes = pchipslopes(x,y[:,r],real.(del[:,r]))
      imagslopes = pchipslopes(x,y[:,r],imag(del[:,r]))
      slopes[:,r] = complex.(realslopes, imagslopes)
    end
  end
  # Save data in a Polynomial struct
  store(x,y,slopes,h,del,dim.y)
end #function pchip


"""
    interpolate(pc::Polynomial{T}, v) where {T} -> output

Interpolate `pc` to `v`, where `v` is a `Real`, `Vector{<:Real}` or `AbstractRange{<:Real}`.
Only values within the original data range are allowed.

Returns the interpolated `output` in the following forms:
- if original `y` data is a vector and `v <: Real`, a `Real` of type `v` is returned
- if `y` is a vector and `v` is a vector or range, a `Vector` of type `v` is returned
- if `y` is an N×M matrix, a vector of length M is returned for a `v` of type `Real`
  and a length(v)×M matrix, if `v` is a vector or range
- for N-dimensional input with `N ≥ 3`, output is of the dimensions in `x`-direction
  for a `Real` `v` or a vector of those output arrays for vectors or ranges of `v`
"""
function interpolate(pc::Polynomial{T}, v::Real, eps::Real=1e-4) where {T}
  # Check that v is within data rangee
  prevfloat(pc.breaks[1]) ≤ v ≤ nextfloat(pc.breaks[end]) || throw(RangeError(v, pc))

  # Get index of next lower data point
  i = lindex(pc.breaks, v)
  # Get index range in coefficient matrix for higher dimensional data
  r = (i-1)prod(pc.dim)+1:prod(pc.dim)i

  # Go to local coordinates
  vi = repeat([v - pc.breaks[i]], prod(pc.dim))

  # Apply nested multiplications
  val = pc.coeffs[r,1]
  for i=2:pc.order
    val = val .* vi .+ pc.coeffs[r,i]
  end
  # Reshape output to the correct dimensions
  pc.dim isa Real ? val[1] : reshape(val, pc.dim)

end #function interpolate

# Method for interpolating data ranges
function interpolate(pc::Polynomial, x::Union{Vector{<:Real},AbstractRange{<:Real}})
  vals = interpolate.(pc, x)
  vals = pc.dim isa Tuple && length(pc.dim) == 1 ?
  Matrix(hcat(vals...)') : vals
end

end # module PCHIP
