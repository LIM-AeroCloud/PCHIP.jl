"""
# Module PCHIP

Piecewise Cubic Hermite Interpolating Polynomial (PCHIP) for interpolation of data
with monotonic x data.

This module is based on the reflectance-fitting package
(https://github.com/rsturley/reflectance-fitting.git) by R. Steven Turley
([Turley, R. Steven, "Cubic Interpolation with Irregularly-Spaced Points in Julia 1.4"
(2018). Faculty Publications. 2177.](https://scholarsarchive.byu.edu/facpub/2177)).

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
- `coeffs::AbstractArray{T,N} where N`: matrix with coefficients needed to solve
  the polynomial at any given point in the `x` range; matrix is of size `pieces⋅dim` × `order`
- `pieces::Int`: number of intervals between `breaks`
- `order::Int`: order of the `Polynomial`
- `dim::Int`: dimension of the y data input

Construct `Polynomial` from the `breaks`, `coeffs`, and y-dimensions `d`:

    Polynomial(
      breaks::Vector{<:Real},
      coeffs::AbstractArray{T,N} where T<:Real where N,
      d::Int=1
    )
"""
struct Polynomial{T<:Real}
  breaks::Vector{T}
  coeffs::AbstractArray{T,N} where N
  pieces::Int
  order::Int
  dim::Union{Int,Tuple{Vararg{Int}}}


  """ Modified internal constructor for `Polynomial`"""
  function Polynomial(
    breaks::Vector{T},
    coeffs::AbstractArray{T,N} where N,
    d::Union{Int,Tuple{Vararg{Int}}}=1
  ) where T
    dlk=length(coeffs); l=length(breaks)-1; dl=prod(d)*l
    k = dlk/dl+100*eps(); k = k < 0 ? ceil(Int, k) : floor(Int, k)
    k ≤ 0 || dl*k ≠ dlk &&  @error "MATLAB:mkpp:PPNumberMismatchCoeffs" l d dlk

    new{T}(breaks, reshape(coeffs, dl ,k), l, k, d)
  end #constructor Polynomial
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
    # range = pc.x[1] < pc.x[end] ? (min=pc.x[1], max=pc.x[end]) :
    #   (min=pc.x[end], max=pc.x[1])
    range = (min=pc.x[1], max=pc.x[end])
    new(val, range)
  end
end
# Format Error message
Base.showerror(io::IO, e::RangeError) = print(io, typeof(e), ": $(e.val) not within data range $(e.range.min)..$(e.range.max)")


"""
    DataError(msg, data)

Warn of any misfits or errors in `data` with a message (`msg`).
Warnings can include unsorted data or data that is not monotonic ascending.
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
    pchip(x::Vector{T1}, y::Vector{T2}) where {T1<:Real, T2<:Real} -> PCHIPdata{T}

Create the `Polynomial` needed for piecewise cubic Hermite interpolation

# Arguments
- `x`: a vector of x values at which the function is known
- `y`: an vector or n-dimensional matrix of y values corresonding to these x values
"""
function pchip(x::Vector{<:Real}, y::AbstractArray{T,N} where T<:Real where N)
  # Check input data and attempt to correct faulty data
  x, y, dim, T = checkinput(x, y)
  # Get slopes at each point
  h, m = diff(x), prod(dim.y)
  del = diff(y,dims=1) ./ repeat(h, outer=[1, m])

  # Compute slopes
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
  store(x,y,slopes,h,del,dim.y)
end #function pchip


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
  y::AbstractArray{T,N},
  s::AbstractArray{T,N},
  dx::Vector{T},
  divdif::AbstractArray{T,N},
  dim::Union{Int,Tuple{Vararg{Int}}}
) where T where N

  n = length(x)
  d = size(y,2)
  dxd = repeat(dx, outer=[1,d]);

  dzzdx, dzdxdx = (divdif .- s[1:n-1,:])./dxd, (s[2:n,:] .- divdif) ./ dxd
  dnm1 = (n-1)d
  Polynomial(x, [reshape(((dzdxdx .- dzzdx) ./ dxd)', dnm1, 1);
    reshape((2dzzdx .- dzdxdx)', dnm1, 1); reshape(s[1:n-1,:]', dnm1, 1);
    reshape(y[1:n-1,:]', dnm1, 1)], dim)
end #function store


"""
    interpolate(pc::PCHIPdata{T}, v, eps::Real=1e-4)::T where {T} -> PCHIPdata{T}

Interpolate `pc` to `v`, where `v` is a `Real`, `Vector{<:Real}` or `AbstractRange`.
Only values within the original data range are allowed. To account for floating point
rounding errors, a correction factor of (1±ε) is applied to the lower/upper bounds:

    pc.x[1]*(1-sign(pc.x[1])*eps)
    pc.x[end]*(1+sign(pc.x[end])*eps)


# Examples
```
x = cumsum(rand(10))
y = cos.(x)
p = pchip(x,y)
v = interpolate(p, 1.2)
```
"""
function interpolate(pc::Polynomial{T}, v::Real, eps::Real=1e-4)::T where {T}

    v > pc.x[1]*(1-sign(pc.x[1])*eps) && v < pc.x[end]*(1+sign(pc.x[end])*eps) ||
      throw(RangeError(v, pc))
    i = lindex(pc.x, v)
    phi(t) = 3*t^2 - 2*t^3
    psi(t) = t^3 - t^2
    H1(x) = phi((pc.x[i+1]-v)/pc.h[i])
    H2(x) = phi((v-pc.x[i])/pc.h[i])
    H3(x) = -pc.h[i]*psi((pc.x[i+1]-v)/pc.h[i])
    H4(x) = pc.h[i]*psi((v-pc.x[i])/pc.h[i])
    pc.y[i]*H1(v) + pc.y[i+1]*H2(v) + pc.d[i]*H3(v) + pc.d[i+1]*H4(v)
end #function interpolate

# Methods for interpolating data ranges
interpolate(pc::Polynomial, x::Vector{<:Real}, eps::Real=1e-4) = interpolate.(pc, x, eps)
interpolate(pc::Polynomial, x::AbstractRange{<:Real}, eps::Real=1e-4) = interpolate.(pc, x, eps)



## Private functions

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

end # module PCHIP
