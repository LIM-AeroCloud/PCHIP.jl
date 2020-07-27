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
export PCHIP, pchip, interpolate


## Define Structs
"""
    PCHIP(x,y,d,h)

Concrete type for holding the data needed for a
Piecewise Cubic Hermite Interpolating Polynomial (PCHIP)

- `x::Vector{Float64}`: strictly monotonic x data
- `x::Vector{Float64}`: strictly monotonic x data
- `d::Vector{Float64}`: slopes at the data points
- `h::Vector{Float64}`: spaces between ith and (i-1)th data point
"""
struct PCHIP{T<:Real}
    x::Vector{T}
    y::Vector{T}
    d::Vector{T}
    h::Vector{T}
end #struct PCHIP

# Ensure PCHIP is seen as scaler during broadcasting
Broadcast.broadcastable(p::PCHIP) = Ref(p)

## Exception handling

"""
    RangeError(val::Real, pc::PCHIP{T}) where {T}

Warn when `val` is out of x range in `pc`.
"""
struct RangeError <: Exception
  val::Real
  range::NamedTuple{(:min,:max),Tuple{Real,Real}}

  function RangeError(val::Real, pc::PCHIP{T}) where {T}
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
  data::Vector{<:Real}
end
# Format Error message
Base.showerror(io::IO, e::DataError) = print(io, typeof(e), ": ", e.msg, "\n", e.data)
# Define default alarm and info message
DataError(msg::String) = DataError(msg, Real[])
DataError(data::Vector{<:Real}) = DataError("data not sorted ascendingly", data)
DataError() = DataError("", Real[])

## Public functions

"""
    pchip(x,y)

Create the PCHIP structure needed for piecewise
continuous cubic spline interpolation

# Arguments
- `x`: an array of x values at which the function is known
- `y`: an array of y values corresonding to these x values
"""
function pchip(x::Vector{T1}, y::Vector{T2}) where {T1<:Real, T2<:Real}
    len = size(x,1)
    if len<3
      throw(DataError("PCHIP requires at least three points for interpolation", x))
    elseif x ≠ sort(x)
      throw(DataError(x))
    end
    T = promote_type(T1, T2)
    T<:Integer && (T = float(T))
    h = x[2:len].-x[1:len-1]
    del = (y[2:len].-y[1:len-1])./h
    # Pre-allocate and fill columns and diagonals
    d = zeros(T, len)
    d[1] = del[1]
    for i=2:len-1
        if del[i]*del[i-1] < 0
            d[i] = 0
        else
            d[i] = (del[i]+del[i-1])/2
        end
    end
    d[len] = del[len-1]
    for i=1:len-1
        if del[i] == 0
            d[i] = 0
            d[i+1] = 0
        else
            alpha = d[i]/del[i]
            beta = d[i+1]/del[i]
            if alpha^2+beta^2 > 9
                tau = 3/sqrt(alpha^2+beta^2)
                d[i] = tau*alpha*del[i]
                d[i+1] = tau*beta*del[i]
            end
        end
    end
    PCHIP{T}(x,y,d,h)
end #function pchip


"""
    interpolate(pc::PCHIP{T}, v::Real, eps::Real=1e-4)::T where {T} -> PCHIP{T}

Interpolate `pc` to the value `v` corresonding to v.
Only values within the original data range ± `eps` to account for floating point
inaccuracies are allowed.

# Examples
```
x = cumsum(rand(10))
y = cos.(x);
cs = PCHIP(x,y)
v = interpolate(cs, 1.2)
```
"""
function interpolate(pc::PCHIP{T}, v::Real, eps::Real=1e-4)::T where {T}

    (v*(1+eps)<first(pc.x) || v*(1-eps)>last(pc.x)) && throw(RangeError(v, pc))
    i = lbound(pc.x, v)
    phi(t) = 3*t^2 - 2*t^3
    psi(t) = t^3 - t^2
    H1(x) = phi((pc.x[i+1]-v)/pc.h[i])
    H2(x) = phi((v-pc.x[i])/pc.h[i])
    H3(x) = -pc.h[i]*psi((pc.x[i+1]-v)/pc.h[i])
    H4(x) = pc.h[i]*psi((v-pc.x[i])/pc.h[i])
    pc.y[i]*H1(v) + pc.y[i+1]*H2(v) + pc.d[i]*H3(v) + pc.d[i+1]*H4(v)
end #function interpolate


## Private functions

"""
    lbound(x::Array{Float64,1}, v::Float64)

Find the index of next lower `x` value in the original data compared to value `v`.
If `v` exists as `x` value, return the index of `x`.
"""
function lbound(x::Vector{<:Real}, v::Real)
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
end #function lbound

end # module PCHIP
