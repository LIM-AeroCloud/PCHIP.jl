PCHIP
=====

Package for fitting of data with a Piecewise Cubic Hermite Interpolating Polynomial (PCHIP).


Installation
------------

_PCHIP_ is an unregistered Julia package, but can be installed by the package manager with: 

```julia
julia> ]
pkg> add https://github.com/pb866/PCHIP.jl.git
pkg> instantiate
```

_PCHIP_ is tested against the latest stable release (1.5.2).


Usage
-----

Import the package with

```julia
using PCHIP
```

_PCHIP_ has 2 exported functions and an exported struct:
- `pchip`
- `interpolate`
- `Polynomial`

Load your `x` and `y` data with the `pchip` function to the `Polynomial{T<:Real}` struct, 
which stores the `breaks` at which the second order derivate is allowed to be discontinuous
and corresponds to the `x` coordinates of the original or measured data. Additional
fields of `Polynomial` are `coeffs` with the matrix holding the coefficients of the 
polynomial, the number of `intervals` between breaks, the `order` of the polynomial,
and the dimensions `dim` of the `y` data in `x`-direction.

The `x` data must be a monotonic vector (either ascending or descending); `y` data 
must be of the same length as `x` data in the first dimension. Duplicate `x` values
are allowed; at sites with multiple `x` values, `NaN` is returned as corresponding 
`y` value.

Optionally, instead of the polynomial, the interpolated values can be returned
if `xi` are given to `pchip` as third optional argument; `xi` must be a `Real` number,
a `Vector{<:Real}` or an `AbstractArray{<:Real}`.

```julia
pchip(x::Vector{<:Real}, y::AbstractArray{T,N} where T<:Real where N, xi=nothing)
```

The `interpolate` function evaluates the `Polynomial{T}` struct at a given data point 
or data range (`v` of type `Real`, `Vector{<:Real}` or `AbstractRange{<:Real}`). 
The returned interpolated `y` value gets promoted to the same 
type `T` as the data held in `Polynomial{T}` no matter of the type of the input value.

```julia
interpolate(pc::Polynomial{T}, v::Real) where {T}
```

Output of `interpolate` or `pchip` if `xi` is given is of the following format:

- if original `y` data is a vector and `v <: Real`, a `Real` of type `v` is returned
- if `y` is a vector and `v` is a vector or range, a `Vector` of element type `v` is returned
- if `y` is an N×M matrix, a vector of length M is returned for a `v` of type `Real`
  and a length(v)×M matrix, if `v` is a vector or range
- for N-dimensional input with `N ≥ 3`, output is of the dimensions in `x`-direction
  of the `y` data for a `Real` `v` or a vector of those output arrays for vectors or ranges of `v`
  

Examples
--------

```julia
x = cumsum(rand(10))
y = cos.(x)
p = pchip(x,y)
v = interpolate(p, [1.2, 2.1, 3.0])
```


```julia
x = cumsum(-rand(10))
y = cos.(x)
p = pchip(x,y, -3:0.01:-1)
```
