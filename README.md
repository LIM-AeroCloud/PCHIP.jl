PCHIP.jl
========

Lightweight package for fitting of 2D data with Piecewise Cubic Hermite Interpolating Polynomial (PCHIP)
based on [reflectance-fitting](https://github.com/rsturley/reflectance-fitting.git)
by R. Steven Turley ([Turley, R. Steven, "Cubic Interpolation with Irregularly-Spaced Points in Julia 1.4"
(2018). Faculty Publications. 2177.](https://scholarsarchive.byu.edu/facpub/2177)).


Installation
------------

_PCHIP_ is an unregistered Julia package, but can be installed by the package manager with: 

```julia
julia> ]
pkg> add https://github.com/pb866/PCHIP.jl.git
pkg> instantiate
```

_PCHIP_ is tested against the latest stable release and the long-term support (LTS) version.


Usage
-----

Import the package with

```julia
using PCHIP
```

_PCHIP_ has 2 exported function and a exported struct:
- `pchip`
- `interpolate`
- `PCHIPdata`

Load your `x` and `y` data with the `pchip` function to the `PCHIPdata{T<:Real}` struct, 
which, in addition to the `x` and `y` data, stores the slopes (or derivatives) at 
each data point together with the spaces between the i-th and (i-1)-th data point. 
Vectors of `x` and `y` musst be subtypes of real and all vectors will be promoted 
to a common type.  
Additionally, `x` data musst be monotonic (either ascending or descending).

```julia
pchip(x::Vector{T1}, y::Vector{T2}) where {T1<:Real, T2<:Real} -> PCHIP{T}
```

The `interpolate` function evaluates the `PCHIPdata{T}` struct at a given data point 
(or value(s) `v` of type `Real`, `Vector{<:Real}` or `AbstractRange{<:Real}`). 
The returned interpolated `y` value gets promoted to the same 
type `T` as the data held in `PCHIPdata{T}` no matter of the type of the input value.

```julia
interpolate(pc::PCHIP{T}, v, eps::Real=1e-4)::T where {T} -> y::T
```

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
p = pchip(x,y)
v = interpolate(p, -3:0.01:-1)
```
