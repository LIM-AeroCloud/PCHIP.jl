module PCHIP


"""
    PCHIP(x,y,d,h)

Concrete type for holding the data needed for a
Piecewise Cubic Hermite Interpolating Polynomial (PCHIP)

- `x::Vector{Float64}`: strictly monotonic x data
- `x::Vector{Float64}`: strictly monotonic x data
- `d::Vector{Float64}`: slopes at the data points
- `h::Vector{Float64}`: spaces between ith and (i-1)th data point
"""
struct PCHIP
    x::Union{Array{Float64,1},
        StepRangeLen{Float64,
            Base.TwicePrecision{Float64},
            Base.TwicePrecision{Float64}}}
    y::Array{Float64,1}
    d::Array{Float64,1}
    h::Array{Float64,1}
end


# Real PCHIP
"""
    pchip(x,y)

Create the PCHIP structure needed for piecewise
continuous cubic spline interpolation

# Arguments
- `x`: an array of x values at which the function is known
- `y`: an array of y values corresonding to these x values
"""
function pchip3(x::Array{Float64,1}, y::Array{Float64,1})
    len = size(x,1)
    if len<3
        error("PCHIP requires at least three points for interpolation")
    end
    h = x[2:len].-x[1:len-1]
    del = (y[2:len].-y[1:len-1])./h
    # Pre-allocate and fill columns and diagonals
    d = zeros(len)
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
    PCHIP(x,y,d,h)
end


"""
    interp(cs::CubicSpline, v::Float)

Interpolate to the value corresonding to v

# Examples
```
x = cumsum(rand(10))
y = cos.(x);
cs = CubicSpline(x,y)
v = interp(cs, 1.2)
```
"""
function interp(pc::PCHIP, v::Float64, eps::Float64=1e-4)

    if v*(1+eps)<first(pc.x)
        error("Extrapolation not allowed, $v<$(first(pc.x))")
    end
    if v*(1-eps)>last(pc.x)
        error("Extrapolation not allowed, $v>$(last(pc.x))")
    end
    i = region(pc.x, v)
    phi(t) = 3*t^2 - 2*t^3
    psi(t) = t^3 - t^2
    H1(x) = phi((pc.x[i+1]-v)/pc.h[i])
    H2(x) = phi((v-pc.x[i])/pc.h[i])
    H3(x) = -pc.h[i]*psi((pc.x[i+1]-v)/pc.h[i])
    H4(x) = pc.h[i]*psi((v-pc.x[i])/pc.h[i])
    pc.y[i]*H1(v) + pc.y[i+1]*H2(v) + pc.d[i]*H3(v) + pc.d[i+1]*H4(v)
end


function region(x::Array{Float64,1}, v::Float64)
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
    mi
end

end # module
