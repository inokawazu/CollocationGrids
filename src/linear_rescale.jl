# LinearRescaleGrid

struct LinearRescaleGrid{T} <: Grid{T}
    grid::Grid{T}
    a::T
    b::T

    function LinearRescaleGrid(grid::Grid{T}, a, b) where T
        a < b || DomainError((a,b), "a and b must satisfy a < b.")
        return new{T}(grid, T(a), T(b))
    end
end

domain(lrg::LinearRescaleGrid) = (lrg.a, lrg.b)
Base.eachindex(lrg::LinearRescaleGrid) = eachindex(lrg.grid)

function gridpoint(lrg::LinearRescaleGrid, i)
    xi = gridpoint(lrg.grid, i)
    return linear_rescale(xi, domain(lrg.grid)..., lrg.a, lrg.b)
end

function cardinal(lrg::LinearRescaleGrid, i, y)
    x = linear_rescale(y, lrg.a, lrg.b, domain(lrg.grid)...)
    return cardinal(lrg.grid, i, x)
end

function derivative(lrg::LinearRescaleGrid, i, j, order)
    factor = linear_rescale_factor(lrg.a, lrg.b, domain(lrg.grid)...)^order # dx/dy 
    preder = derivative(lrg.grid, i, j, order) # d/dx
    return factor * preder
end


function linear_rescale(g::Grid{T}, x, a::Real, b::Real) where T
    return LinearRescaleGrid{T}(g, a, b)
end

"""
linear_rescale(x, a, b, c, d)

Linearly rescales x from [a,b] -> [c,d]
"""
linear_rescale(x, a, b, c, d) = (d - c) * (x - a) / (b - a) + c

"""
linear_rescale_factor(a, b, c, d)

returns the slope of the transformation, `(d - c) / (b - a)`.
"""
linear_rescale_factor(a, b, c, d) = (d - c) / (b - a) # df/dx
