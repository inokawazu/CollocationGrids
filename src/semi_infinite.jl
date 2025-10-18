# Semi infinite grid 
# TODO Not an accurate grid yet, so make accurate

struct SemiInfiniteGrid{T} <: Grid{T}
    grid::Grid{T}
    L::T

    function SemiInfiniteGrid(grid::Grid{T}, L) where T
        L > 0 || DomainError(L, "L must be positive.")
        return new{T}(grid, L)
    end
end

function gridpoint(semi::SemiInfiniteGrid{T}, i) where T
    l, u = domain(semi.grid)
    (; L) = semi
    tosemi(x) = L*(x-l)/(u-x)
    return tosemi(gridpoint(semi.grid, i))
end

function cardinal(semi::SemiInfiniteGrid, j::Integer, y)
    l, u = domain(semi.grid)
    (; L) = semi
    x = (l*L + u*y) / (L + y)
    return cardinal(semi.grid, j, x)
end

function derivative(semi::SemiInfiniteGrid{T}, i, j, order) where T
    l, u = domain(semi.grid)
    (; L) = semi
    fromsemi(y) = (l*L + u*y) / (L + y)
    tosemi(x) = L*(x-l)/(u-x)

    if order == 0
        return T(i == j)
    elseif order == 1
        # d/dy = dx/dy × d/dx
        dCdx = derivative(semi.grid, i, j, 1) # dCij/dx 
        yj = gridpoint(semi, j)
        dxjdy = TaylorDiff.derivative(fromsemi, yj, 1) # dx/dy(yj)
        return dxjdy * dCdx
    elseif order == 2
        dCdx = derivative(semi.grid, i, j, 1) # # dCij/dx
        d2Cijdx2 = derivative(semi.grid, i, j, 2) # d²Cij/dx²
        yj = gridpoint(semi, j)
        dxdy = TaylorDiff.derivative(fromsemi, yj, 1) # dx/dy
        d2xdy2 = TaylorDiff.derivative(fromsemi, yj, 2) # d^2x/dy^2
        return dxdy^2 * d2Cijdx2 + dxdy * d2xdy2 * dCdx
    end

    throw(DomainError(order, "No derivative for this order."))
end

Base.eachindex(semi::SemiInfiniteGrid) = eachindex(semi.grid)

function domain(semi::SemiInfiniteGrid{T}) where T 
    l, u = domain(semi.grid)
    (; L) = semi
    tosemi(x) = L*(x-l)/(u-x)
    tosemi.(domain(semi.grid))
end

