# Interior grid

struct ChebyshevInteriorGrid{T} <: Grid{T}
    N::Int

    function ChebyshevInteriorGrid{T}(N::Integer) where T
        N >= 1 || throw(DomainError(N, "Must be positive."))
        return new{T}(N)
    end
end

function interior_gridpoint(T, i::Integer, N::Integer)
    1 <= i <= N || throw(DomainError(i, "i must be in 1 <= i <= N = $N."))
    return cospi(T(2*i-1)/(2*N))
end

gridpoint(ig::ChebyshevInteriorGrid{T}, i) where T = interior_gridpoint(T, i, ig.N)

function interior_cardinal(T, j::Integer, N::Integer, x)
    t = acos(x)
    xj = interior_gridpoint(T, j, N)
    tj = acos(xj)
    return cos(N*t)*sin(tj)/(N*sin(N*tj)*(cos(t) - cos(tj)))
end

cardinal(ig::ChebyshevInteriorGrid{T}, j, x) where T = interior_cardinal(T, j, ig.N, x)

function interior_derivative(T, i::Integer, j::Integer, N::Integer, order::Integer)::T
    # return _interior_gridpoint_derivative(T, i, j, N, Val(order))

    if order == 0
        return T(i == j)
    elseif order == 1
        xi = interior_gridpoint(T, i, N)
        xj = interior_gridpoint(T, j, N)

        return if i == j
            xj/2/(1-xj^2)
        else 
            (-1)^(i+j) * sqrt((1-xj^2)/(1-xi^2)) / (xi - xj)
        end
    elseif order == 2
        xi = interior_gridpoint(T, i, N)
        xj = interior_gridpoint(T, j, N)

        return if i == j
            xj^2/(1-xj^2)^2 - T(N^2 - 1)/(3*(1-xj^2))
        else
            delta1 = interior_derivative(T, i, j, N, 1)
            delta1 * (xi/T(1-xi^2) - 2/(xi-xj))
        end				   
    end

    throw(DomainError(order, "No derivative for this order."))
end

derivative(ig::ChebyshevInteriorGrid{T}, i, j, order) where T = interior_derivative(T, i, j, ig.N, order)

Base.eachindex(ig::ChebyshevInteriorGrid) = 1:ig.N
domain(_::ChebyshevInteriorGrid{T}) where T = (-one(T), one(T))
