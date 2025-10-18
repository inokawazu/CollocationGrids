# ChebyshevRecombinationGrid - Experimental

struct ChebyshevRecombinationGrid{T, V} <: Grid{T}
    N::Int
    bc_left::Symbol
    bc_right::Symbol
    gridpoints::V

    function ChebyshevRecombinationGrid{T}(N::Integer, left, right) where T
        N > 0 || DomainError(N, "N must be greater than zero.")
        bf(x) = recombi_basis_function(big(T), left, right, N, x)
        gridpoints = find_zeros(bf, -eps(T)-1, 1+eps(T))
        filter!(!isone∘abs, gridpoints)
        gridpoints = T.(gridpoints)
        V = typeof(gridpoints)
        return new{T,V}(N, left, right, gridpoints)
    end
end

ChebyshevRecombinationGrid{T}(N, leftright) where T = ChebyshevRecombinationGrid{T}(N, leftright, leftright)

Base.eachindex(csb::ChebyshevRecombinationGrid) = eachindex(csb.gridpoints)

domain(::ChebyshevRecombinationGrid{T}) where T = (-one(T), one(T))

gridpoint(csg::ChebyshevRecombinationGrid{T}, i) where T = csg.gridpoints[i]

function basis_function(csg::ChebyshevRecombinationGrid{T}, i, x) where T
    recombi_basis_function(T, csg.bc_left, csg.bc_right, i, x)
end

function recombi_basis_function(T, left, right, i, x)
    if left == right == :dirichlet
        return iseven(i) ? chebT(i, x) - 1 : chebT(i, x) - x
    elseif left == right == :neumann
        return if i == 0
            one(T)
        else
            chebT(i, x) - (i/(i+2))^2 * chebT(i+2, x)
        end
    end

    throw(DomainError((left, right), "No basis functions for provided BCs"))
end

function cardinal(csg::ChebyshevRecombinationGrid, i, x)
    phiNp1(x) = basis_function(csg, csg.N, x)
    xi = gridpoint(csg, i)
    return phiNp1(x)/(TaylorDiff.derivative(phiNp1, xi, 1)*(x-xi))
end

function derivative(csg::ChebyshevRecombinationGrid{T}, i, j, order) where T
    if order >= 0 && i == j
        phiNp1(x) = basis_function(csg, csg.N, x)
        xi = gridpoint(csg, i)
        cj2(x, order) = TaylorDiff.derivative(phiNp1, x, order+1)/((order + 1)*TaylorDiff.derivative(phiNp1, xi, 1))
        return cj2(xi, order)

    elseif order >= 0 && i != j
        xi = gridpoint(csg, i)
        cj(x) = cardinal(csg, j, x)
        return order == 0 ? cj(xi) : TaylorDiff.derivative(cj, xi, order)
    end
    throw(DomainError(order, "No derivative for this order."))
end

