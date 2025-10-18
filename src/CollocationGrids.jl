module CollocationGrids

export gridpoints, 
       boundaries,
       derivative_matrix, 
       interpolate, 
       derivative,
       derivative_array_matrix,
       cardinal,
       linear_rescale

export MultiGrid, 
       Grid,
       ChebyshevLobattoGrid,
       ChebyshevInteriorGrid,
       SemiInfiniteGrid,
       FourierGrid

export coordinate_operators, coordinate_vectors,
       function_operator, function_vector

import TaylorDiff, LinearAlgebra
using Roots: find_zeros
using LinearAlgebra
using SparseArrays

abstract type Grid{T} end

# File for grids
include("linear_rescale.jl")
include("cheb.jl")
include("one_hots.jl")

gridpoints(g::Grid{T}) where T = ( gridpoint(g, i) for i in eachindex(g) )

derivative_matrix(g::Grid{T}, order) where T = [
                                                derivative(g, i, j, order) 
                                                for i in eachindex(g), j in eachindex(g)
                                               ]

function interpolate(g::Grid{T}, coefs, x) where T 
    sum(coef * cardinal(g, j, x) for (j, coef) in zip(eachindex(g), coefs))
end

Base.length(g::Grid) = length(eachindex(g))
Base.size(g::Grid) = size(eachindex(g))

boundaries(g::Grid) = sort(extrema(gridpoints(g)))

# MultiGrid

struct MultiGrid{T, N, V} <: Grid{T}
    grids::V
end

Base.iterate(mg::MultiGrid) = iterate(mg.grids)
Base.iterate(mg::MultiGrid, state) = iterate(mg.grids, state)

MultiGrid(grids::Grid{T}...) where T = MultiGrid{T, length(grids), typeof(grids)}(grids)

function Base.eachindex(mg::MultiGrid{T, N}) where {T, N}
    return CartesianIndices(eachindex.(mg.grids))
end

function derivative_matrix(mg::MultiGrid{T, N}, order) where {T, N}
    comp_mats = ntuple(i->sparse(derivative_matrix(mg.grids[i], order[i])), N)
    return kron(reverse(comp_mats)...)
end

function derivative(mg::MultiGrid{T}, I, J, order) where T
    mapreduce(derivative, *, mg.grids, Tuple(I), Tuple(J), order, init = one(T))
end

function derivative(mg::MultiGrid{T, 2}, I, J, order) where T
    derivative(mg.grids[1], I[1], J[1], order[1]) *
    derivative(mg.grids[2], I[2], J[2], order[2])
end

function derivative(mg::MultiGrid{T, 3}, I, J, order) where T
    derivative(mg.grids[1], I[1], J[1], order[1]) *
    derivative(mg.grids[2], I[2], J[2], order[2]) *
    derivative(mg.grids[3], I[3], J[3], order[3])
end

function derivative_array_matrix(grid::MultiGrid{T, N}, order::Int) where {T, N}
    derivative_matrix.(Ref(grid), one_hots(N, order))
end

derivative_matrix(mg::MultiGrid, order::Integer, k::Integer) = derivative_matrix(mg.grids[k], order)

function interpolate(mg::MultiGrid{T, N}, coefs, x) where {T, N}
    sum(coef * cardinal(mg, J, x) for (J, coef) in zip(eachindex(mg), coefs))
end

# function interpolate(mg::MultiGrid{T, N}, coefs::AbstractArray{G, N}, x::NTuple{N}) where {T, G, N}
#     interpolate(mg, coefs, x...)
# end

# function interpolate(mg::MultiGrid{T, N}, coefs::AbstractArray{G, N}, x::Vararg{U, N}) where {T, G, U, N}
#     sum(coef * cardinal(mg, J, x...) for (J, coef) in zip(eachindex(mg), coefs))
# end

function cardinal(mg::MultiGrid{T, N}, J, x) where {T, N}
    mapreduce(cardinal, *, mg.grids, Tuple( J ), x, init=one(T))
end

# function cardinal(mg::MultiGrid{T, N}, J, x) where {T, N}
#     cardinal(mg, J, x...)
# end

function gridpoint(mg::MultiGrid{T, N}, I) where {T, N} 
    return ntuple(d->gridpoint(mg.grids[d], I[d]), Val(N)) # gridpoint.(mg.grids, Tuple(I))
end

function gridpoint(mg::MultiGrid{T, 2}, I) where T
    return (
            gridpoint(mg.grids[1], I[1]),
            gridpoint(mg.grids[2], I[2]),
           )
end

function gridpoint(mg::MultiGrid{T, 3}, I) where T
    return (
            gridpoint(mg.grids[1], I[1]),
            gridpoint(mg.grids[2], I[2]),
            gridpoint(mg.grids[3], I[3]),
           )
end

domain(mg::MultiGrid{T, N}) where {T, N} = domain.(mg.grids)
boundaries(g::MultiGrid) = boundaries.(g.grids)

function derivative_matrices(mg::MultiGrid{T, N}, max_order::Integer) where {T, N}
    ders = Iterators.map(Iterators.product(1:N, 0:max_order)) do (dim, order)
        (dim, order) => derivative_matrix(mg.grids[dim], order)
    end

    return Dict(ders)
end

function apply_derivative_matrix(input::AbstractArray, der_mat::AbstractArray, dim::Integer)
    return apply_derivative_matrix!(similar(input), input, der_mat, dim) 
end

function apply_derivative_matrix!(
        out::AbstractArray, input::AbstractArray,
        der_mat::AbstractArray, dim::Integer)

    fill!(out, 0)

    for (out_i, i) in zip(eachslice(out, dims=dim), axes(der_mat, 1))
        for (in_j, j) in zip(eachslice(input, dims=dim), axes(der_mat, 2))
            @. out_i += der_mat[i, j] * in_j
        end
    end

    return out
end

# FourierGrid 

struct FourierGrid{T} <: Grid{T}
    N::Int

    function FourierGrid{T}(N::Integer) where T
        N >= 1 || throw(DomainError(N, "Must be positive."))
        new{T}(Int(N))
    end
end

function fourier_gridpoint(T, i::Integer, N::Integer)::T
    i >= 0 || throw(DomainError(i, "Must be positive."))
    i < 2N || throw(DomainError(i, "Must be less than 2N = $(2N)"))

    return π*i/N
end

function fourier_cardinal(T, j::Integer, N::Integer, x)
    N >= 0 || throw(DomainError(N, "Must be nonnegative."))
    xj = fourier_gridpoint(T, j, N)
    c(k) = abs(k) == N ? 2 : 1

    return sum(2*cos(k*(x-xj))/c(k) for k in 0:N)/(2N) - 1/(2N)
end

function fourier_derivative(T, i::Integer, j::Integer, N::Integer, order::Integer)::T
    if order == 0
        return T(i == j)
    elseif order == 1
        xi = fourier_gridpoint(T, i, N)
        xj = fourier_gridpoint(T, j, N)
        return i == j ? zero(T) : (-1)^(i+j)/2*cot(T(xi-xj)/2)
    elseif order == 2
        xi = fourier_gridpoint(T, i, N)
        xj = fourier_gridpoint(T, j, N)
        return i == j ? T(-(1+2*N^2)/6) : (-1)^(i+j+1)/2/(sin(T(xi-xj)/2))^2
    end

    return throw(DomainError(order, "No derivative for this order."))
end

gridpoint(fg::FourierGrid{T}, i) where T = fourier_gridpoint(T, i, fg.N)
cardinal(fg::FourierGrid{T}, j, x) where T = fourier_cardinal(T, j, fg.N, x)
derivative(fg::FourierGrid{T}, i, j, order) where T = fourier_derivative(T, i, j, fg.N, order)

domain(_::FourierGrid{T}) where T = (zero(T), T(2pi))

Base.eachindex(fg::FourierGrid) = 0:(2*fg.N-1)

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

# Cheb-Lobbo Grid

struct ChebyshevLobattoGrid{T} <: Grid{T}
    N::Int

    function ChebyshevLobattoGrid{T}(N::Integer) where T
        N >= 1 || throw(DomainError(N, "Must be positive."))
        new{T}(Int(N))
    end
end

function cheb_lob_gridpoint(T, i::Integer, N::Integer)
    0 <= i <= N || throw(DomainError(i, "Grid point must be 0 ≤ i ≤ N = $N."))

    return cospi(T(i)/N)
end

function cheb_lob_derivative(T, i::Integer, j::Integer, N::Integer, order::Integer)
    p(k) = k == 0 || k == N ? 2 : 1

    if order == 0
        return T(i == j)
    elseif order == 1
        return if i == j == 0
            T(1 + 2*N^2)/6
        elseif i == j == N
            -T(1 + 2*N^2)/6
        elseif i == j
            xj = cheb_lob_gridpoint(T, j, N)
            -xj/T(2*(1-xj^2))
        else
            xi = cheb_lob_gridpoint(T, i, N)
            xj = cheb_lob_gridpoint(T, j, N)
            (-1)^(i+j)*T(p(i)/(p(j) * (xi - xj)))
        end
    end

    return throw(DomainError(order, "No derivative for this order."))
end

function cheb_lob_cardinal(T, j::Integer, N::Integer, x)
    p(k) = k == 0 || k == N ? 2 : 1
    xj = cheb_lob_gridpoint(T, j, N)

    return 2/(N*p(j)) * sum(chebT(m, xj)*chebT(m, x)/p(m) for m in 0:N)
end

gridpoint(clb::ChebyshevLobattoGrid{T}, i) where T = cheb_lob_gridpoint(T, i, clb.N)
cardinal(clb::ChebyshevLobattoGrid{T}, j::Integer, x) where T = cheb_lob_cardinal(T, j, clb.N, x)
Base.eachindex(clb::ChebyshevLobattoGrid) = 0:clb.N

derivative(clb::ChebyshevLobattoGrid{T}, i, j, order) where T = cheb_lob_derivative(T, i, j, clb.N, order)
derivative_matrix(g::ChebyshevLobattoGrid{T}, order) where T = [
                                                                derivative(g, i, j, 1) 
                                                                for i in eachindex(g), 
                                                                j in eachindex(g)
                                                               ]^order

domain(::ChebyshevLobattoGrid{T}) where T = (-one(T), one(T))

# ChebyshevRecombinationGrid

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

# ChebyshevSpectralGrid

# Semi infinite grid (Not an accurate grid yet)
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

# Operators

coordinate_operators(grid::Grid) = Diagonal(vec(collect(gridpoints(grid))))
coordinate_vectors(grid::Grid) = vec(collect(gridpoints(grid)))

function coordinate_operators(grid::MultiGrid{T,N}) where {T,N}
    pts = vec(collect(gridpoints(grid)))
    ntuple(i->Diagonal(getindex.(pts, i)), N)
end

function coordinate_vectors(grid::MultiGrid{T,N}) where {T,N}
    pts = vec(collect(gridpoints(grid)))
    ntuple(i->vec(getindex.(pts, i)), N)
end

function_operator(f::Function, grid::Grid) = Diagonal(vec(f.(gridpoints(grid))))
function_vector(f::Function, grid::Grid) = vec(f.(gridpoints(grid)))

# Derivative Matrix Data

# struct DerivativeOperator{T}
#     der_mat::Matrix{T}
#     dim::Int
# end

# function make_derivative_operators(mg::MultiGrid{T,N}, max_order::Integer) where {T, N}
#     dmat_dicts = derivative_matrices(mg, max_order)

#     [ 
#      [DerivativeOperator{T}(dmat_dicts[(dim, order)], dim) for dim in 1:N]
#      for order in 0:max_order 
#     ]
# end

# mutating
# function (op::DerivativeOperator{T})(out, input) where T
#     apply_derivative_matrix!(out, input, op.der_mat, op.dim)
# end

# function (op::DerivativeOperator{T})(input) where T
#     apply_derivative_matrix(input, op.der_mat, op.dim)
# end

# struct ScalarOperator{T}
#     func::Function
#     grid::Grid{T}
# end

# Base.getindex(s::ScalarOperator{T}, inds...) where T = s.func(gridpoint(s.grid, inds...))

# for op in [:+, :*, :\, :/]
#     @eval function Base.$op(so1::ScalarOperator{T}, so2::ScalarOperator{T}) where T
#         so1.grid == so2.grid || error("Scalar operators must have the same grid.")
#         let
#             newfunc(x...) = $op(so1.func(x...), so2.func(x...))
#             return ScalarOperator{T}(newfunc, so1.grid)
#         end
#     end
# end

# function get_coordinate_scalar_operators(grid::MultiGrid{T, N}) where {T, N}
#     Tuple(ScalarOperator(x->x[i], grid) for i in 1:N)
# end

# function get_coordinate_scalar_operators(grid::Grid{T}) where {T}
#     ScalarOperator(identity, grid)
# end

end # module CollocationGrids
