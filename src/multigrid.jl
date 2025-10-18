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

function cardinal(mg::MultiGrid{T, N}, J, x) where {T, N}
    mapreduce(cardinal, *, mg.grids, Tuple(J), x, init=one(T))
end

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

