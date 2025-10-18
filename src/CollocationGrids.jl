module CollocationGrids

export gridpoints, 
       eachgridpoint,
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

export coordinate_operators, 
       coordinate_vectors,
       function_operator,
       function_vector

import TaylorDiff, LinearAlgebra
# TODO: replace finding zeros with more rebust routine
using Roots: find_zeros
using LinearAlgebra
using SparseArrays

abstract type Grid{T} end

# File for grids
include("grid.jl")
include("linear_rescale.jl")
include("cheb.jl")
include("one_hots.jl")
include("multigrid.jl")

# TODO: Not quite sure needed in library, so figure out if so
# function apply_derivative_matrix(input::AbstractArray, der_mat::AbstractArray, dim::Integer)
#     return apply_derivative_matrix!(similar(input), input, der_mat, dim) 
# end

# function apply_derivative_matrix!(
#         out::AbstractArray, input::AbstractArray,
#         der_mat::AbstractArray, dim::Integer)

#     fill!(out, 0)

#     for (out_i, i) in zip(eachslice(out, dims=dim), axes(der_mat, 1))
#         for (in_j, j) in zip(eachslice(input, dims=dim), axes(der_mat, 2))
#             @. out_i += der_mat[i, j] * in_j
#         end
#     end

#     return out
# end

include("fourier.jl")
include("interior_grid.jl")
include("cheb_lobbo.jl")
include("recombination_grid.jl")
# TODO ChebyshevSpectralGrid
include("semi_infinite.jl")
include("operators.jl")

end # module CollocationGrids
