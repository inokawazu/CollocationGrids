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

