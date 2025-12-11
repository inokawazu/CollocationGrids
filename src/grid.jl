gridpoints(g::Grid) = collect(eachgridpoint(g))
eachgridpoint(g::Grid{T}) where T = (gridpoint(g, i) for i in eachindex(g))

derivative_matrix(g::Grid{T}, order) where T = [
                                                derivative(g, i, j, order) 
                                                for i in eachindex(g), j in eachindex(g)
                                               ]

function interpolate(g::Grid{T}, coefs, x) where T 
    sum(coef * cardinal(g, j, x) for (j, coef) in zip(eachindex(g), coefs))
end

Base.length(g::Grid) = length(eachindex(g))
Base.size(g::Grid) = size(eachindex(g))
Base.eltype(g::Grid{T}) where T = T

boundaries(g::Grid) = sort(extrema(gridpoints(g)))
