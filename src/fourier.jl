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

