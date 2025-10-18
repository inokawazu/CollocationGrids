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

