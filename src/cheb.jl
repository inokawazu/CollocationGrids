function chebT(n::Integer, x::T) where T
    n<0 && error("n=$n must be nonnegative.")
    # t = acos(x)
    # return cos(t*n)

    n == 0 && return one(T)
    n == 1 && return x
    
    v0, v1 = one(T), x
    
    while n > 1
        v0, v1 = v1, 2*x*v1-v0
        n -= 1
    end

    return v1
end

function chebU(n::Integer, x::T) where T
    n<0 && error("n=$n must be nonnegative.")
    # t = acos(x)
    # return sin((n+1)*t)/sin(t)

    n == 0 && return one(T)
    n == 1 && return 2x
    
    v0, v1 = one(T), x
    
    while n > 1
        v0, v1 = v1, 2*x*v1-v0
        n -= 1
    end

    return v1
end
