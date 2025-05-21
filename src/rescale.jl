"""
linear_rescale(x, a, b, c, d)

Linearly rescales x from [a,b] -> [c,d]
"""
linear_rescale(x, a, b, c, d) = (d - c) * (x - a) / (b - a) + c

"""
linear_rescale_factor(a, b, c, d)

returns the slope of the transformation, `(d - c) / (b - a)`.
"""
linear_rescale_factor(a, b, c, d) = (d - c) / (b - a) # df/dx
