# CollocationGrids.jl

A Julia package for numerical computation using collocation grids and spectral methods.

## Overview

CollocationGrids.jl provides a comprehensive framework for working with various types of collocation grids commonly used in numerical analysis and spectral methods. The package supports multiple grid types, derivative calculations, interpolation, and multi-dimensional operations.

## Features

- **Multiple Grid Types**: Support for Chebyshev, Fourier, semi-infinite, and custom grids
- **Multi-dimensional Grids**: Tensor product grids for multi-dimensional problems
- **Spectral Differentiation**: Accurate computation of derivatives using spectral methods
- **Interpolation**: Cardinal function-based interpolation
- **Grid Transformations**: Linear rescaling and coordinate transformations

## Installation

```julia
using Pkg
Pkg.add(url = "https://github.com/inokawazu/CollocationGrids")
```

## Quick Start

```julia
using CollocationGrids

# Create a Chebyshev-Lobatto grid with 10 points
grid = ChebyshevLobattoGrid{Float64}(10)

# Get grid points
points = collect(gridpoints(grid))

# Compute first derivative matrix
D1 = derivative_matrix(grid, 1)

# Interpolate a function
f(x) = sin(x)
coefs = f.(points)
interpolated_value = interpolate(grid, coefs, 0.5)
```

## Grid Types

### Single-Dimension Grids

- **`ChebyshevLobattoGrid{T}(N)`**: Chebyshev-Lobatto grid with N+1 points including boundaries
- **`ChebyshevInteriorGrid{T}(N)`**: Chebyshev grid with N interior points (excludes boundaries)
- **`FourierGrid{T}(N)`**: Fourier grid with 2N points for periodic problems
- **`SemiInfiniteGrid`**: Grid for semi-infinite domains
- **`LinearRescaleGrid`**: Linearly rescaled version of any grid

### Multi-Dimensional Grids

```julia
# Create a 2D grid (tensor product)
grid1 = ChebyshevLobattoGrid{Float64}(10)
grid2 = FourierGrid{Float64}(8)
multigrid = MultiGrid(grid1, grid2)

# Get 2D derivative matrices
Dx = derivative_matrix(multigrid, (1, 0))  # ∂/∂x
Dy = derivative_matrix(multigrid, (0, 1))  # ∂/∂y
Dxy = derivative_matrix(multigrid, (1, 1)) # ∂²/∂x∂y
```

## Core Functions

### Grid Operations

```julia
# Get grid points
points = gridpoints(grid)

# Get domain boundaries
bounds = boundaries(grid)

# Get specific grid point
point_i = gridpoint(grid, i)
```

### Differentiation

```julia
# Derivative matrix for order-k derivatives
Dk = derivative_matrix(grid, k)

# Specific derivative element
d_ij = derivative(grid, i, j, order)

# For multi-dimensional grids
derivative_matrices = derivative_array_matrix(multigrid, max_order)
```

### Interpolation

```julia
# Cardinal function at point x for basis function j
card_val = cardinal(grid, j, x)

# Interpolate using coefficients
result = interpolate(grid, coefficients, x)
```

## Grid Construction Examples

### Chebyshev Grids

```julia
# Lobatto grid (includes boundaries at ±1)
clob = ChebyshevLobattoGrid{Float64}(16)

# Interior grid (excludes boundaries)
cint = ChebyshevInteriorGrid{Float64}(16)
```

### Fourier Grids

```julia
# For periodic functions on [0, 2π]
fourier = FourierGrid{Float64}(32)
```

### Custom Domain Grids

```julia
# Rescale Chebyshev grid to [0, 5]
base_grid = ChebyshevLobattoGrid{Float64}(20)
custom_grid = LinearRescaleGrid(base_grid, 0.0, 5.0)

# Or use the convenience function
rescaler = linear_grid_rescale(0.0, 5.0)
custom_grid = rescaler(base_grid)
```

### Semi-Infinite Domains

```julia
# Grid for problems on [0, ∞)
base = ChebyshevLobattoGrid{Float64}(25)
semi = SemiInfiniteGrid(base, 1.0)  # L = 1.0 is the mapping parameter
```

## Advanced Usage

### Coordinate Operators

```julia
# Get coordinate operators for a grid
X = coordinate_operators(grid)           # Returns diagonal matrix
x_vec = coordinate_vectors(grid)         # Returns vector

# For multi-dimensional grids
X, Y = coordinate_operators(multigrid)   # Tuple of operators
```

### Function Operators

```julia
# Create operator from function
f(x) = x^2 + sin(x)
F_op = function_operator(f, grid)        # Diagonal matrix
f_vec = function_vector(f, grid)         # Vector evaluation
```

### Multi-Grid Derivative Arrays

```julia
# Get all derivative matrices up to given order
max_order = 3
all_derivatives = derivative_array_matrix(multigrid, max_order)
```

## Type Parameters

All grids are parameterized by element type `T` (typically `Float64` or `BigFloat`):

```julia
# High precision computation
grid_hp = ChebyshevLobattoGrid{BigFloat}(50)

# Standard precision
grid_std = ChebyshevLobattoGrid{Float64}(50)
```
