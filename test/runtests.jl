using Test

using CollocationGrids

const TEST_ATOL = 1e-4
const TEST_RTOL = 1e-4
const GRID_TEST_N = 10

function everyn2len(iter, len::Integer = GRID_TEST_N)
    @assert len > 0
    every = fld(length(iter), len)
    
    filtered_pairs = Iterators.filter(zip(Iterators.countfrom(0), iter)) do (i, _)
        iszero(mod(i, every))
    end

    return Iterators.map(last, filtered_pairs)
end

function value_test(test_grid, tf)
    val_pnts = Iterators.map(eachgridpoint(test_grid.grid)) do pnt 
        return pnt, tf.f(pnt)
    end |> everyn2len

    @testset "Testng value at $(pnt)" for (pnt, fpnt) in val_pnts
        @test isapprox(fpnt, tf.f(pnt), atol = TEST_ATOL, rtol = TEST_RTOL)
    end
end

function derivative_test(test_grid, tf)
    pnts = gridpoints(test_grid.grid)
    dvals = derivative_matrix(test_grid.grid, 1) * tf.f.(pnts)
    dval_pnts = everyn2len(zip(pnts, dvals))

    @testset "Testng value at $(pnt)" for (pnt, dfpnt) in dval_pnts
        @test isapprox(dfpnt, tf.df(pnt), atol = TEST_ATOL, rtol = TEST_RTOL)
    end
end

struct TestFunction
    label
    f
    df
end

struct TestGrid
    label
    grid
end

const TWO_PI_PERIODIC_TEST_GRIDS = [
              TestGrid("Fourier T = Float64 with order 20", FourierGrid{Float64}(20))
              TestGrid("Fourier T = Float64 with order 10", FourierGrid{Float64}(10))
             ]

const TWO_PI_PERIODIC_TEST_FUNCTIONS = [
                        TestFunction("5", x -> one(x)*5, x -> zero(x))
                        TestFunction("sin(x)", x -> sin(x), x -> cos(x))
                        TestFunction("cos(x)^2", x -> cos(x)^2, x -> - 2 *cos(x) * sin(x))
                       ]

const INTERVAL_TEST_FUNCTIONS = [
                                 TWO_PI_PERIODIC_TEST_FUNCTIONS
                                 TestFunction("x", x -> x, x -> one(x))
                                 TestFunction("x^2", x -> x ^ 2, x -> 2x)
                                ]

const INTERVAL_TEST_GRIDS = [
                             TestGrid(
                                      "Chebyshev-Lobatto Grid T = Float64 with order 20",
                                      ChebyshevLobattoGrid{Float64}(20)
                                     )
                             TestGrid(
                                      "Chebyshev-Lobatto Grid T = Float64 with order 30",
                                      ChebyshevLobattoGrid{Float64}(30)
                                     )
                            ]

append!(
        INTERVAL_TEST_GRIDS,
        map(INTERVAL_TEST_GRIDS) do tg
            grid = linear_rescale(tg.grid, 0, 5)
            label = "with linear rescale from $(boundaries(tg.grid)) to (0, 5)"
            TestGrid(label, grid)
        end
       )

@testset "$(test_grid.label)" for test_grid in TWO_PI_PERIODIC_TEST_GRIDS
    @testset "Testing function: $(test_function.label)" for test_function in TWO_PI_PERIODIC_TEST_FUNCTIONS
        @testset "Value Test" begin
            value_test(test_grid, test_function)
        end

        @testset "Derivative Test" begin
            derivative_test(test_grid, test_function)
        end
    end
end

@testset "$(test_grid.label)" for test_grid in INTERVAL_TEST_GRIDS 
    @testset "Testing function: $(test_function.label)" for test_function in INTERVAL_TEST_FUNCTIONS
        @testset "Value Test" begin
            value_test(test_grid, test_function)
        end

        @testset "Derivative Test" begin
            derivative_test(test_grid, test_function)
        end
    end
end
