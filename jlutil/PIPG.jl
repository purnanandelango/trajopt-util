module PIPG
using LinearAlgebra, StaticArrays

function initialize(prb)

    struct solver_data
        x = [MVector{prb.nx}(randn(prb.nx)) for _ in 1:prb.N]
        u = [MVector{prb.nu}(randn(prb.nu)) for _ in 1:prb.N]         
    end

end

end