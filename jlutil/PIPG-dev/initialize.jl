function initialize(prb,flg::Symbol)
    # prb is a dictionary of problem data

    nx = prb.nx
    nu = prb.nu
    N = prb.N    
    nK0 = [size(prb.F0[t],1) for _ in 1:N]  # Number of 0 cones at each time
    nK1 = [size(prb.F1[t],1) for _ in 1:N]  # Number of ℝ₊ cones at each time

    nz = N*(nx+nu)
    nw = (N-1)*nx + sum(nK0)+sum(nK1)

    struct soldat
        x::Vector{MVector{nx,Float64}} 
        u::Vector{MVector{nu,Float64}}
        ϕ::Vector{MVector{nx,Float64}}
        θ::Vector{Mvector}
        ψ::Vector{Mvector}     

        x̃::Vector{MVector{nx,Float64}} 
        ũ::Vector{MVector{nu,Float64}}
        ϕ̃::Vector{MVector{nx,Float64}}
        θ̃::Vector{Mvector}
        ψ̃::Vector{Mvector}         
    end
    # Sizes of θ, θ̃, ψ, ψ̃ are not specified because the number of cones at two different time instants could be different

    struct soldat_vec
        z::Vector{Float64}
        w::Vector{Float64}
        ξ::Vector{Float64}
        η::Vector{Float64}
    end        

    if flg == :devec    
        sol = soldat([MVector{nx}(randn(nx)) for _ in 1:N],
                     [MVector{nu}(randn(nu)) for _ in 1:N],
                     [MVector{nx}(randn(nx)) for _ in 1:N-1],
                     [MVector{nK0[t]}(randn(nK0[t])) for t in 1:N],
                     [MVector{nK1[t]}(randn(nK1[t])) for t in 1:N],
                     [MVector{nx}(randn(nx)) for _ in 1:N],
                     [MVector{nu}(randn(nu)) for _ in 1:N],
                     [MVector{nx}(randn(nx)) for _ in 1:N-1],
                     [MVector{nK0[t]}(randn(nK0[t])) for t in 1:N],
                     [MVector{nK1[t]}(randn(nK1[t])) for t in 1:N])
    elseif flg == vec:
        sol = soldat_vec(randn(nz),randn(nw),randn(nz),randn(nw))
    else 
        error("Second argument of initialize must be either :devec or :vec")
    end

    return sol 

end