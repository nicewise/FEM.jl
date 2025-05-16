struct constitutive_linear_elastic{T<:AbstractDimension} <: AbstractConstitutiveLaw{T}
    E::Float64
    ν::Float64
    D::Matrix{Float64}

    function constitutive_linear_elastic{T}(E, ν) where T
        D = 0
        if T <: Dim3
            D = [1-ν ν   ν   0    0    0
                 ν   1-ν ν   0    0    0
                 ν   ν   1-ν 0    0    0
                 0   0   0   1-2ν 0    0
                 0   0   0   0    1-2ν 0
                 0   0   0   0    0    1-2ν]*E/(1+ν)/(1-2ν)
        elseif T <: PlaneStrain
            D = [1-ν ν   0
                 ν   1-ν 0
                 0   0   1-2ν]*E/(1+ν)/(1-2ν)
        else
            D = [1 ν 0
                 ν 1 0
                 0 0 (1-ν)/2]*E/(1-ν^2) # FIXME
        end
        return new{T}(E, ν, D)
    end
end

function constitutive_law_apply!(F::constitutive_linear_elastic, p::AbstractPoint)
    ε = p.ε + p.dε
    p.σ .= F.D * ε
    p.D .= F.D
    return nothing
end

struct vm_iso{T<:Dim3} <: AbstractConstitutiveLaw{T}
    E::Float64
    ν::Float64
    σy0::Float64
    H::Float64
    tol::Float64
    G::Float64
    G3::Float64
    D::Matrix{Float64}

    function vm_iso{T}(E, ν, σy0, H) where T
        k = E/3/(1-2ν) # 体积模量
        G = E/2/(1+ν)  # 剪切模量
        D = 3k * Mandel.J + 2G * Mandel.K
        return new{T}(E, ν, σy0, H, σy0*1e-6, G, 3G, D)
    end
end

function constitutive_law_apply!(F::vm_iso{T}, p::AbstractPoint{<:AbstractCellType, T}) where T<:Dim3
    εp = p.statev[1:6] # plastic strain
    γp = p.statev[7] # effective plastic strain

    σ = F.D * (p.ε + p.dε - εp) # trial stress
    s = Mandel.K * σ
    J2 = sum(σ.^2)
    q = sqrt(3/2) * J2
    f = q - (F.σy0 + F.H*γp)
    function yield!()
        Δγ = f / (F.G3 + F.H)
        N = s / J2 # 6x1
        Nb = sqrt(2/3) * N # \bar{N}
        σ -= F.G3 * Δγ * Nb
        εp += Δγ * N
        γp += Δγ
        c1 = 6*F.G^2/(F.G3+H)
        c2 = 6*G^2*Δγ/q
        D = F.D - (c1-c2)*Nb*Nb' - c2*Mandel.K
        #D = F.D - c1 * Nb*Nb'
        return D
    end
    D = f >= F.tol ? yield!() : F.D
    p.statev[1:6] .= εp
    p.statev[7] = γp
    p.D .= D
    p.σ .= σ
    return nothing
end

abstract type AbstractRd end

struct R1 <: AbstractRd
    rc::Float64
    dc::Float64
end

struct R2 <: AbstractRd
    n::Real
    rc::Float64
    dc::Float64
end

struct R3 <: AbstractRd
    rc::Float64
    dc::Float64
end

struct R4 <: AbstractRd
    # 线性
    c0::Float64
    c1::Float64
end

const RD_TYPES = Dict(
    1 => (R1, 2),
    2 => (R2, 3),
    3 => (R3, 2),
    4 => (R4, 2)
)

function Rd_gen(filter::Int, p::Real...)
    # 检查filter有效性
    if !haskey(RD_TYPES, filter)
        throw(ArgumentError("Invalid filter: $filter. Valid options: $(keys(RD_TYPES))"))
    end

    # 获取目标类型及参数数量
    type_def, n_params = RD_TYPES[filter]

    # 校验参数数量
    if length(p) != n_params
        throw(ArgumentError("Filter $filter requires $n_params parameters, got $(length(p))"))
    end

    # 转换类型并构造对象
    return type_def([Float64(x) for x in p]...)
end

function _Rd(R::R1, d)
    @. 4 * R.dc * R.rc * d / (d + R.dc)^2
end

function _Rd(R::R2, d)
    n, rc, dc = R.n, R.rc, R.dc
    @. (d * n * rc) / (dc * (n - 1 + (d/dc)^n))
end

function _Rd(R::R3, d)
    @. d * R.rc * exp(1 - d/R.dc) / R.dc
end

function _Rd(R::R4, d)
    @. R.c0 + R.c1 * d
end

# 标量版本
Rd(R::AbstractRd, d::Real) = _Rd(R, d)

# 带权重的数组版本
function Rd(R::AbstractRd, d::AbstractArray, ω::AbstractArray)
    @assert length(d) == length(ω) "维数不匹配"
    return _Rd(R, d) .* ω
end

# === 每种类型的导数 ===
#TODO chatgpt生成代码，未检查导数是否正确
function _dRd(R::R1, d)
    rc, dc = R.rc, R.dc
    @. 4 * dc * rc * (dc - d) / (d + dc)^3
end

function _dRd(R::R2, d)
    n, rc, dc = R.n, R.rc, R.dc
    @. begin
        a = (d / dc)^n
        num = n * rc * (n - 1 + a - n * a)
        den = dc * (n - 1 + a)^2
        num / den
    end
end

function _dRd(R::R3, d)
    rc, dc = R.rc, R.dc
    @. rc * exp(1 - d/dc) * (1 - d/dc) / dc
end

function _dRd(R::R4, d)
    @. R.c1
end

# 标量输入
dRd(R::AbstractRd, d::Real) = _dRd(R, d)

# 向量输入，必须带权重，输出对角矩阵
function dRd(R::AbstractRd, d::AbstractVector, ω::AbstractVector)
    @assert length(d) == length(ω) "维数不匹配"
    v = _dRd(R, d) .* ω
    return Diagonal(v)
end

# 捕捉非法调用（向量输入但无权重）
function dRd(R::AbstractRd, d::AbstractVector)
    throw(ArgumentError("向量输入时必须提供权重向量 ω"))
end


#FIXME delute跑不出正确的本构曲线
# 对于pcw，体积模量和剪切模量不劣化到小于0已经解决了问题
# 但是delute不知道问题在哪
struct elastic_damage{T<:Dim3} <: AbstractConstitutiveLaw{T}
    l::Float64 # length scale
    x::Int
    y::Int
    k3::Float64
    μ2::Float64
    weaken_ratio::Dict{Bool, NTuple{4, Function}}
    isopen::Function
    Rd::AbstractRd
    tol::Float64

    function elastic_damage{T}(S::Symbol, E::Float64, ν::Float64, Rd::AbstractRd; tol::Float64 = 1e-6, l::Float64 = 0.) where T
        β1 = 16(1 - ν^2) / 9(1 - 2ν)
        β2 = 32(1 - ν) / 15(2 - ν)
        ϑ = (5 - ν) / 3
        weaken_ratio = Dict(
            :DELUTE => Dict(
                true => (
                    d -> max(1e-10, 1 - β1 * d),
                    d -> - β1,

                    d -> max(1e-10, 1 - β2 * ϑ * d),
                    d -> - β2 * ϑ
                ),
                false => (
                    d -> 1,
                    d -> 0,

                    d -> max(1e-10, 1 - β2 * d),
                    d -> - β2
                )
            ),
            :MT => Dict(
                true => (
                    d -> 1 / (1 + β1 * d),
                    d -> - β1 / (1 + β1 * d)^2,

                    d -> 1 / (1 + β2 * ϑ * d),
                    d -> - β2 * ϑ / (1 + β2 * ϑ * d)^2
                ),
                false => (
                    d -> 1,
                    d -> 0,

                    d -> 1 / (1 + β2 * d),
                    d -> - β2 / (1 + β2 * d)^2
                )
            ),
            :PCW => Dict(
                true => (
                    d -> max(1e-10, 1 - 16(1 - ν^2)d / ( 9(1 - 2ν) + 16/3*(1 + ν)^2 * d )),
                    d -> 1296(1 - 2ν)*(ν^2 - 1)/(16d*(ν + 1)^2 - 54ν + 27)^2,

                    d -> max(1e-10, 1 - 480(1 - ν)ϑ * d / ( 225(2 - ν) + 64(4 - 5ν)ϑ * d )),
                    d -> 324000(ν - 5)*(ν - 2)*(ν - 1)/(64d*(ν - 5)*(5ν - 4) - 675ν + 1350)^2
                ),
                false => (
                    d -> 1,
                    d -> 0,

                    d -> max(1e-10, 1 - 480(1 - ν)d / ( 225(2 - ν) + 64(4 - 5ν)d )),
                    d -> -108000(ν - 2)*(ν - 1)/(320d*ν - 256d + 225ν - 450)^2
                )
            )
        )

        k3 = E / (1 - 2ν)
        μ2 = E / (1 + ν)
        D = k3 * Mandel.J + μ2 * Mandel.K
        isopen = σ -> σ' * Mandel.δ / 3 >= 0
        return new{T}(l, 1, 2, k3, μ2, weaken_ratio[S], isopen, Rd, tol)
    end
end

function constitutive_law_apply!(F::elastic_damage{T}, p::AbstractPoint{<:AbstractCellType, T}; before_gradient::Union{Nothing, Bool}=nothing) where T<:Dim3
    if isnothing(before_gradient) || before_gradient
        d = p.statev[1]
        ε = p.ε + p.dε
        _, dkd, _, dμd = F.weaken_ratio[F.isopen(p.D * ε)]
        k3 = F.k3
        μ2 = F.μ2
        g(x) = -ε' * (dkd(d+x) * k3 * Mandel.J + dμd(d+x) * μ2 * Mandel.K) * ε / 2 - Rd(F.Rd, d+x)
        if g(0.) > F.tol
            d += Roots.find_zero(g, 0.0)
        end
        p.statev[1] = d
    end
    if isnothing(before_gradient) || !before_gradient
        d = isnothing(before_gradient) ? p.statev[1] : p.statev[2]
        ε = p.ε + p.dε
        kd, _, μd, _ = F.weaken_ratio[F.isopen(p.D * ε)]
        Chom = kd(d) * k3 * Mandel.J + μd(d) * μ2 * Mandel.K
        p.D .= Chom
        p.σ .= Chom * ε
    end
    return nothing
end



struct PDMT end
struct PDPCW end

struct MicroCrack
    direction::SVector{3}
    ω::Float64 # weight
    N::SVector{6, Float64} # 法向方向张量
    E1::SMatrix{6, 6, Float64} # Walpole 张量基
    E2::SMatrix{6, 6, Float64}
    E3::SMatrix{6, 6, Float64}
    E4::SMatrix{6, 6, Float64}
    E5::SMatrix{6, 6, Float64}
    E6::SMatrix{6, 6, Float64}
    function MicroCrack(n::AbstractVector{<:Real}, w::Float64)
        @assert length(n) == 3 "Input vector 'n' must have length 3, got $(length(n))"
        return MicroCrack(SVector{3}(n), w)
    end
    function MicroCrack(n::SVector{3}, w::Float64)
        N = n ⊗ n
        T = Mandel.δ - N
        return new(
            n,
            w,
            N,
            T ⊗ T / 2,
            N ⊗ N,
            T ⋆ T - T ⊗ T / 2, # ⋆ = ⊗s
            N ⋆ T + T ⋆ N,
            N ⊗ T,
            T ⊗ N
        )
    end
end
compute_walpole(c::MicroCrack, x) =
    x[1] * c.E1 + x[2] * c.E2 + x[3] * c.E3 + x[4] * c.E4 + x[5] * c.E5 + x[6] * c.E6
@inline function Sn(crack::MicroCrack, key::Bool, cn, ct) # cn和ct是书P60页所述cn、ct之倒数，为了少做除法
    if key
        return cn * crack.E2 + ct * crack.E4
    else
        return ct * crack.E4
    end
end
struct AEDDELUTE{T<:Dim3} <: AbstractConstitutiveLaw{T}
    C0::SMatrix{6, 6, Float64, 36}
    cn::Float64
    ct::Float64
    n::Int # 代表性裂隙族数量，ie. 球面积分高斯点数目
    cracks::Vector{MicroCrack}
    Chom::Function
    isopen::Function
    Rd::AbstractRd
    tol::Float64
    function AEDDELUTE{T}(E::Float64, ν::Float64, Rd::AbstractRd; tol::Float64 = 1e-6, n::Int) where T
        known_quadrature = [21, 33, 48, 96, 99, 144, 198]
        @assert n ∈ known_quadrature "只知道($i for i in known_quadrature), 不知道$n 积分点怎么办."
        dnw = sphere_quad(n) # direction and weight
        cracks = MicroCrack[]
        for i in axes(dnw, 1)
            push!(cracks, MicroCrack(dnw[i, 1:3], dnw[i, 4]))
        end
        cn = 16(1-ν^2) / 3E
        ct = (2-ν) * cn
        k3 = E / (1 - 2ν)
        μ2 = E / (1 + ν)
        C0 = k3 * Mandel.J + μ2 * Mandel.K
        κe = Dict(
            true => d -> d < 3(1-2ν)/16(1-ν)^2 ?
                (-16*d*(ν - 1)^2 - 6*ν + 3)/(32*d*ν^2*(ν - 1) - 6*ν + 3) : 0,
            false => d -> 1
        )
        κμ(d) = d < 3(2-ν)/16(1-ν) ? 
            1 - 16(1-ν) / 3(2-ν) * d : 0
        Chom(d::Real, key::Bool) = begin
            κ = κe[key]
            r = 1 - ν - 2κ(d) * ν^2
            return [
                E / r
                E * κ(d) * (1-ν) / r
                E / (1 + ν)
                2κμ(d)
                E * κ(d) * ν / r
                E * κ(d) * ν / r
            ]
        end
        isopen = σ -> [σ' * crack.N for crack in cracks] .>= 0 # σ -> vector{Bool}
        return new{T}(C0, cn, ct, n, cracks, Chom, isopen, Rd, tol)
    end
end
Fd(F::AEDDELUTE{Dim3}, ε::Vector{Float64}, ckeys::Vector{Bool}) =
    # keys to determine the specific micro crack is open or not. ckeys = crack_keys
    [- ε' * F.C0 * Sn(F.cracks[i], key, F.cn, F.ct) * F.C0 * ε / 2 for (i, key) in enumerate(ckeys)]

function constitutive_law_apply!(F::AEDDELUTE{T}, p::AbstractPoint{<:AbstractCellType, T}) where T<:Dim3
    ε = p.ε + p.dε
    d = p.statev[1:F.n]
    stats = F.isopen(p.D * ε) # Vector{Bool} indicate is microcrack open or not
    Fd_const = Fd(F, ε, stats)
    g = [x->Fd_const - Rd.(F.Rd, d+x) for d in d]
    for (i, f) in enumerate(g)
        if f(0.) > F.tol
            d[i] += Roots.find_zero(f, 0.0)
        end
    end
    p.statev[1:F.n] .= d
    Chom = sum(compute_walpole(F.cracks[i], F.Chom(d[i], stats[i])) * F.cracks[i].ω for i in 1:F.n)
    p.D .= Chom
    p.σ .= Chom * ε
    return nothing
end

struct LCP
    n::Int
    model
    function LCP(n) # n denotes the number of variables
        model = Model(PATHSolver.Optimizer)
        #set_optimizer_attribute(model, "output", "no")
        set_optimizer_attribute(model, "silent", true)
        # without this option set to "true"
        # PATHSolver will print a lot of information to stdout
        return new(n, model)
    end
end
function LCP_solve(lcp::LCP, M, q)
    @variable(lcp.model, x[1:n] >= 0)
    @constraint(lcp.model, M * x .+ q ⟂ x)
    optimize!(lcp.model)
    return(value.(x)) # obtain solution
end

struct AEDPCW{T<:Dim3} <: AbstractConstitutiveLaw{T}
    l::Float64
    x::Vector{Int}
    y::Vector{Int}
    n::Int # 代表性裂隙族数量，ie. 球面积分高斯点数目
    ω::Vector{Float64}
    cracks::Vector{MicroCrack}
    C0::SMatrix{6, 6, Float64}
    Tr_gen::Function
    Pd::SMatrix{6, 6, Float64}
    isopen::Function
    Rd::AbstractRd
    tol::Float64
    lcp::LCP
    function AEDPCW{T}(E::Float64, ν::Float64, Rd::AbstractRd; tol::Float64 = 1e-6, l::Float64 = 0., n::Int) where T
        x = collect(1:n)
        y = collect(n+1:2n)
        dnw = sphere_quad(n) # direction and weight
        cracks = MicroCrack[]
        for i in axes(dnw, 1)
            push!(cracks, MicroCrack(dnw[i, 1:3], dnw[i, 4]))
        end
        cn = 16(1-ν^2) / 3E
        ct = (2-ν) * cn
        k3 = E / (1 - 2ν)
        μ2 = E / (1 + ν)
        C0 = k3 * Mandel.J + μ2 * Mandel.K
        Tr_gen(crack::MicroCrack, key::Bool) = C0 * Sn(crack, key, cn, ct) * C0 * crack.ω
        αk = 1 / (k3 + 2μ2)
        αμ = (2k3 + 6μ2) / 5(k3 + 2μ2) / μ2
        Pd = αk * Mandel.J + αμ * Mandel.K
        isopen = σ -> [σ' * crack.N for crack in cracks] .>= 0 # σ -> vector{Bool}
        lcp = LCP(n)
        return new{T}(l, x, y, n, dnw[:, 4], cracks, C0, Tr_gen, Pd, isopen, Rd, tol, lcp)
    end
end

function constitutive_law_apply!(F::AEDPCW{T}, p::AbstractPoint{<:AbstractCellType, T}; before_gradient::Union{Nothing, Bool}=nothing) where T<:Dim3
    if isnothing(before_gradient) || before_gradient
        d_old = p.statev[F.x]
        ε = p.ε + p.dε
        d(x) = d_old + x # vector -> vector
        stats = F.isopen(p.D * ε) # Vector{Bool} indicate is microcrack open or not
        Tr = [F.Tr_gen(crack, stats[i]) for (i, crack) in enumerate(F.cracks)]
        Cd(x) = sum(d(x) .* Tr) # vector -> matrix
        B(x) = (Mandel.I + F.Pd * Cd(x)) \ F.Pd # vector -> matrix
        BCd2(x) = 2 * B(x) * Cd(x) - Mandel.I # vector -> matrix
        Fd(x) = [-ε'*Tr[i]*BCd2(x)*ε/2 for i in 1:F.n]
        g(x) = Fd(x) - Rd(F.Rd, d(x)) # vector -> vector

        CdB(x) = Mandel.I - Cd(x) * B(x)
        BCd(x) = Cd(x) * B(x) - Mandel.I
        function dFd(x)
            o = zeros(F.n, F.n)
            for i in 1:F.n
                o[i, i] = dFd_val(Tr[i], Tr[i]) / 2
                for j in i+1:n
                    o[i, j] = dFd_val(Tr[j], Tr[i])
                end
            end
            return (o + o')
        end
        dFd_val(x, Tj, Ti) = -ε'*CdB(x)*(Tj*B(x)*Ti)*BCd(x)*ε #TODO

        g∂d = dFd(d_old) - dRd(F.Rd, d_old)
        g_penetrator = g(zeros(F.n))
        if sum( g_penetrator .> F.tol ) >= 1
            #Δd = NCP_desent(g, F.n)#optim(g) # TODO
            Δd = LCP_solve(F.lcp, g∂d, g_penetrator)
            d_old = d(Δd)
        end
        p.statev[F.x] .= d_old
    end
    if isnothing(before_gradient) || !before_gradient
        d = isnothing(before_gradient) ? p.statev[F.x] : p.statev[F.y]
        d = p.statev[F.x]
        ε = p.ε + p.dε
        stats = F.isopen(p.D * ε) # Vector{Bool} indicate is microcrack open or not
        Tr = [F.Tr_gen(crack, stats[i]) for (i, crack) in enumerate(F.cracks)]
        cd = sum(d .* Tr)
        b = (Mandel.I + F.Pd * cd) \ F.Pd
        Chom = F.C0 - cd + cd * b * cd
        p.D .= Chom
        p.σ .= Chom * ε
    end
    return nothing
end


struct SEDPCW end


struct APDPCW end

struct TEDMT{T<:Dim3} <: AbstractConstitutiveLaw{T}
    T0::Float64
    Chom::Function
    κhom::Function
    Fd::Function
    Rd::AbstractRd
    tol::Float64
    function TEDMT{T}(E::Float64, ν::Float64, α::Float64, T0::Real, Rd::AbstractRd; tol::Float64 = 1e-6) where T
        k3 = E / (1 - 2ν)
        μ2 = E / (1 + ν)
        κ = k3 * α
        η1 = 16(1-ν^2) / 9(1-2ν)
        η2 = 32(1-ν) * (5-ν) / 45(2-ν)
        Chom(d) = k3 / (1 + η1*d) * Mandel.J + μ2  / (1 + η2 * d) * Mandel.K
        κhom(d) = κ / (1 + η1*d) * Mandel.δ
        Ta(t) = (t - T0) / log(t / T0)
        Sd = η1 / k3 * Mandel.J + η2 / μ2 * Mandel.K
        Fd(σ, d, t) = begin
            ΔT = t - T0
            σ' * Sd * σ / 2 - 3η1 * α * κ * ΔT / (1 + η1 * d)^2 * (ΔT/2 + Ta(t) - t)
        end
        return new{T}(
            T0, Chom, κhom, Fd, Rd, tol
        )
    end
end
function constitutive_law_apply!(F::TEDMT{<:Dim3}, p::AbstractPoint{<:AbstractCellType, <:Dim3})
    d, T = p.statev[1:2]
    ΔT = T - F.T0
    ε = p.ε + p.dε
    σ(d) = F.Chom(d) * ε - F.κhom(d) * ΔT
    g(x) = F.Fd(σ(d+x), d+x, T) - Rd(F.Rd, d+x)
    if g(0.) > F.tol
        d += Roots.find_zero(g, 0.0)
    end
    p.statev[1] = d
    p.D .= F.Chom(d)
    p.σ .= σ(d)
end
