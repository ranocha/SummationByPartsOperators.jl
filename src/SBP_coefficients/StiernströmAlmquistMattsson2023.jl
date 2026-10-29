"""
    StiernströmAlmquistMattsson2023()

Coefficients of the boundary-optimized SBP operators for second derivatives
with variable coefficients given in
- Stiernström, Almquist, Mattsson (2023)
  Boundary-optimized summation-by-parts operators for finite difference
    approximations of second derivatives with variable coefficients.
  Journal of Computational Physics 491, 112376.

These operators live on the same nonuniform grid and use the same first
derivative operators and norm matrix as
[`MattssonAlmquistVanDerWeide2018Accurate`](@ref); the second derivative
operators are fully compatible with them, i.e.,
`D₂(b) = D₁ * Diagonal(b) * D₁ - H⁻¹ * R(b)` with a symmetric, positive
semidefinite remainder `R(b)`. Orders of accuracy 4, 6, 8, 10, and 12 are
available in the interior; the boundary closures are `order ÷ 2 - 1`-th order
accurate.

Use [`var_coef_derivative_operator`](@ref) to create the variable-coefficient
second-derivative operators `D₂(b)` and [`derivative_operator`](@ref) for the
first derivative `D₁` as well as for the constant-coefficient second derivative
`D₂(1)`.
"""
struct StiernströmAlmquistMattsson2023 <: SourceOfCoefficients end

function Base.show(io::IO, source::StiernströmAlmquistMattsson2023)
    if get(io, :compact, false)
        summary(io, source)
    else
        print(io,
              "Stiernström, Almquist, Mattsson (2023) \n",
              "  Boundary-optimized summation-by-parts operators for finite difference \n",
              "    approximations of second derivatives with variable coefficients. \n",
              "  Journal of Computational Physics 491, 112376. \n",
              "See also (first derivatives) \n",
              "  Mattsson, Almquist, van der Weide (2018) \n",
              "  Boundary optimized diagonal-norm SBP operators ('Accurate'). \n",
              "  Journal of Computational Physics 374, pp. 1261-1266.")
    end
end

function construct_grid(::StiernströmAlmquistMattsson2023, accuracy_order, xmin, xmax, N)
    construct_grid(MattssonAlmquistVanDerWeide2018Accurate(), accuracy_order,
                   xmin, xmax, N)
end

@inline function first_derivative_coefficients(::StiernströmAlmquistMattsson2023,
                                               order::Int, T = Float64, mode = FastMode())
    first_derivative_coefficients(MattssonAlmquistVanDerWeide2018Accurate(), order, T, mode)
end

function second_derivative_coefficients(source::StiernströmAlmquistMattsson2023,
                                        order::Int, T = Float64, mode = FastMode())
    if order == 4
        _sam2023_second_derivative_coefficients(source, Val(4), T, mode)
    elseif order == 6
        _sam2023_second_derivative_coefficients(source, Val(6), T, mode)
    elseif order == 8
        _sam2023_second_derivative_coefficients(source, Val(8), T, mode)
    elseif order == 10
        _sam2023_second_derivative_coefficients(source, Val(10), T, mode)
    elseif order == 12
        _sam2023_second_derivative_coefficients(source, Val(12), T, mode)
    else
        throw(ArgumentError("Order of accuracy $order not implemented/derived."))
    end
end

# The constant-coefficient second-derivative operator is `D₂(b)` evaluated for
# `b ≡ 1`, as in the reference implementation.
function _sam2023_second_derivative_coefficients(source::StiernströmAlmquistMattsson2023,
                                                 ::Val{order}, ::Type{T},
                                                 mode) where {order, T}
    p = order ÷ 2
    nrows = 3 * p
    ncols = 4 * p
    left, interior = _sam2023_coefficients(T, order)

    boundary = zeros(T, nrows, ncols)
    for i in 1:nrows, j in 1:ncols, k in axes(left, 1)
        boundary[i, j] += left[k, j, i]
    end
    # All rows of the boundary closure are padded with zeros to the same length
    # so that they can be stored in a homogeneous tuple.
    left_boundary = _sam2023_boundary_rows(T, boundary, Val(nrows), Val(ncols))
    right_boundary = left_boundary

    stencil = zeros(T, 2 * p + 1)
    for j in 1:(2 * p + 1), k in 1:(2 * p + 1)
        stencil[j] += interior[k, j]
    end
    central_coef = stencil[p + 1]
    upper_coef = SVector{p, T}(ntuple(s -> stencil[p + 1 + s], Val(p)))
    lower_coef = SVector{p, T}(ntuple(s -> stencil[p + 1 - s], Val(p)))

    coefficients_D1 = first_derivative_coefficients(MattssonAlmquistVanDerWeide2018Accurate(),
                                                    order, T, mode)
    # The norm weights are padded with ones accordingly.
    weights = coefficients_D1.left_weights
    left_weights = SVector{nrows, T}(ntuple(i -> i <= length(weights) ? weights[i] : one(T),
                                            Val(nrows)))
    right_weights = left_weights

    # The operators are fully compatible, so the first derivative at the
    # boundaries is given by the first and last row of `D₁`.
    left_boundary_derivatives = (coefficients_D1.left_boundary[1],)
    right_boundary_derivatives = (-coefficients_D1.left_boundary[1],)

    DerivativeCoefficients(left_boundary, right_boundary,
                           left_boundary_derivatives, right_boundary_derivatives,
                           lower_coef, central_coef, upper_coef,
                           left_weights, right_weights, mode, 2, order, source)
end

function _sam2023_boundary_rows(::Type{T}, boundary::AbstractMatrix, ::Val{nrows},
                                ::Val{ncols}) where {T, nrows, ncols}
    ntuple(i -> DerivativeCoefficientRow{T, 1, ncols}(SVector{ncols, T}(view(boundary, i,
                                                                             :))),
           Val(nrows))
end

function var_coef_derivative_coefficients(source::StiernströmAlmquistMattsson2023,
                                          derivative_order::Int, accuracy_order::Int, grid,
                                          mode = FastMode())
    @argcheck derivative_order == 2
    T = eltype(grid)
    coefficient_cache = StiernströmAlmquistMattsson2023Cache(T, accuracy_order)
    weights = first_derivative_coefficients(MattssonAlmquistVanDerWeide2018Accurate(),
                                            accuracy_order, T, mode).left_weights

    VarCoefDerivativeCoefficients(coefficient_cache, weights, weights,
                                  mode, derivative_order, accuracy_order, source)
end

"""
    StiernströmAlmquistMattsson2023Cache

Coefficient cache of the variable-coefficient second-derivative operators of
[`StiernströmAlmquistMattsson2023`](@ref).

Since `D₂(b)` depends linearly on the variable coefficients `b`, the action of
the operator can be written as

    (D₂(b) * u)[i] = Σⱼ Σₖ d[k, j, i] * b[k] * u[j] / Δx²

Inside the domain, `d` is translation invariant and stored in `interior`,
where `interior[k, j]` is the coefficient of `b[i + k - p - 1] * u[i + j - p - 1]`
for `p = order ÷ 2`. The boundary closure occupying the first (and, mirrored,
the last) `3p` rows is stored in `left`.
"""
@auto_hash_equals struct StiernströmAlmquistMattsson2023Cache{T, Width, Width2} <:
                         AbstractCoefficientCache{T}
    left::Array{T, 3}
    interior::SMatrix{Width, Width, T, Width2}
end

function StiernströmAlmquistMattsson2023Cache(::Type{T}, order::Int) where {T}
    if order == 4
        _sam2023_cache(T, Val(4))
    elseif order == 6
        _sam2023_cache(T, Val(6))
    elseif order == 8
        _sam2023_cache(T, Val(8))
    elseif order == 10
        _sam2023_cache(T, Val(10))
    elseif order == 12
        _sam2023_cache(T, Val(12))
    else
        throw(ArgumentError("Order of accuracy $order not implemented/derived."))
    end
end

function _sam2023_cache(::Type{T}, ::Val{order}) where {T, order}
    p = order ÷ 2
    Width = 2 * p + 1
    left, interior = _sam2023_coefficients(T, order)
    StiernströmAlmquistMattsson2023Cache{T, Width, Width * Width}(left,
                                                                  SMatrix{Width, Width, T}(interior))
end

function lower_bandwidth(cache::StiernströmAlmquistMattsson2023Cache)
    size(cache.left, 3) - 1
end
function upper_bandwidth(cache::StiernströmAlmquistMattsson2023Cache)
    size(cache.left, 3) - 1
end
left_length(cache::StiernströmAlmquistMattsson2023Cache) = size(cache.left, 3)
right_length(cache::StiernströmAlmquistMattsson2023Cache) = size(cache.left, 3)
function Base.checkbounds(::Type{Bool}, u::AbstractVector,
                          cache::StiernströmAlmquistMattsson2023Cache)
    # The boundary closures at the two ends must not overlap.
    length(u) >= 2 * size(cache.left, 3)
end

function convolve_boundary_coefficients!(dest::AbstractVector,
                                         cache::StiernströmAlmquistMattsson2023Cache,
                                         u::AbstractVector, b::AbstractVector, α,
                                         ::SafeMode)
    _sam2023_convolve_boundary!(dest, cache, u, b, α, false)
end
function convolve_boundary_coefficients!(dest::AbstractVector,
                                         cache::StiernströmAlmquistMattsson2023Cache,
                                         u::AbstractVector, b::AbstractVector, α,
                                         ::Union{FastMode, ThreadedMode})
    _sam2023_convolve_boundary!(dest, cache, u, b, α, false)
end

function convolve_boundary_coefficients!(dest::AbstractVector,
                                         cache::StiernströmAlmquistMattsson2023Cache,
                                         u::AbstractVector, b::AbstractVector, α, β,
                                         ::SafeMode)
    _sam2023_convolve_boundary!(dest, cache, u, b, α, β)
end
function convolve_boundary_coefficients!(dest::AbstractVector,
                                         cache::StiernströmAlmquistMattsson2023Cache,
                                         u::AbstractVector, b::AbstractVector, α, β,
                                         ::Union{FastMode, ThreadedMode})
    _sam2023_convolve_boundary!(dest, cache, u, b, α, β)
end

# `β = false` is used to indicate that `dest` should be overwritten; Julia's
# strong zero `false * x == zero(x)` takes care of uninitialized values.
@inline function _sam2023_convolve_boundary!(dest::AbstractVector,
                                             cache::StiernströmAlmquistMattsson2023Cache{T},
                                             u::AbstractVector, b::AbstractVector,
                                             α, β) where {T}
    @unpack left = cache
    nb, nc, nr = size(left)
    N = length(dest)
    Tacc = typeof(zero(T) * zero(eltype(b)) * zero(eltype(u)))

    @inbounds for i in 1:nr
        acc_left = zero(Tacc)
        acc_right = zero(Tacc)
        for j in 1:nc
            coef_left = zero(typeof(zero(T) * zero(eltype(b))))
            coef_right = zero(typeof(zero(T) * zero(eltype(b))))
            for k in 1:nb
                d = left[k, j, i]
                iszero(d) && continue
                coef_left = muladd(d, b[k], coef_left)
                coef_right = muladd(d, b[N + 1 - k], coef_right)
            end
            acc_left = muladd(coef_left, u[j], acc_left)
            acc_right = muladd(coef_right, u[N + 1 - j], acc_right)
        end
        dest[i] = α * acc_left + β * dest[i]
        dest[N + 1 - i] = α * acc_right + β * dest[N + 1 - i]
    end

    return nothing
end

function _sam2023_stencil_index(i, offset)
    offset > 0 ? :($i + $offset) :
    offset < 0 ? :($i - $(-offset)) : i
end

# The interior stencil is unrolled completely, as for the other operators with
# a repeated interior stencil.
@generated function convolve_interior_coefficients_loopbody(i,
                                                            cache::StiernströmAlmquistMattsson2023Cache{T,
                                                                                                        Width},
                                                            u, b) where {T, Width}
    p = (Width - 1) ÷ 2
    total = nothing
    for j in 1:Width
        coef = nothing
        for k in 1:Width
            term = :(interior[$k, $j] * b[$(_sam2023_stencil_index(:i, k - p - 1))])
            coef = coef === nothing ? term : :($coef + $term)
        end
        term = :(($coef) * u[$(_sam2023_stencil_index(:i, j - p - 1))])
        total = total === nothing ? term : :($total + $term)
    end

    quote
        Base.@_inline_meta
        @unpack interior = cache
        @inbounds @muladd retval = $total
        return retval
    end
end

# The variable-coefficient second-derivative operators of
# Stiernström, Almquist, Mattsson (2023) are constructed as
#
#   D₂(c) = D₁ * C̃₁ * D₁ - H⁻¹ * R(c),
#   R(c)  = Σ_{j=1}^{p} 1 / (αⱼ * h) * D̃ₚ₊ⱼᵀ * C̃ₚ₋ⱼ₊₁ * D̃ₚ₊ⱼ,
#
# cf. Section 3.3 and equation (3.3) of the paper. Here, `D₁` and `H` are the
# boundary-optimized first-derivative operator and norm matrix of
# `MattssonAlmquistVanDerWeide2018Accurate`, the `D̃ₖ` are undivided difference
# approximations of the `k`-th derivative that vanish in the first rows, and the
# `C̃ᵢ` are diagonal matrices containing convex combinations of the variable
# coefficients. Subtracting `H⁻¹ R(c)` from the wide-stencil operator
# `D₁ * C̃₁ * D₁` cancels the outer entries of its stencils, which yields a
# narrow-stencil operator of minimal width `2p + 1`.
#
# All of this is linear in the variable coefficients, so the coefficient tensor
# `d[k, j, i]` described above is assembled directly, once, when building the
# coefficient cache.

"""
    _sam2023_alpha(order)

The positive weights `αⱼ`, `j = 1, …, p`, of the remainder term `R(c)`, where
`αⱼ` belongs to `D̃ₚ₊ⱼ` and `C̃ₚ₋ⱼ₊₁`. See equation (3.3) and Appendix A of
Stiernström, Almquist, Mattsson (2023).
"""
function _sam2023_alpha(order::Integer)
    if order == 4
        return (18, 144)
    elseif order == 6
        return (80, 600, 3600)
    elseif order == 8
        return (350, 2520, 14700, 78400)
    elseif order == 10
        return (1512, 10584, 60480, 317520, 1587600)
    elseif order == 12
        return (6468, 44352, 249480, 1293600, 6403320, 30735936)
    else
        throw(ArgumentError("Order of accuracy $order not implemented/derived."))
    end
end

"""
    _sam2023_averaging(order, i)

The stencil of the diagonal averaging matrix `C̃ᵢ` as a tuple
`(offset, weights, boundary)`. In the interior,

    (C̃ᵢ)ₖₖ = Σₙ weights[n] * c[k + offset + n - 1],

which is a convex combination of the variable coefficients `c` around the `k`-th
grid point. Close to the boundary, this stencil would reach outside of the grid;
`boundary` contains the additional entries `(row, column, weight)` skewing it
into the domain, as in the reference implementation. The right boundary is
handled by mirroring the resulting operator.
"""
function _sam2023_averaging(order::Integer, i::Integer)
    i == 1 && return (0, (1 // 1,), ())
    i == 2 && return (-1, (1 // 2, 1 // 2), ((1, 2, 1 // 2),))
    if order == 6
        i == 3 && return (-1, (1 // 3, 1 // 3, 1 // 3), ((1, 3, 1 // 3),))
    elseif order == 8
        i == 3 && return (-1, (3 // 10, 2 // 5, 3 // 10), ((1, 3, 3 // 10),))
        i == 4 && return (-2, (1 // 4, 1 // 4, 1 // 4, 1 // 4),
                ((1, 3, 1 // 4), (1, 4, 1 // 4), (2, 4, 1 // 4)))
    elseif order == 10
        i == 3 && return (-1, (2 // 7, 3 // 7, 2 // 7), ((1, 3, 2 // 7),))
        i == 4 && return (-2, (1 // 5, 3 // 10, 3 // 10, 1 // 5),
                ((1, 3, 1 // 5), (1, 4, 3 // 10), (2, 4, 1 // 5)))
        i == 5 && return (-2, (1 // 5, 1 // 5, 1 // 5, 1 // 5, 1 // 5),
                ((1, 4, 1 // 5), (1, 5, 1 // 5), (2, 4, 1 // 5),
                 (2, 5, 1 // 5)))
    elseif order == 12
        i == 3 && return (-1, (5 // 18, 4 // 9, 5 // 18), ((1, 3, 5 // 18),))
        i == 4 && return (-2, (5 // 28, 9 // 28, 9 // 28, 5 // 28),
                ((1, 3, 5 // 28), (1, 4, 9 // 28), (2, 4, 5 // 28)))
        i == 5 && return (-2, (1 // 7, 8 // 35, 9 // 35, 8 // 35, 1 // 7),
                ((1, 3, 1 // 7), (1, 4, 1 // 7), (1, 5, 8 // 35),
                 (2, 4, 1 // 7), (2, 5, 1 // 7)))
        i == 6 && return (-3, (1 // 6, 1 // 6, 1 // 6, 1 // 6, 1 // 6, 1 // 6), ())
    end
    throw(ArgumentError("Averaging matrix $i of order $order not implemented/derived."))
end

function _sam2023_rows(::Type{T}, rows::Vararg{NTuple{M, Float64}, N}) where {T, M, N}
    A = Matrix{T}(undef, N, M)
    for i in 1:N, j in 1:M
        A[i, j] = T(rows[i][j])
    end
    return A
end

"""
    _sam2023_undivided_differences(T, order)

The nonzero rows of the boundary closures of the undivided difference operators
`D̃ₚ₊₁, …, D̃₂ₚ` used in the remainder term `R(c)`, where `p = order ÷ 2`.
Each operator `D̃ₖ` has `cld(k, 2)` vanishing rows at the boundary, followed by
`p` rows with `k + 1` nonzero entries; the `r`-th of these rows starts in
column `r`. All remaining rows use the (undivided) central stencil of `D̃ₖ`.

The coefficients are taken from the reference implementation accompanying
Stiernström, Almquist, Mattsson (2023),
https://doi.org/10.5281/zenodo.8018799. They are stored as `Float64` literals
since they are the result of a numerical optimization and are only available
as truncated decimals.
"""
function _sam2023_undivided_differences(::Type{T}, order::Integer) where {T}
    if order == 4
        d3 = _sam2023_rows(T,
                           (-1.727746398798954, 3.7021976718569105, -2.9870306597013294,
                            1.012579386643373),
                           (-0.8173849542405729, 2.691630521668, -2.8374616146508247,
                            0.9632160472233978))
        d4 = _sam2023_rows(T,
                           (1.8176226052481526, -4.7546882767009055, 5.974061319402659,
                            -4.050317546573492, 1.0133218986235861),
                           (0.7946256729910773, -3.588840695557333, 5.674923229301649,
                            -3.8528641888935913, 0.9721559821581974))
        return (d3, d4)
    elseif order == 6
        d4 = _sam2023_rows(T,
                           (5.730211159355065, -12.521994384708053, 11.419402572582198,
                            -5.944279710708975, 1.3166603634797653),
                           (1.4441513881249919, -4.9292485821432015, 6.72861373220115,
                            -4.297441697303727, 1.053925159120787),
                           (1.046607535776914, -4.088738042770866, 6.065823402094353,
                            -4.029148254415782, 1.0054553593153805))
        d5 = _sam2023_rows(T,
                           (-6.719455601453137, 16.377214352871473, -19.171027475746104,
                            14.860699276772438, -6.5833018173988265, 1.2358712649541552),
                           (-1.497152763395936, 6.195174255392029, -11.21435622033525,
                            10.743604243259316, -5.269625795603935, 1.042356280683774),
                           (-1.0511702536596916, 5.110922553463583, -10.109705670157256,
                            10.072870636039454, -5.027276796576904, 1.0043595308908133))
        d6 = _sam2023_rows(T,
                           (7.659106152843694, -20.37392361500009, 28.913418478607,
                            -29.721398553544876, 19.74990545219648, -7.415227589724932,
                            1.1881196746227272),
                           (1.542663188569347, -7.466618770718859, 16.821534330502875,
                            -21.487208486518632, 15.808877386811805, -6.254137684102644,
                            1.0348900354561115),
                           (1.054986321942043, -6.1331070641563, 15.164558505235883,
                            -20.145741272078908, 15.08183038973071, -6.02615718534488,
                            1.0036303046714514))
        return (d4, d5, d6)
    elseif order == 8
        d5 = _sam2023_rows(T,
                           (-8.20980880524119, 18.02402412065559, -17.69209667476945,
                            12.181944917152947, -5.39227640086538, 1.0882128430674802),
                           (-1.391324033305295, 5.098357412264372, -8.913925575302148,
                            8.949251437985165, -4.700301551685109, 0.9579423100430162),
                           (-0.8038859598163582, 4.086138692297563, -8.76459743980955,
                            9.255127895723392, -4.73041456471966, 0.9576313763246136),
                           (-0.8521491233621866, 4.665680161913272, -9.654116593466352,
                            9.766721490121977, -4.912007389167359, 0.9858714539606465))
        d6 = _sam2023_rows(T,
                           (9.26042722370095, -21.899951980342042, 25.706973995952183,
                            -23.795532642767977, 16.17682920259614, -6.529277058404881,
                            1.0805312592656302),
                           (1.4058275749850313, -5.963769968834439, 13.135580487711616,
                            -17.89850287597033, 14.100904655055327, -5.747653860258097,
                            0.9676139873108923),
                           (0.7869238113590604, -4.834088992392305, 13.146896159714325,
                            -18.510255791446784, 14.191243694158981, -5.745788257947681,
                            0.9650693765544026),
                           (0.8420924188251628, -5.598816194295927, 14.481174890199526,
                            -19.533442980243954, 14.736022167502076, -5.915228723763879,
                            0.9881984217769952))
        d7 = _sam2023_rows(T,
                           (-10.257962606384162, 25.816283570514674, -35.082323899232605,
                            40.90934125965728, -37.745934806057654, 22.852469704417086,
                            -7.563718814859411, 1.0718455919447922),
                           (-1.4183700945143243, 6.810922153915118, -18.129991367652288,
                            31.322380032948075, -32.90211086179576, 20.11678851090334,
                            -6.773297911176247, 0.9736795373720869),
                           (-0.7726485722940987, 5.573212298513558, -18.405654623600057,
                            32.392947635031874, -33.11290195303762, 20.110258902816884,
                            -6.755485635880818, 0.970271948450281),
                           (-0.8335597307542618, 6.531952226678581, -20.273644846279335,
                            34.18352521542692, -34.384051724171506, 20.703300533173575,
                            -6.917388952438967, 0.9898672783649978))
        d8 = _sam2023_rows(T,
                           (11.211983054345147, -29.767555907566205, 45.789440151310046,
                            -64.53016279716894, 75.49186961211531, -60.93991921177889,
                            30.254875259437643, -8.574764735558338, 1.0642345748642341),
                           (1.429430378564815, -7.642706523469226, 23.888038475126052,
                            -50.11580805271692, 65.80422172359152, -53.64476936240891,
                            27.093191644704987, -7.789436298976695, 0.9778380155843726),
                           (0.7603564556104127, -6.304846273919941, 24.540872831466743,
                            -51.82871621605099, 66.22580390607524, -53.62735707417836,
                            27.021942543523274, -7.762175587602248, 0.9741194150758726),
                           (0.8261599080832072, -7.465088259061236, 27.031526461705784,
                            -54.693640344683075, 68.76810344834301, -55.2088014217962,
                            27.669555809755867, -7.918938226919982, 0.9911226245726161))
        return (d5, d6, d7, d8)
    elseif order == 10
        d6 = _sam2023_rows(T,
                           (9.769049819830567, -22.16408730199574, 23.898275711233303,
                            -19.89385116004464, 12.62539685474329, -5.15371046484142,
                            0.9189265410746379),
                           (1.3105269062480722, -5.169297753096432, 10.44822539937647,
                            -13.88509253207172, 11.417965409312986, -4.990661471186857,
                            0.8683340414174798),
                           (0.6451169887119983, -3.71636078878867, 10.284728893981919,
                            -15.497774080936004, 12.656179173415861, -5.272476167267739,
                            0.9005859808826359),
                           (0.6401641345069657, -4.5192323253005275, 12.736469756452673,
                            -18.086958257437857, 13.938589493649035, -5.6622569275199846,
                            0.9532241256496956),
                           (0.8600500508856053, -5.6460243554222975, 14.544081782192615,
                            -19.59059143251695, 14.768523072826289, -5.925698060843626,
                            0.989658942878363))
        d7 = _sam2023_rows(T,
                           (-10.633444157744515, 25.55171730225508, -31.6395580874873,
                            33.02683297506298, -28.85621566971037, 18.037986626944967,
                            -6.432485787522465, 0.9451667982016216),
                           (-1.2971946051930945, 5.755263321546313, -14.020486493241325,
                            23.923936099960073, -26.641919288396966, 17.467315149154,
                            -6.0783382899223595, 0.8914241060933615),
                           (-0.6196831570104732, 4.184768290558644, -14.22031288145325,
                            27.12110464163801, -29.531084737970343, 18.45366658543709,
                            -6.304101866178451, 0.9156431249787751),
                           (-0.620960546711953, 5.217915133291754, -17.831057659033743,
                            31.65217695051625, -32.523375485181084, 19.817899246319943,
                            -6.672568879547869, 0.959971240346697),
                           (-0.8524154923673503, 6.587028414659348, -20.361714495069663,
                            34.28353500690466, -34.45988716992801, 20.739943212952692,
                            -6.9276126001485405, 0.9911231229968609))
        d8 = _sam2023_rows(T,
                           (11.447706798618567, -28.904884139025746, 40.25835326040996,
                            -50.64999739918013, 56.821824921675876, -48.101297671853246,
                            25.72994315008986, -7.561334385612973, 0.9596854648778265),
                           (1.2856328177474934, -6.318127111696111, 18.04299286460861,
                            -37.80427246807473, 53.28383857679393, -46.57950706441066,
                            24.313353159689438, -7.131392848746892, 0.9074820740889168),
                           (0.5982000735280487, -4.639124545001227, 18.764346418211453,
                            -43.39376742662081, 59.062169475940685, -49.20977756116557,
                            25.216407464713804, -7.325144999830201, 0.9266911002238212),
                           (0.6046001191624826, -5.910395819932234, 23.77474354537832,
                            -50.643483120826, 65.04675097036217, -52.84773132351985,
                            26.690275518191477, -7.679769922773576, 0.9650100339572174),
                           (0.8457871985025271, -7.528032473896397, 27.14895266009288,
                            -54.853656011047455, 68.91977433985602, -55.30651523454051,
                            27.710450400594162, -7.928984983974887, 0.9922241044136647))
        d9 = _sam2023_rows(T,
                           (-12.220346479758703, 32.228164479135714, -49.72006515955683,
                            73.32928389251629, -101.01269332559123, 108.22791976166981,
                            -77.18982945026958, 34.02600473525838, -8.637169183900438,
                            0.9687307304965969),
                           (-1.2754371756934053, 6.861477623719457, -22.502238094968227,
                            56.12000449056065, -95.91090943822908, 104.803890894924,
                            -72.94005947906831, 32.091267819361015, -8.167338666800251,
                            0.9193420261941577),
                           (-0.5796947369300383, 5.081509889823517, -23.9114283724449,
                            65.09065113993123, -106.31190505669323, 110.72199951262253,
                            -75.6492223941414, 32.9631524992359, -8.340219902014391,
                            0.9351574206107923),
                           (-0.590399097719824, 6.597491849057844, -30.567527415486413,
                            75.965224681239, -117.0841517466519, 118.90739547791966,
                            -80.07082655457444, 34.558964652481095, -8.685090305614956,
                            0.9689184593499157),
                           (-0.8399361406396906, 8.469036533133448, -34.905796277262276,
                            82.28048401657118, -124.05559381174082, 124.43965927771615,
                            -83.1313512017825, 35.68043242788699, -8.930016939722982,
                            0.9930821158404906))
        d10 = _sam2023_rows(T,
                            (12.9576786881247, -35.525089722887365, 59.99547167328367,
                             -101.61365491269612, 166.6135254898348, -216.45583952333962,
                             192.97457362567394, -113.4200157841946, 43.18584591950219,
                             -9.687307304965968, 0.974811851664343),
                            (1.2663266431786644, -7.388019571918344, 27.386715439651102,
                             -79.45976301897477, 159.8515157303818, -209.607781789848,
                             182.35014869767076, -106.97089273120339, 40.83669333400126,
                             -9.193420261941577, 0.928477529002455),
                            (0.5635050680153612, -5.513504360465296, 29.65686950262579,
                             -92.98664448561603, 177.18650842782205, -221.44399902524506,
                             189.12305598535352, -109.87717499745301, 41.701099510071955,
                             -9.351574206107923, 0.9418585809986529),
                            (0.5778889962428166, -7.279834627447505, 38.20940926935802,
                             -108.52174954462714, 195.1402529110865, -237.81479095583933,
                             200.17706638643608, -115.19654884160364, 43.42545152807478,
                             -9.689184593499157, 0.9720394718185964),
                            (0.8347029975818464, -9.410040592370496, 43.63224534657785,
                             -117.54354859510168, 206.75932301956803, -248.8793185554323,
                             207.82837800445623, -118.93477475962331, 44.65008469861491,
                             -9.930821158404905, 0.9937695941338366))
        return (d6, d7, d8, d9, d10)
    elseif order == 12
        d7 = _sam2023_rows(T,
                           (-7.050453821762458, 16.323556769979337, -18.517285691402662,
                            16.903069990805047, -12.900521495071139, 7.824717668026596,
                            -3.1352892049258743, 0.5522057843511536),
                           (-0.7713663600819909, 3.1512297676528003, -6.810512973200659,
                            10.585680229488409, -12.315008394369016, 9.405867614777607,
                            -3.8654404904580697, 0.6195506061909192),
                           (-0.3226806207130543, 1.9678458985349117, -6.367547799908323,
                            13.582054493165302, -17.679957518285956, 12.83136773736317,
                            -4.732759271317664, 0.7216770811616137),
                           (-0.2897599498422067, 2.432123333596445, -9.913384147957748,
                            21.377580445278298, -24.96371798861976, 16.17791596313519,
                            -5.6554717249372475, 0.8347140693470241),
                           (-0.43028213606032106, 4.2103703770684735, -15.829711232892153,
                            29.347405933458543, -30.751033700284403, 18.952643607153902,
                            -6.42930954772397, 0.9299166992799308),
                           (-0.7617781365608933, 6.316706493528616, -19.922469032198602,
                            33.78190955682214, -34.078413772067975, 20.555296711011785,
                            -6.876033783342392, 0.9847819628073188))
        d8 = _sam2023_rows(T,
                           (7.050453821762467, -17.09492313006135, 21.66851545905549,
                            -23.713582964005738, 23.486201724559578, -20.13972606239564,
                            12.541156819703497, -4.417646274809229, 0.6195506061909198),
                           (0.7143091591050187, -3.2169486987428932, 8.129032413761225,
                            -15.699214646211457, 23.981482952559876, -25.082313639406955,
                            15.461761961832279, -4.956404849527353, 0.6682953466302607),
                           (0.29213206789956747, -2.0438756549159063, 7.9665959467484075,
                            -21.27109790873475, 35.35991503657191, -34.21698063296845,
                            18.931037085270656, -5.773416649292909, 0.7556907094214725),
                           (0.26637215914130374, -2.6313682263734606, 12.983764630211125,
                            -34.204128712445275, 49.92743597723952, -43.14110923502718,
                            22.62188689974899, -6.677712554776193, 0.8548590622811767),
                           (0.4100733609469823, -4.738624918943179, 21.106281643856207,
                            -46.95584949353367, 61.502067400568805, -50.540382952410404,
                            25.71723819089588, -7.4393335942394465, 0.9385303628588249),
                           (0.7516151319251657, -7.219093135461276, 26.56329204293147,
                            -54.05105529091543, 68.15682754413595, -54.81412456269809,
                            27.50413513336957, -7.87825570245855, 0.9866588391711931))
        d9 = _sam2023_rows(T,
                           (-7.050453821762475, 17.80923228916639, -24.885464157798413,
                            31.842615377766997, -39.18541637077108, 44.12120901495556,
                            -37.62347045911049, 19.87940823664153, -5.575955455718279,
                            0.6682953466302614),
                           (-0.6669539691445976, 3.276445941549335, -9.498494213133554,
                            22.09688355077953, -42.2525569422805, 56.43520568866565,
                            -46.385285885496835, 22.30382182287309, -6.014658119672346,
                            0.7055921258602363),
                           (-0.2672871814033793, 2.1137687176984663, -9.696641634774577,
                            31.341597391980823, -63.647847065829446, 76.98820642417903,
                            -56.79311125581197, 25.98037492181809, -6.8012163847932525,
                            0.782156066936224),
                           (-0.24708804973573378, 2.82125531669212, -16.439370707786853,
                            51.306193068667916, -89.86938475903113, 97.06749577881116,
                            -67.86566069924696, 30.049706496492867, -7.6937315605305905,
                            0.8705851156672145),
                           (-0.39286387086790425, 5.259831932016188, -27.13664782781512,
                            70.4337742403005, -110.70372132102385, 113.71586164292341,
                            -77.15171457268764, 33.47700117407751, -8.446773265729425,
                            0.9452518688063304),
                           (-0.7426886389663625, 8.121479777393935, -34.1528040551976,
                            81.07658293637314, -122.68228957944471, 123.33178026607071,
                            -82.5124054001087, 35.45215066106348, -8.87992955254074,
                            0.988123585356858))
        d10 = _sam2023_rows(T,
                            (7.050453821762482, -18.476186258311003, 28.161910099347775,
                             -41.341109590900594, 61.282299921550674, -86.37376595723615,
                             94.05867614777624, -66.26469412213844, 27.879777278591394,
                             -6.682953466302614, 0.705592125860237),
                            (0.6268941964775018, -3.3308831364980036, 10.914786591113558,
                             -29.8838860648028, 69.17381165898068, -112.8704113773313,
                             115.96321471374209, -74.3460727429103, 30.07329059836173,
                             -7.055921258602363, 0.7351768214691926),
                            (0.2466529757561427, -2.178601846763726, 11.551532389532584,
                             -44.0923425674166, 106.07974510971573, -153.97641284835805,
                             141.98277813952993, -86.60124973939364, 34.00608192396626,
                             -7.82156066936224, 0.8033771327935791),
                            (0.23087142251463535, -3.003173442654164, 20.27506302406237,
                             -73.29456152666845, 149.78230793171855, -194.13499155762233,
                             169.6641517481174, -100.16568832164289, 38.46865780265295,
                             -8.705851156672145, 0.8832140761940476),
                            (0.37796280007363076, -5.774848874487027, 33.9208097847689,
                             -100.61967748614357, 184.50620220170643, -227.43172328584683,
                             192.8792864317191, -111.5900039135917, 42.23386632864712,
                             -9.452518688063304, 0.9506447012172556),
                            (0.7347407693424404, -9.023866419326595, 42.691005068997,
                             -115.82368990910449, 204.47048263240785, -246.66356053214142,
                             206.28101350027177, -118.17383553687826, 44.399647762703694,
                             -9.88123585356858, 0.9892985172965839))
        d11 = _sam2023_rows(T,
                            (-7.050453821762488, 19.103080454788525, -31.492793235845806,
                             52.2558961820142, -91.16618598635355, 155.54757761621698,
                             -206.92908752510772, 182.2279088358807, -102.22585002150178,
                             36.75624406466438, -7.761513384462607, 0.7351768214691933),
                            (-0.5924756854857472, 3.3811178542850966, -12.374519230919224,
                             39.16048049266345, -107.04747746062039, 206.92908752510738,
                             -255.1190723702326, 204.4517000430033, -110.268732193993,
                             38.807566922313, -8.086945036161119, 0.7592691400398571),
                            (-0.2292203845239877, 2.2391799149842653, -13.526028855965613,
                             59.81813685912054, -166.69674231526758, 282.2900902219898,
                             -312.3621119069659, 238.1534367833325, -124.68896705454297,
                             43.018583681492316, -8.83714846072937, 0.8207915170760159),
                            (-0.21701391114437615, 3.17819153256114, -24.486327517266712,
                             100.78002209916912, -235.37219817841486, 355.91415118897424,
                             -373.2611338458583, 275.45564288451794, -141.05174527639417,
                             47.8821813616968, -9.715354838134523, 0.8935845002936889),
                            (-0.3648850857087133, 6.284354372056049, -41.45876751471754,
                             138.35205654344742, -289.93831774553865, 416.95815935738585,
                             -424.334430149782, 306.8725107623772, -154.8575098717061,
                             51.988852784348175, -10.457091713389811, 0.9550682612282072),
                            (-0.7275857943253935, 9.926253061259255, -52.17789508432967,
                             159.25757362501867, -321.3107584223552, 452.2165276422593,
                             -453.8182297005979, 324.9780477264152, -162.79870846324687,
                             54.346797194627186, -10.882283690262422, 0.9902619055378535))
        d12 = _sam2023_rows(T,
                            (7.050453821762495, -19.695556140274288, 34.87391109013093,
                             -64.63041541293347, 130.32666647901712, -262.5950550768376,
                             413.85817505021544, -437.34698120611364, 306.67755006450534,
                             -147.02497625865752, 46.569080306775646, -8.82212185763032,
                             0.7592691400398578),
                            (0.5625205441379242, -3.427802153520942, 13.874841106233553,
                             -50.02271761282533, 158.42900977125026, -354.7355786144698,
                             510.2381447404652, -490.6840801032079, 330.806196581979,
                             -155.230267689252, 48.52167021696671, -9.111229680478285,
                             0.7792928927215863),
                            (0.2142819290642992, -2.2961219278627834, 15.615594467315484,
                             -78.81028248346854, 250.0451134729014, -483.9258689519824,
                             624.7242238139318, -571.568248279998, 374.0669011636289,
                             -172.07433472596927, 53.02289076437622, -9.849498204912193,
                             0.835348962975199),
                            (0.2050136189556189, -3.3471539032596365, 29.069145008244604,
                             -134.37336279889215, 353.0582972676223, -610.1385448953845,
                             746.5222676917166, -661.0935429228431, 423.15523582918246,
                             -191.5287254467872, 58.29212902880714, -10.723014003524266,
                             0.9022555261620117),
                            (0.35327850260839005, -6.788898256960761, 49.75052101766106,
                             -184.46940872459655, 434.90747661830795, -714.7854160412329,
                             848.668860299564, -736.4940258297053, 464.5725296151183,
                             -207.9554111373927, 62.742550280338875, -11.460819134738486,
                             0.9587627910279094),
                            (0.7210856618217802, -10.828639703191913, 62.61347410119561,
                             -212.3434315000249, 481.9661376335328, -775.2283331010159,
                             907.6364594011958, -779.9473145433965, 488.39612538974063,
                             -217.38718877850874, 65.29370214157454, -11.883142866454243,
                             0.9910661635310787))
        return (d7, d8, d9, d10, d11, d12)
    else
        throw(ArgumentError("Order of accuracy $order not implemented/derived."))
    end
end

# Dense `n × n` matrix of the boundary-optimized first-derivative operator with
# `Δx = 1`. Only the boundary closure at the left end is applied; all other rows
# use the central stencil (truncated at the last column).
function _sam2023_first_derivative_matrix(::Type{T}, order::Integer, n::Integer) where {T}
    coefficients = first_derivative_coefficients(MattssonAlmquistVanDerWeide2018Accurate(),
                                                 order, T)
    D1 = zeros(T, n, n)
    for (i, row) in enumerate(coefficients.left_boundary)
        start = offset(row)
        for (k, value) in enumerate(row.coef)
            D1[i, start + k - 1] = value
        end
    end
    for i in (length(coefficients.left_boundary) + 1):n
        D1[i, i] = coefficients.central_coef
        for s in eachindex(coefficients.lower_coef)
            i - s >= 1 && (D1[i, i - s] = coefficients.lower_coef[s])
        end
        for s in eachindex(coefficients.upper_coef)
            i + s <= n && (D1[i, i + s] = coefficients.upper_coef[s])
        end
    end
    return D1, coefficients.left_weights
end

# Dense `n × n` matrix of the undivided difference operator `D̃ₖ`, see
# `_sam2023_undivided_differences`.
function _sam2023_undivided_difference_matrix(::Type{T}, k::Integer,
                                              block::AbstractMatrix, n::Integer) where {T}
    D = zeros(T, n, n)
    upper = cld(k, 2)
    for i in 1:n, m in 0:k
        j = i + m - upper
        1 <= j <= n || continue
        D[i, j] = T((-1)^(k - m) * binomial(k, m))
    end
    for i in 1:(upper + size(block, 1)), j in 1:n
        D[i, j] = zero(T)
    end
    for r in 1:size(block, 1), c in 1:size(block, 2)
        D[upper + r, r + c - 1] = block[r, c]
    end
    return D
end

# Dense `n × n` matrix mapping the variable coefficients to the diagonal
# of the averaging matrix `C̃ᵢ`, see `_sam2023_averaging`.
function _sam2023_averaging_matrix(::Type{T}, order::Integer, i::Integer,
                                   n::Integer) where {T}
    off, weights, boundary = _sam2023_averaging(order, i)
    W = zeros(T, n, n)
    for r in 1:n, c in eachindex(weights)
        j = r + off + c - 1
        1 <= j <= n || continue
        W[r, j] = T(weights[c])
    end
    for (r, c, weight) in boundary
        W[r, c] = T(weight)
    end
    return W
end

"""
    _sam2023_coefficients(T, order)

Compute the boundary closure `left` and the interior stencil `interior` of the
variable-coefficient second-derivative operator of
[`StiernströmAlmquistMattsson2023`](@ref) with the given order of accuracy, see
[`StiernströmAlmquistMattsson2023Cache`](@ref) for their meaning.
"""
function _sam2023_coefficients(::Type{T}, order::Integer) where {T}
    p = order ÷ 2
    # Size of the boundary closure of `D₂`: number of rows, number of columns
    # they occupy, and number of variable coefficients they depend on.
    nrows = 3 * p
    ncols = 4 * p
    ncoefficients = 4 * p
    # Size of the matrices assembled below. Only the boundary closure at the
    # left end is applied to them; `nwork` is large enough that every entry
    # used below is unaffected by the missing closure at the right end.
    nwork = 6 * p
    # A row of `D₂` in the uniform interior part of the matrices
    interior_row = 5 * p

    D1, weights = _sam2023_first_derivative_matrix(T, order, nwork)
    inv_weights = ones(T, nwork)
    for i in eachindex(weights)
        inv_weights[i] = inv(weights[i])
    end
    blocks = _sam2023_undivided_differences(T, order)
    alphas = _sam2023_alpha(order)

    # `rows[i]` is the index at which row `i` of `D₂` is stored in `d` below;
    # rows that are not needed are marked by `0`.
    rows = zeros(Int, nwork)
    for i in 1:nrows
        rows[i] = i
    end
    rows[interior_row] = nrows + 1
    # `d[k, j, rows[i]]` is the coefficient of `b[k] * u[j]` in `(D₂(b) u)[i]`.
    d = zeros(T, nwork, nwork, nrows + 1)

    # wide-stencil operator `D₁ * C̃₁ * D₁`
    for i in 1:nwork
        r = rows[i]
        iszero(r) && continue
        for k in 1:nwork
            iszero(D1[i, k]) && continue
            for j in 1:nwork
                d[k, j, r] = muladd(D1[i, k], D1[k, j], d[k, j, r])
            end
        end
    end

    # remainder term `-H⁻¹ * R(c)`
    for (jj, α) in enumerate(alphas)
        DD = _sam2023_undivided_difference_matrix(T, p + jj, blocks[jj], nwork)
        W = _sam2023_averaging_matrix(T, order, p - jj + 1, nwork)
        for m in 1:nwork, i in 1:nwork
            r = rows[i]
            iszero(r) && continue
            iszero(DD[m, i]) && continue
            factor_i = DD[m, i] * inv_weights[i] / T(α)
            for j in 1:nwork
                iszero(DD[m, j]) && continue
                factor_ij = factor_i * DD[m, j]
                for k in 1:nwork
                    iszero(W[m, k]) && continue
                    d[k, j, r] -= factor_ij * W[m, k]
                end
            end
        end
    end

    # Restrict the boundary closure to the sparsity pattern of the
    # narrow-stencil operator. The discarded entries vanish analytically; they
    # are of the order of the round-off error of the truncated decimal
    # coefficients of `D₁` and `D̃ₖ`, cf. the reference implementation.
    left = zeros(T, ncoefficients, ncols, nrows)
    for i in 1:nrows
        jmax = i <= 2 * p ? 3 * p : p + i
        for j in 1:jmax, k in 1:ncoefficients
            left[k, j, i] = d[k, j, i]
        end
    end

    width = 2 * p + 1
    interior = zeros(T, width, width)
    for k in 1:width, j in 1:width
        interior[k, j] = d[interior_row - p - 1 + k, interior_row - p - 1 + j, nrows + 1]
    end

    return left, interior
end
