"""
    add_term!(lhs, d, key, coef)

Add `coef * d[key]` to the affine expression `lhs`, but only when `key` is
present in the sparse expression dictionary `d`; absent keys contribute nothing
(i.e. they read as 0). Mutates and returns `nothing`.
"""
function add_term!(lhs::AffExpr, d::Dict, key, coef::Float64)
    e = get(d, key, nothing)
    e === nothing || add_to_expression!(lhs, coef, e)
    return nothing
end

"""
    as_range(ids) -> AbstractVector

Return `1:n` (a range) when `ids` is exactly the contiguous set `1:n`, so that
JuMP's `DenseAxisArray` indexes that dimension with O(1) arithmetic instead of
building a `Dict` and hashing on every lookup. Falls back to `ids` unchanged if
it is not contiguous-from-one.
"""
function as_range(ids::AbstractVector)
    return ids == eachindex(ids) ? Base.OneTo(length(ids)) : ids
end
