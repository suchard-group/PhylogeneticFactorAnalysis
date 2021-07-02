using DataFrames, CSV, Statistics

const BASE_SAMPLER = "old"
const OTHER_SAMPLERS = ["gibbs", "hmc", "ortho"]
const TIME_UNIT = "minute"

const TIME_CONVERSION = Dict("second" => 1.0, "minute" => 60.0, "hour" => 3600.0)

timing_results = "results.csv"
df = DataFrame(CSV.File(timing_results))

df.N = Int.(df.N)
df.K = Int.(df.K)
df.P = Int.(df.P)
df.run = Int.(df.run)

df.min_ept = TIME_CONVERSION[TIME_UNIT] * df.load_min_ess ./ df.time
df.med_ept = TIME_CONVERSION[TIME_UNIT] * df.load_med_ess ./ df.time

function compute_means_and_folds(df::DataFrame, n::Int, k::Int, p::Int)
    sub_df = @view df[(df.N .== n) .* (df.K .== k) .* (df.P .== p), :]
    base_inds = findall(sub_df.sampler .== BASE_SAMPLER)
    base_min = mean(sub_df.min_ept[base_inds])
    base_med = mean(sub_df.med_ept[base_inds])

    s = length(OTHER_SAMPLERS)
    other_mins = zeros(s)
    other_meds = zeros(s)
    for i = 1:s
        other_inds = findall(sub_df.sampler .== OTHER_SAMPLERS[i])
        other_mins[i] = mean(sub_df.min_ept[other_inds])
        other_meds[i] = mean(sub_df.med_ept[other_inds])
    end

    d = s + 1

    new_df = DataFrame(sampler = [BASE_SAMPLER; OTHER_SAMPLERS],
                       N = fill(n, d),
                       P = fill(p, d),
                       K = fill(k, d),
                       min_ept = [base_min; other_mins],
                       med_ept = [base_med; other_meds],
                       min_fold = [1.0; other_mins ./ base_min],
                       med_fold = [1.0; other_meds ./ base_med]
                      )
    return new_df
end

ns = sort(unique(df.N))
ks = sort(unique(df.K))
ps = sort(unique(df.P))

mean_df = DataFrame()

for n in ns
    for p in ps
        for k in ks
            mean_df = [mean_df; compute_means_and_folds(df, n, k, p)]
        end
    end
end

mean_df

function reshape_table(df::DataFrame)

    ns = sort!(unique(df.N))
    ps = sort!(unique(df.P))
    ks = sort!(unique(df.K))

    samplers = unique(df.sampler)
    rows = length(ns) * length(ps) * length(ks)

    @assert size(df, 1) == rows * length(samplers)

    min_name = "minimum ESS per $TIME_UNIT"
    fold_name = "speed increase"
    sampled_min = "$(min_name)_sampled"
    gibbs_min = "$(min_name)_Gibbs"
    hmc_min = "$(min_name)_HMC"
    ortho_min = "$(min_name)_ortho"
    gibbs_fold = "$(fold_name)_Gibbs"
    hmc_fold = "$(fold_name)_HMC"
    ortho_fold = "$(fold_name)_ortho"


    new_df = DataFrame(
                [Int, Int, Int, Float64, Float64, Float64, Float64, Float64, Float64, Float64],
                ["N", "P", "K", sampled_min, gibbs_min, hmc_min, ortho_min, gibbs_fold, hmc_fold, ortho_fold],
                rows)

    row = 0
    for n in ns
        n_inds = findall(isequal(n), df.N)
        ndf = @view df[n_inds, :]
        for p in ps
            p_inds = findall(isequal(p), ndf.P)
            pdf = @view ndf[p_inds, :]
            for k in ks
                k_inds = findall(isequal(k), pdf.K)
                kdf = @view pdf[k_inds, :]

                row += 1
                new_df.N[row], new_df.P[row], new_df.K[row] = n, p, k

                sub_rows =
                        [findfirst(isequal(sampler), kdf.sampler) for
                            sampler in ["old", "gibbs", "hmc", "ortho"]]

                new_df[row, [sampled_min, gibbs_min, hmc_min, ortho_min]] .=
                        kdf.min_ept[sub_rows]

                new_df[row, [gibbs_fold, hmc_fold, ortho_fold]] .=
                        kdf.min_fold[sub_rows[2:end]]
            end
        end
    end
    return new_df
end



table = reshape_table(mean_df)
for i = 1:size(table, 2)
    if eltype(table[!, i]) == Float64
        table[!, i] .= round.(table[!, i], sigdigits=2)
    end
end

function multirow(rows::Int, val)
    return "\\multirow{$rows}{*}{$val}"
end

function my_latexify(df::DataFrame; exclude_header::Bool = true,
                     concatenate_identical::AbstractArray{Int} = Int[],
                     sigdigits::Int = 2,
                     special_formatting::Function = x -> x)
    n, p = size(df)
    table = fill("", n, p)
    rules = Dict{Int, String}()

    for j = 1:p
        if j in concatenate_identical
            concat_inds = find_concat_inds(df[!, j])
            for k = 1:length(concat_inds)
                val, rng = concat_inds[k]
                table[first(rng), j] =
                    multirow(length(rng), parse_value(val, sigdigits = sigdigits))

                if !(last(rng) in keys(rules))
                    rules[last(rng)] = "\\cmidrule{$j-$p}"
                end
            end
        else
            for i = 1:n
                table[i, j] = parse_value(df[i, j], sigdigits = sigdigits)
            end
        end
    end

    table = special_formatting(table)

    rows = [join(table[i, :], " & ") for i = 1:n]
    rows = [row * " \\\\" for row in rows]
    for v in keys(rules)
        rows[v] = rows[v] * rules[v]
    end
    return join(rows, '\n')
end

function make_times(X::Matrix{String}, cols::AbstractVector{Int})
    for col in cols
        for i = 1:size(X, 1)
            X[i, col] = X[i, col] * "\$\\times\$"
        end
    end
    return X
end

times_f(X) = make_times(X, 8:10)

function parse_value(x::Int; args...)
    return string(x)
end

function parse_value(x::AbstractFloat; sigdigits::Int = 2)
    s = string(round(x, sigdigits = sigdigits))
    if endswith(s, ".0") && length(s) >= sigdigits + 2
        s = s[1:(end - 2)]
    end
    return s
end





function find_concat_inds(x::Vector)
    T = eltype(x)
    ps = Pair{T, UnitRange{Int64}}[]
    current = x[1]
    start = 1
    stop = -1
    for i = 2:length(x)
        if x[i] != current
            stop = i - 1
            push!(ps, current => start:stop)
            current = x[i]
            start = i
        end
    end
    push!(ps, current => start:length(x))
    return ps
end



using Latexify

clipboard(
        my_latexify(table,
                    concatenate_identical = 1:2,
                    special_formatting = times_f)
)
