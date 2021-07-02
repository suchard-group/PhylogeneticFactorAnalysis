################################################################################
## Timing how fast inference is under different samplers
################################################################################

using UnPack, DataFrames, TickTock, Statistics, CSV, LinearAlgebra,
      MCMCDiagnostics
using BEASTXMLConstructor, BEASTSimulation, BeastUtils.RunBeast,
      BeastUtils.ESS, BeastUtils.Logs, FactorPostProcessing

const OVERWRITE_XML = false
const GIBBS = "gibbs"
const HMC = "hmc"
const OLD = "old"
const ORTHO = "ortho"
const INFERENCE_REGIMES = [ORTHO, GIBBS, HMC, OLD]
const LOAD_START = "L"
const PREC_START = "factorPrecision"
const DEBUG = true
const TMP_CSV = "tmp.csv"
const RESUME = true
const INITIALIZE_WITH_GIBBS = true

const TIME_MULT = Dict("seconds" => 1.0,
                        "minutes" => 60.0,
                        "hours" => 3600.0,
                        "days" => 86400.0)

const DEFAULT_OPTIONS = MCMCOptions(chain_length = 10_000,
                                    file_log_every = 5,
                                    screen_log_every = 1000,
                                    likelihood_check_count = 0)


struct SimRun
    N::Int
    K::Int
    P::Int
end

function run_all(runs::Array{SimRun}, n_sims::Int, dir::String; resume::Bool = false)
    n = length(runs) * n_sims * length(INFERENCE_REGIMES)
    if resume && isfile(TMP_CSV)
        df = DataFrame(CSV.File(TMP_CSV))
    else
        df = DataFrame()
        df.N = zeros(n)
        df.K = zeros(n)
        df.P = zeros(n)
        df.run = zeros(n)
        df.sampler = fill("NA", n)
        df.time = zeros(n)
        df.load_min_ess = zeros(n)
        df.load_mean_ess = zeros(n)
        df.load_med_ess = zeros(n)
        df.prec_min_ess = zeros(n)
        df.prec_mean_ess = zeros(n)
        df.prec_med_ess = zeros(n)
    end

    current_ind = 1
    for i = 1:n_sims
        for run in runs
            results = do_simulation(run, dir, i, resume = resume)
            for result in results
                df.N[current_ind] = run.N
                df.K[current_ind] = run.K
                df.P[current_ind] = run.P
                df.run[current_ind] = i
                for name in fieldnames(typeof(result))
                    value = getfield(result, name)
                    setindex!(df, value, current_ind, name)
                end
                current_ind += 1
                CSV.write(TMP_CSV, df)
            end
        end
    end
    return df
end

function parse_time(s::String)
    s_split = split(s)
    @assert length(s_split) == 2
    t = parse(Float64, s_split[1])
    unit = s_split[2]
    return t * TIME_MULT[unit]
end

function do_simulation(run::SimRun, dir::String, rep::Int; resume::Bool = false)
    @unpack N, K, P = run

    # check that files exist if resume
    id = OVERWRITE_XML ? "" : "_N$(N)_K$(K)_P$(P)_r$rep"
    xml_paths = Dict{String, String}()
    log_paths = Dict{String, String}()
    time_paths = Dict{String, String}()
    svd_paths = Dict{String, String}()


    for regime in INFERENCE_REGIMES
        bn = regime * id
        path = joinpath(dir, bn * ".xml")
        xml_paths[regime] = path
        log_paths[regime] = joinpath(dir, bn * ".log")
        svd_paths[regime] = joinpath(dir, bn * "_svd.log")
        time_paths[regime] = joinpath(dir, bn * "_timer.txt")
    end

    all_finished = false
    if resume
        finished = Dict{String, Bool}(regime => isfile(svd_paths[regime]) for
                                    regime in INFERENCE_REGIMES)
        fv = values(finished)
        all_finished = all(fv)
    end



    if !all_finished

        # simulate data
        sim_model = BEASTSimulation.randomFactorSimulationModel(N, K, P)
        L = sim_model.extensionModel.L
        L .= Diagonal([sqrt(P) * 0.5^(k) for k = 1:K]) * svd(L).Vt
        # @show diag(sim_model.extensionModel.Λ)[1:20]
        # @show L[:, 1:20]
        # error()
        taxa = sim_model.taxa
        data = simulate(sim_model)
        data_std = data * inv(Diagonal(vec(std(data, dims=1))))
        data_svd = svd(data_std)

        F_init = zeros(N, K)
        L_init = zeros(K, P)
        for i = 1:K
            L_init[i, i] = 1.0
        end
        newick = get_newick(sim_model)

        # make xml


        if INITIALIZE_WITH_GIBBS
            gibbs_init = make_pfa_xml(data, taxa, newick, K, useHMC = false, timing=false, log_factors=true)
            xml_init = "init.xml"
            log_init = "init.log"
            svd_init = "inig_svd.log"
            set_chain_length!(gibbs_init, 1000)
            set_file_logEvery(gibbs_init, 10)
            save_xml(xml_init, gibbs_init)
            run_beast(xml_init, overwrite=true)
            svd_logs(log_init, svd_init, rotate_factors=true)

            _, L = get_log_match(svd_init, "L", burnin=0.5)
            L = reshape(L[end, :], P, K)'
            L_svd = svd(L)
            L_init .= Diagonal(L_svd.S) * L_svd.Vt

            f_cols, F = get_log_match(svd_init, "factors.", burnin=0.5)
            F_cols = reshape(f_cols, K, N)
            F = reshape(F[end, :], K, N)
            for i = 1:N
                taxon = split(F_cols[1, i], '.')[2]
                for j = 2:K
                    @assert split(F_cols[j, i], '.')[2] == taxon
                end

                ind = findfirst(isequal(taxon), taxa)
                F_init[ind, :] .= F[:, i]
            end
        end

        gibbs_bx = make_pfa_xml(data, taxa, newick, K, useHMC=false, timing=true)
        hmc_bx = make_pfa_xml(data, taxa, newick, K, useHMC = true, timing=true)
        old_bx = BEASTXMLConstructor.make_old_pfa_xml(data, taxa, newick, K, timing=true, factors = F_init)
        ortho_bx = BEASTXMLConstructor.make_orthogonal_pfa_xml(data, taxa, newick, K, timing=true, force_ordered=true, shrinkage=20.0)


        ops = get_operators(old_bx)
        factor_op = ops[findfirst(
                            x -> typeof(x) <: BEASTXMLConstructor.FactorTreeGibbsOperatorXMLElement,
                            ops)]
        factor_op_weight = Int(round(N / 4))
        factor_op.weight = factor_op_weight


        all_bx = Dict(GIBBS => gibbs_bx, HMC => hmc_bx, OLD => old_bx, ORTHO => ortho_bx)
        for regime in INFERENCE_REGIMES
            set_options!(all_bx[regime], DEFAULT_OPTIONS)
            set_loadings(all_bx[regime], L_init)
        end

        fle = factor_op_weight + 2

        set_chain_length!(old_bx, div(DEFAULT_OPTIONS.chain_length * fle, 2))
        set_file_logEvery(old_bx, fle)

        ortho_mult = 4
        set_chain_length!(ortho_bx, DEFAULT_OPTIONS.chain_length * ortho_mult)
        set_file_logEvery(ortho_bx, DEFAULT_OPTIONS.file_log_every)

        for regime in INFERENCE_REGIMES
            save_xml(xml_paths[regime], all_bx[regime])
        end

        # run xml & post-processing
        for regime in INFERENCE_REGIMES
            run_beast(xml_paths[regime], overwrite=true)
            svd_logs(log_paths[regime], svd_paths[regime])
        end
    end

    # parse times
    times = Dict{String, Float64}()
    for regime in INFERENCE_REGIMES
        time_string = read(time_paths[regime], String)
        times[regime] = parse_time(time_string)
    end

    # calculate ESSs

    mean_precs = Dict{String, Float64}()
    min_precs = Dict{String, Float64}()
    med_precs = Dict{String, Float64}()

    mean_load = Dict{String, Float64}()
    min_load = Dict{String, Float64}()
    med_load = Dict{String, Float64}()


    for regime in INFERENCE_REGIMES
        svd_path = svd_paths[regime]
        log_path = log_paths[regime]
        load_cols, load_data = get_log_match(svd_path, LOAD_START)
        prec_cols, prec_data = get_log_match(log_path, PREC_START)

        @assert length(load_cols) == K * P
        @assert length(prec_cols) == P

        load_ess = [effective_sample_size(load_data[:, i]) for i = 1:(K * P)]
        prec_ess = [effective_sample_size(prec_data[:, i]) for i = 1:P]

        if DEBUG
            @show regime
            @show load_ess
            @show prec_ess

            print("\n\n")

        end

        mean_precs[regime] = mean(prec_ess)
        min_precs[regime] = minimum(prec_ess)
        med_precs[regime] = median(prec_ess)

        mean_load[regime] = mean(load_ess)
        min_load[regime] = minimum(load_ess)
        med_load[regime] = median(load_ess)
    end


    results = [(sampler = regime, time = times[regime],
                load_min_ess = min_load[regime],
                load_mean_ess = mean_load[regime],
                load_med_ess = med_load[regime],
                prec_min_ess = min_precs[regime],
                prec_mean_ess = mean_precs[regime],
                prec_med_ess = med_precs[regime])
               for regime in INFERENCE_REGIMES]

    return results
end

dir = joinpath(pwd(), "timing")

if length(ARGS) > 0
    dir = ARGS[1]
end

mkpath(dir)
if length(readdir(dir)) > 0 && !RESUME
    error("Directory $dir is not empty. Please delete contents or choose a " *
            "differnt directory.")
end



old_path = pwd()
cd(dir)



# Ns = [20, 30]
# Ks = [1, 2]
# Ps = [10, 20]
# n_sims = 2
Ns = [50, 100, 500, 1000]
Ks = [1, 2, 4]
Ps = [10, 100, 1000]
n_sims = 3

timing_worked = false

try
    runs = [SimRun(n, k, p) for n in Ns, k in Ks, p in Ps]
    results = run_all(runs, n_sims, ".", resume = RESUME)
    CSV.write("results.csv", results)
    global timing_worked = true
catch e
    @error "Something went wrong" exception=(e, catch_backtrace())
end

cd(old_path)
if timing_worked
    println("completed successfully")
else
    println("timing failed")
end

