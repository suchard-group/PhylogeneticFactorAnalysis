using BeastUtils.Logs, DataFrames, CSV, UnPack, RCall, StatsBase, Random,
      PhylogeneticFactorAnalysis


Random.seed!(666)


aquilegia_dir = @__DIR__

cd(aquilegia_dir)

pipe_name = "aquilegiaBinary"
pipe_dir = joinpath(aquilegia_dir, pipe_name)
log_path = joinpath(pipe_dir, pipe_name * "_processed.log")
load_stats_path = joinpath(pipe_dir, pipe_name * "_loadingsStatistics.csv")
trait_labels_path = joinpath(aquilegia_dir, "aquilegia_labels_binary.csv")



cols, data = get_log_match(log_path, "factors.")

pollinators = CSV.read(joinpath(aquilegia_dir, "aquilegia_pollinators.csv"), DataFrame)

n = div(length(cols), 2)
states = size(data, 1)

cols = reshape(cols, 2, n)
data = reshape(data, states, 2, n)

all_df = DataFrame([String, String, Float64, Float64],
                   [:taxon, :pollinator, :f1, :f2],
                   n * states)



function parse_colname(s::String)
    s_split = split(s, '.')
    factor = parse(Int, s_split[end])
    taxon = join(s_split[2:(end - 1)], '.')
    return (taxon = taxon, factor = factor)
end

taxa = [parse_colname.(s).taxon for s in cols[1, :]]
@assert [parse_colname.(s).taxon for s in cols[2, :]] == taxa

all_df.taxon = repeat(taxa, inner=states)

all_df.f1 = vcat([data[:, 1, i] for i = 1:n]...)
all_df.f2 = vcat([data[:, 2, i] for i = 1:n]...)

pollinator_dict = Dict{String, String}([pollinators.taxon[i] => pollinators.pollinator[i] for i = 1:n])
all_df.pollinator = [pollinator_dict[taxon] for taxon in all_df.taxon]

n_total = sample(1:size(all_df, 1), 1000, replace=false)

sampled_df = all_df[n_total, :]

@rput sampled_df
R_PLOT_SCRIPT = PhylogeneticFactorAnalysis.R_PLOT_SCRIPT
@rput R_PLOT_SCRIPT
@rput load_stats_path
@rput trait_labels_path

R"""
library(ggplot2)
library(wesanderson)
library(ggforce)
library(scales)

pollinators_path = "aquilegia_pollinators.csv"
factors_path = file.path("aquilegiaBinary", "aquilegiaBinary_factorMeans.csv")
# factors_path = file.path("aqui_factors.csv")

factors = read.csv(factors_path)
pollinators = read.csv(pollinators_path)
df <- merge(factors,pollinators, by="taxon")


plt <- ggplot(sampled_df, aes(x=f1, y=f2, label=taxon)) +
  geom_point(aes(color=pollinator), size=1, alpha=0.4, shape=16) +
  geom_point(size=3, data=df, aes(color=pollinator, x=f1, y=f2)) +
  scale_color_manual(values=wes_palette("Moonrise2")) +
  scale_y_continuous(breaks = seq(-2, 4, 1)) +
  theme_minimal() +
  xlab("factor 1") +
  ylab("factor 2") +
  coord_fixed(ratio=1)

svg("aquilegia_factors.svg", width=6, height=5)
print(plt)
dev.off()

gc()


source(R_PLOT_SCRIPT)


plot_loadings(load_stats_path, "aquilegiaLoadings.pdf",
              labels_path = trait_labels_path,
              height_scale = 1.2)
"""

