using PhylogeneticFactorAnalysis, RCall


mammals_dir = @__DIR__
cd(mammals_dir)


trait_labels_path = joinpath(mammals_dir, "mammals_labels.csv")
taxon_labels_path = joinpath(mammals_dir, "mammals_classification.csv")
newick = joinpath(mammals_dir, "mammals_newick_processed.txt")

pipe_name = "mammals"

pipe_dir = joinpath(mammals_dir, pipe_name)
load_stats_path = joinpath(pipe_dir, pipe_name * "_loadingsStatistics.csv")
fac_stats_path = joinpath(pipe_dir, pipe_name * "_factorMeans.csv")

R_PLOT_SCRIPT = PhylogeneticFactorAnalysis.R_PLOT_SCRIPT

@rput R_PLOT_SCRIPT
@rput trait_labels_path
@rput taxon_labels_path
@rput load_stats_path
@rput fac_stats_path
@rput newick

R"""
source(R_PLOT_SCRIPT)


plot_loadings(load_stats_path, "mammalsLoadings.pdf",
              labels_path = trait_labels_path,
              height_scale = 1.2)
"""

R"""

plot_factor_tree("mammals", newick, fac_stats_path,
                 class_path = taxon_labels_path,
                 tip_labels=FALSE,
                 layout="circular",
                 line_width=0.2)
"""