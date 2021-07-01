using BeastUtils.DataStorage, BeastUtils.MatrixUtils, BEASTDataPrep

using PhyloNetworks, CSV, DataFrames, LinearAlgebra
const MISSING_VAL = -999.0


cd(@__DIR__)
mammals_dir = @__DIR__


newick_path = joinpath(mammals_dir, "mammals_newick.txt")
panth_path = joinpath(mammals_dir, "PanTHERIA_1-0_WR05_Aug2008.txt")
new_newick_path = joinpath(mammals_dir, "mammals_newick_processed.txt")
csv_path = joinpath(mammals_dir, "mammals.csv")
order_path = joinpath(mammals_dir, "mammals_classification.csv")
# standardized_csv_path = "mammals_standardized.csv"


original_newick = read(newick_path, String)

panth_df = CSV.read(panth_path, DataFrame)
panth_taxa = ["$(panth_df.MSW05_Genus[i])_$(panth_df.MSW05_Species[i])" for i = 1:size(panth_df, 1)]

keep_traits = ["5-1_AdultBodyMass_g",
               "3-1_AgeatFirstBirth_d",
               "7-1_DispersalAge_d",
               "9-1_GestationLen_d",
               "14-1_InterbirthInterval_d",
               "15-1_LitterSize",
               "16-1_LittersPerYear",
               "17-1_MaxLongevity_m",
               "23-1_SexualMaturityAge_d",
               "24-1_TeatNumber",
               "25-1_WeaningAge_d"]

df = DataFrame(taxon = panth_taxa)
order_df = DataFrame(taxon = panth_taxa, order = panth_df.MSW05_Order)

df = [df DataFrame(panth_df[!, keep_traits])]
for i = 2:size(df, 2)
    df[!, i] = Float64.(df[!, i])
end

X = Matrix(df[!, 2:end])
keep_taxa = [!all(isequal(MISSING_VAL), @view X[i, :]) for i = 1:size(X, 1)]

df = df[keep_taxa, :]
order_df = order_df[keep_taxa, :]

for i = 2:size(df, 2)
    for j = 1:size(df, 1)
        if df[j, i] == MISSING_VAL
            df[j, i] = NaN
        end
    end
end

df, newick = conform_tree_and_data(df, original_newick)

net = readTopology(newick)
v1 = vcv(net)

PhyloNetworks.cleanAfterRead!(net, true) #this will not work if root is multifurcating
directEdges!(net)

# Check that the network is bifurcating and unchanged
for node in net.node
    if node === net.node[net.root]
        @assert length(node.edge) == 2
    elseif node.leaf
        @assert length(node.edge) == 1
    else
        @assert length(node.edge) == 3
    end
end

v2 = vcv(net)

@assert maximum(Matrix(v1) - Matrix(v2)) == 0.0
@assert names(v1) == names(v2)

newick = writeTopology(net)


log_traits = ["5-1_AdultBodyMass_g",
                 "3-1_AgeatFirstBirth_d",
                 "7-1_DispersalAge_d",
                 "9-1_GestationLen_d",
                 "14-1_InterbirthInterval_d",
                 "15-1_LitterSize",
                 "16-1_LittersPerYear",
                 "17-1_MaxLongevity_m",
                 "23-1_SexualMaturityAge_d",
                 "24-1_TeatNumber",
                 "25-1_WeaningAge_d"]

for trait in log_traits
    df[!, trait] .= log.(df[!, trait])
end


CSV.write(csv_path, df)
write(new_newick_path, newick)


# mammals_classification

new_classes = Dict("Microbiotheria" => "other Australidelphia",
                   "Notoryctemorphia" => "other Australidelphia",
                   "Dermoptera" => "other Euarchonta",
                   "Scandentia" => "other Euarchonta",
                   "Tubulidentata" => "Afrotheria",
                   "Proboscidea" => "Afrotheria",
                   "Hyracoidea" => "Afrotheria",
                   "Sirenia" => "Afrotheria",
                   "Paucituberculata" => "Ameridelphia",
                   "Pilosa" => "Xenarthra",
                   "Erinaceomorpha" => "Eulipotyphla",
                   "Macroscelidea" => "Afrotheria",
                   "Cingulata" => "Xenarthra",
                   "Peramelemorphia" => "other Australidelphia",
                   "Afrosoricida" => "Afrotheria",
                   "Dasyuromorphia" => "other Australidelphia",
                   "Didelphimorphia" => "Ameridelphia",
                   "Soricomorpha" => "Eulipotyphla"
                   )
ks = keys(new_classes)

n = size(order_df, 1)

order_df.class = fill("", n)

for i = 1:n
    order = order_df.order[i]
    if order in ks
        order_df.class[i] = new_classes[order]
    else
        order_df.class[i] = order
    end
end

CSV.write(order_path, DataFrame(taxon = order_df.taxon, classification = order_df.class))


# classes = unique(order_df.class)



# counts = [count(isequal(x), order_df.class) for x in classes]
# perm = sortperm(counts)

# for i = 1:length(classes)
#     println("$(classes[perm[i]]):\t$(counts[perm[i]])")
# end

# @show length(classes)
# @show length(unique(order_df.order))