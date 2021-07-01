using BEASTDataPrep, CSV

cd(@__DIR__)
nexus_path = "Aquilegia.nex"


df = parse_continuous_data(nexus_path)

# genera = Dict(
#                  "BA" => "Aquilegia_barnebyi",
#                  "JO" => "Aquilegia_jonesii",
#                  "CA.VA.1" => "Aquilegia_canadensis",
#                  "CHAP.NM.3" => "Aquilegia_chaplinei",
#                  "CH.CHI.2" => "Aquilegia_chrysantha",
#                  "COOC.CO.2" => "Aquilegia_coerulea",
#                  "CH.NM.3" => "Aquilegia_",
#                  "COAL.WY.1" => "Aquilegia_",
#                  "COCO.WY.3" => "Aquilegia_",
#                  "COOC.UT.4" => "Aquilegia_",
#                  "COPI.AZ.1" => "Aquilegia_",
#                  "DE.AZ.1" => "Aquilegia_",
#                  "EL.NM.4" => "Aquilegia_",
#                  "EX.CA.2" => "Aquilegia_",
#                  "FL.UT.5" => "Aquilegia_",
#                  "FO.CA.BAL.3" => "Aquilegia_",
#                  "FO.CA.BTR.1" => "Aquilegia_",
#                  "HI.TX.4" => "Aquilegia_",
#                  "BR.ALB.5" => "Aquilegia_",
#                  "LA.WY.2" => "Aquilegia_",
#                  "LO.AZ.2" => "Aquilegia_",
#                  "LO.TX.5" => "Aquilegia_",
#                  "MI.CO.1" => "Aquilegia_",
#                  "PU.CA.PP.3" => "Aquilegia_",
#                  "PU.NV.3" => "Aquilegia_",
#                  "SA.CO.1" => "Aquilegia_",
#                  "SC.NV.4" => "Aquilegia_",
#                  "SH.NV.4" => "Aquilegia_",
#                  "SK.SON.2" => "Aquilegia_",
#                  "TR.UT.3" => "Aquilegia_",
#                  "OW" => "Aquilegia_",
# )

# from Table S1 of https://doi.org/10.1038/nature05857
genera = Dict(
              "BA" => "barnebyi",
              "BR" => "brevistyla",
              "CA" => "canadensis",
              "CHAP" => "chaplinei",
              "COAL" => "coerulea (var. alpina)",
              "COCO" => "coerulea (var. coerulea)",
              "COOC" => "coerulea (var. ochrolueca)",
              "COPI" => "coerulea (var. pinetorum)", # (likely coded as PI in original)
              "DE" => "desertorum",
              "EL" => "elegantula",
              "EX" => "eximia",
              "FL" => "flavescens",
              "FO" => "formosa",
              "HI" => "hinckleyana",
              "JO" => "jonesii",
              "LA" => "laramiensis",
              "LO" => "longissma",
              "MI" => "micrantha",
              "OW" => "NOT USED", # not in original analysis (and not used in this analysis)
              "PU" => "pubescens",
              "SA" => "saximontana",
              "SC" => "scopulorum",
              "SH" => "shockleyi",
              "SK" => "skinneri",
              "TR" => "triternata"
)

descriptive_labels = Dict("taxon" => "taxon",
                          "OR." => "orientation",
                          "SP_L" => "spur_length",
                          "LOG_SP_L" => "log_spur_length",
                          "BL_L" => "blade_length",
                          "SE_L" => "sepal_length",
                          "SP-chroma" => "spur_chroma",
                          "SP-hue" => "spur_hue",
                          "SP-brightness" => "spur_brightness",
                          "BL-chroma" => "blade_chroma",
                          "BL-hue" => "blade_hue",
                          "BL-brightness" => "blade_brightness",
                          "PC1" => "NA",
                          "PC2" => "NA",
                          "Syndrome" => "syndrome",
                          "CHS" => "NA",
                          "CHI" => "NA",
                          "F3H" => "NA",
                          "DFR" => "NA",
                          "ANS" => "NA",
                          "UF3GT" => "NA",
                          "CHI_transform" => "NA",
                          "anthocyanins" => "anthocyanins"
)


# plot_transformed_data(df)

keep_columns = ["taxon",
                "orientation",
                "log_spur_length",
                "blade_length",
                "sepal_length",
                "spur_chroma",
                "spur_hue",
                "spur_brightness",
                "blade_chroma",
                "blade_hue",
                "blade_brightness",
                "anthocyanins",
                "syndrome"]

labels = names(df)

keep_inds = findall([descriptive_labels[x] in keep_columns for x in labels])
df = df[!, keep_inds]

# p = plot_transformed_data(df)

using PhyloNetworks

nexus = read(nexus_path, String)
sans_comments = BEASTDataPrep.remove_comments_nexus(nexus)
tmp_nexus = "tmp.nex"
@assert !isfile(tmp_nexus)
write(tmp_nexus, sans_comments)
trees = readNexusTrees(tmp_nexus)
tree = trees[1]
rm(tmp_nexus)

taxa_codes =
"1 BA,
2 JO,
3 CA.VA.1,
4 CHAP.NM.3,
5 CH.CHI.2,
6 COOC.CO.2,
7 CH.NM.3,
8 COAL.WY.1,
9 COCO.WY.3,
10 COOC.UT.4,
11 COPI.AZ.1,
12 DE.AZ.1,
13 EL.NM.4,
14 EX.CA.2,
15 FL.UT.5,
16 FO.CA.BAL.3,
17 FO.CA.BTR.1,
18 HI.TX.4,
19 BR.ALB.5,
20 LA.WY.2,
21 LO.AZ.2,
22 LO.TX.5,
23 MI.CO.1,
24 PU.CA.PP.3,
25 PU.NV.3,
26 SA.CO.1,
27 SC.NV.4,
28 SH.NV.4,
29 SK.SON.2,
30 TR.UT.3,
31 OW"

taxa_codes = split(taxa_codes, ",\n")
nodes = tree.node
node_names = [node.name for node in nodes]
for code in taxa_codes
    original_label, new_label = split(code)
    tree_ind = findfirst(isequal(original_label), tree.names)
    node_ind = findfirst(isequal(original_label), node_names)
    if isnothing(tree_ind)
        @warn "Could not find taxon $new_label in tree."
    else
        tree.names[tree_ind] = new_label
        nodes[node_ind].name = new_label
    end
end

newick = writeTopology(tree)


df, newick = conform_tree_and_data(df, newick)

CSV.write("aquilegia.csv", df)
write("aquilegia_newick.txt", newick)
