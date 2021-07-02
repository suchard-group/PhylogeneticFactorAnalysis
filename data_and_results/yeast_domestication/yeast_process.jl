using CSV, Plots, DataFrames, Statistics

cd(@__DIR__)
yeast_dir = @__DIR__
# data_dir = raw"C:\Users\gabeh\OneDrive\SD_storage\many_traits\public\integrated_factors\data"
yeast_path = joinpath(yeast_dir, "gallone_s5_values.csv")
meta_path = joinpath(yeast_dir, "gallone_s5_metadata.csv")

df = CSV.read(yeast_path, DataFrame)

n, p = size(df)
eltypes = [eltype(df[!, i]) for i = 1:p]

# reorder phenotypes
#NOTE: YOU HAVE TO DO THIS OR THE WRONG VALUES WILL BE LABELED MISSING
new_order = [14, 15, 16, 11, 12, 13, 1, 2, 3, 4,
          5, 6, 7, 8, 9, 10, 61, 62, 63, 43, 46,
         44, 47, 45, 48, 39, 40, 41, 42, 27, 28,
         29, 49, 50, 51, 52, 53, 54, 55, 56, 57,
         36, 37, 38, 33, 35, 32, 31, 34, 64, 65,
         25, 19, 21, 23, 17, 26, 20, 22, 24, 18,
         30, 68, 71, 76, 77, 78, 69, 70, 67, 73,
         79, 80, 74, 75, 72, 82, 66, 81, 58, 59,
         60]

new_order .+= 2

@assert length(unique(new_order)) == p - 2

new_order = [[1, 2]; new_order]
@assert sort(new_order) == collect(1:length(new_order))

select!(df, new_order)

# missing values simply stored as zeros (and there are real zeros in the data set)
# the missing values were determined by close inspection of Table 3A in Gallone et al. (2016) https://doi.org/10.1016/j.cell.2016.08.020
missing_dict = Dict(1:6 => ["SP007", "WI018", "SA002", "BE100", "BE099"],
                    17:19 => ["BE051"],
                    50:51 => ["BE019", "BE056"],
                    63:76 => ["BE051", "BE036", "BE037", "BE101", "BE035", "SA007", "BE032", "WI017", "BE033", "WI018", "BE034", "SP007"],
                    78:78 => ["BE051", "BE036", "BE037", "BE101", "BE035", "SA007", "BE032", "WI017", "BE033", "WI018", "BE034", "SP007"],
                    79:79 => ["BE036", "BE037", "BE035", "SA007", "BE032", "WI017", "BE033", "WI018", "BE034", "SP007"],
                    82:82 => ["WI004"])

ks = keys(missing_dict)
for k in ks
    taxa = missing_dict[k]
    rows = findall(x -> x in taxa, df[!, 1])
    @assert length(rows) == length(taxa)

    cols = (first(k) + 2):(last(k) + 2)

    for row in rows
        for col in cols
            df[row, col] = NaN
        end
    end
end

mdf = CSV.read(meta_path, DataFrame)
mdf = mdf[!, new_order]
@assert names(mdf) == names(df)

# Need to manually modify two flocculation measurements that are outside of the theoretical range
six_standardized = 2.531
i1 = findfirst(x -> x == "BE050", df[!, 1])
i2 = findfirst(x -> x == "BE062", df[!, 1])
df.Flocculation[i1] = six_standardized
df.Flocculation[i2] = six_standardized



nms = names(df)
select!(df, [[1]; 3:length(nms)])
rename!(df, 1 => "taxon")

CSV.write(joinpath(yeast_dir, "yeast_continuous.csv"), df, quotestrings = true)

# for labelling the plot

pretty_names = ["Sample name",
                "Origin",
                "KCl (500 mM)",
                "KCl (1000 mM)",
                "KCl (1500 mM)",
                "NaCl (250 mM)",
                "NaCl (500 mM)",
                "NaCl (1000 mM)",
                "4°C",
                "16°C",
                "30°C",
                "40°C",
                "Ethanol (5%)",
                "Ethanol (7%)",
                "Ethanol (9%)",
                "Ethanol (10%)",
                "Ethanol (11%)",
                "Ethanol (12%)",
                "Sulfite (1.5 mM)",
                "Sulfite (2.25 mM)",
                "Sulfite (3.0 mM)",
                nms[22],
                nms[23],
                nms[24],
                nms[25],
                nms[26],
                nms[27],
                "Copper (0.075 mM)",
                "Cadmium (0.3 mM)",
                "Cadmium (0.4 mM)",
                "Cadmium (0.5 mM)",
                "Actidione (0.2 mg/L)",
                "Actidione (0.4 mg/L)",
                "Actidione (0.5 mg/L)",
                "Acetic acid (50 mM)",
                "Acetic acid (75 mM)",
                "Acetic acid (100 mM)",
                "Levulinic acid (25 mM)",
                "Levulinic acid (50 mM)",
                "Levulinic acid (75 mM)",
                "Formic acid (25 mM)",
                "Formic acid (50 mM)",
                "Formic acid (75 mM)",
                "Glycerol (20°C)",
                "Melibiose (20°C)",
                "Sorbitol (20°C)",
                "Ethanol (20°C)",
                "Galactose (20°C)",
                "Fructose (20°C)",
                "Sucrose (20°C)",
                "Maltose (20°C)",
                nms[52],
                nms[53],
                "Glucose (10°C)",
                "Ethanol (10°C)",
                "Fructose (10°C)",
                "Sucrose (10°C)",
                "Maltose (10°C)",
                "Glucose (39°C)",
                "Ethanol (39°C)",
                "Fructose (39°C)",
                "Sucrose (39°C)",
                "Maltose (39°C)",
                nms[64],
                nms[65],
                nms[66],
                nms[67],
                nms[68],
                nms[69],
                nms[70],
                nms[71],
                nms[72],
                nms[73],
                nms[74],
                nms[75],
                nms[76],
                nms[77],
                nms[78],
                "4-VG",
                nms[80],
                nms[81],
                "Ethanol production",
                nms[83],
                nms[84]
                ]

pretty_names = pretty_names[3:end]

# CSV.write("yeast_continuous.txt", df, quotestrings = false, delim='\t')

## Save labels

categories = Vector{String}(undef, length(pretty_names))
categories[1:41] .= "environmental stress"
categories[42:62] .= "nutrient stress"
categories[63:78] .= "aroma"
categories[79:end] .= "reproduction"

labels_df = DataFrame(label = string.(names(df))[2:end], pretty = pretty_names, cat = categories)
CSV.write(joinpath(yeast_dir, "yeast_labels.csv"), labels_df, quotestrings = true)
