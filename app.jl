module App

using PlotlyBase
using DataFrames
using MultivariateStats, CSV
import Colors
import ColorSchemes

using GenieFramework
@genietools

db_multi = DataFrame()

db_multi[!,"System name"] = ["White Matter", "Telencephalon", "Diencephalon", "Mesencephalon","Metencephalon",
                          "Mylencephalon", "Spinal Cord", "Olfactory Epithelium",  "Hypothalamus",
                         "E. coli", "P. aeruginosa", "S. enterica", "V. cholerae", "Zebrafish Embryo",
                        "D. melongaster", "C. elegans",  "P. mammillata", "Cancer organoid", 
                        "Sphere packing", "1:1:4 ellipsoids", "1:4:4: ellipsoids", "1:2:3: ellipsoids", 
                        "Polydisperse packing","Poisson-Voronoi", "Glassy material","Star positions",  
                        "Diffusion limited aggregation", "A. thaliana","Polyurethane foam", "Fluid foam", "Snowflake yeast"]
                        db_multi[!,"Number of samples"] = vcat(ones(9),5*ones(4),10,14,ones(3),5*ones(6),3,1,5,6,1,1,1)
db_multi[!,"Data source"] = vcat(repeat(["Ref [1]"],9),repeat(["Ref [2]"],4), "Ref [3]","Ref [4]","Ref [5]", "Ref [6]", "Ref [7]", 
    repeat(["Simulated with code from Ref [8]"],5), "Generated for this project" , "Ref [9]", 
    "Ref [10]",  "Generated with code from https://github.com/fogleman/dlaf", "Ref [11]", "Ref [12]", "Ref [13]", "Ref [14]")
db_multi[!,"Number of points per sample (approximate)"] = vcat(5_000,20_000,35_000,12_500,90_000,7_500,150_000,25_000,7_500,18_000*ones(4),14_000,20_000,43_500,3_200,2_100,10_000*ones(6),4096,110_000,10_000,2_000,18_000,21_000,2_000)

#TODO fix zebrafish region numbering
db_multi[!,"Search word"] = ["zebrafish_region_2","zebrafish_region_3", "zebrafish_region_4","zebrafish_region_5","zebrafish_region_6",
                        "zebrafish_region_7","zebrafish_region_8", "zebrafish_region_9","zebrafish_region_10",
                        "ecoli", "pseudomonas", "salmonella", "vibrio",  "zebrafish_embryo",
                        "fly_embryo","worm",  "ascidian", "Guo_organoid", 
                         "PackedSpheres", "PackedEllipses", "PackedMandM", "PackedIreg", "PolySpheres", "PV/PV", 
                         "Glassy", "HYGStarDatabase",  "DLA",  "Plant","foam", "bubbles", "yeast"]

d_mat = CSV.read(joinpath("data", "total_distance_compute.txt"), DataFrame)
                                        
color_dict = Dict("cheng_zebrafish_region_2" =>"#A6CEE2", 
      "cheng_zebrafish_region_3"=> "#2179B4",
      "cheng_zebrafish_region_4"=> "#B4D88B",
      "cheng_zebrafish_region_5"=> "#36A047",
      "cheng_zebrafish_region_6"=> "#F6999A",
      "cheng_zebrafish_region_7"=> "#E21F26",
      "cheng_zebrafish_region_8"=> "#FDBF6F",
      "cheng_zebrafish_region_9"=> "#F57F20",
      "cheng_zebrafish_region_10"=> "#CAB3D6",
      "data/motif/Biofilm/ecoli_Pos11-kde2011.h5"=>"#EC523F", 
      "data/motif/Biofilm/ecoli_Pos14_3-kde2011.h5"=>"#EC523F", 
      "data/motif/Biofilm/ecoli_Pos19-kde2011.h5"=>"#EC523F", 
      "data/motif/Biofilm/ecoli_Pos32-kde2011.h5"=>"#EC523F", 
      "data/motif/Biofilm/ecoli_Pos6_3-kde2011.h5"=>"#EC523F", 
      "data/motif/Biofilm/pseudomonas_Pos11.h5"=>"#40A44A", 
      "data/motif/Biofilm/pseudomonas_Pos28.h5"=>"#40A44A", 
      "data/motif/Biofilm/pseudomonas_Pos41.h5"=>"#40A44A", 
      "data/motif/Biofilm/pseudomonas_Pos47.h5"=>"#40A44A", 
      "data/motif/Biofilm/pseudomonas_Pos5.h5"=>"#40A44A", 
      "data/motif/Biofilm/salmonella_Pos10.h5"=>"#B276B2", 
      "data/motif/Biofilm/salmonella_Pos13.h5"=>"#B276B2", 
      "data/motif/Biofilm/salmonella_Pos16.h5"=>"#B276B2", 
      "data/motif/Biofilm/salmonella_Pos19.h5"=>"#B276B2", 
      "data/motif/Biofilm/salmonella_Pos21.h5"=>"#B276B2", 
      "data/motif/Biofilm/vibrio_Pos11-kdv615.h5"=>"#AB8E30", 
      "data/motif/Biofilm/vibrio_Pos16-kdv615.h5"=>"#AB8E30", 
      "data/motif/Biofilm/vibrio_Pos23-kdv615.h5"=>"#AB8E30", 
      "data/motif/Biofilm/vibrio_Pos27-kdv615.h5"=>"#AB8E30", 
      "data/motif/Biofilm/vibrio_Pos5-kdv615.h5"=>"#AB8E30", 
      "data/motif/Keller2008/zebrafish_embryo_150.h5"=>"#9BD2D9",
      "data/motif/Keller2008/zebrafish_embryo_230.h5"=>"#7ABDD0", 
      "data/motif/Keller2008/zebrafish_embryo_310.h5"=>"#60A5C6", 
      "data/motif/Keller2008/zebrafish_embryo_390.h5"=>"#4D8EBF", 
      "data/motif/Keller2008/zebrafish_embryo_470.h5"=>"#4375B5", 
      "data/motif/Keller2008/zebrafish_embryo_550.h5"=>"#405CA7", 
      "data/motif/Keller2008/zebrafish_embryo_630.h5"=>"#3D448A", 
      "data/motif/Keller2008/zebrafish_embryo_710.h5"=>"#333265", 
      "data/motif/Keller2008/zebrafish_embryo_790.h5"=>"#21203F", 
      "data/motif/Keller2008/zebrafish_embryo_870.h5"=>"#111321", 
      "data/motif/Keller2010/fly_embryo_51.h5"=>"#FCDEB3", 
      "data/motif/Keller2010/fly_embryo_61.h5"=>"#FDD49E", 
      "data/motif/Keller2010/fly_embryo_71.h5"=>"#FDC892", 
      "data/motif/Keller2010/fly_embryo_81.h5"=>"#FBBA84", 
      "data/motif/Keller2010/fly_embryo_91.h5"=>"#F9A470", 
      "data/motif/Keller2010/fly_embryo_101.h5"=>"#F68D5C", 
      "data/motif/Keller2010/fly_embryo_111.h5"=>"#F47A51", 
      "data/motif/Keller2010/fly_embryo_121.h5"=>"#ED6648", 
      "data/motif/Keller2010/fly_embryo_131.h5"=>"#E34C34",
      "data/motif/Keller2010/fly_embryo_141.h5"=>"#D83327", 
      "data/motif/Keller2010/fly_embryo_151.h5"=>"#C42126", 
      "data/motif/Keller2010/fly_embryo_161.h5"=>"#B21F24", 
      "data/motif/Keller2010/fly_embryo_171.h5"=>"#991D20", 
      "data/motif/Keller2010/fly_embryo_181.h5"=>"#7E1416", 
      "worm"=>"#FFDE17", 
      "ascidian"=>"#FDBF6D", 
      "Guo_organoid"=>"#2279B5", 
      "data/motif/Packing/PackedEllipses1.h5"=>"#E21F26", 
      "data/motif/Packing/PackedMandM1.h5"=>"#36A047", 
      "data/motif/Packing/PackedIreg1.h5"=>"#F57F20", 
      "data/motif/Packing/PackedSpheres1.h5"=>"#2179B4", 
      "data/motif/Packing/PolySpheres1.h5"=>"#942768", 
      "data/motif/Packing/PackedEllipses2.h5"=>"#E21F26", 
      "data/motif/Packing/PackedMandM2.h5"=>"#36A047", 
      "data/motif/Packing/PackedIreg2.h5"=>"#F57F20", 
      "data/motif/Packing/PackedSpheres2.h5"=>"#2179B4", 
      "data/motif/Packing/PolySpheres2.h5"=>"#942768", 
      "data/motif/Packing/PackedEllipses3.h5"=>"#E21F26", 
      "data/motif/Packing/PackedMandM3.h5"=>"#36A047", 
      "data/motif/Packing/PackedIreg3.h5"=>"#F57F20", 
      "data/motif/Packing/PackedSpheres3.h5"=>"#2179B4", 
      "data/motif/Packing/PolySpheres3.h5"=>"#942768", 
      "data/motif/Packing/PackedEllipses4.h5"=>"#E21F26", 
      "data/motif/Packing/PackedMandM4.h5"=>"#36A047",
      "data/motif/Packing/PackedIreg4.h5"=>"#F57F20", 
      "data/motif/Packing/PackedSpheres4.h5"=>"#2179B4", 
      "data/motif/Packing/PolySpheres4.h5"=>"#942768", 
      "data/motif/Packing/PackedEllipses5.h5"=>"#E21F26", 
      "data/motif/Packing/PackedMandM5.h5"=>"#36A047", 
      "data/motif/Packing/PackedIreg5.h5"=>"#F57F20", 
      "data/motif/Packing/PackedSpheres5.h5"=>"#2179B4", 
      "data/motif/Packing/PolySpheres5.h5"=>"#942768", 
      "data/motif/PV/PV_1.h5"=>"#942768", 
      "data/motif/PV/PV_2.h5"=>"#942768", 
      "data/motif/PV/PV_3.h5"=>"#942768", 
      "data/motif/PV/PV_4.h5"=>"#942768", 
      "data/motif/PV/PV_5.h5"=>"#942768", 
      "data/motif/Glassy/glass_T_044_1.h5"=>"#E12028", 
      "data/motif/Glassy/glass_T_044_2.h5"=>"#E12028", 
      "data/motif/Glassy/glass_T_044_3.h5"=>"#E12028", 
      "data/motif/HYGStarDatabase/star_positions.h5"=>"#8E9738", 
      "data/motif/DLA/DLA_1.h5"=>"#8E9838", 
      "data/motif/DLA/DLA_2.h5"=>"#8E9838", 
      "data/motif/DLA/DLA_3.h5"=>"#8E9838", 
      "data/motif/DLA/DLA_4.h5"=>"#8E9838", 
      "data/motif/DLA/DLA_5.h5"=>"#8E9838",
      "Plant_1"=>"#ED513F",
      "Plant_2"=>"#ED513F",
      "Plant_4"=>"#ED513F",
      "Plant_13"=>"#ED513F",
      "Plant_15"=>"#ED513F",
      "Plant_18"=>"#ED513F",
      "foam"=>"#e21f26",
      "bubbles"=>"#2179B4",
      "yeast"=>"#a6cee2")


symbol_dict = Dict("White Matter"=>"diamond",
                  "Telencephalon"=>"diamond",
                  "Diencephalon"=>"diamond",
                  "Mesencephalon"=>"diamond",
                  "Metencephalon"=>"diamond",
                  "Mylencephalon"=>"diamond",
                  "Spinal Cord"=>"diamond",
                  "Olfactory Epithelium"=>"diamond",
                  "Hypothalamus"=>"diamond",
                  "E. coli"=>"square", 
                  "P. aeruginosa"=>"square", 
                  "S. enterica"=>"square", 
                  "V. cholerae"=>"square", 
                  "Zebrafish Embryo"=>"circle",
                  "D. melongaster"=>"hexagon", 
                  "C. elegans"=>"triangle-down",  
                  "P. mammillata"=>"triangle-down", 
                  "Cancer organoid"=>"triangle-down", 
                  "Sphere packing"=>"triangle-up-open",
                  "1:1:4 ellipsoids"=>"triangle-up-open",
                  "1:4:4: ellipsoids"=>"triangle-up-open",
                  "1:2:3: ellipsoids"=>"triangle-up-open",
                  "Polydisperse packing"=>"triangle-up-open", 
                  "Glassy material"=>"square-open",  
                  "Diffusion limited aggregation"=>"circle-open",
                  "Star positions"=>"star-open", 
                  "Poisson-Voronoi"=>"square-open",
                  "A. thaliana"=>"triangle-down",
                  "Polyurethane foam"=>"diamond-open", 
                  "Fluid foam"=>"diamond-open", 
                  "Snowflake yeast"=>"star")


const multi_table_options = DataTableOptions(columns = Column(["System name", "Number of samples","Number of points per sample (approximate)","Data source"]))

function replace_names(text_name)
  idx = findfirst(occursin.(db_multi[:,"Search word"],text_name))
  return db_multi[idx,"System name"]
end

function restricted_distance_matrix(ii)
  key_words  = db_multi[ii,"Search word"]
  idx_keep = [any(occursin.(key_words, n)) for n in names(d_mat)]
  text_names = names(d_mat)[idx_keep]
  return d_mat[idx_keep,idx_keep], text_names
end
                                    
function plot_data()
  PlotData( x = (1:10),
            y = (1:10),
            plot = StipplePlotly.Charts.PLOT_TYPE_SCATTER
          )
end

function plot_data_MDS(mds_coord,text_names)
  PlotData(
      x = mds_coord[:,1],
      y = mds_coord[:,2],
      name = "number of casts",
      mode = "markers",
      text = replace_names.(text_names),
      marker = Dict(:color => [color_dict[t] for t in text_names],
            :symbol=>[symbol_dict[replace_names(t)] for t in text_names]),
      plot = StipplePlotly.Charts.PLOT_TYPE_SCATTER
    )
end

function plot_layout(xtitle, ytitle)
  PlotLayout(
    xaxis = [PlotLayoutAxis(xy = "x",title = xtitle,range=[-9, 9])],
    yaxis = [PlotLayoutAxis(xy = "y", title = ytitle, scaleanchor="x",scaleratio=1,range=[-4, 4])]
  )
end

@handlers begin

  @in multisystems_selection = DataTableSelection()

  @out multisystems = DataTable(db_multi,multi_table_options) 
  @out multisystems_pagination = DataTablePagination(rows_per_page=50)
  
  @out MDS12_data = [plot_data()]
  @out MDS12_layout = PlotLayout(plot_bgcolor = "#fff",xaxis = [PlotLayoutAxis(title = "plot loading...")])

  @out MDS23_data = [plot_data()]
  @out MDS23_layout = PlotLayout(plot_bgcolor = "#fff",xaxis = [PlotLayoutAxis(title = "plot loading...")])

  @onchangeany isready, multisystems_selection begin
    ii = union(getindex.(multisystems_selection, "__id"))
    if length(ii) == 0
      ii = 1:size(db_multi)[1]
    end
    d_mat_r, text_names  = restricted_distance_matrix(ii)
    MDS_coords = permutedims(MultivariateStats.transform(MultivariateStats.fit(MDS,
        Matrix(d_mat_r), maxoutdim=3, distances=true)))
    
    MDS12_data = [plot_data_MDS(MDS_coords[:,1:2],text_names)]
    MDS23_data = [plot_data_MDS(MDS_coords[:,3:-1:2],text_names)]
    
    MDS12_layout = plot_layout("MDS PC1", "MDS PC2")
    MDS23_layout = plot_layout("MDS PC3", "MDS PC2")
    
    @show "App is loaded"

  end
end

@page("/", "app.jl.html")

end
