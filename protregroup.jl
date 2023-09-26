proj_info = (id = "sdenolly_canxapms",
             data_ver = "20211125",
             msfolder = "mq20211122")
using Pkg
Pkg.activate(@__DIR__)

using Revise
using CSV, RData, DataFrames

@info "Project '$(proj_info.id)' dataset version=$(proj_info.data_ver)"

const misc_scripts_path = joinpath(base_scripts_path, "misc_jl");
const analysis_path = joinpath(base_analysis_path, proj_info.id)
const data_path = joinpath(analysis_path, "data")
const results_path = joinpath(analysis_path, "results")
const scratch_path = joinpath(analysis_path, "scratch")
const plots_path = joinpath(analysis_path, "plots")

# https://github.com/innatelab/analysis_utils_jl/blob/main/protgroup_assembly.jl
includet(joinpath(misc_scripts_path, "protgroup_assembly.jl"));

peptides_rdata = load(joinpath(data_path, proj_info.msfolder,
                               "$(proj_info.id)_$(proj_info.msfolder)_$(proj_info.data_ver)_peptides.RData"))
peptides_df = peptides_rdata["peptides.df"]
proteins_df = peptides_rdata["proteins.df"]
proteins_df.is_reverse = occursin.(Ref(r"^REV__"), proteins_df.protein_ac)
proteins_df.is_contaminant = occursin.(Ref(r"^CON__"), proteins_df.protein_ac)
# strip (some) isoform annotations from the description
proteins_df.protein_description_orig = copy(proteins_df.protein_description)
proteins_df.protein_description = replace.(replace.(proteins_df.protein_description_orig,
                                           Ref(r"Isoform \S+ of " => "")),
                                           Ref(r", isoforms? \S+$" => ""))

peptide2acs = Dict(r.peptide_id => (Set{String}(split(coalesce(r.protein_acs, r.lead_razor_protein_ac), ';')), r.peptide_rank)
                   for r in eachrow(peptides_df))

protregroups = ProtgroupAssembly.assemble_protgroups(peptide2acs, verbose=true,
                                                     nspec_peptides=2, rspec_peptides=0.25)
protregroups_df = ProtgroupAssembly.dataframe(protregroups,
                            protein_ranks=Dict(r.protein_ac => ProtgroupXMatch.rank_uniprot(r) for r in eachrow(proteins_df)),
                            protgroup_col=:protregroup_id, protein_col=:protein_ac,
                            proteins_info=proteins_df,
                            extra_cols = [:organism, :gene_name, :protein_name, :protein_description])
rename!(protregroups_df, :gene_name => :gene_names, :protein_name => :protein_names,
                         :protein_description => :protein_descriptions)
protregroups_df.is_reverse = occursin.(Ref(r"(^|;)REV__"), protregroups_df.protein_acs)
protregroups_df.is_contaminant = occursin.(Ref(r"(^|;)CON__"), protregroups_df.majority_protein_acs)

CSV.write(joinpath(data_path, proj_info.msfolder,
                   "$(proj_info.id)_$(proj_info.msfolder)_$(proj_info.data_ver)_protregroups.txt"),
          protregroups_df, delim='\t', missingstring="")
