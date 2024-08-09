import pandas as pd
from data_manipulation import data_manipulation

# data = data_manipulation("./PI_DataSet.txt")
# # data.filter_muations() # filtering out single mutations
# data.remove_indels()
# # data.filterna_column(["ATV", "IDV"])
# data.count_mutations(5)
# data.explain_input_file
# data.count_drugs()
# data.create_mutfile("ATV", "./drugs/ATV/ATV_mutfile.txt")

# Using aligned pdb files, as these have ligand at the end and all water and stuff like that removed
drug_classes = {
                "ATV": "2fxe_aligned.pdb", 
                # "ATV": "2fxe_aligned_fixed.pdb", 
                "IDV": "1hsg_aligned.pdb", 
                # "SQV": "2nmw.pdb", # Nicholas Raiken investigates this one
                # "APV", # APV not found in columns for PI # Nicholas Raiken investigates this one
                "LPV": "1mui_aligned.pdb", 
                "NFV": "1ohr_aligned.pdb", 
                # "RTV",# RTV not found in columns for PI # Nicholas Raiken investigates this one
                # "TPV": "1d4y.pdb" # Nicholas Raiken investigates this one
}

for drug in drug_classes.keys():
    pdb = drug_classes[drug]
    pdb_path = f"./drugs/{drug}/{pdb}"
    mutfile = f"{drug}_mutfile.txt"
    data = data_manipulation("./PI_DataSet.txt")
    data.filter_muations() # filtering out single mutations
    data.remove_indels()
    data.filterna_column([drug])
    data.count_mutations(5)
    # data.explain_input_file 
    
    output_file = f"./drugs/{drug}/prism_rosetta_XXX_{pdb[:-4]}.txt"
    data.add_output(filename=output_file, pdb_path=pdb_path)
    data.join_data()
    # one plot for all mutations present:
    data.plot_pred_exp(drug, outpath=f"./drugs/{drug}/{drug}_exp_pred_scatter.png")
    # Create plots for each n_muts
    for mut in [2,3,4,5]:
         data.plot_pred_exp(drug, outpath=f"./drugs/{drug}/{drug}_exp_pred_scatter_{mut}_muts.png", n_muts=mut)
    


    data.relax_structure_script(drug=drug,pdb=pdb, bash_location=f"./drugs/{drug}/relax.sh", mutfile=mutfile)
    data.plot_heatmap(drug=drug, n_muts=100,outpath=f"./drugs/{drug}/{drug}_heatmap_absolute.png")
    data.plot_heatmap(drug=drug, n_muts=100,outpath=f"./drugs/{drug}/{drug}_heatmap_percent.png", percent=True)
    data.plot_n_mutation_histogram(drug=drug, outpath=f"./drugs/{drug}/{drug}_histogram_maxmuts.png")
    # This project calls for 5 mutations to be relied upon when running rosetta
    data.filter_n_mutations(5)
    data.create_mutfile(drug, f"./drugs/{drug}/" + mutfile, pdb_path=pdb_path, n_muts=5) # Maximum 5 mutations in each variant
    data.plot_heatmap(drug=drug, n_muts=5,outpath=f"./drugs/{drug}/{drug}_heatmap_absolute_5muts.png")
    data.plot_heatmap(drug=drug, n_muts=5,outpath=f"./drugs/{drug}/{drug}_heatmap_percent_5muts.png", percent=True)
    data.plot_n_mutation_histogram(drug=drug, n_muts=5,outpath=f"./drugs/{drug}/{drug}_histogram_5muts.png")