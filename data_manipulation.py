import pandas as pd
import numpy as np
import math
from Bio.PDB.PDBParser import PDBParser # To parse pdb files
import plotly.express as px
from correlation import correlation # statistics.correlation updated for 3.12.4

AA_MAP = {'ALA': 'A', 
          'ARG': 'R', 
          'ASN': 'N', 
          'ASP': 'D', 
          'CYS': 'C', 
          'GLN': 'Q', 
          'GLU': 'E', 
          'GLY': 'G', 
          'HIS': 'H', 
          'ILE': 'I', 
          'LEU': 'L', 
          'LYS': 'K', 
          'MET': 'M', 
          'PHE': 'F', 
          'PRO': 'P', 
          'SER': 'S', 
          'THR': 'T', 
          'TRP': 'W', 
          'TYR': 'Y', 
          'VAL': 'V'
}


class data_manipulation():
    def __init__(self, txt: str, mut_col: str = "CompMutList", rosetta_path: str = "/groups/sbinlab/software/PRISM_tools/rosetta_stability-v0.2.6/software/rosetta_ddG_pipeline/") -> None:
        assert txt.endswith(".txt"), "Data filepath must be a txt!"
        self.data = pd.read_csv(txt, sep="\t")
        self.mut_col = mut_col
        # Filling in mutations NA with empty string - opting to remove NaN instead
        self.data = self.data.dropna(subset = [self.mut_col]) # Relied on in further analysis
        # self.data[self.mut_col] = self.data[self.mut_col].fillna("")
        # Removal of space between mutations
        self.data[self.mut_col] = self.data[self.mut_col].apply(lambda x: x.replace(" ", ""))
        self.rosetta_path = rosetta_path
        self.out_data = pd.DataFrame()

    @property
    def shape(self) -> tuple[int,int]:
        return self.data.shape
    
    @property
    def columns(self) -> list[str]:
        return self.data.columns

    def sample(self, n = 5) -> pd.DataFrame:
        return self.data.sample(n)
    
    def head(self, n = 5) -> pd.DataFrame:
        return self.data.head(n)
        
    @property
    def explain_input_file(self) -> None:
        """
        Printing how input data should look for reference:
        mutfile.txt
        R A57 G R B57 G L A63 P L B63 P
        L A19 I L B19 I I A93 L I B93 L
        E A35 D E B35 D L A63 P L B63 P
        R A57 K R B57 K V A77 I V B77 I
        L A63 P L B63 P I A72 T I B72 T
        So for example, the first line:
        Mutation R57G on chain A, hence why it says A57.
        After that it's the same mutation, but on chain B - hence why B57.
        """
        print("""######################################################
                Input Data Explanation
mutfile.txt
R A57 G R B57 G L A63 P L B63 P
L A19 I L B19 I I A93 L I B93 L
E A35 D E B35 D L A63 P L B63 P
R A57 K R B57 K V A77 I V B77 I
L A63 P L B63 P I A72 T I B72 T
So for example, the first line:
Mutation R57G on chain A, hence why it says A57.
After that it's the same mutation, but on chain B - hence why B57.
######################################################
        """)

    def get_data(self) -> pd.DataFrame:
        """Return the clinical data"""
        return self.data

    def get_out_data(self) -> pd.DataFrame:
        """Return the output datafframe"""
        return self.out_data

    def filter_muations(self, multiple_muts: bool = True) -> None:
        """
        Filters mutations in the data
        @type multiple_muts: bool
        @param multiple_muts: Whether multiple mutations should be kept (True) or whether to only have single mutations (False)
        @rtype: None
        @returns: Updates data without single/multiple mutations. And updates user on amount of mutations filtered out
        """
        mut_before = self.data.shape[0]
        if multiple_muts:
            self.data = self.data[self.data[self.mut_col].str.contains(",")]
            print(mut_before-self.data.shape[0], "single mutations filtered out")
        else:
            self.data = self.data[~self.data[self.mut_col].str.contains(",")]
            print(mut_before-self.data.shape[0], "multiple mutations filtered out")

    def filterna_column(self, col: str | list[str]) -> None:
        """
        Removes rows with Na in specified columns
        @type col: str | list[str]
        @param col: Name of the columns to remove rows with NaN
        @rtype: None
        @returns: Updates self.data with rows with only valid entries (Not NaN!)
        """
        # Typecast to list to work with it
        rows = self.data.shape[0]
        if isinstance(col, str):
            col = [col]
        self.data = self.data.dropna(subset = col)
        print(rows-self.data.shape[0], "rows removed, as they contained NaN in columns:", col)
        
    def remove_stop_codon_mutations(self) -> None:
        """
        Removing all rows with stop codon mutations
        @rtype: None
        @returns: Updated data without stop codon mutations
        """
        mut_before = self.data.shape[0]
        self.data = self.data[~self.data[self.mut_col].str.contains("\*")]
        print(mut_before-self.data.shape[0], "stop codon mutations filtered out")

    def remove_deletions(self) -> None:
        """
        Removing all rows with deletions
        @rtype: None
        @returns: Updated data without deletions
        """
        mut_before = self.data.shape[0]
        self.data = self.data[~self.data[self.mut_col].str.contains("\~")]
        print(mut_before-self.data.shape[0], "deletions filtered out")

    def remove_insertions(self) -> None:
        """
        Removing all rows with insertions
        @rtype: None
        @returns: Updated data without insertions
        """
        mut_before = self.data.shape[0]
        self.data = self.data[~self.data[self.mut_col].str.contains("\#")]
        print(mut_before-self.data.shape[0], "insertions filtered out")

    def remove_indels(self) -> None:
        """
        Removing all rows with deletions, insertions or stop codons - in that order
        @rtype: None
        @returns: Updated data without deletions, insertions and stop codons
        """
        self.remove_deletions()
        self.remove_insertions()
        self.remove_stop_codon_mutations()


    def __create_n_mutations(self) -> None: 
        """
        Create the n_mutations columns with the amount of mutations present in each column
        """
        def count_commas(x: str) -> int:
            """Count the amount of commas in a string"""
            try:
                l_comma = list(set(x.replace(" ", "").split(","))) # Only unique elements with list -> set -> list                
                return len(l_comma)
            except:
                return np.nan # error
        self.data["n_mutations"] = self.data[self.mut_col].apply(lambda x: count_commas(x))  
    def count_mutations(self, n = 5) -> dict:
        """
        Count amount of the number of mutations - so how many single mutations, how many double mutations etc.
        @type n: int
        @param n: Amount of mutations to print amount of
        @rtype: dict
        @returns: Print out the number of each number of mutations - aswell as a dictionary with keys equal mutation count and values equal amount
        """
        self.__create_n_mutations()
        # print(self.data["n_mutations"].unique())
        # Print amount of 1..5
        assert n > 1, "n must be larger than one"
        # No need in continuing beyond the maximum mutation amount
        if self.data["n_mutations"].max(axis = 0) < n:
            n = int(self.data["n_mutations"].max(axis=0))
        else:
            pass
        amount = {}
        for i in range(1, n + 1):
            amount[i] = self.data[self.data["n_mutations"] == i].shape[0]
            print(f"Amount of rows with {i} mutations: {amount[i]}")
        return amount
    
    def count_drugs(self, drugs: list[str] = ["FPV","ATV","IDV","LPV","NFV","SQV","TPV","DRV"]) -> dict:
        """
        Count the amount of rows with data for each drug
        @type drugs: list[str]
        @param drugs: Name of columns with drug fold to count valid data for
        @rtype: dict
        @returns: Dictionary with the amount of rows for each drug. Also print it to screen
        """
        drug_dict = {}
        for drug in drugs:
            assert drug in self.data.columns, f"{drug} cannot be found in the list of available drugs within the data"
            drug_dict[drug] = self.data.dropna(subset = [drug]).shape[0]
            print(f"{drug_dict[drug]} valid rows for drug: {drug}")
        return drug_dict
    
    def pdb_to_seq(self, pdb_path: str, warnings: bool = False) -> dict:
        """
        Construct a sequence from a pdb file
        @type pdb_path: str
        @param pdb_path: Path to a pdb file
        @type warnings: bool
        @param warnings: Whether warnings should be printed or not. Default is False: They will not be printed
        @rtype: dict
        @returns: A dictionary with the keys being the chains present and the values being the sequence belonging to each chain.
        """
        seq_dict = {}
        pdb = PDBParser(QUIET=True)
        structure = pdb.get_structure('struct', pdb_path)  
        # Getting to each residue in each chain  
        for model in structure: 
            for chain in model:
                print(chain.id)
                seq = ""
                for residue in chain:
                    try:
                        seq += str(AA_MAP[residue.resname])
                    except KeyError:
                        if warnings:
                            print(f"Invalid amino acid: {residue.resname} found in pdb: {pdb_path}")
                        else:
                            pass
                seq_dict[str(chain.id)] = seq
        return seq_dict

    def filter_n_mutations(self, n_muts: int)  -> None:
        """
        Filters the amount of mutations present by given integer
        @type n_muts: int
        @param n_muts: maximum amount of mutations to be present. Filtered with <=
        @returns: Overwrites the dataframe with rows with n_muts amount of mutations or lower
        """    
        self.__create_n_mutations()
        self.data = self.data[self.data["n_mutations"] <= n_muts]

    def create_mutfile(self, drug: str, outpath: str, pdb_path: str,chains: list[str] = ["A", "B"], n_muts: int = 5) -> str:
        """
        Creating the mutation file for the given drug
        @type drug: str
        @param drug: Name of the column containing mutation data to transform. e.g.FPV, ATV, IDV etc.
        @type outpath: str
        @param outpath: Path to store the output .txt file in
        @type pdb_path: str
        @param pdb_path: Path to pdb file to verify mutations line up with sequence
        @type chains: list[str]
        @param chains: List of chain ID to create each mutation at
        @type n_muts: int
        @param n_muts: Amount of mutations maximum allowed to be present for each row
        @rtype: str
        @returns: Create a .txt file of the valid mutations at the output path. Returns the raw filedata as a string.
        """
        assert drug in self.data.columns, f"Drug: {drug} can't be found in the columns of the data"
        assert outpath.endswith(".txt"), f"Output filepath: {outpath} must be written to a .txt file"
        # Amount of mutations in each column:
        self.__create_n_mutations()
        # Sequence to verify position match up with AA
        sequence_dict = self.pdb_to_seq(pdb_path=pdb_path)
        # Valid data for drug
        drug_data = self.data.dropna(subset = [drug])
        # Maximum of n_muts present for the drug
        drug_data = drug_data[drug_data["n_mutations"] <= n_muts]
        print("Printing", drug_data.shape[0], "with maximum", n_muts, "present", "rows of data to", outpath)
        filedata = ""
        errors = 0
        error_pos = {}
        for muts in drug_data[self.mut_col]:
            # Print to file here - or maybe stack data first to avoid too much IO operation
            line = ""
            for mut in muts.strip().split(","):
                # Mutation information
                pos = mut[1:-1]
                wt = mut[0]
                mut_to = mut[-1]
                # print("Pos:", pos, "wt", wt, "mut", mut)
                #if str(mut[-1]) == "X" or str(wt) == "X": # X is unsure what the Amino Acid is! - should probably implement check for checking for AMINO_ACIDS instead
                # More robust check
                if str(mut_to) not in AA_MAP.values():
                    print(f"{mut_to} is not a valid amino acid")
                    break
                for chain in chains:
                    # Chain information
                    muterror = False
                    line += wt + " " + chain + pos + " " + mut_to + " "
                    try:
                        if sequence_dict[chain][int(pos)-1] != wt:
                            print("mut:", mut, "wt:", wt)
                            print(f"Amino acid on chain {chain} on position {pos} is {sequence_dict[chain][int(pos)-1]} and not the supplied: {wt}\nThe neighbouring amino acids -3/+3 is: {sequence_dict[chain][int(pos)-3: int(pos)+3]}")
                            muterror = True
                            if int(pos) in error_pos.keys():
                                error_pos[int(pos)] = error_pos[int(pos)] + [wt]
                            else:
                                error_pos[int(pos)] = [wt]
                            break
                        else:
                            pass
                    except ValueError: # In case typecasting to int not possible
                        # A ValueError could be mutation M36MI - here it is uncertain whether the mutation is to M or I
                        # print(mut)
                        muterror = True
                        break
                if muterror: # Error in mutation
                    errors += 1
                    break
                else:
                    pass
            else: # Only adding data, when there isn't a break!
                line += "\n"
                filedata += line
        with open(outpath, "w") as text_file:
            text_file.write(filedata)
        print(f"Mutation errors encountered:", errors, "out of", len(drug_data[self.mut_col]), "rows")
        print("Positions where errors were encountered:", error_pos)
        return filedata
    
    def relax_structure_script(self, drug: str, pdb: str, bash_location: str, mutfile: str):
        """
        Creating a script for relaxing a pdb structure
        How the relax script should look (Example for ATV):
#!/bin/sh
#SBATCH --job-name=hiv
#SBATCH --time=48:00:00
#SBATCH --mem 5000
#SBATCH --partition=sbinlab_ib
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=fvt270@alumni.ku.dk

# Run this command first before slurm submission
conda activate /groups/sbinlab/software/PRISM_tools/py3_ros_ddG_env

dir_py=/groups/sbinlab/fvt270/PRISM/software/rosetta_ddG_pipeline
dir_run=/groups/sbinlab/fvt270/summer_project/drugs/ATV

#/sbinlab/software/Rosetta_2021_Aug_c7009b3/source/bin/score_jd2.linuxgccrelease -s $dir_run/2fxe_aligned.pdb

### ATV
python3 $dir_py/run_pipeline.py -s $dir_run/2fxe_aligned.pdb -o $dir_run/output_atv -i create -mm mut_file -m $dir_run/ATV_mutfile_two_muts.txt --chainid AB --run_struc AB --overwrite_path True --slurm_partition sbinlab_ib --ligand True --dump_pdb True

### Stability_calculations -i proceed instead of -i create
python3 $dir_py/run_pipeline.py -s $dir_run/2fxe_aligned.pdb -o $dir_run/output_atv -i proceed -mm mut_file -m $dir_run/ATV_mutfile_two_muts.txt --chainid AB --run_struc AB --overwrite_path True --slurm_partition sbinlab_ib --ligand True --dump_pdb True

        """
        assert bash_location.endswith(".sh"), "Bash location must be a path to a bash script"
        # Newline characters are wierd when submitting to slurm...
        slurm_options = f"""#!/bin/sh
#SBATCH --job-name=hiv_{drug}
#SBATCH --time=48:00:00
#SBATCH --mem 5000
#SBATCH --partition=sbinlab_ib

conda activate /groups/sbinlab/software/PRISM_tools/py3_ros_ddG_env

dir_py=/groups/sbinlab/fvt270/PRISM/software/rosetta_ddG_pipeline
dir_run=/groups/sbinlab/fvt270/summer_project/drugs/{drug}
"""
        relax_setup = f"""
### {drug}
python3 $dir_py/run_pipeline.py -s $dir_run/{pdb} -o $dir_run/output_{drug} -i create -mm mut_file -m $dir_run/{mutfile} --chainid AB --run_struc AB --overwrite_path True --slurm_partition sbinlab_ib --ligand True --dump_pdb True

### Stability_calculations -i proceed instead of -i create
python3 $dir_py/run_pipeline.py -s $dir_run/{pdb} -o $dir_run/output_{drug} -i proceed -mm mut_file -m $dir_run/{mutfile} --chainid AB --run_struc AB --overwrite_path True --slurm_partition sbinlab_ib --ligand True --dump_pdb True

### ONLY ddg_calculation
#python3 $dir_py/run_pipeline.py -s $dir_run/{pdb} -o $dir_run/output_{drug} -i ddg_calculation -mm mut_file -m $dir_run/{mutfile} --chainid AB --run_struc AB --overwrite_path True --slurm_partition sbinlab_ib --ligand True --dump_pdb True
"""
        with open(bash_location, "w") as bash_file:
            bash_file.write(slurm_options)
            bash_file.write(relax_setup)

    def plot_heatmap(self, drug: str, n_muts: int = 5, outpath: str | None = None, percent: bool = False) -> None:
        """
        Plots a heatmap with the mutational landscape of the current data.
        @type drug: str
        @param drug: Drug to investigate the heatmap of. The drug must be in the columns of the data.
        @type n_muts: int
        @param n_muts: How many mutations to maximally investigate mutational amount for. Default is 5.
        @type outpath: str | None
        @param outpath: The outptfilepath. If None is given, the heatmap will simplt be shown with fig.show()
        @type percent: bool
        @param percent: if to use absolute counts or percentages. True for percentage, False for absolute numbers.
        @rtype: None
        @returns: Creates a heatmap of the number of mutations to and from each amino acid for the given drug.
        Either shows the plot to the user or saves it to a specified path - given user input.
        """
        # Setting up correct columns and data:
        # Amount of mutations in each column:
        self.__create_n_mutations()
        # Valid data for drugs
        assert drug in self.data.columns, f"{drug} not found in the available columns. Available columns:\n{self.data.columns}"
        drug_data = self.data.dropna(subset = [drug])
        # Maximum of n_muts present for the drug
        drug_data = drug_data[drug_data["n_mutations"] <= n_muts]
        print("Creating heatmap from shape:", drug_data.shape[0], "with a maximum of", n_muts, "mutations present")


        # Constructing heatmap matrix
        # Constructing matrix for each amino acid
        heatmap_matrix = []
        AA_to_NUMBER = {} # Dictionary for assigning number to each amino acid
        for i, AA in enumerate(list(AA_MAP.values())):
            heatmap_matrix.append([0 for j in AA_MAP.values()])
            AA_to_NUMBER[AA] = i
        # Creating number system for Amino Acids instead of letters for indexing

        # Iterating mutations
        for muts in drug_data[self.mut_col]:
            for mut in muts.strip().split(","):
                # Mutation information
                pos = mut[1:-1]
                wt = mut[0]
                mut_to = mut[-1]
                try: # If invalid mutation e.g. V82VA - where it is not certain whether it is V or A
                    i = AA_to_NUMBER[wt]
                    j = AA_to_NUMBER[mut_to]
                    heatmap_matrix[i][j] = heatmap_matrix[i][j] + 1 
                except:
                    continue
        if percent:
            count = 0
            for row in heatmap_matrix:
                count += sum(row)
            new_matrix = heatmap_matrix.copy()
            for i, row in enumerate(heatmap_matrix):
                for j, col in enumerate(row):
                    new_matrix[i][j] = (heatmap_matrix[i][j] / count) * 100
            heatmap_matrix = new_matrix
            fig = px.imshow(heatmap_matrix, 
                            text_auto=True,
                            labels=dict(x="Wild-Type amino acid", y="Mutated amino acid", color="Percent"),
                            x=list(AA_MAP.values()),
                            y=list(AA_MAP.values()),
                            width=1200, 
                            height=1200
                            )
        else:
            fig = px.imshow(heatmap_matrix, 
                            text_auto=True,
                            labels=dict(x="Wild-Type amino acid", y="Mutated amino acid", color="Amount"),
                            x=list(AA_MAP.values()),
                            y=list(AA_MAP.values()),
                            width=1200, 
                            height=1200
                            )
        if outpath is None: # Simply show plot
            fig.show()
        else:
            assert isinstance(outpath, str), "Output filepath must be a valid string"
            fig.write_image(outpath)
            print("Heatmap saved to:", outpath)


    def plot_n_mutation_histogram(self, drug: str, outpath: str | None= None, n_muts: int | None = None) -> None:
        """
        Plots a histogram with the amount of mutations present for mutations with data for column {drug}.
        @type drug: str
        @param drug: Which drug to build histogram from
        @type outpath: str | None
        @param outpath: If None, then call .show() on figure. Else save to path stated by outpath
        @type n_muts: int | None
        @param n_muts: Maximum number of mutations to plot. If None, simply take the maximum number of mutations present.
        @returns: Either saves plot to outpath or call .show() on the figure.
        """
        # Data for drug
        assert drug in self.columns, f"{drug} not found in available columns for data"
        drug_data = self.data.dropna(subset = [drug])
        self.__create_n_mutations()
        # Maximum of n_muts present for the drug
        if n_muts is not None:
            drug_data = drug_data[drug_data["n_mutations"] <= n_muts]
            muts = n_muts
        else:
            muts = drug_data["n_mutations"].max() # don't touch maximum number of muts
        title = f"Mutational landscape for drug: {drug}<br>Maximum number of mutations present: {muts}"
        x_title = f"Number of mutations for drug: {drug}"
        fig = px.histogram(drug_data, x = "n_mutations", title=title, labels={'x':x_title, 'y':'Count'})
        if outpath is None: # Simply show plot
            fig.show()
        else:
            assert isinstance(outpath, str), "Output filepath must be a valid string"
            fig.write_image(outpath)
            print(f"Histogram for drug: {drug} saved to:", outpath)
        

    def __beautify_mutation_data(self, variant: str) -> str | None:
        if "=" in variant:
            return None
        else:
            return variant.replace(":", ",")

    def __standardize_mutation_data(self, variant: str, sequence: str) -> str:
        if "," not in variant:
            return variant
        seqlen = len(sequence)
        data = ""
        for mut in variant.split(","):
            wt = mut[0]
            pos = int(mut[1:-1])
            mut_to = mut[-1]
            if pos-1 > seqlen:
                # Because of the way mutations are handled: 
                # Mutating both chain A and B in our case so position gets assigned uniquely
                continue
            else:
                data += "," + mut
        return data.strip(",")

    def add_output(self, filename: str, pdb_path: str) -> None:
        """
        Add output from filename "prism_rosetta_XXX_PDB.txt" to stats
        @type filename: str
        @param filename: Filename of the output file following the prism_rosetta_XXX_PDB.txt format
        @type pdb_path: str
        @param pdb_path: Path to pdb structure file to get the sequence.
        @rtype: None
        @returns: Overwrites self.out_data with your output data
        """
        output = pd.read_csv(filename, sep="\s", comment='#', engine="python")
        output["muts"] = output["variant"].apply(lambda x: self.__beautify_mutation_data(x))
        sequence = self.pdb_to_seq(pdb_path)["A"]
        output = output[output["muts"].astype(str).ne('None')]
        output["muts"] = output["muts"].apply(lambda x: self.__standardize_mutation_data(variant=x,sequence=sequence))

        self.out_data = output


    def join_data(self):
        """
        Left joins the data stored in self.data and self.out_data into self.out_data
        @rtype: None 
        @returns: Overwrites self.out_data with the joined data
        """
        assert not self.out_data.empty, "No output data has been assigned.\nYou can assign output data with self.add_output(filename: str, pdb_path: str)"
        self.out_data = pd.merge(self.out_data, self.data, left_on="muts", right_on="CompMutList", how="inner")
        def count_commas(x) -> int:
            return int(len(x.split(",")))
        self.out_data["n_mutation"] = self.out_data["CompMutList"].apply(lambda x: count_commas(x))
        print(self.out_data.head())

    def plot_pred_exp(self, drug: str, outpath: str | None = None, n_muts: None | list[int] | int = None, color_dict: dict = {
              '1': 'orange',
              '2': 'gold',
              '3': 'green',
              '4':'red',
              '5':'blue'}, x_range: list = [-5, 20], y_range: list = [-3,5]) -> None:
        """
        Plot log(experimental drug fold resistance) data against predicted ddG
        @type drug: str
        @param drug: Name of drug to plot (column name)

        @type outpath: str | None
        @param outpath: If None, then call .show() on figure. Else save to path stated by outpath

        @type n_muts: None | list[int] | int
        @param n_muts: Amount of or which mutations to plot. If None (default) then plot everything present in self.out_data. 
        If int, use exactly n_mutations in n_muts.
        If list[int], use only number of mutations present in n_muts list.

        @type color_dict: dict
        @param color_dict: Map coloring of n_mutations - to make them discernible

        @type x_range: list[float] | list[int]
        @param x_range: Set the x_range for the plot

        @type y_range: list[float] | list[int]
        @param y_range: Set the y_range for the plot
        
        @rtype: None
        @returns: Either saves plot to outpath or call .show() on the figure.

        """
        assert not self.out_data.empty, "No output data has been assigned.\nYou can assign output data with self.add_output(filename: str, pdb_path: str)"
        assert drug in self.out_data.columns, "Output data and input data have not been joined\nYou can do that with self.add_output(filename: str, pdb_path: str) and then self.join_data()"
        logdrug = "log(" + drug + ")"
        self.out_data[logdrug] = self.out_data[drug].apply(lambda x: math.log(x))
        self.out_data["str_n_mutations"] = self.out_data["n_mutation"].apply(lambda x: str(x))
        self.out_data = self.out_data.sort_values(by=['n_mutation']) # To have consistent ordering
        plot_data = self.out_data.copy()
        if n_muts is None:
            pass
        elif isinstance(n_muts, int):
            plot_data = plot_data[plot_data["n_mutation"] == n_muts]
        elif isinstance(n_muts, list):
            plot_data = plot_data[plot_data["n_mutation"].isin(n_muts)]
        else:
            print("n_muts is neither: None, list[int] or int. Will treat n_muts as if it were None.")    

        fig = px.scatter(data_frame = plot_data, 
                         y = logdrug, 
                         x = "mean_ddG", 
                         color = "str_n_mutations", 
                        color_discrete_map=color_dict, 
                        trendline="ols"
        )
        spearman = round(correlation(plot_data["mean_ddG"], plot_data[logdrug], method="ranked"),4)
        pearson = round(correlation(plot_data["mean_ddG"], plot_data[logdrug], method="linear"),4)
        print(f"Spearman rank correlation: {spearman}\nPearson rank correlation: {pearson}")
        full_fig = fig.full_figure_for_development(warn=False) # To access axis ranges
        x0, x1 = full_fig.layout.xaxis.range
        y0, y1 = full_fig.layout.yaxis.range
        fig.update_layout(
                  yaxis_title='Experimental ' + logdrug,
                  xaxis_title="Rosetta predicted ddG",
                  margin=dict(l=5, r=5, t=40, b=5),
                  legend_title_text='N mutations',
                  title_x=0.5, title_y=0.96,
                  xaxis_range=x_range,
                  yaxis_range=y_range
                  )
        fig.add_shape(
            type='line',
            line=dict(
        color="black",
        width=2,
        dash="dot",
    ),
            x0=x0,
            x1=x1,
            y0=y0,
            y1=y1
        )
        x_sum = abs(x_range[0]) + x_range[1]
        y_sum = abs(y_range[0]) + y_range[1]
        fig.add_annotation(x=x_range[0]+(x_sum/50), # 10% of x-axis placement of text
                           y=y_range[0]+(y_sum*0.95), # 10% from the top y-axis placement of text
            text=f"Pearson Rank Correlation: {pearson}",
            showarrow=False,
            xanchor="left"
            # yshift=10
            )
        fig.add_annotation(x=x_range[0]+(x_sum/50), # 10% of x-axis placement of text
                           y=y_range[0]+(y_sum*0.9), # 15% from the top y-axis placement of text
            text=f"Spearman Rank Correlation: {spearman}",
            showarrow=False,
            xanchor="left"
            # yshift=10
            )
        if outpath is None: # Simply show plot
            fig.show()
        else:
            assert isinstance(outpath, str), "Output filepath must be a valid string"
            fig.write_image(outpath)
            print(f"Experimental and predicted scatter plot for drug: {drug} saved to:", outpath)


# Standard deviations!
# Look into extreme values
# Axis ranges should be similar across plots
# Diagonals doesn't necessarily make sense
# Rosetta can't really handle mutations to Proline!
# Spearman rank correlations