import os

def write_input_file(input_file_path):
    """ Path to input file (single protein structure or multiple protein structures within a single file) or name of input file if in current working directory. 
    
    Input: 
    - input_file_path: path to input protein file for IFD
    
    Return: string containing input file information for IFD job """

    return f"""INPUT_FILE  {input_file_path}"""


def write_trim_sidechains(binding_site):
    """ Writing the TRIM_SIDECHAINS Stage of the IFD Input File. This stage specifices which side chains should be temporarily mutated to alanine for Glide docking step.
    Includes the default settings which mutate side chains based on Bfactor cutoff. 

    Input: 
    - binding_site: the position number of the co-crystallized ligand in the protein input file (e.g. C:502)
     
    Returns: string containing settings / info for the trim sidechains stage of the IFD job """

    # defining default settings for trimming
    RESIDUES = 'AUTO'
    METHOD = 'BFACTOR'
    DISTANCE_CUTOFF = 5.0
    BFACTOR_CUTOFF = 40.0
    MAX_RESIDUES = 3

    # writing string in correct format
    trim_sidechain_to_write = f"""STAGE TRIM_SIDECHAINS
    RESIDUES {RESIDUES}
    METHOD {METHOD}
    BINDING_SITE {binding_site}
    DISTANCE_CUTOFF {DISTANCE_CUTOFF}
    BFACTOR_CUTOFF {BFACTOR_CUTOFF}
    MAX_RESIDUES {MAX_RESIDUES}"""

    return trim_sidechain_to_write

def write_glide_docking2_one(binding_site, ligand_file, hydrogen_bond_constraints):
    """ Writing the first GLIDE_DOCKING2 Stage. Performs both grid generation and ligand docking with GlideDocking2. 
    Includes many of the default settings for docking and gridgen. 
        
    Input: 
    - binding_site: the position number of the co-crystallized ligand in the protein input file (e.g. C:502)
    - ligand_file: the file path to .mae file containing ligand to dock 
    - hydrogen_bond_constraints: list of tuples containing constraints information on the CRBN protein when docking with ligand; each tuple contains the atom number and 'acceptor' or 'donor' 
                                 (e.g. [(8440, 'aacceptor'), (8479, 'donor'), (8451, ''donor')])
                                 
    Return: string containing info/setting for glide docking to write to .inp file """

    # glide docking strings to join (all the separate strings to then connect at the end)
    glide_docking = []

    # defining default settings for grid generation and docking
    INNERBOX = 10.0
    OUTERBOX = "auto"
    LIGAND_FILE = ligand_file
    LIGANDS_TO_DOCK = "all"
    GRIDGEN_RECEP_CCUT = 0.25
    GRIDGEN_RECEP_VSCALE = 0.70
    GRIDGEN_FORCEFIELD = "OPLS_2005"
    DOCKING_PRECISION = "SP"
    DOCKING_LIG_CCUT = 0.15
    DOCKING_CV_CUTOFF = 100.0
    DOCKING_LIG_VSCALE = 0.50
    DOCKING_POSES_PER_LIG = 20
    DOCKING_FORCEFIELD = "OPLS_2005"
    DOCKING_RINGCONFCUT = 2.5
    DOCKING_AMIDE_MODE = "penal"

    binding_site = f"ligand {binding_site}"

    # writing the default settings to correctly formatted string 
    glide_docking_settings = f"""STAGE GLIDE_DOCKING2
    BINDING_SITE {binding_site}
    INNERBOX {INNERBOX}
    OUTERBOX {OUTERBOX}
    LIGAND_FILE  {LIGAND_FILE}
    LIGANDS_TO_DOCK {LIGANDS_TO_DOCK}
    GRIDGEN_RECEP_CCUT {GRIDGEN_RECEP_CCUT}
    GRIDGEN_RECEP_VSCALE {GRIDGEN_RECEP_VSCALE}
    GRIDGEN_FORCEFIELD {GRIDGEN_FORCEFIELD}
    DOCKING PRECISION {DOCKING_PRECISION}
    DOCKING_LIG_CCUT  {DOCKING_LIG_CCUT}
    DOCKING_CV_CUTOFF  {DOCKING_CV_CUTOFF}
    DOCKING_LIG_VSCALE {DOCKING_LIG_VSCALE}
    DOCKING_POSES_PER_LIG {DOCKING_POSES_PER_LIG}
    DOCKING_FORCEFIELD {DOCKING_FORCEFIELD}
    DOCKING_RINGCONFCUT {DOCKING_RINGCONFCUT}
    DOCKING_AMIDE_MODE {DOCKING_AMIDE_MODE}"""
    glide_docking.append(glide_docking_settings)

    # defining the hbond constraints
    gridgen_constraints = [] # defines constrants for gridgen 
    use_cons = [] # defines constraints to use
    hbond_features = [] #defines chemical patterns to reward/penalize
    num_constraints = len(hydrogen_bond_constraints) #num of constraints

        # iterating through all constraints
    for i, constraint in enumerate(hydrogen_bond_constraints):

        # defining atom description and number for gridgen
        atom_description = f'"atom{i} {constraint[0]}"'
        gridgen_constraints.append(atom_description)

        # getting h bond docking patterns (based on whether it is acceptor or donor)
        feature = get_h_bond_constraint_feature(i+1, constraint[1]) # numbering for this section is 1-indexing
        hbond_features.append(feature)

        # adding info to use_cons
        use_cons.append(f'atom{i}:{i+1}') # recognizing switch to 1-indexing
    
    # formatting strings correctly and adding to glide_docking (to combine later)
    glide_docking.extend(hbond_features)
    glide_docking.append(f'    GRIDGEN_HBOND_CONSTRAINTS {",".join(gridgen_constraints)}')
    glide_docking.append('    DOCKING_[[CONSTRAINT_GROUP:1]]')
    glide_docking.append(f'    DOCKING_USE_CONS {",".join(use_cons)}')
    glide_docking.append(f'    DOCKING_NREQUIRED_CONS {num_constraints}')

    # combining with newline separators and returning
    return '\n'.join(glide_docking)

def write_glide_docking2_two(ligand_file, hydrogen_bond_constraints):
    """ Writing the second GLIDE_DOCKING2 Stage. Performs both grid generation and ligand docking with GlideDocking2. 
    Includes many of the default settings for docking and gridgen. 
    
    NO INPUT OF BINDING_SITE AS DEFAULT NUMBERING SCHEME NOW LABELS LIGAND AS Z:999

    Input: 
    - ligand_file: the file path to .mae file containing ligand to dock 
    - hydrogen_bond_constraints: list of tuples containing constraints information on the CRBN protein when docking with ligand; each tuple contains the atom number and 'acceptor' or 'donor' 
                                 (e.g. [(8440, 'aacceptor'), (8479, 'donor'), (8451, ''donor')])
                                 
    Return: string containing info/setting for glide docking to write to .inp file """

    # glide docking strings to join (all the separate strings to then connect at the end)
    glide_docking = []

    # defining default settings for grid generation and docking
    BINDING_SITE = "ligand Z:999"
    INNERBOX = 10.0
    OUTERBOX = "auto"
    LIGAND_FILE = ligand_file
    LIGANDS_TO_DOCK = "self"
    GRIDGEN_RECEP_CCUT = 0.25
    GRIDGEN_RECEP_VSCALE = 1.00
    GRIDGEN_FORCEFIELD = "OPLS_2005"
    DOCKING_PRECISION = "SP"
    DOCKING_LIG_CCUT = 0.15
    DOCKING_CV_CUTOFF = 0.0
    DOCKING_LIG_VSCALE = 0.80
    DOCKING_POSES_PER_LIG = 1
    DOCKING_FORCEFIELD = "OPLS_2005"
    DOCKING_RINGCONFCUT = 2.5
    DOCKING_AMIDE_MODE = "penal"

    # writing the default settings to correctly formatted string 
    glide_docking_settings = f"""STAGE GLIDE_DOCKING2
    BINDING_SITE {BINDING_SITE}
    INNERBOX {INNERBOX}
    OUTERBOX {OUTERBOX}
    LIGAND_FILE  {LIGAND_FILE}
    LIGANDS_TO_DOCK {LIGANDS_TO_DOCK}
    GRIDGEN_RECEP_CCUT {GRIDGEN_RECEP_CCUT}
    GRIDGEN_RECEP_VSCALE {GRIDGEN_RECEP_VSCALE}
    GRIDGEN_FORCEFIELD {GRIDGEN_FORCEFIELD}
    DOCKING PRECISION {DOCKING_PRECISION}
    DOCKING_LIG_CCUT  {DOCKING_LIG_CCUT}
    DOCKING_CV_CUTOFF  {DOCKING_CV_CUTOFF}
    DOCKING_LIG_VSCALE {DOCKING_LIG_VSCALE}
    DOCKING_POSES_PER_LIG {DOCKING_POSES_PER_LIG}
    DOCKING_FORCEFIELD {DOCKING_FORCEFIELD}
    DOCKING_RINGCONFCUT {DOCKING_RINGCONFCUT}
    DOCKING_AMIDE_MODE {DOCKING_AMIDE_MODE}"""
    glide_docking.append(glide_docking_settings)

    # defining the hbond constraints
    gridgen_constraints = [] # defines constrants for gridgen 
    use_cons = [] # defines constraints to use
    hbond_features = [] #defines chemical patterns to reward/penalize
    num_constraints = len(hydrogen_bond_constraints) #num of constraints

        # iterating through all constraints
    for i, constraint in enumerate(hydrogen_bond_constraints):

        # defining atom description and number for gridgen
        atom_description = f'"atom{i} {constraint[0]}"'
        gridgen_constraints.append(atom_description)

        # getting h bond docking patterns (based on whether it is acceptor or donor)
        feature = get_h_bond_constraint_feature(i+1, constraint[1]) # numbering for this section is 1-indexing
        hbond_features.append(feature)

        # adding info to use_cons
        use_cons.append(f'atom{i}:{i+1}') # recognizing switch to 1-indexing
        
    # formatting strings correctly and adding to glide_docking (to combine later)
    glide_docking.extend(hbond_features)
    glide_docking.append(f'    GRIDGEN_HBOND_CONSTRAINTS {",".join(gridgen_constraints)}')
    glide_docking.append('    DOCKING_[[CONSTRAINT_GROUP:1]]')
    glide_docking.append(f'    DOCKING_USE_CONS {",".join(use_cons)}')
    glide_docking.append(f'    DOCKING_NREQUIRED_CONS {num_constraints}')


    # combining with newline separators and returning
    return '\n'.join(glide_docking)

def get_h_bond_constraint_feature(i, type):
    """ Type refers to the atom type on the receptor protein interacting in H-bond restraint (i.e. carbonyl group on CRBN is acceptor, amino group on CRBN is donor)"""
    assert type in ['acceptor', 'donor']

    if type == 'donor':
        docking_patterns = f"""
    DOCKING_[[FEATURE:{i}]]
        DOCKING_PATTERN1   "[N]#C 1 include"
        DOCKING_PATTERN2   "[n] 1 include"
        DOCKING_PATTERN3   "N(=N=N) 1 include"
        DOCKING_PATTERN4   "N(=N)=N 1 include"
        DOCKING_PATTERN5   "[N;X2]=C[N;X3] 1 include"
        DOCKING_PATTERN6   "[NX2] 1 include"
        DOCKING_PATTERN7   "[N;X1] 1 include"
        DOCKING_PATTERN8   "[N;X2] 1 include"
        DOCKING_PATTERN9   "[N;X3] 1 include"
        DOCKING_PATTERN10   "[#7] 1 include"
        DOCKING_PATTERN11   "[O;X1]~[N;X3]~[O;X1] 1,3 include"
        DOCKING_PATTERN12   "[O;X1;-]C=[O;X1] 1,3 include"
        DOCKING_PATTERN13   "[O-][P,S] 1 include"
        DOCKING_PATTERN14   "[O-] 1 include"
        DOCKING_PATTERN15   "[o] 1 include"
        DOCKING_PATTERN16   "[O;X1] 1 include"
        DOCKING_PATTERN17   "[O;X2] 1 include"
        DOCKING_PATTERN18   "[#8] 1 include"
        DOCKING_PATTERN19   "[S-;X1] 1 include"
        DOCKING_PATTERN20   "[F-] 1 include"
        DOCKING_PATTERN21   "[Cl-] 1 include"
        DOCKING_PATTERN22   "[n;X3] 1 exclude"
        DOCKING_PATTERN23   "[N;X3][*]=[*] 1 exclude"
        DOCKING_PATTERN24   "[N;X3][c,n]~[c,n,o] 1 exclude"
        DOCKING_PATTERN25   "[NX3;+] 1 exclude"
        DOCKING_PATTERN26   "[NH3] 1 exclude"
        DOCKING_PATTERN27   "[NX4;+] 1 exclude"
        DOCKING_PATTERN28   "[s] 1 exclude"
        DOCKING_PATTERN29   "[SX1]=[*] 1 exclude"
        DOCKING_PATTERN30   "[S-]C 1 exclude"
        DOCKING_PATTERN31   "[SX2;H1] 1 exclude"
        DOCKING_PATTERN32   "[SX2] 1 exclude"
        DOCKING_PATTERN33   "S(=O)(=O)([C,N])([C,N]) 1 exclude"
        DOCKING_PATTERN34   "S(=O)([#6])([#6]) 1 exclude"
        DOCKING_PATTERN35   "[SX4] 1 exclude" """

    elif type == 'acceptor':
        docking_patterns = f"""
    DOCKING_[[FEATURE:{i}]]
        DOCKING_PATTERN1   "[#1][#7] 1 include"
        DOCKING_PATTERN2   "[#1][S;X2] 1 include"
        DOCKING_PATTERN3   "[#1][O-] 1 include"
        DOCKING_PATTERN4   "[#1][O;X2] 1 include" """
    
    return docking_patterns

def write_residues_refinement():
    """ Writing the Compile Residue and Prime Refinement stage. Compiles residues that have any atoms that are 5 Angstroms from the ligand which is then optimized 
    and minimized through Prime refinement. 
    
    Return: string containing info/setting for residues refinement to write to .inp file"""

    # writing the stage to compile a list of residues for refinement based on distance cutoff
    DISTANCE_CUTOFF = 5.0
    compile_residues = f"""STAGE COMPILE_RESIDUE_LIST
    DISTANCE_CUTOFF {DISTANCE_CUTOFF}"""
    
    # writing the prime refinement stage
    NUMBER_OF_PASSES = 1
    USE_MEMBRANE = 'no'
    OPLS_VERSION = 'OPLS_2005'
    prime_refinement = f"""STAGE PRIME_REFINEMENT
    NUMBER OF PASSES {NUMBER_OF_PASSES}
    USE_MEMBRANE {USE_MEMBRANE}
    OPLS_VERSION {OPLS_VERSION}"""
    
    # combining and returning 
    return ('\n \n').join([compile_residues, prime_refinement])

def write_sort_and_filter():
    """ Writes the two sort consecutive sort and filter stages. Groups all structures and poses and then sorts by r_psp_Prime_Energy. First stage keeps all poses 
    with r_psp_Prime_Energy value within 30.0 of lowest value and second stage keeps the 20 poses with the lowest r_psp_Prime_Energy value. """

    filter_property = "r_psp_Prime_Energy"
    pose_keep_one = 30.0
    pose_keep_two = 20

    sort_and_filter_one = f"""STAGE SORT_AND_FILTER
    POSE_FILTER {filter_property}
    POSE_KEEP {pose_keep_one}"""
    
    sort_and_filter_two = f"""STAGE SORT_AND_FILTER
    POSE_FILTER {filter_property}
    POSE_KEEP {pose_keep_two}#"""
    
    return ('\n \n'.join([sort_and_filter_one, sort_and_filter_two]))

def write_scoring():
    """ Writes the scoring stage which scores the structure by property. Default is to score by 1 x r_i_glide_gscore and 0.05 x r_psp_Prime_Energy 
    
    Returns: string containg settings / info for the scoring stage """

    SCORE_NAME = "r_psp_IFDScore"
    TERM_ONE_PROP = "r_i_glide_gscore"
    TERM_ONE_WEIGHT = 1.0
    TERM_TWO_PROP = "r_psp_Prime_Energy"
    TERM_TWO_WEIGHT = 0.05
    REPORT_FILE = "report.csv"

    return f"""STAGE SCORING
    SCORE_NAME {SCORE_NAME}
    TERM {TERM_ONE_WEIGHT},{TERM_ONE_PROP},0
    TERM {TERM_TWO_WEIGHT},{TERM_TWO_PROP},1
    REPORT_FILE {REPORT_FILE}"""

def main(input_file_path, binding_site, ligand_file, hydrogen_bond_constraints, directory = os.getcwd()):
    """ Writes an .inp file containing all the necessary stages of protein-protein docking to the specified directory (default is cwd).
    
    Input:
    - input_file_path: path to input protein file for IFD
    - binding_site: the position number of the co-crystallized ligand in the protein input file (e.g. C:502)
    - ligand_file: the file path to .mae file containing ligand to dock
    - hydrogen_bond_constraints: list of tuples containing constraints information on the CRBN protein when docking with ligand; each tuple contains the atom number and 'acceptor' or 'donor' 
                                 (e.g. [(8440, 'acceptor'), (8479, 'donor'), (8451, ''donor')])
    - directory: directory to write .inp file to (default is cwd)"""

    IFD_input = [] # stores all input file information to combine

    # writing all stages in correct order
    IFD_input.append(write_input_file(input_file_path))
    IFD_input.append(write_trim_sidechains(binding_site))
    IFD_input.append(write_glide_docking2_one(binding_site, ligand_file, hydrogen_bond_constraints))
    IFD_input.append(write_residues_refinement())
    IFD_input.append(write_sort_and_filter())
    IFD_input.append(write_glide_docking2_two(ligand_file, hydrogen_bond_constraints))
    IFD_input.append(write_scoring())

    # combine 
    inp = "\n \n".join(IFD_input)

    # write to directory/InducedFit.inp
    with open(f'{directory}/InducedFit.inp', 'w') as file:
        file.write(inp)

    # return path to file
    return os.path.join(directory, 'InducedFit.inp')
if __name__ == "__main__":
    main('InducedFit_5HXB-2_03Angstrom_rec.mae','C:502', 'InducedFit_5HXB-2_03Angstrom.maegz', [(8479, 'donor'), (8440, 'acceptor'), (8451, 'donor')])