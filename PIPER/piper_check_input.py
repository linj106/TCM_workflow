#Import Python modules
import logging
import re

# #Import Schrodinger modules 
from schrodinger.structure import StructureReader, StructureWriter
from schrodinger.application.prepwizard2.diagnostics import get_problems

###Initiate logger###
logger = logging.getLogger(__name__)

# checking protein file inputs
def file_type_error(file, allowed_file_type = ['mae', 'maegz', 'pdb']):
    """ Given a specific file, checks whether a file_type_error arises where the file is not in the allowed file type inputs to PIPER
    
    Return: boolean (True if file-type-error, False if no error) """

    file_no_path = file.split('/')[-1]  #split file input by '/' in case user inputs in file path and retrieves the file name (last item in split list)
    file_type = file_no_path.split('.')[-1] # splits by '.'; file type is the last thing in split list
    if file_type not in allowed_file_type:
        return True
    else:
        return False

def protein_file_type_error(protein_file_path):
    """ Given the user input of a protein file path (or name if in current working directory), a protein_file_type_error 
    arises when the file ext of the input file is not in the allowed file type inputs
    
    Return: boolean (True if file-type-error, False if no error) """

    if file_type_error(protein_file_path, allowed_file_type = ['mae', 'maegz', 'pdb']): # calls file-type-error with allowed file types for inputs
        logger.critical('File type error, protein files must end in .mae, .maegz, or .pdb')
        return True
    else:
        return False 

def protein_loading_error(protein_file, chain):
    """ Given specific protein file and chain to act as ligand or receptor, checks whether the protein file is loadable and the selected protein chain is present.
    Logs more specific error messages related to loading the protein
    
    Return: boolean (True if error occurs, False if no error)"""

    # tries loading in the protein and catches specific errors
    try:
        protein_chains = set() 
        for structure in StructureReader(protein_file): # reads all the structures in file 
            structure_exists = True # structures do exist in the file
            for atom in structure.atom: # adds in all chains in the file 
                protein_chains.add(atom.chain)

        if not structure_exists: # no structure found
            logger.critical(f'No protein structure exists within the file.')
            return True
        
        if (chain is not None) and (chain not in protein_chains): # user-inputted chain doesn't exist in the file 
            logger.critical(f'Selected chain does not exist in protein file. Please select another chain from {protein_chains} or upload another file')
            return True
    except FileNotFoundError as fnf: # file not found / invalid path
        logger.critical(f'FileNotFoundError: The protein file path is invalid and file was not found: {fnf}. Please check the path {protein_file}')
        return True
    except Exception as e: # other exceptions / error arise in trying to load and check protein file
        logger.critical(f'Unknown error occured in loading protein: {e}. ')
        return True
    
    # returns false if no errors found 
    return False 

def protein_diagnostics_error(protein_file_path):
    """ Given the user input of a protein file path (or name if in current working directory), checks that this protein file was 
    prepared with no issues being found in the protein file. 
    
    Return: boolean (True if issues are found, False if no issues) """

    error = False 
    structures = StructureReader(protein_file_path) # reading in the protein file 

    for st in structures: # iterating through the protein file
        problems = get_problems(st) # getting problems
        valences = problems.invalid_types # checking valence error and raising error if not empty list
        if valences:
            logger.critical('Protein not correctly prepared, issues with valence identified.')
            error = True

        missing = problems.missing # checking missing atom error and raising error if not empty list 
        if missing:
            logger.critical('Protein not correctly prepared, issues with missing atoms identified')
            error = True

        overlapping = problems.overlapping # checking overlapping positions error and raising error if not empty list
        if overlapping:
            logger.critical('Protein not correctly prepared, issues with overlapping positions')
            error = True

        alternates = problems.alternates # checking alternate conformations error and raising error if not empty list 
        if alternates:
            logger.critical('Protein not correctly prepared, issues with alternate conformations')
            error = True 

    return error

def invalid_protein_error(protein_file, chain):
    """ Given specific protein file and chain to act as ligand or receptor, checks that this protein file is loadable and the chain is actually present. Also checks that
    the protein is already prepared by the user. 
    
    Return: boolean (True if protein error arises, False if no error)"""

    # checks valid filetype of protein_file
    if protein_file_type_error(protein_file):
        return True

    # checks protein file loadability and existence of chain
    if protein_loading_error(protein_file, chain):
        return True

    #v checks protein preparation
    if protein_diagnostics_error(protein_file):
        return True
    
    #returns false if no errors found
    return False 

# checking the constraints .txt file inputs

def regex_mismatch(pattern, input):
    """ Given a specific regex pattern and input to check, checks whether a regex_mismatch_error occurs.
    
    Input:
    - pattern: regex pattern
    - input: string to check
    
    Return: boolean (True if pattern does not match, False if it does match) """

    if re.fullmatch(pattern, input):
        return False
    
    else:
        return True

def residue_error(residue_info, protein_type, args):
    """ Checks that the residue is the correct pattern (e.g. HIS375) and checks that the residue exists within the protein file (including part of the specific 
    chain if parsed as argument).
    
    Input:
    - residue_info: str representing residue (AAA37)
    - protein_type: type of protein (either receptor or ligand)
    - args: user parsed arguments
    
    Return: True (if error), False (if no error)""" \
    
    # checking regex pattern of residue_info - string of 3 letters or numbers for residue (either AA or ligand) followed by res num

    if regex_mismatch('^[A-Za-z0-9]{3}\d+$', residue_info):
        logger.critical(f'Residue information of {residue_info} incorrect; must be three letter/num code followed by residue number')
        return True 
    
    if protein_type == 'receptor':
        protein_file = args.receptor_prot
        protein_chain = args.receptor_chain

    elif protein_type == 'ligand':
        protein_file = args.ligand_prot
        protein_chain = args.ligand_chain

    else:
        raise AssertionError("protein type can only be receptor or ligand")
    
    for st in StructureReader(protein_file):
        for chain in st.chain:
            for residue in chain.residue:
                if (residue.pdbres.replace(" ", "") == residue_info[:3]) and (residue.resnum == int(residue_info[3:])): # sometimes pdbres has empty spaces, so remove those
                    if (protein_chain is None) or (residue.chain == protein_chain):
                        return False # residue correctly found
    
    return True # residue not found 
                
def distance_constraint_error(line_num, args, constraint_line_split):
    """ Checks line in constraint file input if it starts with distance. Makes sure that it follows the format:
    "distance" | dmin (float) | dmax (float) | REC_RESIDUE (e.g. "HIS375") | LIG_RESIDUE (e.g. "HIS375")
    
    Input:
    - line_num: line number in txt file
    - args: user-parsed arguments
    - constraint_line_split: split line of constraint file as list (result of file.readline().split())
    
    Return: boolean (True if error, False if no error)"""

    # check dmin and dmax are floats
    try: 
        dmin = float(constraint_line_split[1])
        dmax = float(constraint_line_split[2])
    except ValueError:
        logger.critical(f'Check the distance constraint in line {line_num}. dmin and dmax must be floats')
        return True 
    
    # check that the correct number of arguments exist
    if len(constraint_line_split) != 5:
        logger.critical(f'Check the distance constraint in line {line_num}. Must have exactly five arguments per line.')
        return True

    # check the receptor residue input
    if residue_error(constraint_line_split[3], "receptor", args):
        logger.critical(f'Check the distance constraint in line {line_num}. Receptor residue does not exist in receptor protein + chain')
        return True
    
    # check the ligand residue input
    if residue_error(constraint_line_split[4], "ligand", args):
        logger.critical(f'Check the distance constraint in line {line_num}. Ligand residue does not exist in ligand protein + chain')
        return True
    
    return False

def attraction_constraint_error(line_num, args, constraint_line_split):
    """ Checks line in constraint file input if it starts with attraction. Makes sure that it follows the format:
    "attraction" | attraction_bonus (float btw 0.11 and 0.99) | protein_type (either "receptor" or "ligand") | Residues Involved (separated by space)
    
    Input:
    - line_num: line number in file
    - args: user-parsed arguments
    - constraint_line_split: split line of constraint file as list (result of file.readline().split())
    
    Return: boolean (True if error, False if no error)"""

    # check attraction bonus float btw 0.11 and 0.99
    try: 
        bonus = float(constraint_line_split[1])
        assert (bonus >= 0.11 and bonus <= 0.99)
    except ValueError:
        logger.critical(f'Check the attraction constraint in line {line_num}. bonus must be a float')
        return True 
    except AssertionError:
        logger.critical(f'Check the attraction consrtraint in line {line_num}. bonus must be between 0.11 and 0.99')

    # check that the correct number of arguments exist
    if len(constraint_line_split) < 4:
        logger.critical(f'Check the attraction constraint in line {line_num}. Must have at least 4 arguments.')
        return True
    
    protein_type = constraint_line_split[2]
    
    # check the residues inputs if receptor protein 
    if protein_type == 'receptor':
        protein_file = args.receptor_prot
        for residue in constraint_line_split[3:]:
            if residue_error(residue, protein_file, args):
                logger.critical(f'Check the attraction constraint in line {line_num}. Residue {residue} does not exist in receptor protein + chain')
                return True
    
    # check the ligand residue input
    elif protein_type == 'ligand':
        protein_file = args.ligand_prot
        for residue in constraint_line_split[3:]:
            if residue_error(constraint_line_split[4], protein_file, args):
                logger.critical(f'Check the attraction constraint in line {line_num}. Residue {residue} does not exist in ligand protein + chain')
                return True
    
    # otherwise invalid protein_type
    else:
        logger.critical(f'Check the attraction constraint in line {line_num}. Protein type must be receptor or ligand')
        return True
    
    return False

def repulsion_constraint_error(line_num, args, constraint_line_split):
    """ Checks line in constraint file input if it starts with repulsion. Makes sure that it follows the format:
    "repulsion" | protein_type (either "receptor" or "ligand") | Residues Involved (separated by space)
    
    Input:
    - line_num: line number in txt file 
    - args: user-parsed arguments
    - constraint_line_split: split line of constraint file as list (result of file.readline().split())
    
    Return: boolean (True if error, False if no error)"""

    # check that the correct number of arguments exist
    if len(constraint_line_split) < 3:
        logger.critical(f'Check the repulsion constraint in line {line_num}. Must have at least 3 arguments.')
        return True
    
    protein_type = constraint_line_split[1]
    
    # check the residues inputs if receptor protein 
    if protein_type == 'receptor':
        protein_file = args.receptor_prot
        for residue in constraint_line_split[2:]:
            if residue_error(residue, protein_file, args):
                logger.critical(f'Check the repulsion constraint in line {line_num}. Residue {residue} does not exist in receptor protein + chain')
                return True
    
    # check the ligand residue input
    elif protein_type == 'ligand':
        protein_file = args.ligand_prot
        for residue in constraint_line_split[2:]:
            if residue_error(constraint_line_split[4], protein_file, args):
                logger.critical(f'Check the repulsion constraint in line {line_num}. Residue {residue} does not exist in ligand protein + chain')
                return True
    
    # otherwise invalid protein_type
    else:
        logger.critical(f'Check the repulsion constraint in line {line_num}. Protein type must be receptor or ligand')
        return True
    
    return False

def invalid_constraints_error(args, input_constraints_txt):
    """ Checks whether the input constraints txt file is valid by performing the specific checks above per line. 
    
    Input:
    - args: user-parsed arguments
    - input_constraints_txt: file path to constraints txt file
    
    Return: boolean (True if error, False if no error) """

    # iterating through lines 
    with open(input_constraints_txt, 'r') as f:
        for line_num, line in enumerate(f, start = 1):
            line_split = line.strip().split()
            # if line starts with distance, check against distance
            if line_split[0] == 'distance':
                if distance_constraint_error(line_num, args, line_split):
                    return True
            
            # if line starts with attraction, check against attraction
            elif line_split[0] == 'attraction':
                if attraction_constraint_error(line_num, args, line_split):
                    return True
            
            # if line starts with repulsion, check against repulsion
            elif line_split[0] == 'repulsion':
                if repulsion_constraint_error(line_num, args, line_split):
                    return True

            # else check if comment
            elif line_split[0].startswith('#'):
                continue

            # finally invalid line
            else:
                logger.critical(f'Line {line_num} is invalid; must be comment or start with [distance, attraction, repulsion]')
                return True 

def check_parsed_args(parser, system_arg, args, unknowns):
    """ Given the result of parsing user arguments, checks whether the arguments 
    are valid. Logs specific error messages and returns false if certain fatal issues are found. 
    Otherwise, log warnings and proceeds with runs.
    
    Return: boolean (True if no fatal errors, False if fatal errors)"""
    
    # checks if no inputs are provided to the script (no args except calling script)
    if len(system_arg) == 1: 
        # print error to logger
        logger.critical('PIPER requires protein inputs (one receptor and one ligand) to dock for basic functionality')

        # print help to console / terminal window
        parser.print_help()

        return False # fatal error (no proteins to dock)
    
    # checks if unknown inputs exist
    if unknowns:
        #Document warning and ignore variables
        logger.warning('ignoring unrecognized arguments: %s'%unknowns)
    
    # checks if receptor protein info is valid
    if invalid_protein_error(args.receptor_prot, args.receptor_chain) is True:
        logger.critical('Receptor protein failed checks. See above for more detail.')
        return True

    # checks if ligand protein info is valid
    if invalid_protein_error(args.ligand_prot, args.ligand_chain) is True:
        logger.critical('Ligand protein failed checks. See above for more detail.')
        return True 
    
    # check constraints
    if args.constraint is not None:
        if invalid_constraints_error(args, args.constraint):
            logger.critical('Constraints file is invalid. See above for more details.')
            return True 

    return False 

def check_inputted_args(args):
    """ Checks inputted args (such as parsed in TCM workflow but specific piper args are passed into piper). 
    Same as check_parsed_args but removes checks on system args and unknowns (these should be handled globally
    under TCM checks). """

    # checks if receptor protein info is valid
    if invalid_protein_error(args.receptor_prot, args.receptor_chain) is True:
        logger.critical('Receptor protein failed checks. See above for more detail.')
        return True

    # checks if ligand protein info is valid
    if invalid_protein_error(args.ligand_prot, args.ligand_chain) is True:
        logger.critical('Ligand protein failed checks. See above for more detail.')
        return True 
    
    # check constraints
    if invalid_constraints_error(args, args.constraint):
        logger.critical('Constraints file is invalid. See above for more details.')
        return True 

    return False 

