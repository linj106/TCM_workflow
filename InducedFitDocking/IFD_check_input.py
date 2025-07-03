#Import Python modules
import logging

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
        logger.critical("File type error, protein files must end in .mae, .maegz, or .pdb")
        return True 
    else:
        return False 
        

def protein_loading_error(protein_file, chain):
    """ Given specific protein file and chain to act as ligand or receptor, checks whether the protein file is loadable and the selected protein chain is present.
    Logs more specific error messages related to loading the protein
    
    Return: boolean (True if error occurs, False if no error)"""

    error = False
    # tries loading in the protein and catches specific errors
    try:
        protein_chains = set() 
        for structure in StructureReader(protein_file): # reads all the structures in file 
            structure_exists = True # structures do exist in the file
            for atom in structure.atom: # adds in all chains in the file 
                protein_chains.add(atom.chain)

        if not structure_exists: # no structure found
            logger.critical(f'No protein structure exists within the file.')
            return True # critical error and can't check other conditions
        
        if (chain is not None) and (chain not in protein_chains): # user-inputted chain doesn't exist in the file 
            logger.critical(f'Selected chain does not exist in protein file. Please select another chain from {protein_chains} or upload another file')
            error = True

    except Exception as e: # other exceptions / error arise in trying to load and check protein file
        logger.critical(f'Error occured in loading protein: {e}. ')
        error = True
    
    # returns error boolean
    return error 


def invalid_protein_error(protein_file):
    """ Given specific protein file as receptors in induced fit docking, checks that this protein file is loadable with structures. 
    
    Return: boolean (True if protein error arises, False if no error)"""

    # checks valid filetype of protein_file
    if protein_file_type_error(protein_file):
        return True # critical error and can't check other things

    # checks protein file loadability
    if protein_loading_error(protein_file, chain = None):
        return True
    
# check the input ligand 
def ligand_file_type_error(ligand_file_path):
    """ Given the user input of a ligand file path (or name if in current working directory), a ligand_file_type_error 
    arises when the file ext of the input file is not in the allowed file type inputs
    
    Return: boolean (True if file-type-error, False if no error) """

    if file_type_error(ligand_file_path, allowed_file_type = ['mae', 'sdf', 'mol2', 'smi', 'maegz']): # calls file-type-error with ligand allowed file-type-errors
        logger.critical("File type error, ligand files must end in .mae, .maegz, .sdf, .mol2, or .smi")
        return True 
    else:
        return False 

def invalid_ligand_error(ligand_file_path):
    """ Given inputted ligand file, checks that the file extension is correct.
    
    Return: boolean (True if ligand error arises, False if no error)"""

    if ligand_file_type_error(ligand_file_path):
        return True
    else: 
        return False

def check_parsed_args(parser, system_arg, args, unknowns):
    """ Given the result of parsing user arguments, checks whether the arguments 
    are valid. Logs specific error messages and returns false if certain fatal issues are found. 
    Otherwise, log warnings and proceeds with runs.
    
    Return: boolean (True if no fatal errors, False if fatal errors)"""
    
    # checks if no inputs are provided to the script (no args except calling script)
    if len(system_arg) == 1: 
        # print error to logger
        logger.critical('IFD requires at least two inputs (ligand and protein) to dock for basic functionality')

        # print help to console / terminal window
        parser.print_help()

        return False # fatal error (no proteins to dock)
    
    # checks if unknown inputs exist
    if unknowns:
        #Document warning and ignore variables
        logger.warning('ignoring unrecognized arguments: %s'%unknowns)
    
    error = False

    # checks if receptor protein info is valid
    if invalid_protein_error(args.proteins) is True:
        logger.critical('Receptor protein input failed checks. See above for more detail.')
        error = True
    
    # check if ligand input is valid
    if invalid_ligand_error(args.ligand) is True:
        logger.critical('Ligand input failed checks, See above for more details.')

    return error

def check_inputted_args(args):
    """ Checks inputted args (such as parsed in TCM workflow but specific piper args are passed into piper). 
    Same as check_parsed_args but removes checks on system args and unknowns (these should be handled globally
    under TCM checks).
    
    Return: boolean (True if no fatal errors, False if fatal errors)"""
    
    error = False

    # checks if receptor protein info is valid
    if invalid_protein_error(args.proteins) is True:
        logger.critical('Receptor protein input failed checks. See above for more detail.')
        error = True
    
    # check if ligand input is valid
    if invalid_ligand_error(args.ligand) is True:
        logger.critical('Ligand input failed checks, See above for more details.')

    return error
