#Import Python modules
import logging
import textwrap
import argparse
import sys
import os

#Import Schrodinger modules 
from schrodinger.structure import StructureReader, StructureWriter
from schrodinger.application.prepwizard2.diagnostics import get_problems

###Initiate logger###
logger = logging.getLogger(__name__)

def file_type_error(file_path, allowed_file_type):
    """ Given the user input of a file path (or name if in current working directory), a file_type_error 
    arises when the file ext of the input file is not in the allowed file type inputs
    
    Return: boolean (True if file-type-error, False if no error) """

    file_no_path = file_path.split('/')[-1]  #split file input by '/' in case user inputs in file path and retrieves the file name (last item in split list)
    *file_name, file_type = file_no_path.split('.') # splits by '.'; file type is the last thing in split list and file name collects everything else
    if file_type in allowed_file_type:
        return True
    else:
        return False

#functions to check proteins
def protein_file_type_error(protein_file_path):
    """ Given the user input of a protein file path (or name if in current working directory), a protein_file_type_error 
    arises when the file ext of the input file is not in the allowed file type inputs
    
    Return: boolean (True if file-type-error, False if no error) """

    return file_type_error(protein_file_path, allowed_file_type = ['.mae', '.maegz', '.pdb']) # calls file-type-error with allowed file types for inputs

def protein_loading_error(protein_file_path, chain = None):
    """ Given the user input of a protein file path (or name if in current working directory) and optional chain, checks whether the 
    protein file is loadable and the selected protein chain is present (if argument passed). Logs more specific error messages related to loading the protein
    
    Return: boolean (True if error occurs, False if no error)"""

    # tries loading in the protein and catches specific errors
    try:
        protein_chains = set() 
        for structure in StructureReader(protein_file_path): # reads all the structures in file 
            structure_exists = True # structures do exist in the file and are loaded
            for atom in structure.atom: # adds in all chains in the file 
                protein_chains.add(atom.chain)

        if not structure_exists: # no structure found
            logger.critical(f'No structure exists within the protein file.')
            return True
        
        if chain is not None and chain not in protein_chains: # user-inputted chain doesn't exist in the file 
            logger.critical(f'Selected chain does not exist in protein file. Please select another chain from {protein_chains} or upload another file')
            return True
        
    except FileNotFoundError as fnf: # file not found / invalid path
        logger.critical(f'FileNotFoundError: The protein file path is invalid and file was not found: {fnf}. Please check the path {protein_file_path}')
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

def protein_invalid_error(protein_file_path, chain = None):
    """ Given the user input of a protein file path (or name if in current working directory) and optional chains, performs multiple checks to ensure the 
    input is valid. 
    
    Return: boolean (True if protein error arises, False if no error)"""

    # checks valid filetype of protein_file
    if protein_file_type_error(protein_file_path):
        return True

    # checks protein file loadability and existence of chain
    if protein_loading_error(protein_file_path, chain):
        return True

    # checks protein preparation
    if protein_diagnostics_error(protein_file_path):
        return True
    
    #returns false if no errors found
    return False 

#functions to check ligand
def ligand_file_type_error(ligand_file_path):
    """ Given the user input of a ligand file path (or name if in current working directory), a ligand_file_type_error 
    arises when the file ext of the input file is not in the allowed file type inputs
    
    Return: boolean (True if file-type-error, False if no error) """

    return file_type_error(ligand_file_path, allowed_file_type = ['.mae', '.sdf', '.mol2', '.smi', '.maegz']) # calls file-type-error with ligand allowed file-type-errors
    
def ligand_loading_error(ligand_file_path):
    """ Given the user input of a ligand file path (or name if in current working directory), cehsks whether the ligand file is loadable
    
    Return: boolean (True if file-type-error, False if no error) """

    # tries loading in the ligand and catches specific errors
    try:
        for structure in StructureReader(ligand_file_path): # reads all the structures in file 
            structure_exists = True # structures do exist in the file and are loaded

        if not structure_exists: # no structure found
            logger.critical(f'No structure exists within the ligand file.')
            return True
        
    except FileNotFoundError as fnf: # file not found / invalid path
        logger.critical(f'FileNotFoundError: The ligand file path is invalid and file was not found: {fnf}. Please check the path {ligand_file_path}')
        return True
    
    except Exception as e: # other exceptions / error arise in trying to load and check protein file
        logger.critical(f'Unknown error occured in loading ligand: {e}. ')
        return True
    
    # returns false if no errors found 
    return False 

def ligand_diagnostics_error(ligand_file_path):
    """ Given the user input of a ligand file path (or name if in cwd), checks whether the 
    ligand was successfully run through LigPrep and doesn't have any errors.
    
    Return: boolean (True if ligand diagnostic error, False if no error)"""

    return False 
    #MUST ADD 

def ligand_invalid_error(ligand_file_path):
    """ Given the user input of a ligand file path (or name if in cwd), performs multiple checks to ensure the input is valid 
    
    Return: boolean (True if ligand error arises, False if no error) """

    # checks file type
    if ligand_file_type_error(ligand_file_path):
        return True
    
    # checks loadability
    if ligand_loading_error(ligand_file_path):
        return True
    
    # checks ligand preparation
    if ligand_diagnostics_error(ligand_file_path):
        return True

    # returns false if no errors found
    return False

# check parsed arguments 
def check_args(parser, system_arg, args, unknowns):
    """ Given the result of parsing user arguments, checks whether the arguments 
    are valid. Logs specific error messages and returns false if certain fatal issues are found. 
    Otherwise, log warnings and proceeds with runs.
    
    Return: boolean (True if no fatal errors, False if fatal errors)"""
    
    error = False
    
    # checks if no inputs are provided to the script (no args except calling script)
    if len(system_arg) == 1: 
        # print error to logger
        logger.critical('TCM requires protein and ligand inputs for basic functionality')

        # print help to console / terminal window
        parser.print_help()

        return True # fatal error (no arguments to check)
    
    # checks if cereblon protein info is valid
    if protein_invalid_error(args.cereblon) is True:
        logger.critical('Cereblon protein failed checks. See above for more detail.')
        error = True 

    # checks if protein of interest info is valid
    if protein_invalid_error(args.protein) is True:
        logger.critical('Target protein of interest failed checks. See above for more detail.')
        error = True
    
    # checks if ligand is valid
    if ligand_invalid_error(args.ligand):
        logger.critical('Ligand failed checks. See above for more detail.')
        error = True
    
    # checks if unknown inputs exist
    if unknowns:
        #Document warning and ignore variables
        logger.warning('Ignoring unrecognized arguments: %s'%unknowns)
        
    return error 
