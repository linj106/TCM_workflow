#Import Python modules
import logging
import textwrap
import argparse
import sys
import os
import re

#Import Schrodinger modules 
from schrodinger.structure import StructureReader, StructureWriter
from schrodinger.application.prepwizard2.diagnostics import get_problems

###Initiate logger###
logger = logging.getLogger(__name__)

# checking protein file inputs
def file_type_error(file, allowed_file_type = ['.mae', '.maegz', '.pdb']):
    """ Given a specific file, checks whether a file_type_error arises where the file is not in the allowed file type inputs to PIPER
    
    Return: boolean (True if file-type-error, False if no error) """

    file_no_path = file.split('/')[-1]  #split file input by '/' in case user inputs in file path and retrieves the file name (last item in split list)
    *file_name, file_type = file_no_path.split('.') # splits by '.'; file type is the last thing in split list and file name collects everything else
    if file_type in allowed_file_type:
        return True
    else:
        return False

def protein_file_type_error(protein_file_path):
    """ Given the user input of a protein file path (or name if in current working directory), a protein_file_type_error 
    arises when the file ext of the input file is not in the allowed file type inputs
    
    Return: boolean (True if file-type-error, False if no error) """

    return file_type_error(protein_file_path, allowed_file_type = ['.mae', '.maegz', '.pdb']) # calls file-type-error with allowed file types for inputs

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
        
        if chain not in protein_chains: # user-inputted chain doesn't exist in the file 
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
    print(protein_file)
    if file_type_error(protein_file):
        return True

    # checks protein file loadability and existence of chain
    if protein_loading_error(protein_file, chain):
        return True

    # checks protein preparation
    # if protein_preparation_error(protein_file):
    #     return True
    
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
        return True 
    else:
        return False
    
# def residue_not_found_error(residues_list):

# def check_distance_constraints(distance_constraint_file):
#     """ Checks distance constraints file and makes sure it follows basic format. 
    
#     Input:
#     - distance_constraint_file: .txt file with distance constraint information
    
#     Return: True or False """

#     # building regex with specific required pattern of single distance_constraint file 

#     residue_regex = r'^[A-Z]{3}\d+$' # requires three capital letters (residue type) followed at least one digit (residue num) 
#     distance_regex = r'^(?:0|0?\.\d+|[1-9]\d*(?:\.\d*)?)$' # requires a non-negative number for distance 

#     # iterating and checking through constraint file 
#     with open(distance_constraint_file, 'r') as file: 
#         for constraint in file: # iterating through lines (each one a separate distance constraint)
#             constraint = constraint.strip()
#             information = constraint.split() 
            
#             # checking that each line has exactly 4 arguments separated by space
#             if len(information) != 4:
#                 logger.critical(f'Check distance constraints txt file; each line must have 4 inputs separated by space')
#                 return True 

#             # defining specific info from line 
#             receptor_residue = information[0]
#             ligand_residue = information[1]
#             dmin = information[2]
#             dmax = information[3]
            
#             # match regex for residues
#             if regex_mismatch(residue_regex, receptor_residue) or regex_mismatch(residue_regex, ligand_residue):
#                 logger.critical(f'Check distance constrants txt file; residues are represented by capital 3 letter amino acid code followed by residue number (e.g. HIS375)')
#                 return True
            
#             # match regex for distances
#             if regex_mismatch(distance_regex, dmin) or regex_mismatch(distance_regex, dmax):
#                 logger.critical(f'Check distanct constraints txt file; dmin and dmax must be non-negative valid numbers')

            
            

def check_args(parser, system_arg, args, unknowns):
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
    
    return False 
