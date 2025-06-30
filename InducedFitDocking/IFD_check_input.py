#Import Python modules
import logging

# #Import Schrodinger modules 
from schrodinger.structure import StructureReader, StructureWriter
from schrodinger.application.prepwizard2.diagnostics import get_problems

###Initiate logger###
logger = logging.getLogger(__name__)

# checking protein file inputs
def file_type_error(file, allowed_file_type):
    """ Given a specific file, checks whether a file_type_error arises where the file is not in the allowed file type inputs to IFD
    
    Return: boolean (True if file-type-error, False if no error) """

    file_no_path = file.split('/')[-1]  #split file input by '/' in case user inputs in file path and retrieves the file name (last item in split list)
    *file_name, file_type = file_no_path.split('.') # splits by '.'; file type is the last thing in split list and file name collects everything else
    if file_type in allowed_file_type:
        return True
    else:
        return False
    
# check the input protein poses
def protein_file_type_error(protein_file_path):
    """ Given the user input of a protein file path (or name if in current working directory), a protein_file_type_error 
    arises when the file ext of the input file is not in the allowed file type inputs
    
    Return: boolean (True if file-type-error, False if no error) """

    return file_type_error(protein_file_path, allowed_file_type = ['.mae', '.maegz', '.pdb']) # calls file-type-error with allowed file types for inputs

def protein_loading_error(protein_file):
    """ Given specific protein file, checks whether the protein file is loadable.
    Logs more specific error messages related to loading the protein
    
    Return: boolean (True if error occurs, False if no error)"""

    # tries loading in the protein and catches specific errors
    try:
        protein_chains = set() 
        for structure in StructureReader(protein_file): # reads all the structures in file 
            structure_exists = True # structures do exist in the file

        if not structure_exists: # no structure found
            logger.critical(f'No protein structure exists within the file.')
            return True
        
    except FileNotFoundError as fnf: # file not found / invalid path
        logger.critical(f'FileNotFoundError: The protein file path is invalid and file was not found: {fnf}. Please check the path {protein_file}')
        return True
    except Exception as e: # other exceptions / error arise in trying to load and check protein file
        logger.critical(f'Unknown error occured in loading protein: {e}. ')
        return True
    
    # returns false if no errors found 
    return False 

# check the input ligand 
def ligand_file_type_error(protein_file_path):
    """ Given the user input of a ligand file path (or name if in current working directory), a ligand_file_type_error 
    arises when the file ext of the input file is not in the allowed file type inputs
    
    Return: boolean (True if file-type-error, False if no error) """

    return file_type_error(protein_file_path, allowed_file_type = ['.mae', '.maegz', '.pdb']) # calls file-type-error with allowed file types for inputs

