#Import Python modules
import logging
import sys
import os
import json

###Initiate logger###
logger = logging.getLogger(__name__)

def read_user_json(user_input_json):
    """ Reads in custom user json file that contains the parameters for PIPER job. This overrides institution default settings but will still be updated 
    with user commandline arguments. Essentially new user default parameters. 
    
    Input: path to json file to read in
    Returns: Python dictionary containing the parameters of PIPER path """

    #Open institution-based parameters json file for reading
    with open(user_input_json, "r") as injson:
        #Put parameters in a list
        parameters = json.load(injson)

    return parameters

def read_institution_json(institutional_parameters):
    """ Reads in json file in PPW path that contains the default parameters for PIPER job.
    
    Input: institutional_parameters - default_parameters.json file to read for PIPER job
    Returns: Python dictionary containing the parameters of PIPER path 
    """
    #Open institution-based parameters json file for reading
    with open(institutional_parameters, "r") as injson:
        #Put parameters in a list
        parameters = json.load(injson)

    #Capture current step
    logger.info("json file found. Institution parameters read in.")
    
    #Write parameters for debugging
    logger.debug("Parameters: %s"%parameters)

    #Return parameters list
    return parameters

def write_default_PIPER(user_dir_path):
    """ Writes default PIPER settings as an example to user's designated path which can be modified with custom 
    settings. 
    
    Input: user directory path to folder to write PIPER_default_params.json to 
    Returns: str statement that confirms the writing of file to specific path 
    """

    #Getting parameters dictionary
    parameters = {
    'poses': 30,
    'rotations': 70000,
    'refinement_protocol':'None',
    'raw': False,
    'OMPI': 10,
    'JOBID': True,
    'use_nonstandard_residue': 'yes',
    'HOST': 'js-hopi-sge-xl-rhel8:10',
    'TMPLAUNCHDIR': True,
    'DEBUG': False
    }

    #Convert dictionary to json object
    json_object = json.dumps(parameters, indent=4)

    #Open user path json file and write 
    with open(os.path.join(user_dir_path, "PIPER_default_params.json"), 'w') as outfile:
        outfile.write(json_object)

    return f'Default PIPER json file written to {user_dir_path}/PIPER_default_params.json'

def get_default_PIPER():
    """ Retrieves the default PIPER settings if no default parameters json file are found. Default settings
    are explicitly written into this Python script

    Returns: dict with default PIPER settings
    """

    # defines default parameters for PIPER job
    parameters = {
    'poses': 30,
    'rotations': 70000,
    'refinement_protocol':'None',
    'raw': False,
    'OMPI': 10,
    'JOBID': True,
    'use_nonstandard_residue': 'yes',
    'HOST': 'js-hopi-sge-xl-rhel8:10',
    'TMPLAUNCHDIR': True,
    'DEBUG': False
    }

    return parameters

def main(PIPER_module_path, PIPER_default_path, user_json_path = None):
    """ Retrieves the correct default settings based on user inputs. 
    
    Input:
    - PIPER_module_path: path to PIPER
    - PIPER_default_path: path to location of default PIPER json
    - user_json_path: path to user's custom json file to use (optional)
    Output: dict with correct PIPER default arguments """

    # checks if user inputted custom default arguments and uses these arguments if it exists
    if user_json_path is not None:
        return read_user_json(user_json_path)

    # otherwise tries to read in institutional arguments
    try:
        return read_institution_json(PIPER_default_path)
    
    # if institutional file not found, then gets default arguments defined in script 
    except FileNotFoundError:
        logger.warning(f'Institutional json file not found. Reading in default arguments from {PIPER_module_path}/piper_default.py.')
    
    return get_default_PIPER()