##Import Python modules
import logging
import sys
import os
import json

###Initiate logger###
logger = logging.getLogger(__name__)

def read_user_json(user_input_json):
    """ Reads in custom user json file that contains the parameters for IFD job. This overrides institution default settings but will still be updated 
    with user commandline arguments. Essentially new user default parameters. 
    
    Input: path to json file to read in
    Returns: Python dictionary containing the parameters of IFD path """

    #Open institution-based parameters json file for reading
    with open(user_input_json, "r") as injson:
        #Put parameters in a list
        parameters = json.load(injson)

    return parameters

def read_institution_json(insitutional_parameters):
    """ Reads in json file in IFD path that contains the default parameters for IFD job.
    
    Input: path to institution default paramaters .json file for IFD
    Returns: Python dictionary containing the parameters of IFD path 
    """
    #Open institution-based parameters json file for reading
    with open(insitutional_parameters, "r") as injson:
        #Put parameters in a list
        parameters = json.load(injson)

    #Capture current step
    logger.info("json file found. Institution parameters read in.")
    
    #Write parameters for debugging
    logger.debug("Parameters: %s"%parameters)

    #Return parameters list
    return parameters

def write_default_IFD(user_dir_path):
    """ Writes default IFD settings as an example to user's designated path which can be modified with custom 
    settings. 
    
    Input: user directory path to folder to write IFD_default_params.json to 
    Returns: str statement that confirms the writing of file to specific path 
    """

    #Getting parameters dictionary
    parameters = {
    'NGLIDDECPU': 20,
    'NPRIMECPU': 20, 
    'NOLOCAL': True, 
    'HOST': 'js-hopi-sge-xxl-rhel8',
    'SUBHOST': 'js-hpi-sge-xxl-rhel8',
    'TMPLAUNCHDIR': True}

    #Convert dictionary to json object
    json_object = json.dumps(parameters, indent=4)

    #Open user path json file and write 
    with open(os.path.join(user_dir_path, "IFD_default_params.json"), 'w') as outfile:
        outfile.write(json_object)

    return f'Default IFD json file written to {user_dir_path}/IFD_default_params.json'

def get_default_IFD():
    """ Retrieves the default IFD settings if no default parameters json file are found. Default settings
    are explicitly written into this Python script

    Returns: dict with default IFD settings
    """

    # defines default parameters for PIPER job
    parameters = {
    'NGLIDDECPU': 20,
    'NPRIMECPU': 20, 
    'NOLOCAL': True, 
    'HOST': 'js-hopi-sge-xxl-rhel8',
    'SUBHOST': 'js-hpi-sge-xxl-rhel8',
    'TMPLAUNCHDIR': True}

    return parameters

def main(IFD_module_path, IFD_default_path, user_json_path = None):
    """ Retrieves the correct default settings based on user inputs. 
    
    Input:
    - IFD_module_path: path to IFD scripts installation path
    - IFD_default_path: path to location of IFD default institutional arguments 
    - user_json_path: path to user's custom json file to use (optional)
    Output: dict with correct PIPER default arguments """

    # checks if user inputted custom default arguments and uses these arguments if it exists
    if user_json_path is not None:
        return read_user_json(user_json_path)

    # otherwise tries to read in institutional arguments
    try:
        parameters = read_institution_json(IFD_default_path)
    
    # if institutional file not found, then gets default arguments defined in script 
    except FileNotFoundError:
        logger.warning(f'Institutional json file not found. Reading in default arguments from {IFD_module_path}/IFD_default.py.')
        parameters = get_default_IFD()
    
    return parameters
