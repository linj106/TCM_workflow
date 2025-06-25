#Import Python modules
import logging
import sys
import os

#Import PIPER modules
import piper_parseargs as piper_parseargs
import piper_check_input
import piper_default    
import piper_run
import piper_constraints

###Initiate logger###
logger = logging.getLogger(__name__)

#Get path to user directory
master_dir = os.getcwd()

#Get Schrodinger environmental variable
SCHRODINGER = os.getenv('SCHRODINGER')

#Get Schrodinger release version. Assumes pathname has version
#E.g., "/schrodinger/2023-2/"
schrodinger_version = os.path.basename(SCHRODINGER)

#Get home path environmental variable
homepath = os.getenv('HOME')

#Initiate logger
logger = logging.getLogger()

#Building piper directory and .log file
piper_dir = os.path.join(master_dir, 'PIPER')
os.makedirs(piper_dir, exist_ok=False)
fh = logging.FileHandler(os.path.join(piper_dir, 'PIPER.log'))

#Logger settings
fh.setFormatter(logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')) #Create writing format for logging
logger.addHandler(fh) #adding handler 
logger.setLevel(logging.DEBUG) #Set logger default to debug (all warning lvls allowed)

#Getting PIPER installation path 
PIPER_path = os.path.dirname(__file__)

def update_default_w_args(default, args):
    """ Updates the default arguments with user-inputted arguments to get
    final PIPER settings
    
    Inputs:
    - default: dict containing default args
    - args: user arguments, output of argparse 
    
    Returns: dictionary with final PIPER settings """
    user_input = {k:v for k,v in vars(args).items() if v is not None}
    default.update(user_input)

    # delete input files from parameters; only want to display PIPER settings
    del default['receptor_prot']
    del default['receptor_chain']
    del default['ligand_prot']
    del default['ligand_chain']

    return default 

def piper():
    """ Runs PIPER job and associated tasks. 
    
    Requires:
    - piper_dir: directory of PIPER job (to store info and results of current PIPER job)
    - SCHRODINGER: location of SCHRODINGER installation"""

    #Logging the start of PIPER module
    logger.info(f'PIPER started. Results and information in {piper_dir}')

    #Parsing the arguments and checking for errors
    args = piper_parseargs.parse_and_check_args()

    #If there are constraints, processing those and writing to constraints.json in piper_dir
    args = piper_constraints.main(args, piper_dir)

    #Getting default settings
    default = piper_default.main(PIPER_path)

    #Updating default with arguments to get final input
    params = update_default_w_args(default, args)
    logger.info(f'Final settings of PIPER job: {params}')

    #Getting run command and running PIPER job
    piper_run.piper(args, params, SCHRODINGER)

    #Logging run submission
    logger.info(f'PIPER protein-protein docking started. Results and more information found in {piper_dir}')

if __name__ == '__main__':

    ## Setting up paths and loggers
    piper()
