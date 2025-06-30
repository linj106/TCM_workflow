#Import Python modules
import logging
import sys
import os

#Import PIPER modulestcm_dir
import piper_parseargs
import piper_check_input
import piper_default    
import piper_run
import piper_constraints

###Initiate logger###
logger = logging.getLogger()

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

def run_piper(piper_dir, SCHRODINGER, piper_args = None):
    """ Runs PIPER job and associated tasks. 
    
    Requires:
    - piper_dir: directory of PIPER job (to store info and results of current PIPER job)
    - SCHRODINGER: location of SCHRODINGER installation
    - piper_args: piper specific args (used in args passed into module)"""

    #Logging the start of PIPER module
    logger.info(f'PIPER started. Results and information in {piper_dir}')

    #Parsing and checking the arguments
    if piper_args is None: # if no arguments are passed in, parse and check
        args = piper_parseargs.parse_and_check_args()
    else: # o.w. check args and exit if error
        args = piper_args
        if piper_check_input.check_inputted_args(args):
            sys.exit(0)

    #Building constraints and incorporating file into args
    args = piper_constraints.main(args, piper_dir)

    #Getting default settings (from user if argument passed in)
    default = piper_default.main(PIPER_path, os.path.join(os.getcwd(), args.default))

    #Updating default with arguments to get final input
    params = update_default_w_args(default, args)
    logger.info(f'Final settings of PIPER job: {params}')

    #Getting run command and running PIPER job
    piper_run.piper(args, params, SCHRODINGER, piper_dir)

    #Logging run submission
    logger.info(f'PIPER protein-protein docking started. Results and more information found in {piper_dir}')

def set_up(jobname = 'prot_prot_docking'):
    """ Initial setup of defining paths and loggers if the script is called as a standalaone function. 
    
    Input:
    - output_directory: user-inputted output directory to store .log and results (optional); default is None in which the log and results are stored in a folder 
    called PIPER in the cwd """

    #Get path to user directory
    master_dir = os.getcwd()

    #Get Schrodinger environmental variable
    SCHRODINGER = os.getenv('SCHRODINGER')

    #Get Schrodinger release version. Assumes pathname has version
    #E.g., "/schrodinger/2023-2/"
    schrodinger_version = os.path.basename(SCHRODINGER)

    #Get home path environmental variable
    homepath = os.getenv('HOME')

    #Creating directory for IFD workflow and results in cwd 
    piper_dir = os.path.join(master_dir, jobname)
    os.makedirs(piper_dir, exist_ok=True)
    fh = logging.FileHandler(os.path.join(piper_dir, f'{jobname}.log'))

    #Logger settings
    fh.setFormatter(logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')) #Create writing format for logging
    logger.addHandler(fh) #adding handler 
    logger.setLevel(logging.DEBUG) #Set logger default to debug (all warning lvls allowed)

    return master_dir, piper_dir, SCHRODINGER

if __name__ == '__main__':
    """ If the script is called in by name (as a standalone module), it will define necessary paths and logger info and run the job. """

    # Initially parsing user arguments to find whether an output directory was inputted
    jobname = 'prot_prot_docking'
    for i, arg in enumerate(sys.argv):
        if (arg == "--jobname") and i + 1 < len(sys.argv): # identifying the flag of name
            jobname = sys.argv[i+1] # saving jobname
            break
    
    # Setting up paths and loggers with specific output directory or None
    master_dir, piper_dir, SCHRODINGER = set_up(jobname)

    # Running the piper job 
    run_piper(piper_dir, SCHRODINGER)
