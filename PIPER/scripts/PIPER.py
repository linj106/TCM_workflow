#Import Python modules
import logging
import sys
import os

#Import PIPER modulestcm_dir
import piper_parseargs as piper_parseargs
import piper_check_input
import piper_default    
import piper_run

###Initiate logger###
logger = logging.getLogger(__name__)

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

def run_piper(piper_dir, SCHRODINGER, args = None):
    """ Runs PIPER job and associated tasks. 
    
    Requires:
    - piper_dir: directory of PIPER job (to store info and results of current PIPER job)
    - SCHRODINGER: location of SCHRODINGER installation
    - args: """

    #Logging the start of PIPER module
    logger.info(f'PIPER started. Results and information in {piper_dir}')

    #Parsing and checking the arguments
    if args is None: 
        args = piper_parseargs.parse_and_check_args()
    
    #Getting default settings
    default = piper_default.main(PIPER_path)

    #Updating default with arguments to get final input
    params = update_default_w_args(default, args)
    logger.info(f'Final settings of PIPER job: {params}')

    #Getting run command and running PIPER job
    piper_run.piper(args, params, SCHRODINGER)

    #Logging run submission
    logger.info(f'PIPER protein-protein docking started. Results and more information found in {piper_dir}')

def set_up(output_directory = None):
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

    if output_directory is None:  #if user did not define an output directory

        #Creating directory for TCM workflow and results
        piper_dir = os.path.join(master_dir, 'PIPER')
        os.makedirs(piper_dir, exist_ok=False)
        fh = logging.FileHandler(os.path.join(piper_dir, 'PIPER.log'))

    else: #user did define output directory
        
        # joins the paths: output directory is either a full path in which os.path.join() discards the master dir
        # or output directory is a file path starting at cwd in which this is appended
        piper_dir = os.path.join(master_dir, output_directory)

        fh = logging.FileHandler(os.path.join(piper_dir, 'PIPER.log'))

    #Logger settings
    fh.setFormatter(logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')) #Create writing format for logging
    logger.addHandler(fh) #adding handler 
    logger.setLevel(logging.DEBUG) #Set logger default to debug (all warning lvls allowed)
    logger.info(f'PIPER started. Results and info can be found in {piper_dir}')

    return master_dir, piper_dir, SCHRODINGER

if __name__ == '__main__':
    """ If the script is called in by name (as a standalone module), it will define necessary paths and logger info and run the job. """

    # Initially parsing user arguments to find whether an output directory was inputted
    output_directory = None
    for i, arg in enumerate(sys.argv):
        if (arg == "--output" or arg == "-o") and i + 1 < len(sys.argv): # identifying the flag of output
            output_directory = sys.argv[i+1] # saving output directory
            break
    
    # Setting up paths and loggers with specific output directory or None
    master_dir, piper_dir, SCHRODINGER = set_up(output_directory)

    # Running the piper job 
    run_piper(piper_dir, SCHRODINGER)
