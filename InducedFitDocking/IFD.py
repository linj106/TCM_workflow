#Import Python modules
import logging
import sys
import os

#Import IFD modules
import IFD_find_info
import IFD_parseargs
import IFD_run
import IFD_write_input_file
import IFD_default

###Initiate logger###
logger = logging.getLogger()

#Getting PIPER installation path 
IFD_path = os.path.dirname(__file__)

def update_default_w_args(default, args):
    """ Updates the default arguments with user-inputted arguments to get
    final PIPER settings
    
    Inputs:
    - default: dict containing default args
    - args: user arguments, output of argparse 
    
    Returns: dictionary with final PIPER settings """
    user_input = {k:v for k,v in vars(args).items() if v is not None}
    default.update(user_input)

    return default 

def set_up(jobname = 'induced_fit_docking', output_directory = None):
    """ Initial setup of defining paths and loggers if the script is called as a standalaone function. Builds
    directory and log file based on jobname in current working directory. 
    
    Input:
    - jobname: user-inputted jobname to name directory, files, and hpc job; default is 'induced_fit_docking' """

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
    if output_directory is None:
        ifd_dir = os.path.join(master_dir, jobname)
        os.makedirs(ifd_dir, exist_ok=True)
    else:
        ifd_dir = output_directory
    fh = logging.FileHandler(os.path.join(ifd_dir, f'{jobname}_CLI.log'))

    #Logger settings
    fh.setFormatter(logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')) #Create writing format for logging
    logger.addHandler(fh) #adding handler 
    logger.setLevel(logging.DEBUG) #Set logger default to debug (all warning lvls allowed)
    logger.info(f'IFD started. Results and info can be found in {ifd_dir}')

    return master_dir, ifd_dir, SCHRODINGER

def run_ifd(ifd_dir, SCHRODINGER, ifd_args = None):
    """ Runs IFD job and associated tasks.
    
    Requires:
    - piper_dir: directory of IFD job (to store info and results of current IFD job)
    - SCHRODINGER: location of SCHRODINGER installation
    - ifd_args: ifd specific args (used in args passed into module) """

    # Logging the start of PIPER module
    logger.info(f'Induced Fit Docking started. Results and information in {ifd_dir}')

    # Parsing and checking the arguments
    if ifd_args is None: # no arguments were passed in, parse and check 
        args = IFD_parseargs.parse_and_check_args()
    else: # arguments were passed in, just check the args
        args = ifd_args
        if IFD_parseargs.check_inputted_args(args):
            sys.exit(0)
    
    # Building input file and deleting used arguments from args Namespace
    args, input_file_name = IFD_write_input_file.make_input_file_from_args(args, ifd_dir)

    # Getting default settings (from user if argument passed in which overrides other settings)
    institutional_default_setings = "ADD_HERE"
    user_default_settings = args.default
    default = IFD_default.main(IFD_path, institutional_default_setings, user_default_settings)
    
    # Deleting the default arguments from args to only get cmdline linux ones
    del args.default

    # Update default with user parsed args
    params = update_default_w_args(default, args)
    logger.info(f'Final settings of PIPER job: {params}')

    # IFD run
    IFD_run.ifd(args, params, SCHRODINGER, ifd_dir, input_file_name)

if __name__ == '__main__':
    """ If the script is called in by name (as a standalone module), it will define necessary paths and logger info and run the job. """

    # Initially parsing user arguments to find whether an output directory and jobname was inputted
    jobname = 'prot_prot_docking'
    output_directory = None

    for i, arg in enumerate(sys.argv):
        if (arg == "--jobname") and i + 1 < len(sys.argv): # identifying the flag of name
            jobname = sys.argv[i+1] # saving jobname
        if ((arg == "--output") or (arg == "-o")) and i + 1 < len(sys.argv): # identifying the flag of output directory
            output_directory = os.path.join(os.getcwd(), sys.argv[i+1])
    
    # Setting up paths and loggers with specific output directory or None
    master_dir, ifd_dir, SCHRODINGER = set_up(jobname, output_directory)

    # Running the piper job 
    run_ifd(ifd_dir, SCHRODINGER)