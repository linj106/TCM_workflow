#Import Python modules
import logging
import sys
import os

#Import IFD modules
import IFD_find_info
import IFD_parseargs
import IFD_run
import IFD_write_input_file

###Initiate logger###
logger = logging.getLogger()

#Getting PIPER installation path 
IFD_path = os.path.dirname(__file__)

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
        ifd_dir = os.path.join(master_dir, 'IFD')
        os.makedirs(ifd_dir, exist_ok=False)
        fh = logging.FileHandler(os.path.join(ifd_dir, 'IFD.log'))

    else: #user did define output directory
        
        # joins the paths: output directory is either a full path in which os.path.join() discards the master dir
        # or output directory is a file path starting at cwd in which this is appended
        ifd_dir = os.path.join(master_dir, output_directory)

        fh = logging.FileHandler(os.path.join(ifd_dir, 'IFD.log'))

    #Logger settings
    fh.setFormatter(logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')) #Create writing format for logging
    logger.addHandler(fh) #adding handler 
    logger.setLevel(logging.DEBUG) #Set logger default to debug (all warning lvls allowed)
    logger.info(f'IFD started. Results and info can be found in {ifd_dir}')

    return master_dir, ifd_dir, SCHRODINGER

def run_ifd(ifd_dir, SCHRODINGER):
    
    # Parsing the arguments
    args = IFD_parseargs.parse_args()
    params = vars(args)

    # IFD run
    IFD_run.ifd(args, params, SCHRODINGER)

if __name__ == '__main__':
    """ If the script is called in by name (as a standalone module), it will define necessary paths and logger info and run the job. """

    # Initially parsing user arguments to find whether an output directory was inputted
    output_directory = None
    for i, arg in enumerate(sys.argv):
        if (arg == "--output" or arg == "-o") and i + 1 < len(sys.argv): # identifying the flag of output
            output_directory = sys.argv[i+1] # saving output directory
            break
    
    # Setting up paths and loggers with specific output directory or None
    master_dir, ifd_dir, SCHRODINGER = set_up(output_directory)

    # Running the piper job 
    run_ifd(ifd_dir, SCHRODINGER)