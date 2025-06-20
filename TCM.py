#Import Python modules
import logging
import sys
import os
import shutil
import subprocess
import time
import tarfile
import glob

#Import necessary modules for TCM functionality
import tcm_parseargs
import tcm_check_input

#Import necessary modules for parts of workflow (PIPER, IFD, etc.)
from PIPER import PIPER

#Get TCM installation path 
TCM_path = os.path.dirname(__file__)

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

#Creating directory for TCM workflow and results
tcm_dir = os.path.join(master_dir, 'TCM')
os.makedirs(tcm_dir, exist_ok=False)
fh = logging.FileHandler(os.path.join(master_dir, 'TCM', 'TCM.log'))
logger.info(f'TCM Workflow started. Results and info can be found in {tcm_dir}')

#Logger settings
fh.setFormatter(logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')) #Create writing format for logging
logger.addHandler(fh) #adding handler 
logger.setLevel(logging.DEBUG) #Set logger default to debug (all warning lvls allowed)

def run_job(command):
    """ Runs specific linux command """
    #Run provided command, joining list with space. Pipe stdout and sdterror to log file
    process = subprocess.run(' '.join(command), stdout=subprocess.PIPE, \
        stderr=subprocess.STDOUT, shell=True, text=True)
    
    #Iterate over sdtout and sdterror
    for line in process.stdout.split('\n'):
        #Ignore blank lines
        if line != "":
            #Ignore ExitStatus
            if "ExitStatus" not in line:
                #Write to log file for debugging
                logger.debug(line)

def parsing_and_checking():
    """ Parsing user arguments to TCM Workflow and Checking Protein / Ligand Inputs """

    #Parsing Arguments and Checking User Inputs 
    logger.info('Parsing Arguments and Checking User Inputs')

    args = tcm_parseargs.parse_args() # calling parse_args function which has checks wrapped into the function
    
    logger.info(f'Parsing complete with no fatal errors. The following arguments were recognized: {vars(args)}')

def run_piper(SCHRODINGER, tcm_dir, args):
    """ Runs PIPER protein-protein docking with user-parsed arguments.
    
    Input: 
    args - results of parsing user arguments
    SCHRODINGER - path to schrodinger installation
    piper dir - directory to store piper job results """

    # making piper directory
    piper_dir = os.path.join(tcm_dir, 'PIPER')
    os.makedirs(tcm_dir, exist_ok=False)

    # entering piper directory
    os.chdir(piper_dir)

    # logging the start
    logger.info(f"Initiating PIPER Protein-Protein Docking. Results will be found in {piper_dir}")

    # calling and running PIPER 
    PIPER.run_piper(piper_dir, SCHRODINGER, args)

    logger.info("Completed PIPER Protein-Protein Docking.")

def main():

    # parsing arguments 
    args = tcm_parseargs.parse_args()

    # running PIPER 
    run_piper(SCHRODINGER, tcm_dir, args)