#Import Python modules
import logging
import sys
import os
import shutil
import subprocess
import time
import tarfile
import glob
import argparse

#Import necessary modules for TCM functionality
import tcm_parseargs
import tcm_check_input

#Import necessary modules for parts of workflow (PIPER, IFD, etc.)
from PIPER import PIPER
from InducedFitDocking import IFD

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

# Initially looking at sys.argv to use user input of name to create custom directory
for i,item in enumerate(sys.argv):
    if item in ['-n', '--name']:
        name = sys.argv[i+1]
        break

#Creating directory for TCM workflow and results
tcm_dir = os.path.join(master_dir, f'TernaryComplexModeling_{name}')
os.makedirs(tcm_dir, exist_ok=True)
fh = logging.FileHandler(os.path.join(tcm_dir, f'TernaryComplexModeling_{name}.log'))
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

    args, args_by_group = tcm_parseargs.parse_args() # calling parse_args function which has checks wrapped into the function
    
    logger.info(f'Parsing complete with no fatal errors. The following arguments were recognized: {vars(args)}')

def run_piper(SCHRODINGER, tcm_dir, args, args_by_group):
    """ Runs PIPER protein-protein docking with user-parsed arguments.
    
    Input: 
    SCHRODINGER - path to schrodinger installation
    tcm dir - directory of tcm job
    args - results of user-parsed arugments
    args_by_group - dictionary mapping arguments by group
     
    Return: piper_dir - directory of results from PIPER """

    # making piper directory
    piper_dir = os.path.join(tcm_dir, f'prot_prot_docking_{args.name}')
    os.makedirs(piper_dir, exist_ok=False)

    # entering piper directory
    os.chdir(piper_dir)

    # adding necessary arguments for PIPER (defining the input protein files as receptors and ligands, etc)
    args.receptor_prot = args.cereblon
    args.ligand_prot = args.protein

    # getting all arguments for PIPER + adding custom job name with input naming scheme
    args_piper = argparse.Namespace(**{k: getattr(args,k) for k in args_by_group['piper']})
    args_piper.jobname = f'prot_prot_docking_{args.name}'
    
    # logging the start
    logger.info(f"Initiating PIPER Protein-Protein Docking. Results will be found in {piper_dir}")

    # calling and running PIPER 
    PIPER.run_piper(piper_dir, SCHRODINGER, args_piper)

    logger.info(f"Completed PIPER Protein-Protein Docking. Results found in {piper_dir}")

    return piper_dir

def run_IFD(SCHRODINGER, tcm_dir, piper_dir, args, args_by_group):
    """ Runs Induced Fit Docking with user-parsed arguments.
    
    Input: 
    SCHRODINGER - path to schrodinger installation
    tcm_dir - directory of tcm job
    piper dir - directory of piper job results
    args - user-parsed arguments 
    args_by_group - dictionary mapping arguments by group """

    # making IFD directory
    ifd_dir = os.path.join(tcm_dir, f'InducedFitDocking_{args.name}')
    os.makedirs(ifd_dir, exist_ok=False)

    # entering ifd directory
    os.chdir(ifd_dir)

    # getting all arguments for IFD + adding custom job name with input naming scheme + input files
    args_ifd = argparse.Namespace(**{k: getattr(args,k) for k in args_by_group['ifd']})
    args_ifd.jobname = f'InducedFitDocking_{args.name}'
    args_ifd.proteins = os.path.join(piper_dir, f'prot_prot_docking_{args.name}-out.maegz')
    args_ifd.ligand = args.ligand
    
    # logging the start
    logger.info(f"Initiating Induced Fit Docking. Results will be found in {ifd_dir}")

    # calling and running PIPER 
    IFD.run_ifd(ifd_dir, SCHRODINGER, args_ifd)

    logger.info(f"Completed Induced Fit Docking. Results found in {ifd_dir}")

def main():

    # parsing arguments 
    args, args_by_group = tcm_parseargs.parse_args()

    # running PIPER
    piper_dir = run_piper(SCHRODINGER, tcm_dir, args, args_by_group)

    # running IFD 
    run_IFD(SCHRODINGER, tcm_dir, piper_dir, args, args_by_group)