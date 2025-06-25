#Import Python modules
import logging
import subprocess
import textwrap
import argparse
import sys
import os

#Import IFD modules
import IFD_write_input_file
import IFD_find_info

###Initiate logger###
logger = logging.getLogger(__name__)

def run_job(command):
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

# building run command from params dictionary
def build_params_command(params):
    """ Builds terminal commands from params dictionary.

    Input: params dictionary  
    Returns: list of run command"""

    cmd = []

    # iterating over all parameters and adding to cmd in correct format
    for k,v in params.items():
        print(k,v)
        print(type(v))
        #True or False parameters convert to the format where flag is included (T) or not included
        #e.g. 'raw':True goes to -raw while 'raw':False has no flag in cmd
        if v is True:
            cmd.append(f'-{k}')

        elif v is False:
            continue

        #int converts to flag of key with int value as following argument
        #e.g. 'poses':30 goes to '-poses 30' in cmd
        elif isinstance(v, int):
            cmd.append(f'-{k} {v}')
            print(k, 'int')
        #str converts to flag of key with str value as following argument
        #e.g. 'ligand_chain': 'A' goes to '-ligand_chain A' in cmd
        #special handling for host information
        elif isinstance(v, str):
            cmd.append(f'-{k} {v}')
            print(k, 'string')

    print(cmd)
    return cmd

def get_inputs(args):
    """ From user-parsed arguments, identifies the input protein poses and ligand file to dock and returns these as full paths. Deletes these files from the args
    
    Input: user-parsed args with .ligand, .proteins
    Returns: tuple of (full_input_protein_path, full_ligand_file_path) """

    cwd = os.getcwd() # getting cwd
    protein = args.proteins
    ligand = args.ligands

    # getting full path to protein file
    protein_path = os.path.join(cwd, protein)

    # getting full path to ligand file
    ligand_path = os.path.join(cwd, ligand)

    # parsed input args so deleting
    del args.proteins
    del args.ligands

    return (protein_path, ligand_path)

def ifd(args, params, SCHRODINGER):
    """ Runs Induced Fit Docking job 

    Input:
    - args: original user-parsed arguments
    - params: dict of parsed args that correspond to job settings
    - SCHRODINGER: directory of schrodinger installation"""

    # get run command
    run_cmd = [f'{SCHRODINGER}/ifd']

    # find input 
    input_protein_file, input_ligand_file = get_inputs(args)

    # parsing input protein structure to find W380, H380, and cocrystallized ligand
    results_dictionary = IFD_find_info.parse_structure(input_protein_file)

    # defining binding site
    ligand_binding_site = results_dictionary['lig_id']

    # defining h_bond_constraints (either default based on W380 and H378 or custom user constraints)
    if args.h_bond_constraints is None:
        h_bond_constraints = [(results_dictionary['H378_o'], 'acceptor'), (results_dictionary['H378_h'], 'donor'), (results_dictionary['W380_h'], 'donor')] # defining known h_bond_constraints
    else:
        h_bond_constraints = args.h_bond_constraints
    del args.h_bond_constraints
    
    # writing input file
    input_file_path = IFD_write_input_file.main(input_protein_file, ligand_binding_site, input_ligand_file, h_bond_constraints)
    input_file_path = input_file_path.split('/')[-1]
    # building command
    params_cmd = build_params_command(params)
    print(params_cmd)

    # jobname 
    job_name = [f'-jobname IFD_job']

    # combine cmd
    command = run_cmd + [input_file_path] + params_cmd
    print(command)
    # log the cmd
    logger.info("Running IFD: %s"%' '.join(command))

    run_job(command)
