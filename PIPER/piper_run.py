#Import Python modules
import logging
import subprocess
import textwrap
import argparse
import sys
import os

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

        #str converts to flag of key with str value as following argument
        #e.g. 'ligand_chain': 'A' goes to '-ligand_chain A' in cmd
        #special handling for host information to add OMPI after host 
        elif isinstance(v, str):
            if k == 'HOST':
                cmd.append(f'-{k} {v}:{params['OMPI'] if 'OMPI' in params else 10}')
            else:
                cmd.append(f'-{k} {v}')
    
    return cmd

def get_inputs(args):
    """ From user-parsed arguments, identifies the input receptor and ligand structure and chain information and 
    builds list of correct commands for PIPER job submission.
    
    Input: user-parsed args with .receptor, .receptor_chain, .ligand, .ligand_chain
    Returns: list of correct cmds for the inputs to PIPER job """

    # getting the information
    rec = args.receptor_prot
    rec_chain = args.receptor_chain
    lig = args.ligand_prot
    lig_chain = args.ligand_chain

    return [f'-receptor {rec}', f'-receptor_chain {rec_chain}', f'-ligand {lig}', f'-ligand_chain {lig_chain}']

def piper(args, params, SCHRODINGER, piper_dir):
    """ Runs piper job.
    
    Input:
    - args: original user-parsed arguments
    - params: dict of parsed args that correspond to job settings
    - SCHRODINGER: directory of schrodinger installation """

    #get run command
    run_cmd = [f'{SCHRODINGER}/run -FROM psp piper.py']
    
    #input receptor and ligand information
    input = get_inputs(args)

    #build cmd of params
    params_cmd = build_params_command(params)

    #job name
    job_name = [f'-jobname {args.jobname}' if args.jobname is not None else f'-jobname prot_prot_docking']

    #combine cmd
    command = run_cmd + job_name + input + params_cmd

    #log the cmd
    logger.info("Running PIPER protein-protein docking: %s"%' '.join(command))

    # enter into piper_dir (to store results there) and run 
    os.chdir(piper_dir)
    run_job(command)





            