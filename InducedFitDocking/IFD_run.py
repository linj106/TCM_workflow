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
import IFD_default

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
def build_params_command(params, cmd_line = ['NGLIDECPU','NPRIMECPU', 'NOLOCAL', 'HOST', 'SUBHOST', 'TMPLAUNCHDIR', 'DEBUG']):
    """ Builds terminal commands from params dictionary.

    Input: params dictionary  
    - cmd_line: list of arguments that are related to command line
    Returns: list of run command"""

    cmd = []

    # iterating over all parameters and adding to cmd in correct format
    for k,v in params.items():
        if k not in cmd_line:
            continue
        #True or False parameters convert to the format where flag is included (T) or not included
        #e.g. 'raw':True goes to -raw while 'raw':False has no flag in cmd
        elif v is True:
            cmd.append(f'-{k}')

        elif v is False:
            continue

        #int converts to flag of key with int value as following argument
        #e.g. 'poses':30 goes to '-poses 30' in cmd
        elif isinstance(v, int):
            cmd.append(f'-{k} {v}')

        #str converts to flag of key with str value as following argument
        #e.g. 'ligand_chain': 'A' goes to '-ligand_chain A' in cmd
        #special handling for host information
        elif isinstance(v, str):
            cmd.append(f'-{k} {v}')


    return cmd

def get_inputs(args):
    """ From user-parsed arguments, identifies the input protein poses and ligand file to dock and returns these as full paths. Deletes these files from the args
    
    Input: user-parsed args with .ligand, .proteins
    Returns: tuple of (full_input_protein_path, full_ligand_file_path) """

    protein = args.proteins
    ligand = args.ligand

    return protein, ligand

def make_input_file(args, ifd_dir):
    """Uses arguments in args to create .inp file for IFD. Importantly, certain
    parsed arguments are required for the .inp file and others are cmdline arguments
    that go to linux command. Uses specific arguments for .inp file and deletes the
    ones that are used for input file to only have arguments for linux commands.
    
    Input:
    - args: user-parsed args
    - ifd_dir: directory of ifd results 
    
    Return: 
    - args: modified args with input file arguments deleted
    - input_file_path: name of input file (must enter ifd dir later) """
    
    # getting input protein and ligand from args
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

    # writing input file
    file_out_name = f'{args.jobname}.inp' if args.jobname is not None else 'induced_fit_docking.inp' #filename of .inp defines job name
    if args.template is None: # no user parsed template
        input_file_path = IFD_write_input_file.write_input_file_from_scratch(input_protein_file, ligand_binding_site, input_ligand_file, h_bond_constraints, file_out_path = os.path.join(ifd_dir, file_out_name))
    else: #user parsed template
        template_path = os.path.join(os.getcwd(), args.template)
        input_file_path = IFD_write_input_file.write_input_file_from_template(template_path, input_protein_file, ligand_binding_site, input_ligand_file, h_bond_constraints,  file_out_path = os.path.join(ifd_dir, file_out_name))

    # schrodinger does not allow full file paths for input_file_path, therefore must enter ifd_dir and use name of file
    input_file_path = input_file_path.split('/')[-1]

    # deleting parsed arguments for input file so only cmdline schrodinger IFD-recognized args are left
    del args.proteins
    del args.ligand
    del args.template
    del args.h_bond_constraints
    del args.jobname

    return args, input_file_path


def ifd(args, params, SCHRODINGER, ifd_dir, input_file_name):
    """ Runs Induced Fit Docking job 

    Input:
    - args: original user-parsed arguments
    - params: dict of parsed args that correspond to job settings
    - SCHRODINGER: directory of schrodinger installation
    - ifd_dir: directory containing job results and log 
    - input_file_name: name of .inp file in ifd_dir """

    # get run command
    run_cmd = [f'{SCHRODINGER}/ifd']

    # building command
    params_cmd = build_params_command(params)

    # combine cmd
    command = run_cmd + [input_file_name] + params_cmd

    # log the cmd
    logger.info("Running IFD: %s"%' '.join(command))
    
    # entering ifd_dir to run job 
    os.chdir(ifd_dir)
    run_job(command)
