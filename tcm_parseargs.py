#Import Python modules
import logging
import textwrap
import argparse
import sys
import os

#Import TCM functionality
import tcm_check_input

#Import parsers from other components of workflow
from PIPER import piper_parseargs

###Initiate logger###
logger = logging.getLogger(__name__)

class ArgumentParser(argparse.ArgumentParser):
    def error(self, message):
        logger.critical(message)
        print("An error has occurred. Check log file for information.")
        run_command(('kill','0'))

    """Disables prefix matching in ArgumentParser."""
    def _get_option_tuples(self, option_string):
        """Prevent argument parsing from looking for prefix matches."""
        return []

# Building Parser and Parsing Arguments
def build_parser():
    """ Builds our custom parser based on argparse.ArgumentParser. Defines the specific user inputs and
    potential flags with associated help. Parser arguments are largely derived from available settings to PIPER in Schrodinger
    as seen in https://learn.schrodinger.com/private/edu/release/current/Documentation/html/utilities/program_utility_usage/piper.html?Highlight=PIPER%20command.
    
    Returns: object of class argparse.ArgumentParser with defined user inputs"""

    # creating object of class ArgumentParser with program name and description
    parser = argparse.ArgumentParser(
        prog = 'Protein Preparation Workflow', 
        formatter_class=argparse.RawDescriptionHelpFormatter,
        usage='%(prog)s [options]',
        description=textwrap.dedent('''\
        ----------------------------------------------
        Protein Preparationm Workflow:
            a) Selects and deletes all glycerol molecules --optional
            b) Configures the global settings (e.g. pH, small molecules to process, etc.)
            c) Preprocess (e.g. fill in missing side chains, cap termini, etc.)
            d) Optimize H-bond Assignments 
            e) Minimize and Delete Waters
         '''))
    
    # organizing the parser arguments by groups
    input = parser.add_argument_group('PROTEIN AND LIGAND INPUTS TO FORM TERNARY COMPLEX') # 2 protein (CRBN + POI) + ligand (CeLMod)
    piper = parser.add_argument_group('PIPER protein-protein docking custom settings') # inputs to change settings of PIPER protein-protein docking
    ifd = parser.add_argument_group('Induced Fit Docking (IFD) custom settings') # inputs to change settings of Induced-Fit Docking
    mdfit = parser.add_argument_group('MDFit custom settings') # inputs to change settings of MDFit
    fep = parser.add_argument_group('Free Energy Pertubation (FEP) custom settings') # inputs to change settings of FEP  
    job_control = parser.add_argument_group('JOB CONTROL INFORMATION') #arguments related to job control / server submission (e.g. HOST and JOBNAME)

    # defining subparsers and 
    subparsers = parser.add_subparsers(help = 'Organizing the ')
    # adding specific arguments to our input group
    input.add_argument('-c', '--cereblon', dest = 'cereblon', required = True, help = 'prepared protein file of cereblon; must be .mae or .pdb')
    input.add_argument('-p', '--poi', '--protein', dest = 'protein', required = True, help = 'prepared protein file of target protein / protein of interest; must be .mae or .pdb')
    input.add_argument('-l', '--ligand', dest = 'ligand', required = True, help = 'file of ligand / CELMod; must be .mae, .sdf, .mol2, or .smi')

    return parser

def parse_args(master_dir):
    """ Builds and parse user's inputs"""

    # building parser with defined user inputs
    parser = build_parser() 

    # using parser to parse user inputs
    # collecting known and unknown arguments
    args, unknowns = parser.parse_known_args()

    # changing file paths to full paths (in case we change working directories)
    args.cereblon = os.path.join(master_dir, args.cereblon)
    args.protein = os.path.join(master_dir, args.protein)
    args.ligand = os.path.join(master_dir, args.ligand)

    # adding necessary arguments for PIPER (defining the specific protein files as receptors and ligands, etc)
    args.receptor_prot = args.cereblon
    args.ligand_prot = args.protein
    
    # making sure that the arguments are valid before proceeding
    if tcm_check_input.check_args(parser ,sys.argv, args, unknowns) is True:
        sys.exit(0) # exits because fatal error found
    
    return args
