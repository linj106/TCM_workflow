#Import Python modules
import logging
import textwrap
import argparse
import sys
import os
import subprocess

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
        subprocess.run(['kill','0'])

    """Disables prefix matching in ArgumentParser."""
    def _get_option_tuples(self, option_string):
        """Prevent argument parsing from looking for prefix matches."""
        return []

# Building Parser and Parsing Arguments
def build_parser():
    """ Builds our custom parser based on argparse.ArgumentParser. Defines the specific user inputs and
    potential flags with associated help.
    
    Returns: 
    - object of class argparse.ArgumentParser with defined user inputs
    - dictionary of which dest correspond to which stage of the workflow (e.g. input, PIPER, etc.) """

    # creating object of class ArgumentParser with program name and description
    parser = argparse.ArgumentParser(
        prog = 'Ternary Complex Modeling Workflow', 
        formatter_class=argparse.RawDescriptionHelpFormatter,
        usage='%(prog)s [options]',
        description=textwrap.dedent('''\
        ----------------------------------------------
        Ternary Complex Modelling Workflow:
                                    
        Predicts ternary complexes using Schrodinger tools. Requires input 
        of prepared cereblon (cocrystallized with IMiD ligand) and protein
        of interest (POI). Workflow consists of the following modules.
                                    
            a) PIPER protein-protein docking between receptor (cereblon) and ligand (POI) proteins
            b) Induced Fit Docking of generated PIPER poses with ligand
            c) MDFit analysis of stability of protein complexes and ligand
            d) FEP+ calculations using library of CELMoDs to evaluate predicted ternary complex
        ----------------------------------------------
         '''))
    
    # argument group dictionary 
    args_by_group = {}

    # organizing the parser arguments by groups
    input = parser.add_argument_group('PROTEIN AND LIGAND INPUTS TO FORM TERNARY COMPLEX') # 2 protein (CRBN + POI) + ligand (CeLMod)
    piper = parser.add_argument_group('PIPER protein-protein docking custom settings') # inputs to change settings of PIPER protein-protein docking
    ifd = parser.add_argument_group('Induced Fit Docking (IFD) custom settings') # inputs to change settings of Induced-Fit Docking
    mdfit = parser.add_argument_group('MDFit custom settings') # inputs to change settings of MDFit
    fep = parser.add_argument_group('Free Energy Pertubation (FEP) custom settings') # inputs to change settings of FEP  
    job_control = parser.add_argument_group('JOB CONTROL INFORMATION') #arguments related to job control / server submission (e.g. HOST and JOBNAME)

    # adding specific arguments to our input group
    input.add_argument('-c', '--cereblon', dest = 'cereblon', type = str, required = True, help = 'prepared protein file of cereblon; must be .mae or .pdb')
    input.add_argument('-p', '--poi', '--protein', dest = 'protein', type = str, required = True, help = 'prepared protein file of target protein / protein of interest; must be .mae or .pdb')
    input.add_argument('-l', '--ligand', dest = 'ligand', type = str, required = True, help = 'file of ligand / CELMod; must be .mae, .sdf, .mol2, or .smi')
    input.add_argument('-n', '--name', dest = 'name', type = str, required = True, help = 'naming scheme to identify job by (e.g. by protein of interest like 5HXB-2)')

    # adding specific argument into piper group
    piper.add_argument('--piper_settings', dest = 'piper_settings', type = str, required = True, help = 'path to json file containing settings to apply to piper job')
    
    # adding specific argument into IFD group
    ifd.add_argument('--ifd_settings', dest = 'ifd_settings', type = str, required = True, help = 'path to json file containing settings to apply to ifd job')
    
    # building args by group list to separate Namespace args
    args_by_group['ifd'] = ['ligand', 'ifd_settings']
    args_by_group['piper'] = ['receptor_prot', 'ligand_prot', 'piper_settings']
    
    return parser, args_by_group

def parse_args(master_dir):
    """ Builds and parse user's inputs. 
    
    Returns: 
    -args: Namespace object containing user parsed args
    -args_by_group: dictionary grouping arguments by module in workflow """

    # building parser with defined user inputs
    parser, args_by_group = build_parser() 

    # using parser to parse user inputs
    # collecting known and unknown arguments
    args, unknowns = parser.parse_known_args()

    # changing file paths to full paths (in case we change working directories)
    args.cereblon = os.path.join(master_dir, args.cereblon)
    args.protein = os.path.join(master_dir, args.protein)
    args.ligand = os.path.join(master_dir, args.ligand)
    
    # making sure that the arguments are valid before proceeding
    if tcm_check_input.check_args(parser ,sys.argv, args, unknowns) is True:
        sys.exit(0) # exits because fatal error found
    
    return args, args_by_group
