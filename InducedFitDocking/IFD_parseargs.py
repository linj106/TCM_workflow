#Import Python modules
import logging
import textwrap
import argparse
import subprocess
import os
import sys

# Import IFD functionality script
import IFD_check_input

###Initiate logger###
logger = logging.getLogger(__name__)

# Initiate class of argument parser with some custom functions
class ArgumentParser(argparse.ArgumentParser):
    
    def error(self, message):
        """ADD COMMENT"""
        logger.critical(message)
        print("An error has occurred. Check log file for information.")
        subprocess.run(['kill', '0'])

    """Disables prefix matching in ArgumentParser."""
    def _get_option_tuples(self, option_string):
        """Prevent argument parsing from looking for prefix matches."""
        return []
    
class ParseKeyValuePairs(argparse.Action):
    """Custom argpars.Action that takes in the arguments after a specific flag and groups
    them into paired tuples. For example, --hbond 8479 donor 8440 acceptor becomes [(8479, donor),(8440, acceptor)] """
    def __call__(self, parser, namespace, values, option_string=None):

        #check that even number of arguments passed in
        if len(values) % 2 != 0:
            raise argparse.ArgumentError(self, "Every atom must be paired with either donor or acceptor")
        
        #parsing and organizing into tuples
        results = []
        for i in range(0, len(values), 2):
            atom_num = values[i]
            h_type = values[i+1]
            try:
                int(atom_num)
                assert h_type in ['donor', 'acceptor']
                results.append((atom_num, h_type))

            except ValueError or AssertionError:
                raise argparse.ArgumentError(self, "Invalid arguments, must have atom number followed by donor or acceptor")
        
        setattr(namespace, self.dest, results) # setting our parsed results 

# Building Parser and Parsing Arguments
def str2bool(v: str) -> bool:
    """Convert string to boolean."""
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError(f"Invalid boolean value: {v}")

def full_path(path):
    """ Given path, returns the entire file path"""
    return os.path.join(os.getcwd(), path)

def build_parser():
    """ Builds our custom parser based on argparse.ArgumentParser. Defines the specific user inputs and
    potential flags with associated help. Parser arguments are largely derived from available settings to PIPER in Schrodinger
    as seen in https://learn.schrodinger.com/private/edu/release/current/Documentation/html/utilities/program_utility_usage/piper.html?Highlight=PIPER%20command.
    
    Returns: object of class argparse.ArgumentParser with defined user inputs"""

    # creating object of class ArgumentParser with program name and description
    parser = argparse.ArgumentParser(
        prog = 'Induced Fit Docking', 
        formatter_class=argparse.RawDescriptionHelpFormatter,
        usage='%(prog)s [options]',
        description=textwrap.dedent('''\
        ----------------------------------------------
        Induced Fit Docking:
        Command Line Implementation of Induced Fit Docking from Schrodinger. Parses 
        user arguments for generated poses between cereblon and protein of interest
        and ligand to dock to create proper input files for system call to IFD.
        ----------------------------------------------              
        '''))
    
    # organizing the parser arguments by groups
    input = parser.add_argument_group('PROTEIN/LIGAND INPUTS') # inputs to IFD (protein and ligand files)
    job_control = parser.add_argument_group('JOB CONTROL INFORMATION') #arguments related to job control / server submission (e.g. HOST and JOBNAME)
    h_bond_constraints = parser.add_argument_group('HYDROGEN BONDING CONSTRAINTS BETWEEN LIGAND AND CRBN IN IFD') #inputs related to h bond constraints while docking ligand 
    default = parser.add_argument_group('DEFAULT SETTINGS') # arguments related to changing default IFD settings (impt for module in TCM)

    # adding specific arguments to our input group
    input.add_argument('-l', '--ligands', '--ligand', type = full_path, dest = 'ligand', required = True, help = 'input file of ligand to dock; input must be .mae (either compressed or uncompressed)')
    input.add_argument('-p', '--protein', '--proteins', type = full_path, dest = 'proteins', required = True, help = 'input file of protein poses; input must be .mae (either compressed or uncompressed)')
    input.add_argument('--template', dest = 'template', type = full_path, help = 'template .inp file to use for IFD jobs')

    # adding specific arguments to change server/job info group
    job_control.add_argument('--NGLIDECPU', dest = 'NGLIDECPU', type = int, help = 'maximum number of Glide jobs to run simultaneously')
    job_control.add_argument('--NPRIMECPU', dest = 'NPRIMECPU', type = int, help = 'maximum number of Prime jobs to run simultaneously')
    job_control.add_argument('--NOLOCAL', dest = 'NOLOCAL', type = str2bool, help = 'forces driver to run in local temp directory; requires bool; external HOST requires NOLOCAL')
    job_control.add_argument('--HOST', dest = 'HOST', type = str, help = 'specific host on BMS RHEL8 cluster to submit job to')
    job_control.add_argument('--SUBHOST', dest = 'SUBHOST', type = str, help = 'specific host to run subjobs on')
    job_control.add_argument('--TMPLAUNCHDIR', dest = 'TMPLAUNCHDIR', type = str2bool, help = 'launches temporary directory to store the data used by system; requires bool')
    job_control.add_argument('--jobname', dest = 'jobname', type = str, help = 'custom name for job to display on BMS RHEL8 cluster')
    job_control.add_argument('-d, --debug', dest = 'DEBUG', type = str2bool, help = 'shows details of job control to help with debugging; requires bool')
    job_control.add_argument('-o', '--output', dest = 'output', type = full_path, help = 'directory to place results and loggers in; must already exist')
    
    # adding specific arguments to add constraints
    h_bond_constraints.add_argument('--hbond', '--constraints', nargs = '+', dest = 'h_bond_constraints', type = full_path, action = ParseKeyValuePairs, 
                                    help = """hydrogen bond restraints are represented by the atom number on CRBN and acceptor or donor (amine hydrogens are donors and carbonyl oxygens are acceptors);
                                    ex: --hbond 8479 donor 8440 acceptor 8451 donor -> means that atom 8479 amd 8451 are h-bond donors and 8440 is h-bond acceptor""")
    
    # adding specific arguments to change default settings (also for use in modules in which TCM workflow requires default json files to change settings of jobs)
    default.add_argument('--default', dest = 'default', type = full_path, help = 'json file containing the default settings for IFD job')

    return parser


def parse_and_check_args():
    """ Builds parser and parse user's inputs. 
    
    Returns:
    args - recognized user-parsed arugments """

    # building parser with defined user inputs
    parser = build_parser() 

    # using parser to parse user inputs
    # collecting known and unknown arguments
    args, unknowns = parser.parse_known_args()

    # making sure that the arguments are valid before proceeding
    if IFD_check_input.check_parsed_args(parser ,sys.argv, args, unknowns) is True:
        sys.exit(0) # exits because fatal error found
    
    return args
