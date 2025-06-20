#Import Python modules
import logging
import textwrap
import argparse
import sys
import os

#Import PIPER functionality scripts
import piper_check_input
import piper_default
import piper_run

###Initiate logger###
logger = logging.getLogger(__name__)

# Initiate class of argument parser with some custom functions
class ArgumentParser(argparse.ArgumentParser):
    
    def error(self, message):
        """ADD COMMENT"""
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
    input = parser.add_argument_group('PROTEIN INPUTS TO DOCK') # protein inputs to dock (one as receptor and other as ligand)
    constraints = parser.add_argument_group('DISTANCE AND REPULSION CONSTRAINTS') # arguments related to user inputted distance constraints (e.g. HDX-MS data)
    job_control = parser.add_argument_group('JOB CONTROL INFORMATION') #arguments related to job control / server submission (e.g. HOST and JOBNAME)
    options = parser.add_argument_group('PIPER DOCKING SETTINGS') #arguments related to specific PIPER settings (e.g. number of poses)

    # adding specific arguments to our input group
    input.add_argument('-r', '--receptor', dest = 'receptor_prot', required = True, help = 'protein file acting as receptor in PIPER docking; input file must be .mae or .pdb')
    input.add_argument('--r_chain','--receptor_chain', dest = 'receptor_chain', help = 'specific chain in receptor protein to use as receptor')
    input.add_argument('-l', '--ligand', dest = 'ligand_prot', required = True, help = 'protein file acting as ligand in PIPER docking')
    input.add_argument('--l_chain', '--ligand_chain', dest = 'ligand_chain', help = 'specific chain in ligand protein to use as ligand')
    
    # adding specific arguments to change job settings / options 
    options.add_argument('--use_nonstandard_residue', dest = 'use_nonstandard_residue', help = 'whether or not to use nonstandard residues in docking calculations')
    options.add_argument('--rotations', dest = 'rotations', help = 'number of rotation matrices to use from rotation file')
    options.add_argument('--poses', dest = 'poses', help = 'max number of different poses to return from docking')
    options.add_argument('--raw', dest = 'raw', help = 'store all poses in pose-viwer format without refinement')

    # adding specific arguments to change server/job info group
    job_control.add_argument('--host', dest = 'HOST', help = 'specific host on BMS RHEL8 cluster to submit job to')
    job_control.add_argument('--jobname', dest = 'jobname', help = 'custom name for job to display on BMS RHEL8 cluster')
    job_control.add_argument('--ompi', dest = 'OMPI', help = 'number of processors to run job on')
    job_control.add_argument('-d, --debug', dest = 'DEBUG', help = 'shows details of job control to help with debugging')
    job_control.add_argument('--job_id', dest = 'JOBID', help = 'runs the job through job control layer')
    job_control.add_argument('-o', '--output', dest = 'output', help = 'output directory to save results to; directory must already exist')
    
    # adding specific arguments to add constraints
    constraints.add_argument('--dist', '--distance_constraint', nargs = '+', dest = 'distance_constraint', 
        help = """
        Distance Constraints Argument
        - input: .txt file containing distance constraints info (example below)
        - use new line for different distance constraint pairs
        - REC_RESIDUE (3 Letter AA + Residue Num) LIG_RESIDUE (3 Letter AA + Residue Num) DMIN DMAX

        ### example_distance.txt
        HIS354 LYS672 2.0 5.0
        GLY354 HIS765 2.0 7.0 
        
        """)
    
    constraints.add_argument('--attract', '--attraction', nargs = '+', dest = 'attraction', 
        help = """
        Attraction Constraints Argument
        - input: one or more .txt files each containing separate attraction constraints info 

        Guidance on .txt file
        - LINE1: residues involved in attraction in the form of 3 letter AA + residue num (e.g. HIS375)
        - LINE2: attraction bonus 
        - LINE3: protein ('receptor' or 'ligand')

        ### example_attraction.txt
        HIS354 LYS672 GLY354 HIS765
        0.11
        receptor
        
        """)
    
    constraints.add_argument('--repel', '--repulsion', nargs = '+', dest = 'repulsion', 
        help = """
        Repulsion Constraints Argument
        - input: one or more .txt files each containing separate repulsion constraints info 

        Guidance on .txt file
        - LINE1: residues involved in repulsion in the form of 3 letter AA + residue num (e.g. HIS375)
        - LINE2: protein in which the residues are on ('receptor' or 'ligand')

        ### example_repulsion.txt
        HIS354 LYS672 GLY354 HIS765
        receptor
        
        """)

    return parser

def parse_and_check_args():
    """ Builds parser and parse user's inputs. Then, runs checks on user inputs and exits if fatal error found.
    
    Return: 
    args - recognized user-parsed arguments """

    # building parser with defined user inputs
    parser = build_parser() 

    # using parser to parse user inputs
    # collecting known and unknown arguments
    args, unknowns = parser.parse_known_args()

    # making sure that the arguments are valid before proceeding
    if piper_check_input.check_args(parser ,sys.argv, args, unknowns) is True:
        sys.exit(0) # exits because fatal error found
    
    return args



    

