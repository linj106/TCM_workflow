#Import Python modules
import logging
import textwrap
import argparse
import sys
import os

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
    input = parser.add_argument_group('PROTEIN/LIGAND INPUTS') # inputs to IFD (protein and ligand files)
    job_control = parser.add_argument_group('JOB CONTROL INFORMATION') #arguments related to job control / server submission (e.g. HOST and JOBNAME)
    h_bond_constraints = parser.add_argument_group('HYDROGEN BONDING CONSTRAINTS BETWEEN LIGAND AND CRBN IN IFD') #inputs related to h bond constraints while docking ligand 

    # adding specific arguments to our input group
    input.add_argument('-l', '--ligands', '--ligand', dest = 'ligands', required = True, help = 'input file of ligand to dock; input must be .mae (either compressed or uncompressed)')
    input.add_argument('-p', '--protein', '--proteins', dest = 'proteins', required = True, help = 'input file of protein poses; input must be .mae (either compressed or uncompressed)')

    # adding specific arguments to change server/job info group
    job_control.add_argument('--NGLIDECPU', dest = 'NGLIDECPU', help = 'maximum number of Glide jobs to run simultaneously')
    job_control.add_argument('--NPRIMECPU', dest = 'NPRIMECPU', help = 'maximum number of Prime jobs to run simultaneously')
    job_control.add_argument('--NOLOCAL', dest = 'NOLOCAL', action='store_true', help = 'forces driver to run in local temp directory, should be used with remote HOST arguments')
    job_control.add_argument('--HOST', dest = 'HOST', help = 'specific host on BMS RHEL8 cluster to submit job to')
    job_control.add_argument('--SUBHOST', dest = 'SUBHOST', help = 'specific host to run subjobs on')
    job_control.add_argument('--TMPLAUNCHDIR', dest = 'TMPLAUNCHDIR', action = 'store_true', help = 'launches temporary directory to store the data used by system')
    job_control.add_argument('--jobname', dest = 'jobname', help = 'custom name for job to display on BMS RHEL8 cluster')
    job_control.add_argument('-d, --debug', dest = 'DEBUG', help = 'shows details of job control to help with debugging')
    job_control.add_argument('-o', '--output', dest = 'output', help = 'output directory to save results to; directory must already exist')
    
    # adding specific arguments to add constraints
    h_bond_constraints.add_argument('--hbond', '--constraints', nargs = '+', dest = 'h_bond_constraints', action = ParseKeyValuePairs, 
                                    help = """hydrogen bond restraints are represented by the atom number on CRBN and acceptor or donor (amine hydrogens are donors and carbonyl oxygens are acceptors);
                                    ex: --hbond 8479 donor 8440 acceptor 8451 donor -> means that atom 8479 amd 8451 are h-bond donors and 8440 is h-bond acceptor""")
    return parser

def parse_args():
    """ Builds parser and parse user's inputs. 
    
    Returns:
    args - recognized user-parsed arugments """

    # building parser with defined user inputs
    parser = build_parser()

    # using parse to parse user inputs
    # collecting known and unknown arguments
    args, unknowns = parser.parse_known_args()

    logger.info('This works')
    return args