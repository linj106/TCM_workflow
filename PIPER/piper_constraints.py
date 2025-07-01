#Import Python modules
import logging
import sys
import os
import json

###Initiate logger###
logger = logging.getLogger(__name__)

def write_constraint_json(constraints_list, path = None):
    """ Writes the constraint dictionary to a constraints.json file in the path or cwd if path = none
    
    Input:
    - constraints_list: list containing dictionary representation of all constraints 
    - path: path of folder to write json file to; default is None """
    if path is None:
        with open("constraints.json", "w") as f:
            json.dump(constraints_list, f, indent = 2)
    else:
        with open(os.path.join(path, "constraints.json"), "w") as f:
            json.dump(constraints_list, f, indent = 2)
        
def build_asl_string(residue_number, residue_type):
    """ Builds asl string representation to select the specific residue number and type. 
    
    Input:
    - residue number: specific residue number in protein (int); e.g. 355
    - residue type: specific residue identity (str); e.g. HIS 
    
    Returns: asl formatted string that contains selection information """
    return f'(res.num {residue_number} and res.ptype "{residue_type}")'

def build_repulsion_constraint(selected_residues, protein_type):
    """ Builds dictionary representation of a single repulsion constraint.  
    
    Input:
    - selected residues: residues in specific protein involved in repulsion; list of tuples (residue number: type = int, residue type: type = str)
    - protein_type: whether the selected proteins are on the receptor or ligand protein 
    
    Return: dictionary representing repulsion of the selected residues in protein-protein interaction """
    
    # initializing dictionary representation of repulsion info
    repulsion_dictionary = {}
    repulsion_dictionary["constraint_type"] = "repulsion" # defining constraint as repulsion
    
    # getting protein_type involved in repulsion (the protein in which the selected residues are a part of)
    assert protein_type in ["receptor", "ligand"]
    repulsion_dictionary["protein_type"] = protein_type
    
    # building atom selection language (asl) string to select specific atoms
    residues_str = []
    for residue in selected_residues:
        residues_str.append(build_asl_string(residue[0], residue[1]))
    asl = ' OR '.join(residues_str)

    repulsion_dictionary["asl"] = asl

    return repulsion_dictionary

def build_attraction_constraint(selected_residues, protein_type, bonus = 0.11):
    """ Builds dictionary representation of a single attraction constraint. 
    
    Input:
    - selected residues: residues in specific protein involved in repulsion; list of tuples (residue number - int, residue type - str)
    - protein_type: whether the selected proteins are on the receptor or ligand protein 
    - bonus: value added to the default of 1 to define scaling factor for attractive potential 
    
    Return: dictionary representing attraction of the selected residues in protein-protein interaction """

    # building dictionary representation of attraction info
    attraction_dictionary = {}
    attraction_dictionary["constraint_type"] = "attraction" # defining constraint as attraction

    # getting protein_type involved in repulsion (the protein in which the selected residues are a part of)
    assert protein_type in ["receptor", "ligand"]
    attraction_dictionary["protein_type"] = protein_type
    
    # building atom selection language (asl) string to select specific atoms
    residues_str = []
    for residue in selected_residues:
        residues_str.append(build_asl_string(residue[0], residue[1]))
    asl = ' OR '.join(residues_str)

    attraction_dictionary["asl"] = asl

    # setting attraction term (with value 10 x bonus)
    attraction_dictionary["attraction"] = bonus * 10

    return attraction_dictionary

def build_distance_pair(receptor_residue, ligand_residue, dmax, dmin = 2):
    """ Builds dictionary representation of a single distance pair constraint (integrated within broader distance constraint). 
    
    Input: 
    - receptor_residue: tuple of (residue number, residue type) identifying the receptor residue involved in distance constraint
    - ligand_residue: tuple of (residue number, residue type) identifying the ligand residue involved in distance constraint 
    - dmax: maximum distance in Angstroms allowed in protein-protein docking; int 
    - dmin: minimum distance allowed in protein-protein docking; int; default is 2 Angstroms 
    
    Return: dictionary representing distance constraint between receptor and ligand pair in protein-protein docking """

    distance_pair = {}

    # identifying receptor asl
    distance_pair["rec_asl"] = build_asl_string(receptor_residue[0], receptor_residue[1])

    # identifying ligand asl
    distance_pair["lig_asl"] = build_asl_string(ligand_residue[0], ligand_residue[1])

    # setting dmin and dmax
    distance_pair["dmin"] = dmin
    distance_pair["dmax"] = dmax

    return distance_pair

def build_distance_constraint(distance_pairs_list):
    """ Builds dictionary representation of all distance pair constraints.
    
    Input:
    - distance_pairs_list: list of distance pairs information represented by tuple of (rec residue num, rec residue type, lig residue num, lig residue type, dmin, dmax) 
    
    Return:
    - dictionary representing all distance constraints; all inputted distance constraints must be fulfilled - required = len(distance_pairs_list) """

    # defining distance constraints and base term 
    distance_constraint = {}
    distance_constraint["constraint_type"] = "distance"
    distance_constraint["required"] = len(distance_pairs_list)

    # adding distance pairs
    distance_pairs = []
    for pair in distance_pairs_list:
        receptor_residue = (pair[0], pair[1])
        ligand_residue = (pair[2], pair[3])
        dmin = pair[4]
        dmax = pair[5]
        distance_pairs.append(build_distance_pair(receptor_residue, ligand_residue, dmax, dmin))
    
    distance_constraint["distance_pairs"] = distance_pairs

    return distance_constraint

def parse_constraint_file(constraint_file):
    """ Builds dictionary representation by reading through constraint file
    
    Input:
    - constraint_file: path to .txt file that contains constraints
     
    Return:
    - dictionary rep of constraints """
    all_constraints = []
    distance_pair_list = []

    with open(constraint_file, 'r') as f:
        for line_number, line in enumerate(f, start=1):
            split_line = line.strip().split()

            # if distance constraint, add info to distance_pair
            if split_line[0] == "distance":
                dmin = float(split_line[1])
                dmax = float(split_line[2])
                receptor_residue = split_line[3]
                ligand_residue = split_line[4]

                # splitting residue info into res num and res type
                # e.g. 'HIS238' -> 'HIS', '238'
                receptor_residue_type = receptor_residue[:3]
                receptor_residue_num = int(receptor_residue[3:])

                ligand_residue_type = ligand_residue[:3]
                ligand_residue_num = int(ligand_residue[3:])

                # appending distance pair
                distance_pair = (receptor_residue_num, receptor_residue_type, ligand_residue_num, ligand_residue_type, dmin, dmax)
                distance_pair_list.append(distance_pair)
            
            # if attraction constraint
            elif split_line[0] == "attraction":
                bonus = round(float(split_line[1]), 2)
                protein_type = split_line[2]
                selected_residues = [(int(residue[3:]), residue[:3]) for residue in split_line[3:]] # iterating through remaining list and processing residue
                all_constraints.append(build_attraction_constraint(selected_residues, protein_type, bonus)) # build attraction constraint dict and add
            
            # if repulsion constraint
            elif split_line[0] == "repulsion":
                protein_type = split_line[1]
                selected_residues = [(int(residue[3:]), residue[:3]) for residue in split_line[2:]] # iterating through remaining list and processing residue
                all_constraints.append(build_repulsion_constraint(selected_residues, protein_type)) # build repulsion constraint dict and add

            # else is comment and ignore
            else:
                continue

    # build distance constraint and add
    all_constraints.append(build_distance_constraint(distance_pair_list))
    
    return all_constraints

def main(args, path_to_dir):
    """ Writes constraints list from args and saves to constraints.json file in path_to_dir. Modifies 
    args to delete parsed files arguments and replace with proper argument to correct constraints_file. 
    
    Returns:
    - new args with constraint_file argument """

    # builds and writes constraints file 
    constraints_list = parse_constraint_file(args.constraint_file)
    write_constraint_json(constraints_list, path_to_dir)

    # adding new argument with path to constraint file
    args.constraints_file = os.path.join(path_to_dir, "constraints.json")
    
    return args

# DEPRECATED
def parse_distance_constraint_file(distance_constraint_file):
    """ Builds dictionary representation by reading through distance constraint file.
    
    Input:
    - distance_constraint_file: .txt file that contains distance restraints
    
    Return:
    - dictionary representation of distance constraint """

    distance_pair_list = []

    with open(distance_constraint_file, 'r') as file: # opening distance_constraint file
        for constraint in file: # iterating through lines
            constraint = constraint.strip()
            information = constraint.split() #splitting line by space
            
            # defining specific information 
            receptor_residue = information[0]
            ligand_residue = information[1]
            dmin = float(information[2])
            dmax = float(information[3])

            # splitting residue info into res num and res type
            # e.g. 'HIS238' -> 'HIS', '238'
            receptor_residue_type = receptor_residue[:3]
            receptor_residue_num = int(receptor_residue[3:])

            ligand_residue_type = ligand_residue[:3]
            ligand_residue_num = int(ligand_residue[3:])

            # appending distance pair
            distance_pair = (receptor_residue_num, receptor_residue_type, ligand_residue_num, ligand_residue_type, dmin, dmax)
            distance_pair_list.append(distance_pair)

    return build_distance_constraint(distance_pair_list)

def parse_attraction_file(attraction_constraint_file):
    """ Builds dictionary representation by reading through attraction constraint file
    
    Input:
    - attraction_constraint_file: .txt file that contains attraction constraint info (for a single attraction set)
    
    Return:
    - dictionary representation of attraction constraint """

    with open(attraction_constraint_file, 'r') as file:

        # reading in line 1 which contains all the residues
        residues = file.readline().strip().split()

        # reading in line 2 which contains the attraction bonus
        bonus = round(float(file.readline().strip()), 2)

        # reading in line 3 which contains either receptor or ligand
        receptor_or_lig = file.readline().strip()

    # building parsed residues into selected residue list with proper formatting
    selected_residues = []
    for residue in residues:
        residue_num = int(residue[3:])
        residue_type = residue[:3]
        selected_residues.append((residue_num, residue_type))
    
    return build_attraction_constraint(selected_residues, receptor_or_lig, bonus)

def parse_repulsion_file(repulsion_constraint_file):
    """ Builds dictionary representation by reading through repulsion constraint file
    
    Input:
    - attraction_constraint_file: .txt file that contains repulsion constraint info (for a single attraction set)
    
    Return:
    - dictionary representation of repulsion constraint """

    with open(repulsion_constraint_file, 'r') as file:

        # reading in line 1 which contains all the residues
        residues = file.readline().strip().split()

        # reading in line 2 which contains either receptor or ligand
        receptor_or_lig = file.readline().strip()

    # building parsed residues into selected residue list with proper formatting
    selected_residues = []
    for residue in residues:
        residue_num = int(residue[3:])
        residue_type = residue[:3]
        selected_residues.append((residue_num, residue_type))
    
    return build_repulsion_constraint(selected_residues, receptor_or_lig)

def build_constraints(args):
    """ Parses user inputs of args and builds dictionary representation of user-inputted restraints.
    
    Input: 
    - user-parsed arguments with attributes of .distance_constraint, .attraction, and .repulsion (if user inputted)
    
    Returns: list containing dictionary representation of all constraints ready to write to json """

    all_constraints = []

    if args.distance_constraint is not None: # if distance constraints exist from user
        print(args.distance_constraint, args.attraction)
        for file in args.distance_constraint: # go through list of files and parse each separately
            print(file)
            all_constraints.append(parse_distance_constraint_file(file))

    if args.attraction is not None: # if attraction constraints exist from user
        for file in args.attraction: # go through list of files and parse each separately
            all_constraints.append(parse_attraction_file(file))

    if args.repulsion is not None: # if repulsion constraints exist from user
        for file in args.repulsion: # go through list of files and parse each separately
            all_constraints.append(parse_repulsion_file(file))

    return all_constraints
