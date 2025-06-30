#Import Python modules
import logging

#Import Schrodinger modules
from schrodinger.structure import StructureReader
from schrodinger.structutils.analyze import find_ligands

###Initiate logger###
logger = logging.getLogger(__name__)

def get_ligand_name(input_protein_file):
    """ Given input protein file, finds the name of the ligand and makes sure all the proteins in the 
    input_protein file have exactly one ligand with this same name. Raises assertion errors if
    this is not the case.
    
    Input:
    -input_protein_file: file path to input protein (either multiple structures or single structure)
    
    Returns: name of ligand in protein files which should match the residue.pdbres attribute """
    final_ligand = None

    # finds the name of the ligand (all protein files should have the same ligand name and only one ligand per structure)
    for st in StructureReader(input_protein_file): 
        ligand = find_ligands(st) #finding ligand in every structure (returns as list of all ligands)
        assert len(ligand) == 1, "Please make sure that the input protein file has only one ligand per protein structures"
        ligand_name = ligand[0].pdbres # iterating through ligand name and checks that the ligand is the same across all structures 
        if final_ligand is None:
            final_ligand = ligand_name
        else:
            assert ligand_name == final_ligand, "Please ensure that the same ligand exists in all protein files"
    
    return final_ligand

def find_ligand(residue, ligand):
    """Checks if the residue is our desired ligand. If it is the residue, it processes the residue information
    to return our ligand id to be used for our .inp file. Otherwise, it returns None. 
    
    Input: 
    -residue: class Object of the current residue in the structure (created through iteration w. schrodinger.StructureReader)
    -chain: class Object of the current chain in the structure (created through iteration w. schrodinger.StructureReader)
    -ligand: name of the ligand to match
    
    Returns: ligand_id (str) or None """

    if residue.pdbres == ligand: # if this is the correct residue
        # builds the correct ligand identifier w. residue name and chian
        lig_res_num = residue.resnum 
        lig_chain = residue.chain
        ligand_id = f'{lig_chain}:{lig_res_num}'
        return ligand_id 
        
    else:
        return None

# def find_W380_backbone(structure, residue = 'TRP380'):""
def W380_backbone(residue):
    """Glutarimide ring of IMiDs and CELMoDs are known to make a hydrogen bond with the hydrogen bound 
    to the backbone nitrogen of W380 (this hydrogen acts as an h bond donor). This function checks if 
    the residue is TRP380 / W380. If it is, it returns the atom number corresponding to that specific hydrogen. Otherwise,
    returns None
    
    Input: 
    -residue: class Object of the current residue in the structure (created through iteration w. schrodinger.StructureReader)
    
    Returns: atom_num(int) or None  """

    if residue.resnum == 380 and residue.pdbres.strip() == 'TRP': # identified TRP380
        atom = residue.getAtomByPdbName(' H  ') # hydrogen corresponding to backbone nitrogen of W380 (pdb name of H)
        return int(atom.atom_name[1:])
    
def H378_backbone(residue):
    """Glutarimide ring of IMiDs and CELMoDs are known to make a hydrogen bond with the carbonyl
    oxygen if the backbone of H378. This function checks if the residue is HIS378 / H378. """

    if residue.resnum == 378 and residue.pdbres.strip() == 'HIS': # identified HIS378
        atom = residue.getBackboneOxygen() #getting carbonyl backbone oxygen
        return int(atom.atom_name[1:])
    
def H378_side_chain(residue):
    """Glutarimide ring of IMiDs and CELMoDs are known to make a hydrogen bond with the hydrogen bound 
    to the ND1 nitrogen in sidechain of W380. This function checks if the residue is HIS378 / H378. """

    if residue.resnum == 378 and residue.pdbres.strip() == 'HIS': # identified HIS378
        atom = residue.getAtomByPdbName(' HD1')
        return int(atom.atom_name[1:])
    
def parse_structure(input_protein_file):
    """ Given the input proteins file, parses the first structure of the protein file
    to find the numbers for the ligand and specific hydrogen bonding sites on CRBN (backbone of 
    HIS378, sidechain of HIS378, and backbone of TRP380). 
    
    Input:
    - input_protein_file: path to protein file (.mae or .maegz)
    
    Return: results dictionary such as {'lig_id':C:502, 'H378_o': 'O8440', 'H378_h': 'H8451', 'W380_h':'H8479'}
    """
    
    results = {}

    # finding the name of the ligand
    lig = get_ligand_name(input_protein_file)

    # iterating through the residues in the structure
    with StructureReader(input_protein_file) as reader:
        for structure in reader:
            for residue in structure.residue:
                if 'lig_id' not in results: #ligand not yet found
                    lig_id = find_ligand(residue, lig) # checks if this is ligand and returns correct info if it is (o.w. None)
                    if lig_id is not None:
                        results['lig_id'] = lig_id
                if 'W380_h' not in results: #hydrogen on W380 backbone not yet found
                    backbone_h = W380_backbone(residue) # checks if this is W380 and returns correct info if it is (o.w. None)
                    if backbone_h is not None:
                        results['W380_h'] = backbone_h
                if 'H378_o' not in results: #oxygen on H378 backbone not yet found
                    backbone_o = H378_backbone(residue) # checks if this is H378 and returns correct info if it is (o.w. None)
                    if backbone_o is not None:
                        results['H378_o'] = backbone_o
                if 'H378_h' not in results: #hydrogen on H378 side chain not yet found
                    side_chain_h = H378_side_chain(residue) # checks if this is H378 and returns correct info if it is (o.w. None)
                    if side_chain_h is not None:
                        results['H378_h'] = side_chain_h
    
    return results 
                    
if __name__ == '__main__':
    structure = 'InducedFit_5HXB-2_03Angstrom_rec.mae'
    res = parse_structure(structure)
    print(res)
