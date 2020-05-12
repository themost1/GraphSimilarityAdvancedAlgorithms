
import numpy as np
import matplotlib.pyplot as plt
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.DSSP import DSSP
from Bio.PDB.Polypeptide import PPBuilder
from PPBuilderPlus import PPBuilderPlus
# create parser
parser = PDBParser()
AA_CODES = 'ACDEFGHIKLMNPQRSTVWXY'
Q8_CODES = 'GHITEBS-'
AA_DIMS = 21
Q8_DIMS = 8


def int_encoding(seq, code):
  # Input: sequence(str), CODE
  # Output: (n,1) vector

  encoding = np.zeros((len(seq),), dtype=np.int64)
  for idx, seq_chr in enumerate(seq):
    encoding[idx] = code.find(seq_chr)

  return encoding
def get_aa_encoded(protein_file):
    structure = Bio.PDB.PDBParser(QUIET=True).get_structure(pdb_path[:-4], pdb_path)
    ppb = PPBuilder()
    pdb_aas = []
    for pp in ppb.build_peptides(structure): 
        pdb_aa = str(pp.get_sequence())
        pdb_aas.append(pdb_aa)
    encoded = int_encoding(pdb_aas, AA_CODES)
    return encoded
def get_calpha_coords(protein_file): 
	structure = parser.get_structure('PHA-L', protein_file)

	model = structure[0]
	chain = model['A']

	cmat = np.zeros((len(chain), 3))
	# this example uses only the first residue of a single chain.

	# it is easy to extend this to multiple chains and residues.
	for row,r1 in enumerate(chain):  
		v = r1['CA'].get_vector()
		cmat[row,:] = v
	return cmat		

def get_dcalpha(protein_file): 
# read structure from file
	structure = parser.get_structure('PHA-L', protein_file)

	model = structure[0]
	chain = model['A']

	dmat = np.zeros((len(chain), len(chain)))
	# this example uses only the first residue of a single chain.

	# it is easy to extend this to multiple chains and residues.
	for row,r1 in enumerate(chain): 
		for col,r2 in enumerate(chain): 
			distance = r1['CA'] - r2['CA']
			dmat[row, col] = distance
	return dmat		
def dcalpha_to_adj(dcalpha):
    #normalized = dcalpha / np.max(dcalpha) # TODO: May need to use global scaling factor
    #dcalpha_adj = 1 - normalized

    dcalpha_adj = dcalpha + np.eye(dcalpha.shape[0]) # Add self-loop of spurious distance 1 for inverse distance calculation
    dcalpha_adj = 1 / dcalpha_adj
    assert(np.sum((dcalpha_adj == 0).astype('int32')) == 0) # Check underflow

    return dcalpha_adj

# Converts the C-alpha distance matrix into a protein contact map
# Note: Studies suggest defining threshold to be in 10-18 Angstrom range
def dcalpha_to_contact(dcalpha, threshold=10):
    adj_contact_map = (dcalpha < threshold).astype('float32')
    degree_contact_map = np.diag(np.sum(adj_contact_map, axis=1))

    return adj_contact_map, degree_contact_map
def makeProteinGraph(protein_file):
	aa = get_aa_encoded(protein_file)
	print(aa)
	#coords = get_calpha_coords(protein_file)
	dmat = get_dcalpha(protein_file)
	adj_contact_map = dcalpha_to_adj(dmat)
	dcalpha_to_contact = dcalpha_to_contact(dmat)
makeProteinGraph('/home/arjunsrivatsa/GraphSimilarityAdvancedAlgorithms/Proteins_to_compare/5jo9_model_4_100k_ground_truth.pdb' 
)