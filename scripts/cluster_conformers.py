import argparse
import os
from rdkit import Chem
from rdkit.Chem import rdMolAlign
from rdkit.ML.Cluster import Butina

def main():
    # Setup argument parsing
    parser = argparse.ArgumentParser(description="Process conformers and save cluster centroids.")
    parser.add_argument('conformers_path', type=str, help='Path to the conformers .xyz file')
    parser.add_argument('save_base_path', type=str, help='Base directory to save the cluster centroids')
    parser.add_argument('n_atoms', type=int, help='Number of atoms that compose a structure')

    # Parse arguments
    args = parser.parse_args()

    conformers_path = args.conformers_path
    save_base_path = args.save_base_path
    n_atoms = args.n_atoms

    # Read the conformers file
    with open(conformers_path) as file:
        all_lines = file.readlines()

    xyz_blocks = []
    block = None
    comparison_string = f" {n_atoms}\n"
    # Split the conformers into separate blocks
    for line in all_lines:
        if line == comparison_string:
            if block is not None:
                xyz_blocks.append(block)
            block = line
            continue
        block = block + line
    xyz_blocks.append(block)

    # Convert XYZ blocks to RDKit molecules and remove hydrogens
    mols = []
    for block in xyz_blocks:
        mol = Chem.rdmolfiles.MolFromXYZBlock(block)
        if mol is not None:
            mol = Chem.RemoveAllHs(mol)
            mols.append(mol)

    # Calculate pairwise distances between conformers
    dists = []
    for i in range(len(mols)):
        for j in range(i):
            dists.append(rdMolAlign.AlignMol(mols[i], mols[j]))

    # Cluster conformers using Butina clustering algorithm
    clusts = Butina.ClusterData(dists, len(mols), 1.2, isDistData=True, reordering=True)

    # Save the cluster centroids
    cluster_id = 0
    for clust in clusts:
        mol_id = clust[0]
        xyz_string = xyz_blocks[mol_id]

        # Construct the save path
        xyz_path = os.path.join(save_base_path, f"cluster_centroid_{cluster_id:04d}.xyz")
        with open(xyz_path, "w") as file:
            file.write(xyz_string)

        cluster_id += 1

if __name__ == "__main__":
    main()
