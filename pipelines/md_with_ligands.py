import os
import subprocess
import argparse
from Bio import PDB

def run_command(command, workdir=None):
    """Utility function to run shell commands."""
    result = subprocess.run(command, shell=True, cwd=workdir, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    return result.stdout.decode('utf-8')

def detect_ligands(input_pdb):
    """Detect non-standard residues (ligands) from the PDB file using Biopython."""
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure('protein', input_pdb)

    ligands = []
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.id[0] != ' ' and PDB.is_aa(residue) is False and residue.resname != "HOH":
                    ligands.append(residue.resname)
    
    return list(set(ligands))  # Return unique ligands

def split_protein_and_ligand(input_pdb, ligand_name, target_dir):
    """Split the protein and ligand into separate PDB files."""
    protein_pdb = os.path.join(target_dir, "protein_only.pdb")
    ligand_pdb = os.path.join(target_dir, f"{ligand_name}.pdb")

    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure('protein', input_pdb)

    io = PDB.PDBIO()

    # Save protein without ligand
    class ProteinSelect(PDB.Select):
        def accept_residue(self, residue):
            return residue.resname != ligand_name

    io.set_structure(structure)
    io.save(protein_pdb, ProteinSelect())

    # Save ligand separately
    class LigandSelect(PDB.Select):
        def accept_residue(self, residue):
            return residue.resname == ligand_name

    io.set_structure(structure)
    io.save(ligand_pdb, LigandSelect())

    return protein_pdb, ligand_pdb

def modify_topology_file(protein_top, ligand_itp, topol_top, ligands):
    """Modify the topol.top file to include ligand topology."""
    with open(protein_top, 'r') as top_file:
        top_lines = top_file.readlines()

    with open(topol_top, 'w') as new_top_file:
        for line in top_lines:
            new_top_file.write(line)

        # Include the ligand topology file
        new_top_file.write(f'#include "{ligand_itp}"\n')

        # Add the ligand to the system definition
        new_top_file.write("\n[ system ]\n")
        new_top_file.write("Protein with ligand in water\n")

        new_top_file.write("\n[ molecules ]\n")
        new_top_file.write("Protein_chain_A    1\n")  # Assuming 1 protein chain
        for ligand in ligands:
            new_top_file.write(f"{ligand} 1\n")

def run_gromacs_pipeline(input_pdb, target_dir):
    """Run the GROMACS MD simulation pipeline with GROMOS 53a6 force field."""
    os.makedirs(target_dir, exist_ok=True)

    # Step 1: Detect ligands (ignore water molecules)
    ligands = detect_ligands(input_pdb)
    print(f"Detected ligands: {ligands}")

    # Step 2: Split protein and ligand
    ligand_name = ligands[0]  # Assume one ligand for simplicity
    protein_pdb, ligand_pdb = split_protein_and_ligand(input_pdb, ligand_name, target_dir)

    # Step 3: Convert protein PDB to GRO using pdb2gmx
    processed_gro = os.path.join(target_dir, "processed.gro")
    protein_top = os.path.join(target_dir, "topol.top")
    command = f"gmx pdb2gmx -f {protein_pdb} -o {processed_gro} -water spce -ff gromos53a6 -p {protein_top}"
    print("Running pdb2gmx for protein...")
    run_command(command)

    # Step 4: Get FMN topology from ATB
    print("Assuming FMN topology generated via ATB")
    ligand_itp = os.path.join(target_dir, "fmn.itp")  # Path to FMN .itp from ATB

    # Step 5: Modify the topology to include the FMN ligand
    modify_topology_file(protein_top, ligand_itp, os.path.join(target_dir, "final_topol.top"), ligands)

    # Step 6: Define the simulation box
    box_gro = os.path.join(target_dir, "newbox.gro")
    command = f"gmx editconf -f {processed_gro} -o {box_gro} -c -d 1.0 -bt cubic"
    print("Running editconf...")
    run_command(command)

    # Step 7: Solvate the system
    solvated_gro = os.path.join(target_dir, "solvated.gro")
    command = f"gmx solvate -cp {box_gro} -cs spc216.gro -o {solvated_gro} -p {os.path.join(target_dir, 'final_topol.top')}"
    print("Running solvate...")
    run_command(command)

    # Step 8: Add ions to neutralize
    ions_tpr = os.path.join(target_dir, "ions.tpr")
    command = f"gmx grompp -f ions.mdp -c {solvated_gro} -p {os.path.join(target_dir, 'final_topol.top')} -o {ions_tpr}"
    run_command(command)

    command = f"gmx genion -s {ions_tpr} -o {os.path.join(target_dir, 'solvated_ions.gro')} -p {os.path.join(target_dir, 'final_topol.top')} -pname NA -nname CL -neutral"
    run_command(command)

    # Step 9: Energy minimization (Steepest Descent)
    em_tpr = os.path.join(target_dir, "em.tpr")
    em_gro = os.path.join(target_dir, "em.gro")
    command = f"gmx grompp -f minim.mdp -c {os.path.join(target_dir, 'solvated_ions.gro')} -p {os.path.join(target_dir, 'final_topol.top')} -o {em_tpr}"
    run_command(command)
    command = f"gmx mdrun -v -deffnm em"
    run_command(command)

    # Step 10: NVT Equilibration (Berendsen Thermostat)
    nvt_tpr = os.path.join(target_dir, "nvt.tpr")
    command = f"gmx grompp -f nvt.mdp -c {em_gro} -p {os.path.join(target_dir, 'final_topol.top')} -o {nvt_tpr}"
    run_command(command)
    command = f"gmx mdrun -v -deffnm nvt"
    run_command(command)

    # Step 11: NPT Equilibration (Parrinello-Rahman Pressure Coupling)
    npt_tpr = os.path.join(target_dir, "npt.tpr")
    command = f"gmx grompp -f npt.mdp -c {os.path.join(target_dir, 'nvt.gro')} -p {os.path.join(target_dir, 'final_topol.top')} -o {npt_tpr}"
    run_command(command)
    command = f"gmx mdrun -v -deffnm npt"
    run_command(command)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run GROMACS MD simulation pipeline with GROMOS force field and FMN ligand")
    parser.add_argument("pdb_file", help="Input PDB file")
    parser.add_argument("target_dir", help="Directory to save output files")
    
    args = parser.parse_args()

    run_gromacs_pipeline(args.pdb_file, args.target_dir)
