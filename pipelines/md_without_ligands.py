import os
import subprocess
import argparse

def run_command(command, workdir=None):
    """Utility function to run shell commands and handle errors."""
    result = subprocess.run(command, shell=True, cwd=workdir, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    return result.stdout.decode('utf-8')

def create_mdp_files(target_dir):
    """Create necessary MDP files for ions, energy minimization, NVT, and NPT equilibration."""
    ions_mdp = os.path.join(target_dir, "ions.mdp")
    minim_mdp = os.path.join(target_dir, "minim.mdp")
    nvt_mdp = os.path.join(target_dir, "nvt.mdp")
    npt_mdp = os.path.join(target_dir, "npt.mdp")
    
    # ions.mdp
    with open(ions_mdp, 'w') as f:
        f.write("""
        integrator = steep
        emtol = 1000.0
        nsteps = 50000
        nstenergy = 10
        coulombtype = PME
        rlist = 1.2
        rcoulomb = 1.2
        rvdw = 1.2
        """)
    
    # minim.mdp
    with open(minim_mdp, 'w') as f:
        f.write("""
        integrator = steep
        emtol = 1000.0
        nsteps = 50000
        nstenergy = 10
        coulombtype = PME
        rlist = 1.2
        rcoulomb = 1.2
        rvdw = 1.2
        """)
    
    # nvt.mdp
    with open(nvt_mdp, 'w') as f:
        f.write("""
        integrator = md
        nsteps = 50000
        dt = 0.002
        tcoupl = Berendsen
        tc-grps = Protein Water_and_ions
        tau_t = 0.1
        ref_t = 300
        """)
    
    # npt.mdp
    with open(npt_mdp, 'w') as f:
        f.write("""
        integrator = md
        nsteps = 50000
        dt = 0.002
        tcoupl = Berendsen
        tc-grps = Protein Water_and_ions
        tau_t = 0.1
        ref_t = 300
        pcoupl = Parrinello-Rahman
        ref_p = 1.0
        compressibility = 4.5e-5
        """)
    
    return ions_mdp, minim_mdp, nvt_mdp, npt_mdp

def run_gromacs_pipeline(input_pdb, target_dir):
    """Run the GROMACS MD simulation pipeline with GROMOS 53a6 force field, without ligands."""
    os.makedirs(target_dir, exist_ok=True)

    # Create necessary MDP files
    ions_mdp, minim_mdp, nvt_mdp, npt_mdp = create_mdp_files(target_dir)

    # Step 1: Convert protein PDB to GRO using pdb2gmx
    processed_gro = os.path.join(target_dir, "processed.gro")
    protein_top = os.path.join(target_dir, "topol.top")
    command = f"gmx pdb2gmx -f {input_pdb} -o {processed_gro} -water spce -ff gromos53a6 -p {protein_top}"
    print("Running pdb2gmx for protein...")
    run_command(command)

    # Step 2: Define the simulation box
    box_gro = os.path.join(target_dir, "newbox.gro")
    command = f"gmx editconf -f {processed_gro} -o {box_gro} -c -d 1.0 -bt cubic"
    print("Running editconf to define box...")
    run_command(command)

    # Step 3: Solvate the system
    solvated_gro = os.path.join(target_dir, "solvated.gro")
    command = f"gmx solvate -cp {box_gro} -cs spc216.gro -o {solvated_gro} -p {protein_top}"
    print("Running solvate...")
    run_command(command)

    # Step 4: Add ions to neutralize the system
    ions_tpr = os.path.join(target_dir, "ions.tpr")
    command = f"gmx grompp -f {ions_mdp} -c {solvated_gro} -p {protein_top} -o {ions_tpr} -maxwarn 2"
    run_command(command)

    command = f"gmx genion -s {ions_tpr} -o {os.path.join(target_dir, 'solvated_ions.gro')} -p {protein_top} -pname NA -nname CL -neutral"
    run_command(command)

    # Step 5: Energy minimization (Steepest Descent)
    em_tpr = os.path.join(target_dir, "em.tpr")
    em_gro = os.path.join(target_dir, "em.gro")
    command = f"gmx grompp -f {minim_mdp} -c {os.path.join(target_dir, 'solvated_ions.gro')} -p {protein_top} -o {em_tpr} -maxwarn 2"
    run_command(command)
    command = f"gmx mdrun -v -deffnm {os.path.join(target_dir, 'em')}"
    run_command(command)

    # Step 6: NVT Equilibration (Berendsen Thermostat)
    nvt_tpr = os.path.join(target_dir, "nvt.tpr")
    command = f"gmx grompp -f {nvt_mdp} -c {em_gro} -p {protein_top} -o {nvt_tpr} -maxwarn 2"
    run_command(command)
    command = f"gmx mdrun -v -deffnm {os.path.join(target_dir, 'nvt')}"
    run_command(command)

    # Step 7: NPT Equilibration (Parrinello-Rahman Pressure Coupling)
    npt_tpr = os.path.join(target_dir, "npt.tpr")
    command = f"gmx grompp -f {npt_mdp} -c {os.path.join(target_dir, 'nvt.gro')} -p {protein_top} -o {npt_tpr} -maxwarn 2"
    run_command(command)
    command = f"gmx mdrun -v -deffnm {os.path.join(target_dir, 'npt')}"
    run_command(command)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run GROMACS MD simulation pipeline without ligands")
    parser.add_argument("pdb_file", help="Input PDB file")
    parser.add_argument("target_dir", help="Directory to save output files")
    
    args = parser.parse_args()

    run_gromacs_pipeline(args.pdb_file, args.target_dir)
