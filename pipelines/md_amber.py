import os
import subprocess
import argparse

def run_command(command, workdir=None):
    """Utility function to run shell commands and handle errors."""
    result = subprocess.run(command, shell=True, cwd=workdir, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    
    stdout = result.stdout.decode('utf-8')
    stderr = result.stderr.decode('utf-8')
    
    # Print the command output (stdout) and error output (stderr) for tracking progress
    print(stdout)
    if stderr:
        print(stderr)
    
    # Raise an error if the command failed
    if result.returncode != 0:
        raise subprocess.CalledProcessError(result.returncode, command)
    
    return stdout

def create_mdp_files(target_dir):
    """Create necessary MDP files for energy minimization, NVT, and NPT equilibration with explicit solvation."""
    minim_mdp = os.path.join(target_dir, "minim.mdp")
    nvt_mdp = os.path.join(target_dir, "nvt.mdp")
    npt_mdp = os.path.join(target_dir, "npt.mdp")
    
    # minim.mdp for explicit solvation
    with open(minim_mdp, 'w') as f:
        f.write("""
        integrator = steep
        emtol = 1000.0
        nsteps = 50000
        emstep = 0.01
        cutoff-scheme = Verlet
        coulombtype = PME
        rcoulomb = 1.0
        rvdw = 1.0
        pbc = xyz
        energygrps = System
        """)

    # nvt.mdp for explicit solvation (constant temperature)
    with open(nvt_mdp, 'w') as f:
        f.write("""
        integrator = md
        nsteps = 50000
        dt = 0.002
        tcoupl = V-rescale
        tc-grps = System
        tau_t = 0.1
        ref_t = 300
        coulombtype = PME
        rcoulomb = 1.0
        rvdw = 1.0
        pbc = xyz
        energygrps = System
        """)

    # npt.mdp for explicit solvation (constant temperature and pressure)
    with open(npt_mdp, 'w') as f:
        f.write("""
        integrator = md
        nsteps = 50000
        dt = 0.002
        tcoupl = V-rescale
        tc-grps = System
        tau_t = 0.1
        ref_t = 300
        pcoupl = Parrinello-Rahman
        pcoupltype = isotropic
        tau_p = 2.0
        ref_p = 1.0
        compressibility = 4.5e-5
        coulombtype = PME
        rcoulomb = 1.0
        rvdw = 1.0
        pbc = xyz
        energygrps = System
        """)

    return minim_mdp, nvt_mdp, npt_mdp

def run_gromacs_pipeline(input_pdb, target_dir):
    """Run the GROMACS MD simulation pipeline with AMBER force field and explicit solvation."""
    os.makedirs(target_dir, exist_ok=True)

    # Create necessary MDP files
    minim_mdp, nvt_mdp, npt_mdp = create_mdp_files(target_dir)

    # Step 1: Convert protein PDB to GRO using pdb2gmx with AMBER force field
    processed_gro = os.path.join(target_dir, "processed.gro")
    protein_top = os.path.join(target_dir, "topol.top")
    command = f"gmx pdb2gmx -f {input_pdb} -o {processed_gro} -ff amber99sb -p {protein_top} -water spce"
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
    print("Running solvate and adding water molecules...")
    run_command(command)

    # Step 4: Add ions to neutralize the system
    ions_tpr = os.path.join(target_dir, "ions.tpr")
    command = f"gmx grompp -f {minim_mdp} -c {solvated_gro} -p {protein_top} -o {ions_tpr} -maxwarn 2"
    print("Running grompp to prepare for ion addition...")
    run_command(command)

    ionized_gro = os.path.join(target_dir, "solvated_ions.gro")
    command = f"gmx genion -s {ions_tpr} -o {ionized_gro} -p {protein_top} -pname NA -nname CL -neutral"
    print("Adding ions to neutralize the system...")
    run_command(command)

    # Step 5: Energy minimization (Steepest Descent)
    em_tpr = os.path.join(target_dir, "em.tpr")
    em_gro = os.path.join(target_dir, "em.gro")
    command = f"gmx grompp -f {minim_mdp} -c {ionized_gro} -p {protein_top} -o {em_tpr} -maxwarn 2"
    run_command(command)
    command = f"gmx mdrun -v -deffnm {os.path.join(target_dir, 'em')}"
    run_command(command)

    # Step 6: NVT Equilibration (constant temperature)
    nvt_tpr = os.path.join(target_dir, "nvt.tpr")
    command = f"gmx grompp -f {nvt_mdp} -c {em_gro} -p {protein_top} -o {nvt_tpr} -maxwarn 2"
    run_command(command)
    command = f"gmx mdrun -v -deffnm {os.path.join(target_dir, 'nvt')}"
    run_command(command)

    # Step 7: NPT Equilibration (constant pressure)
    npt_tpr = os.path.join(target_dir, "npt.tpr")
    command = f"gmx grompp -f {npt_mdp} -c {os.path.join(target_dir, 'nvt.gro')} -p {protein_top} -o {npt_tpr} -maxwarn 2"
    run_command(command)
    command = f"gmx mdrun -v -deffnm {os.path.join(target_dir, 'npt')}"
    run_command(command)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run GROMACS MD simulation pipeline with explicit solvation using AMBER force field")
    parser.add_argument("pdb_file", help="Input PDB file")
    parser.add_argument("target_dir", help="Directory to save output files")
    
    args = parser.parse_args()

    run_gromacs_pipeline(args.pdb_file, args.target_dir)
