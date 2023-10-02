import sys, os, shutil, glob
sys.path.append("../HPC_Jobs/")
import numpy as np
from Sbatch import JobScripts
from task import Task
from excel import Excel
import subprocess

class GromacsEquill(Task):
    def __init__(self,job):
        super().__init__(job)
        self.newPath = f"{self.job.location}Gromacs"

    def writeInputFile(self, temp =310):
        with open(f"{self.newPath}/nvt.mdp","w") as file:
            file.writelines([
                "title               = equilibration NVT\n",
                "define              = -DPOSRES\n",
                "; Run parameters\n",
                "integrator          = md\n",
                "nsteps              = 400000\n",
                "dt                  = 0.001\n",
                "; Output control\n",
                "nstxout             = 1000\n",
                "nstvout             = 1000\n",
                "nstenergy           = 1000\n",
                "nstlog              = 1000\n",
                "; Neighborsearching and short-range nonbonded interactions\n",
                "cutoff-scheme       = Verlet\n",
                "nstlist             = 20\n",
                "ns_type             = grid\n",
                "pbc                 = xyz\n",
                "rlist               = 0.5\n",
                "; Electrostatics\n",
                "coulombtype         = PME\n",
                "pme_order           = 4\n",
                "fourierspacing      = 0.16\n",
                "ewald_rtol          = 1e-05\n",
                "; Temperature coupling\n",
                "tcoupl              = V-rescale\n",
                "tc-grps             = system\n",
                "tau_t               = 0.1\n",
                "ref_t               = 353\n",
                "; Pressure coupling\n",
                "Pcoupl              = no\n",
                "; Velocity generation\n",
                "gen_vel             = yes\n",
                "gen_temp            = 280\n",
                "gen_seed            = -1\n",
                "; Bond parameters\n",
                "constraint_algorithm    = lincs         ; holonomic constraints\n",
                "constraints             = h-bonds       ; all bonds (even heavy atom-H bonds) constrained\n",
                "lincs_iter              = 1             ; accuracy of LINCS\n",
                "lincs_order             = 4             ; also related to accuracy\n",
                "; Simulated annealing\n",
                "annealing	= single 	    ; single sequence of points for each T-coupling group\n",
                "annealing_npoints	= 2		        ; two points - start and end temperatures\n",
                "annealing_time 	= 0 400   	    ; time frame of heating - heat over period of 500 ps\n"
                f"annealing_temp	= 0 {temp}\n"
            ])


    def generateJobScript(self):
        with open(f"{self.newPath}/nvt.sh","w") as file:
            file.writelines([
            "#!/usr/local_rwth/bin/zsh\n",
            "module load GCC/11.2.0 OpenMPI/4.1.1 GROMACS/2021.5-PLUMED-2.8.0\n",
            "gmx grompp -f nvt.mdp -c em.gro -p System.top -r em.gro -o nvt.tpr\n",
            "gmx mdrun -v -deffnm nvt -tableb table_d0.xvg\n"
            ])



    def submit(self):
        ret = subprocess.run(f"{self.newPath}/nvt.sh",
                    capture_output = True, 
                    text = True)
        print(ret)
        

if __name__ == "__main__":
    with Excel() as scheduler:
        jobs = scheduler.readJobs()
    
    task = Gromacs(jobs[0])
    #task.moveFiles()
    task.writeInputFile()
    #task.generateJobScript()
    # task.submit()
