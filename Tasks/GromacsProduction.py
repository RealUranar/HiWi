import sys, os, shutil, glob
sys.path.append("../HPC_Jobs/")
import numpy as np
from Sbatch import JobScripts
from task import Task
from excel import Excel
import subprocess

class GromacsProd(Task):
    def __init__(self,job):
        super().__init__(job)
        self.newPath = f"{self.job.location}Gromacs"


    def writeInputFile(self, temp =310):
        with open(f"{self.newPath}/prod.mdp","w") as file:
            file.writelines([
                "; Run parameters\n",
                "integrator      = md            ; leap-frog integrator\n",
                "nsteps          = 40000000       ; 2 * 5000000 = 10000000 ps\n",
                "dt              = 0.002         ; 2 fs\n",
                "define              = -DPOSRES\n",
                "; Output control\n",
                "nstxout         = 10000         ; save coordinates every 1.0 ps\n",
                "nstvout         = 10000         ; save velocities every 1.0 ps\n",
                "nstenergy       = 2000          ; save energies every 1.0 ps\n",
                "nstlog          = 10000         ; update log file every 1.0 ps\n\n",
                "; Output frequency and precision for xtc file\n",
                "nstxtcout       = 1000\n",
                "xtc-precision   = 1000\n",
                "; This selects the subset of atoms for the xtc file. You can select multiple groups. By default all atoms will be written.\n\n",
                "xtc-grps        = F1\n\n",
                "; Selection of energy groups\n",
                "energygrps      = F1\n\n",
                "; Bond parameters\n",
                "continuation            = no            ; Restarting after NVT\n",
                "constraint_algorithm    = lincs         ; holonomic constraints\n",
                "constraints             = h-bonds       ; all bonds (even heavy atom-H bonds) constrained\n",
                "lincs_iter              = 1             ; accuracy of LINCS\n",
                "lincs_order             = 4             ; also related to accuracy\n\n",
                "; Neighborsearching\n",
                "cutoff-scheme   = Verlet\n",
                "ns_type         = grid          ; search neighboring grid cells\n",
                "nstlist         = 10            ; 20 fs, largely irrelevant with Verlet scheme\n",
                "rcoulomb        = 0.9           ; short-range electrostatic cutoff (in nm)\n",
                "rvdw            = 0.9           ; short-range van der Waals cutoff (in nm)\n\n",
                "; Electrostatics\n",
                "coulombtype     = PME           ; Particle Mesh Ewald for long-range electrostatics\n",
                "pme_order       = 4             ; cubic interpolation\n",
                "fourierspacing  = 0.16          ; grid spacing for FFT\n",
                "rlist           = 1.0\n\n",
                "; Temperature coupling is on\n",
                "tcoupl          = V-rescale     ; modified Berendsen thermostat\n",
                "tc-grps         = System        ;\n",
                "tau_t           = 0.1           ; time constant, in ps\n",
                f"ref_t           = {temp}          ; reference temperature, one for each group, in K\n\n",
                "; Periodic boundary conditions\n",
                "pbc             = xyz           ; 3-D PBC\n",
                "; Dispersion correction\n",
                "DispCorr        = EnerPres      ; account for cut-off vdW scheme\n",
                "; Velocity generation\n",
                "gen_vel         = yes           ; Velocity generation is on\n",
                "gen_seed        = -1\n"

            ])

        with open(f"{self.newPath}/plumed.dat", "w") as file:
            file.writelines([
                "UNITS LENGTH=A TIME=0.001\n",
                "t: TORSION ATOMS=279,278,277,276\n",
                "a: ALPHABETA ATOMS1=279,278,277,276 REFERENCE=3.14\n\n",
                "METAD ...\n",
                "LABEL=metad\n",
                "ARG=t\n",
                "PACE=200\n",
                "HEIGHT=1.0\n",
                "SIGMA=0.1\n",
                "GRID_MIN=-pi\n",
                "GRID_MAX=pi\n",
                "GRID_BIN=100\n",
                "CALC_RCT\n",
                "FILE=HILLS\n",
                "BIASFACTOR=60\n",
                f"TEMP={temp:.1f}\n",
                "... METAD\n\n",
                "PRINT FILE=COLVAR ARG=t,a,metad.*\n",
                "FLUSH STRIDE=1\n",
            ])


    def generateJobScript(self):
        with open(f"{self.newPath}/prod.sh","w") as file:
            file.writelines([
            "#!/usr/local_rwth/bin/zsh\n",
            "module load GCC/11.2.0 OpenMPI/4.1.1 GROMACS/2021.5-PLUMED-2.8.0\n",
            "gmx grompp -f prod.mdp -c nvt.gro -r nvt.gro -p System.top -o prod.tpr\n",
            ])
        JobScripts().writeGausianJob(name = self.job.name, location=self.newPath)



    def submit(self):
        ret = subprocess.run(f"{self.newPath}/prod.sh",
                    capture_output = True, 
                    text = True)
        print(ret)
        return super().submit(self.newPath)
        

if __name__ == "__main__":
    with Excel() as scheduler:
        jobs = scheduler.readJobs()
    
    task = GromacsProd(jobs[0])
    #task.moveFiles()
    task.writeInputFile()
    #task.generateJobScript()
    # task.submit()
