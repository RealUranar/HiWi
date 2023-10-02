import sys, os, shutil, glob
sys.path.append("../HPC_Jobs/")
import numpy as np
from Sbatch import JobScripts
from task import Task
from excel import Excel
from scipy.interpolate import CubicSpline
import subprocess

class GromacsEnergy(Task):
    def __init__(self,job):
        super().__init__(job)
        self.newPath = f"{self.job.location}Gromacs"

    def moveFiles(self):
        os.mkdir(self.newPath) # Create new SubFolders
        try:
            groFilePath = glob.glob(f"{self.job.location}Amber/*.gro")[0]
            topFilePath = glob.glob(f"{self.job.location}Amber/*.top")[0]
        except IndexError:
            print("gro/top-File not Found")
            return

        shutil.copy(groFilePath, f"{self.newPath}/System.gro")
        shutil.copy(topFilePath, f"{self.newPath}/System.top")
        self.__getEnergys()

    def writeFiles(self):
        with open(f"{self.newPath}/posre.itp","w") as file:
            file.writelines([
                "[ position_restraints ]\n",
                "; atom  type      fx      fy      fz\n",
                "52      1       1000000    1000000    1000000\n",
                "51      1       1000000    1000000    1000000\n",
                "49      1       1000000    1000000    1000000"
            ])

        with open(f"{self.newPath}/table_fourier.itp","w") as file:
            file.writelines([
                "; ai    aj    ak    al  funct   n   k\n",
                "11   12   13   14       8       0   1 "  
            ])

    def writeInputFile(self):
        with open(f"{self.newPath}/em.mdp","w") as file:
            file.writelines([
                "title = Energy Minimization\n",
                "define = -DPOSRES  ; Define position restraints\n\n",
                "; minim.mdp - used as input into grompp to generate em.tpr\n",
                "integrator	= steep		; Algorithm (steep = steepest descent minimization)\n",
                "emtol		= 100.0  	; Stop minimization when the maximum force < 1000.0 kJ/mol/nm\n",
                "emstep          = 0.01          ; Energy step size\n",
                "nsteps		= 500000	  	; Maximum number of (minimization) steps to perform\n\n",
                "; Parameters describing how to find the neighbors of each atom and how to calculate the interactions\n",
                "nstlist		    = 15                 ; Frequency to update the neighbor list and long range forces\n",
                "cutoff-scheme       = Verlet\n",
                "ns_type		    = grid		; Method to determine neighbor list (simple, grid)\n",
                "coulombtype	    = PME		; Treatment of long range electrostatic interactions\n",
                "rcoulomb	    = 1.0		; Short-range electrostatic cut-off\n",
                "rvdw		    = 1.0		; Short-range Van der Waals cut-off\n",
                "pbc		    = xyz 		; Periodic Boundary Conditions (yes/no)\n\n",
                "; Output control\n",
                "nstxout             = 100\n",
                "nstvout             = 1000\n",
                "nstenergy           = 1000\n",
                "nstlog              = 1000\n"
            ])


    def __getEnergys(self):
        energyFiles = []
        subfolders = ["singlet_left", "singlet_right", "triplet_left", "triplet_right"]
        for subFolder in subfolders:
            energyFiles.append(glob.glob(f"{self.job.location}Orca_Dihedral/{subFolder}/*relaxscanscf.dat")[0])

        sing_left = np.genfromtxt(energyFiles[0])
        sing_right = np.genfromtxt(energyFiles[1])
        trip_left = np.genfromtxt(energyFiles[2])
        trip_right = np.genfromtxt(energyFiles[3])

        sing = np.append(sing_left[::-1].T, sing_right[1:].T, axis=1)
        trip = np.append(trip_left[::-1].T, trip_right[1:].T, axis=1)

        energy_combined = np.where(sing[1]-trip[1] < 0, sing[1], trip[1])
        phi = sing[0]


        phi_renormalized = phi - phi[0]
        phi_ges = np.append(-phi_renormalized[::-1], phi_renormalized[1:])

        E_normal_in_kJ = (energy_combined - energy_combined.min()) *2625.5
        E_ges = np.append(E_normal_in_kJ, E_normal_in_kJ[-2::-1])

        cs = CubicSpline(phi_ges, E_ges, bc_type='periodic')

        y = CubicSpline.__call__(cs, x = phi_ges, nu=1)
        y_minus = y[::-1]

        np.savetxt(f"{self.newPath}/table_d0.xvg", np.column_stack((phi_ges, E_ges, y_minus)), fmt="%12.8f\t %12.8f\t %12.8f")

    def generateJobScript(self):
        with open(f"{self.newPath}/em.sh","w") as file:
            file.writelines([
            "#!/usr/local_rwth/bin/zsh\n",
            "module load GCC/11.2.0 OpenMPI/4.1.1 GROMACS/2021.5-PLUMED-2.8.0\n",
            "gmx grompp -f em.mdp -c System.gro -p System.top -r System.gro -o em.tpr\n",
            "gmx mdrun -v -deffnm em -tableb table_d0.xvg\n"
            ])

    def submit(self):
        ret = subprocess.run(f"{self.newPath}/em.sh",
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
