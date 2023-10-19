import os, shutil, sys, glob
sys.path.append("Modules/Misc")
from Sbatch import JobScripts
from task import Task
import numpy as np

import subprocess

class GromacsProd(Task):
    def __init__(self,job):
        super().__init__(job)
        self.newPath = f"{self.job.location}Gromacs"
        self.executionOrder = [self.writeInputFile,
                               self.generateJobScript,
                               self.submit]

    def writeInputFile(self):
        self._getEnergys()
        self._writeTableFourier()
        shutil.copy("Modules/GromacsScripts/prod.mdp", f"{self.newPath}")
        shutil.copy("Modules/GromacsScripts/plumed.dat", f"{self.newPath}")
        
        with open(f"{self.job.location}Gromacs/System.gro","r") as file:
            structure = file.read()
        dihedral = self._findSubstring(smilesString="CN=NC" ,inStructure=structure, inFormat="gro")[6]
        dihedralString = f"{dihedral[0]+1},{dihedral[1]+1},{dihedral[2]+1},{dihedral[3]+1}"
        with open(f"{self.newPath}/plumed.dat", "r") as file:
            lines = file.readlines()
        
        with open(f"{self.newPath}/plumed.dat", "w") as file:
            for line in lines:
                if "ATOMS=" in line:
                    file.write(f"t: TORSION ATOMS={dihedralString}\n")
                elif "ATOMS1=" in line:
                    file.write(f"a: ALPHABETA ATOMS1={dihedralString} REFERENCE=3.14")
                else:
                    file.write(line)

    def generateJobScript(self):
        with open(f"{self.newPath}/prod.sh","w") as file:
            file.writelines([
            "#!/usr/local_rwth/bin/zsh\n",
            "module load GCC/11.2.0 OpenMPI/4.1.1 GROMACS/2021.5-PLUMED-2.8.0\n",
            "gmx grompp -f prod.mdp -c nvt.gro -r nvt.gro -p System.top -o prod.tpr\n",
            ])
        os.chmod(f"{self.newPath}/prod.sh", 0o755)
        JobScripts().writeGromacsJob(name = self.job.id, location=self.newPath)


    def submit(self):
        ret = subprocess.run(f"./prod.sh",
                    capture_output = True,
                    text = True,
                    cwd=self.newPath)
        print(f"Setup for Gromacs production job {self.job.name} finished with code {ret.returncode}")
        self.job.updateJob(GromacsProduction = 2)
        print(f"Submitted Gromacs Production job {self.job.name}")
        return super().submit(self.newPath)
        
        
    def isFinished(self):
        tail = self._readTail(self.newPath, file = "prod.log")
        hasFinished = "Finished mdrun" in str(tail)
        succesfull = "Constraint error in algorithm" not in str(tail)

        if hasFinished:
            if succesfull:
                self.job.updateJob(GromacsProduction = 1)
                print(f"Gromacs Job {self.job.name} has finished succesfull")
            else:
                self.job.updateJob(GromacsProduction = -1)
                print(f"Gromacs Job {self.job.name} run into a problem")
        else:
            print(f"Gromacs Job {self.job.name} is still running")


    def _writeTableFourier(self):
        with open(f"{self.job.location}Amber/System.gro","r") as file:
            structure = file.read()
        dihedral = self._findSubstring(smilesString="CN=NC" ,inStructure=structure, inFormat="gro")[0]
        
        with open(f"{self.newPath}/table_fourier.itp", "w") as file:
            file.writelines([
                "; ai    aj    ak    al  funct   n   k\n",
                f"{dihedral[0]+1}   {dihedral[1]+1}   {dihedral[2]+1}   {dihedral[3]+1}       8       0   1   \n"  
            ])

    def _changeTOPFile(self):
        with open(f"{self.newPath}/System.top", "r") as file:
            lines = file.readlines()

        temp = []
        for i, line in enumerate(lines):
            if "[ system ]" in line:
                temp.insert(i, '\n#include "table_fourier.itp"\n')
            temp.append(line)

        with open(f"{self.newPath}/System.top", "w") as file:
            file.writelines(temp)

    def _getEnergys(self):
            from scipy.interpolate import CubicSpline
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

if __name__ == "__main__":
    import sys
    sys.path.append("Modules/Misc")
    from job import Job
    job = Job(name = "Test", id = 666, location="Calculations/TESTING/", tasks={"Amber":1})

    task = GromacsProd(job)
    #task.moveFiles()
    task.writeInputFile()
    #task.generateJobScript()
    # task.submit()
