import os, shutil, glob, sys
import subprocess
from task import Task

import numpy as np

class GromacsEnergy(Task):
    def __init__(self,job):
        super().__init__(job)
        self.newPath = f"{self.job.location}Gromacs"
        self.executionOrder = [self.moveFiles,
                               self.writeInputFile,
                               self.generateJobScript,
                               self.submit]

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


    def writeInputFile(self):
        self._getEnergys()
        self._changeGROFile()
        self._changeTOPFile()
        self._writePosRe()
        shutil.copy("Modules/GromacsScripts/em.mdp", f"{self.newPath}")
        shutil.copy("Modules/GromacsScripts/table_fourier.itp", f"{self.newPath}")

    def generateJobScript(self):
        with open(f"{self.newPath}/em.sh","w") as file:
            file.writelines([
            "#!/usr/local_rwth/bin/zsh\n",
            "module load GCC/11.2.0 OpenMPI/4.1.1 GROMACS/2021.5-PLUMED-2.8.0\n",
            "gmx grompp -f em.mdp -c System.gro -p System.top -r System.gro -o em.tpr\n",
            "gmx mdrun -v -deffnm em -tableb table_d0.xvg\n"
            ])
        os.chmod(f"{self.newPath}/em.sh", 0o755)


    def submit(self):
        ret = subprocess.run(f"./em.sh",
                    capture_output = True, 
                    text = True,
                    cwd=self.newPath)
        
        if ret.returncode != 0:
             self.job.updateJob(GromacsEnergy = -1)
        else: 
            self.job.updateJob(GromacsEnergy= 1, GromacsEquil= 3)
        print(f"Gromacs energy minimization returned code: {ret.returncode}")
        
    def _writeTableFourier(self):
        dihedral = self._findNNDihedral(xyzFile="Orca_Opt/orca.xyz")
        with open(f"{self.newPath}/table_fourier.itp", "w") as file:
            file.writelines([
                "; ai    aj    ak    al  funct   n   k\n",
                f"{dihedral[0]+1}   {dihedral[1]+1}   {dihedral[2]+1}   {dihedral[3]+1}       8       0   1   \n"
            ])
#         


    def _writePosRe(self):
        sys.path.append("Modules/Misc")
        from InputFileReader import Reader

        with open(f"{self.job.location}Orca_Opt/orca.xyz", "r") as file:
            atomsInStartAtom = int(file.readline())
        
        input = Reader.readInputFile(f"{self.job.location}Input")
        freezeAtoms = input["freezeFragmentAt"] + atomsInStartAtom #calculate the Atoms to freeze in the combined molecule

        with open(f"{self.newPath}/posre.itp", "w") as file:
            file.writelines([
                "[ position_restraints ]\n",
                "; atom  type      fx      fy      fz\n"])
            for atom in freezeAtoms:
                file.writelines([
                    f"{atom}      1       1000000    1000000    1000000\n"
                ])

    
    def _changeTOPFile(self):
        with open(f"{self.newPath}/System.top", "r") as file:
            lines = file.readlines()

        temp = []
        for i, line in enumerate(lines):
            if "[ bonds ]" in line:
                temp.insert(i-1, '\n#ifdef POSRES\n#include "posre.itp"\n#endif\n')
            elif "[ system ]" in line:
                temp.insert(i, '\n#include "table_fourier.itp"\n')
            temp.append(line)

        with open(f"{self.newPath}/System.top", "w") as file:
            file.writelines(temp)

    def _changeGROFile(self):
        with open(f"{self.newPath}/System.gro", "r") as file:
            lines = file.readlines()

        old = np.array(lines[-1].split(), dtype=float)
        new = old * np.array([2,0.90,0.90])
        new = np.array([x if x>2.15 else 2.15 for x in new])  #THIS IS A HACK TO MAKE GROMACS HAPPY!!! CHANGE THIS!!!
        lines[-1] = f"   {new[0]:.5f}   {new[1]:.5f}   {new[2]:.5f}\n"

        with open(f"{self.newPath}/System.gro", "w") as file:
            file.writelines(lines)


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
    sys.path.append("Modules/Misc")
    from excel import Excel
    with Excel() as scheduler:
        jobs = scheduler.readJobs()
    
    task = GromacsEnergy(jobs[0])
    task._writePosRe()
    task._writeTableFourier()
    #task.moveFiles()
    #task.writeInputFile()
    #task.generateJobScript()
    # task.submit()
