import os, shutil, glob
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
            print(f"gro/top-File not Found for Job {self.job.name}")
            raise FileNotFoundError

        shutil.copy(groFilePath, f"{self.newPath}/System.gro")
        shutil.copy(topFilePath, f"{self.newPath}/System.top")


    def writeInputFile(self):
        self._getEnergys()
        self._changeGROFile()
        self._changeTOPFile()
        self._writePosRe()
        self._writeTableFourier()
        shutil.copy("Modules/GromacsScripts/em.mdp", f"{self.newPath}")

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
        dihedral = self._findSubstring(smilesString="CN=NC" ,inFile="Amber/System.gro")[0]
        with open(f"{self.newPath}/table_fourier.itp", "w") as file:
            file.writelines([
                "; ai    aj    ak    al  funct   n   k\n",
                f"{dihedral[0]+1}   {dihedral[1]+1}   {dihedral[2]+1}   {dihedral[3]+1}       8       0   1   \n"  
            ])
#         


    def _writePosRe(self):
        freezeAtoms = self._findSubstring(smilesString="CSC" ,inFile="Amber/System.gro")[0]

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
        from convertFile import convertFile
        from unitcell import makeUnitcell
        with open(f"{self.newPath}/System.gro", "r") as file:
            filePDB = convertFile(file.read(), inFormat="gro", outFormat="pdb")

        newCell = makeUnitcell(filePDB)
        newCellGRO = convertFile(newCell, inFormat="pdb", outFormat="gro")

        with open(f"{self.newPath}/System.gro", "w") as file:
            file.write(newCellGRO)


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
    
    task = GromacsEnergy(job)
    out = task._changeGROFile()
    print(out)
    # task._writePosRe()
    # task._writeTableFourier()
    #task.moveFiles()
    #task.writeInputFile()
    #task.generateJobScript()
    # task.submit()
