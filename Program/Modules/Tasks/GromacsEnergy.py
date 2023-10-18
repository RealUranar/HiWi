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
        self._changeGROFile()
        self._changeTOPFile()
        self._writePosRe()
        shutil.copy("Modules/GromacsScripts/em.mdp", f"{self.newPath}")

    def generateJobScript(self):
        with open(f"{self.newPath}/em.sh","w") as file:
            file.writelines([
            "#!/usr/local_rwth/bin/zsh\n",
            "module load GCC/11.2.0 OpenMPI/4.1.1 GROMACS/2021.5-PLUMED-2.8.0\n",
            "gmx grompp -f em.mdp -c System.gro -p System.top -r System.gro -o em.tpr\n",
            "gmx mdrun -v -deffnm em -ntmpi 1\n"
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
        


    def _writePosRe(self):
        with open(f"{self.job.location}Amber/System.gro","r") as file:
            structure = file.read()
        freezeAtoms = self._findSubstring(smilesString="CSC" ,inStructure=structure, inFormat="gro")[0]

        with open(f"{self.newPath}/posre.itp", "w") as file:
            file.writelines([
                "[ position_restraints ]\n",
                "; atom  type      fx      fy      fz\n"])
            for atom in freezeAtoms:
                file.writelines([
                    f"{atom+1}      1       1000000    1000000    1000000\n"
                ])

    
    def _changeTOPFile(self):
        with open(f"{self.newPath}/System.top", "r") as file:
            lines = file.readlines()

        temp = []
        for i, line in enumerate(lines):
            if "[ bonds ]" in line:
                temp.insert(i-1, '\n#ifdef POSRES\n#include "posre.itp"\n#endif\n')
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
