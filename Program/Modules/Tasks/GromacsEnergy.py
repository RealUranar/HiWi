import os, shutil, glob
import subprocess
from task import Task
from InputFileReader import Reader
from writePlumed import writePlumed
from writeMDP import writeMdpFile


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

        if Reader(f"{self.job.location}Input").getKeyword("tasks")[0] == "rates":  #Write a plumed file to restrain the CNNC dihedral at 0 degrees
            with open(f"{self.job.location}Amber/System.gro","r") as file:
                structure = file.read()
            dihedral = self._findSubstring(smilesString="*N=N*" ,inStructure=structure, inFormat="gro")[6]

            with open(f"{self.newPath}/plumedRestraint.dat", "w") as file:
                file.write(writePlumed(dihedral=dihedral, RESTRAIN=True, PRINT=False))

        with open(f"{self.newPath}/em.mdp", "w") as file:
            file.write(writeMdpFile(job_type="energy_minimization"))

    def generateJobScript(self):
        with open(f"{self.newPath}/em.sh","w") as file:
            file.writelines([
            "#!/usr/local_rwth/bin/zsh\n",
            "module load GCC/11.2.0 OpenMPI/4.1.1 GROMACS/2021.5-PLUMED-2.8.0\n",
            "gmx grompp -f em.mdp -c System.gro -p System.top -r System.gro -o em.tpr\n",
            "gmx mdrun -v -deffnm em -ntmpi 1"
            ])
            if Reader(f"{self.job.location}Input").getKeyword("tasks")[0] == "rates":
                file.write(" -plumed plumedRestraint.dat")
        os.chmod(f"{self.newPath}/em.sh", 0o755)


    def submit(self):
        ret = subprocess.run(f"./em.sh",
                    capture_output = True, 
                    text = True,
                    cwd=self.newPath)
        
        if ret.returncode != 0:
             self.job.updateJob(failedtasks = self.job.getRunningTask()[0])
        else: 
            self.job.updateJob(finnishedtasks = self.job.getRunningTask()[0])
        print(f"Gromacs energy minimization returned code: {ret.returncode}")
        


    def _writePosRe(self):
        with open(f"{self.job.location}Amber/System.gro","r") as file:
            structure = file.read()
        freezeAtoms = self._findSubstring(smilesString="csc" ,inStructure=structure, inFormat="gro")[0]

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
            if "[ system ]" in line:
                temp.insert(i, '\n#ifdef POTENTIAL\n#include "table_fourier.itp"\n#endif\n')
            temp.append(line)

        with open(f"{self.newPath}/System.top", "w") as file:
            file.writelines(temp)

    def _changeGROFile(self):
        with open(f"{self.newPath}/System.gro", "r") as file:
            groFile = file.readlines()

        groFile[-1] = "5.00000   2.40000   2.40000\n"

        with open(f"{self.newPath}/System.gro", "w") as file:
            file.writelines(groFile)


if __name__ == "__main__":
    import sys
    sys.path.append("Modules/Misc")
    from job import Job
    job = Job(name = "Test", id = 666, location="Calculations/TESTING/", tasks={"Amber":1})
    
    task = GromacsEnergy(job)
    task._changeTOPFile()
    #out = task._changeTOPFile()
    #print(out)
    # task._writePosRe()
    # task._writeTableFourier()
    #task.moveFiles()
    #task.writeInputFile()
    #task.generateJobScript()
    # task.submit()
