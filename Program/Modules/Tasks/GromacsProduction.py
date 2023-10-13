import os, shutil, sys
sys.path.append("Modules/Misc")
from Sbatch import JobScripts
from task import Task

import subprocess

class GromacsProd(Task):
    def __init__(self,job):
        super().__init__(job)
        self.newPath = f"{self.job.location}Gromacs"
        self.executionOrder = [self.writeInputFile,
                               self.generateJobScript,
                               self.submit]

    def writeInputFile(self):
        shutil.copy("Modules/GromacsScripts/prod.mdp", f"{self.newPath}")
        shutil.copy("Modules/GromacsScripts/plumed.dat", f"{self.newPath}")
        
        dihedral = self._findSubstring(smilesString="CN=NC" ,inFile=f"Gromacs/System.gro")[6]
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
