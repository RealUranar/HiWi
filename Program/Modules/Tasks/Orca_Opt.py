import sys, os, shutil, glob
sys.path.append("../HPC_Jobs/")

from Sbatch import JobScripts
from task import Task
from excel import Excel


class Orca_opt(Task):
    def __init__(self,job):
        super().__init__(job)
        self.newPath = f"{self.job.location}Orca_Opt"

    def moveFiles(self):
        os.mkdir(self.newPath) # Create new SubFolders
        xyzFilePath = glob.glob(f"{self.job.location}/*.xyz")[0]  #Get Path to the xyz File
        shutil.copy(xyzFilePath, f"{self.newPath}/startMolecule.xyz") #Copy xyz File to new directory
        os.remove(xyzFilePath)

    def writeInputFile(self):
        dihedral = self._findNNDihedral()
        with open(f"{self.newPath}/orca.inp","w") as file:
            
            file.writelines(
                ["! BP def2-SVP def2/J Opt UKS\n",
                 "%pal\n",
                 "nprocs 16\n",
                 "end\n\n",
                 "%geom\n",
                 "Constraints\n"
                 "{"+ f'D {" ".join(map(str, dihedral))} 270.0 C ' + "}\n",
                 "end\nend\n\n",
                 "* xyzfile 0 1 startMolecule.xyz\n"]
            )

    def generateJobScript(self):
        JobScripts().writeOrcaJob(name = self.job.name, location=self.newPath)

    def submit(self):
        return super().submit(self.newPath)
        
    def isFinished(self):
        hasFinished, succesfull = False, False
        tail = self._readTail(self.newPath)
        hasFinished = "TOTAL RUN TIME:" in tail
        succesfull = "****ORCA TERMINATED NORMALLY****" in tail
        if succesfull == False:
            self.job.updateJob(Orca_opt = -1)
        if hasFinished == False:
            return
        self.job.updateJob(Orca_opt = 1, Orca_Dihedral = 3)

if __name__ == "__main__":
    with Excel() as scheduler:
        jobs = scheduler.readJobs()
    
    task = Orca_opt(jobs[0])
    #task.moveFiles()
    task.writeInputFile()
    task.generateJobScript()
    #task.submit()
