import os, shutil, glob
from Sbatch import JobScripts
from task import Task



class Orca_opt(Task):
    def __init__(self,job):
        super().__init__(job)
        self.newPath = f"{self.job.location}Orca_Opt"
        self.executionOrder = [self.moveFiles,
                               self.writeInputFile,
                               self.generateJobScript,
                               self.submit]

    def moveFiles(self):
        os.mkdir(self.newPath) # Create new SubFolders
        xyzFilePath = glob.glob(f"{self.job.location}/*.xyz")[0]  #Get Path to the xyz File
        shutil.copy(xyzFilePath, f"{self.newPath}/startMolecule.xyz") #Copy xyz File to new directory
        os.remove(xyzFilePath)

    def writeInputFile(self):
        dihedral = self._findNNDihedral()[0]
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
        JobScripts().writeOrcaJob(name = self.job.id, location=self.newPath)

    def submit(self):
        self.job.updateJob(Orca_Opt = 2)
        return super().submit(self.newPath)
        
    def isFinished(self):
        hasFinished, succesfull = False, False
        tail = self._readTail(self.newPath)
        hasFinished = "TOTAL RUN TIME:" in tail or ".... aborting the run" in tail
        succesfull = "****ORCA TERMINATED NORMALLY****" in tail
        if hasFinished == True:
            if succesfull == True:
                self.job.updateJob(Orca_Opt = 1, Orca_Dihedral = 3)
                print(f"Orca opt Job {self.job.name} has finished succesfull")
            else:
                self.job.updateJob(Orca_Opt = -1)
                print(f"Orca opt Job {self.job.name} run into a problem")
        else:
            print(f"Orca opt Job {self.job.name} is still running")

if __name__ == "__main__":
    from excel import Excel
    with Excel() as scheduler:
        jobs = scheduler.readJobs()
    
    task = Orca_opt(jobs[0])
    #task.moveFiles()
    task.writeInputFile()
    task.generateJobScript()
    #task.submit()
