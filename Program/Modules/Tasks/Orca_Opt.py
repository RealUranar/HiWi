import os, shutil, glob, sys
sys.path.append("Modules/Misc")
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
        print(self.newPath)
        os.mkdir(self.newPath) # Create new SubFolders
        xyzFilePath = glob.glob(f"{self.job.location}/*.xyz")[0]  #Get Path to the xyz File
        shutil.copy(xyzFilePath, f"{self.newPath}/startMolecule.xyz") #Copy xyz File to new directory
        os.remove(xyzFilePath)

    def writeInputFile(self):
        with open(f"{self.newPath}/startMolecule.xyz","r") as file:
            structure = file.read()
        dihed = self._findSubstring(smilesString="*N=N*", inStructure=structure)[0]
        with open(f"{self.newPath}/orca.inp","w") as file:
            
            file.writelines(
                ["! BP def2-SVP def2/J Opt UKS\n",
                 "%pal\n",
                 "nprocs 16\n",
                 "end\n\n",
                 "%geom\n",
                 "Constraints\n"
                 "{"+ f'D {dihed[0]-1} {dihed[1]-1} {dihed[2]-1} {dihed[3]-1} 270.0 C ' + "}\n",
                 "end\nend\n\n",
                 "* xyzfile 0 1 startMolecule.xyz\n"]
            )

    def generateJobScript(self):
        JobScripts().writeOrcaJob(name = self.job.id, location=self.newPath)

    def submit(self):
        self.job.updateJob(runningtasks = self.job.getNextTask()[0])
        return super().submit(self.newPath)
        
    def isFinished(self):
        tail = self._readTail(self.newPath)
        hasFinished = "TOTAL RUN TIME:" in tail or ".... aborting the run" in tail
        succesfull = "****ORCA TERMINATED NORMALLY****" in tail
        if hasFinished == True:
            if succesfull == True:
                self.job.updateJob(finnishedtasks = self.job.getNextTask()[0])
                print(f"Orca opt Job {self.job.name} has finished succesfull")
            else:
                self.job.updateJob(failedtasks = self.job.getNextTask()[0])
                print(f"Orca opt Job {self.job.name} run into a problem")
        else:
            print(f"Orca opt Job {self.job.name} is still running")

if __name__ == "__main__":
    import sys
    sys.path.append("Modules/Misc")
    from job import Job
    job = Job(name = "Test", id = 666, location="Calculations/TESTING/", tasks={"Amber":1})
    
    task = Orca_opt(job)
    #task.moveFiles()
    task.writeInputFile()
    #task.generateJobScript()
    #task.submit()
