import os, shutil, glob
from Sbatch import JobScripts
from task import Task



class Orca_Dihedral(Task):
    def __init__(self,job):
        self.subFolders = ["singlet_left", "singlet_right", "triplet_left", "triplet_right"]
        self.spins = [1,3]
        self.angles = ["270, 180", "270, 360"]
        self.steps = 30
        super().__init__(job)
        self.newPath = f"{self.job.location}Orca_Dihedral"
        self.executionOrder = [self.moveFiles,
                        self.writeInputFile,
                        self.generateJobScript,
                        self.submit]

    def moveFiles(self):
        optiFilePath = glob.glob(f"{self.job.location}/Orca_Opt/orca.xyz")[0]  #Get Path to the optimized structure

        for subFolder in self.subFolders:
            os.makedirs(f"{self.newPath}/{subFolder}") # Create new SubFolders
            shutil.copy(optiFilePath, f"{self.newPath}/{subFolder}/orcaOpt.xyz") #Copy optimized structure to new directory
            

    def writeInputFile(self):
        with open(f"{self.job.location}/Orca_Opt/orca.xyz","r") as file:
            structure = file.read()
        dihed = self._findSubstring(smilesString="*N=N*", inStructure=structure)[0]

        folderNr = 0
        for spin in self.spins:
            for angle in self.angles:
                with open(f"{self.newPath}/{self.subFolders[folderNr]}/orca.inp","w") as file:
                    folderNr += 1
                    file.writelines(
                        ["! BP def2-SVP def2/J Opt UKS\n",
                        "%pal\n",
                        "nprocs 16\n",
                        "end\n\n",
                        "%geom Scan\n",
                        f"D {dihed[0]-1} {dihed[1]-1} {dihed[2]-1} {dihed[3]-1} = {angle}, {self.steps}\n",
                        "end\nend\n\n",
                        f"* xyzfile 0 {spin} orcaOpt.xyz\n"]
                    )

    def generateJobScript(self):
        for subFolder in self.subFolders:
            JobScripts().writeOrcaJob(name = self.job.id, location=f"{self.newPath}/{subFolder}")

    def submit(self):
        self.job.updateJob(runningtasks = self.job.getNextTask()[0])
        for subFolder in self.subFolders:
            super().submit(f"{self.newPath}/{subFolder}")
        return        


    def isFinished(self):
        def oracFinished(path):
            tail = self._readTail(path)
            hasFinished = "TOTAL RUN TIME:" in tail or "slurmstepd: error:" in tail
            succesfull = "****ORCA TERMINATED NORMALLY****" in tail
            return hasFinished, succesfull
        
        allDone = True
        for subfolder in self.subFolders:
            hasFinished, succesfull = oracFinished(f"{self.newPath}/{subfolder}/")
            
            if hasFinished:
                if succesfull:
                    print(f"Orca dihedral Job {self.job.name}/{subfolder} has finished succesfull")
                else:
                    self.job.updateJob(failedtasks = self.job.getNextTask()[0])
                    allDone = False
                    print(f"Gromacs Job {self.job.name}/{subfolder} run into a problem")
            else:
                allDone = False
                print(f"Orca dihedral Job {self.job.name}/{subfolder} is still running")
        if allDone:
            self.job.updateJob(finnishedtasks = self.job.getNextTask()[0])
        


if __name__ == "__main__":
    from excel import Excel
    with Excel() as scheduler:
        jobs = scheduler.readJobs()
    
    task = Orca_Dihedral(jobs[0])
    #task.moveFiles()
    task.writeInputFile()
    #task.generateJobScript()
    # task.submit()
