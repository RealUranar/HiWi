import sys, os, shutil, glob
sys.path.append("../HPC_Jobs/")

from Sbatch import JobScripts
from task import Task
from excel import Excel


class Orca_Dihedral(Task):
    def __init__(self,job):
        self.subFolders = ["singlet_left", "singlet_right", "triplet_left", "triplet_right"]
        self.spins = [1,3]
        self.angles = ["270, 180", "270, 360"]
        self.steps = 30
        super().__init__(job)
        self.newPath = f"{self.job.location}Orca_Dihedral"

    def moveFiles(self):
        optiFilePath = glob.glob(f"{self.job.location}/Orca_Opt/orca.xyz")[0]  #Get Path to the optimized structure

        for subFolder in self.subFolders:
            os.makedirs(f"{self.newPath}/{subFolder}") # Create new SubFolders
            shutil.copy(optiFilePath, f"{self.newPath}/{subFolder}/orcaOpt.xyz") #Copy optimized structure to new directory
            

    def writeInputFile(self):
        dihedral = " ".join(map(str,self._findNNDihedral()))

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
                        f"D {dihedral} = {angle}, {self.steps}\n",
                        "end\nend\n\n",
                        f"* xyzfile 0 {spin} orcaOpt.xyz\n"]
                    )

    def generateJobScript(self):
        for subFolder in self.subFolders:
            JobScripts().writeOrcaJob(name = self.job.name, location=f"{self.newPath}/{subFolder}")

    def submit(self):
        for subFolder in self.subFolders:
            super().submit(f"{self.newPath}/{subFolder}")
        return        



if __name__ == "__main__":
    with Excel() as scheduler:
        jobs = scheduler.readJobs()
    
    task = Orca_Dihedral(jobs[0])
    #task.moveFiles()
    task.writeInputFile()
    task.generateJobScript()
    # task.submit()
