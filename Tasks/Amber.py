import sys, os, shutil, glob
sys.path.append("../HPC_Jobs/")
import subprocess
from Sbatch import JobScripts
from task import Task
from excel import Excel




class Amber(Task):
    def __init__(self,job):
        super().__init__(job)
        self.newPath = f"{self.job.location}Amber"

    def moveFiles(self):
        os.mkdir(self.newPath) # Create new SubFolders
        gaussOutFile = glob.glob(f"{self.job.location}Gaussian/*.log")[0]  #Get Path to the log file
        shutil.copy(gaussOutFile, f"{self.newPath}/") #Copy log File to new directory

    def writeInputFile(self):
        with open(f"{self.newPath}/tleap.in", "w") as file:
                file.writelines([
                    "#tleap.in\n",
                    "source leaprc.gaff\n",
                    "#source leaprc.water.tip3p\n",
                    "loadAmberPrep orca.prep\n",
                    "loadamberparams orca.frcmod\n"
                    "SYS = loadpdb NEWPDB4x4.PDB\n",
                    "SaveAmberParm SYS System.prmtop System.inpcrd\n",
                    "quit",
                ])

    def submit(self):
        ret = subprocess.run("/home/seifert/src/HiWi/Tasks/amber.sh",
                            capture_output = True, 
	                        text = True)
        print(ret)


            

if __name__ == "__main__":
    with Excel() as scheduler:
        jobs = scheduler.readJobs()
    
    task = Amber(jobs[0])
    task.moveFiles()