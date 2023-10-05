import os, shutil, glob
import subprocess
from task import Task


class Amber(Task):
    """
    A class representing a task to run Amber simulations.

    Attributes:
    job (Job): The job to run the task for.
    newPath (str): The path to the new directory to create for the task.
    executionOrder (list): The order in which to execute the task's methods.

    Methods:
    moveFiles(): Creates a new directory and moves necessary files to it.
    writeInputFile(): Writes an input file for the Amber simulation.
    submit(): Submits the Amber simulation and updates the job status.
    """

    def __init__(self,job):
        super().__init__(job)
        self.newPath = f"{self.job.location}Amber"
        self.executionOrder = [self.moveFiles,
                        self.writeInputFile,
                        self.submit]

    def moveFiles(self):
        os.mkdir(self.newPath) # Create new SubFolders
        gaussOutFile = glob.glob(f"{self.job.location}Gaussian/*.log")[0]  #Get Path to the log file
        shutil.copy(gaussOutFile, self.newPath) #Copy log File to new directory
        shutil.copy("Modules/AmberScripts/amber.sh", self.newPath)
        shutil.copy("Modules/AmberScripts/unitcell.py", self.newPath)
        shutil.copy("Modules/AmberScripts/convert2Gromax.py", self.newPath)
        os.chmod(f"{self.newPath}/amber.sh", 0o755)

    def writeInputFile(self):
        with open(f"{self.newPath}/tleap.in", "w") as file:
                file.writelines([
                    "#tleap.in\n",
                    "source leaprc.gaff\n",
                    "#source leaprc.water.tip3p\n",
                    "loadAmberPrep amber.prep\n",
                    "loadamberparams amber.frcmod\n"
                    "SYS = loadpdb NEWPDB4x4.PDB\n",
                    "SaveAmberParm SYS System.prmtop System.inpcrd\n",
                    "quit",
                ])

    def submit(self):
        self.job.updateJob(Amber = 2)
        ret = subprocess.run(f"./amber.sh",
                            capture_output = True, 
	                        text = True,
                            cwd=self.newPath)
        if ret.returncode != 0:
             self.job.updateJob(Amber = -1)
             
        self.job.updateJob(Amber = 1, GromacsEnergy= 3)
        print(f"Amber file conversion returned code: {ret.returncode}")


if __name__ == "__main__":
    from excel import Excel
    with Excel() as scheduler:
        jobs = scheduler.readJobs()
    
    task = Amber(jobs[0])
    task.moveFiles()
