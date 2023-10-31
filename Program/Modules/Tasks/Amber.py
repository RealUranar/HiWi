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

    def isFinished(self):
        if len(glob.glob(f"{self.newPath}/SHIFTED.PDB")) == 0:
            return
        
        import parmed as pmd
        try:
            self._runCondaScript(script="PropPDB -p SHIFTED.PDB -o NEWPDB4x4.PDB -ix 1 -iy 4 -iz 4", taskFolder=self.newPath)
            self._runCondaScript(script="tleap -f tleap.in", taskFolder=self.newPath)

            # Save a GROMACS topology and GRO file
            amber = pmd.load_file(f"{self.newPath}/System.prmtop", f"{self.newPath}/System.inpcrd")
            amber.save(f"{self.newPath}/System.top")
            amber.save(f"{self.newPath}/System.gro")

            self.job.updateJob(Amber = 1, GromacsEnergy= 3)
            print(f"Amber Job for {self.job.name} succesfull!")

        except Exception as e:
            print(f"Amber Job Error: {e}")
            self.job.updateJob(Amber = -1)
            print(f"Amber Job for {self.job.name} failed!")

    def submit(self):
        from unitcell import makeUnitcell
        self.job.updateJob(Amber = 2)

        try:
            self._runCondaScript(script="antechamber -fi gout -fo prepi -c resp -i gauss.log -o amber.prep -rn F1 -at gaff2", taskFolder=self.newPath)
            self._runCondaScript(script="parmchk2 -i amber.prep -f prepi -o amber.frcmod", taskFolder=self.newPath)
            
            with open(f"{self.newPath}/NEWPDB.PDB", "r") as file:
                outFile = makeUnitcell(inFile=file.read(), z_ySideLengh=6)
            with open(f"{self.newPath}/CHANGEME.PDB", "w") as file:
                file.write(outFile)
            print(f"Amber job {self.job.name} Part 1 succesfull!\n Shift the cell in 'CHANGEME.PDB' and save as 'SHIFTED.PDB' to continue.")

        except Exception as e:
            print(f"Amber Job Error: {e}")
            self.job.updateJob(Amber = -1)
            print(f"Amber Job for {self.job.name} failed!")



if __name__ == "__main__":
    from excel import Excel
    with Excel() as scheduler:
        jobs = scheduler.readJobs()
    
    task = Amber(jobs[0])
    task.moveFiles()
