import os, shutil
from task import Task
import subprocess

class GromacsEquill(Task):
    def __init__(self,job):
        super().__init__(job)
        self.newPath = f"{self.job.location}Gromacs"
        self.executionOrder = [self.writeInputFile,
                               self.generateJobScript,
                               self.submit]
        
    def writeInputFile(self, temp =310):
        shutil.copy("Modules/GromacsScripts/nvt.mdp", f"{self.newPath}")

    def generateJobScript(self):
        with open(f"{self.newPath}/nvt.sh","w") as file:
            file.writelines([
            "#!/usr/local_rwth/bin/zsh\n",
            "module load GCC/11.2.0 OpenMPI/4.1.1 GROMACS/2021.5-PLUMED-2.8.0\n",
            "gmx grompp -f nvt.mdp -c em.gro -p System.top -r em.gro -o nvt.tpr\n",
            "gmx mdrun -v -deffnm nvt -tableb table_d0.xvg\n"
            ])
        os.chmod(f"{self.newPath}/nvt.sh", 0o755)


    def submit(self):
        print("Running equillibration.... this may take a while (;")
        ret = subprocess.run(f"./nvt.sh",
                    capture_output = True, 
                    text = True,
                    cwd=self.newPath)
        if ret.returncode != 0:
             self.job.updateJob(GromacsEquil = -1)
        else:
            self.job.updateJob(GromacsEquil = 1, GromacsProduction= 3)
        print(f"Gromacs equillibration returned code: {ret.returncode}")
        

if __name__ == "__main__":
    from excel import Excel
    with Excel() as scheduler:
        jobs = scheduler.readJobs()
    
    task = Gromacs(jobs[0])
    #task.moveFiles()
    task.writeInputFile()
    #task.generateJobScript()
    # task.submit()
