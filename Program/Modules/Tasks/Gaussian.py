import os, shutil, glob, sys

sys.path.append("Modules/Misc")
from InputFileReader import Reader
from Sbatch import JobScripts
from task import Task
from molecule import Molecule

import numpy as np

class Gaussian_opt(Task):
    def __init__(self,job):
        super().__init__(job)
        self.newPath = f"{self.job.location}Gaussian"
        self.inputFile = ""
        self.executionOrder = [self.moveFiles,
                               self.writeInputFile,
                               self.generateJobScript,
                               self.submit]

    def moveFiles(self):
        os.mkdir(f"{self.job.location}/Gaussian") # Create new SubFolders
        optiFilePath = glob.glob(f"{self.job.location}/Orca_Opt/orca.xyz")[0]  #Get Path to the optimized structure
        shutil.copy(optiFilePath, f"{self.newPath}/orca_opt.xyz") #Copy xyz File to new directory


    def writeInputFile(self):
        comFile = self._setupMolecule()
        with open(f"{self.newPath}/combined.com", "w") as file:
                    file.writelines([
                        "%Chk=checkpoint.chk\n",
                        "#P RHF/6-31G* Opt",
                        "\n\nCOMMENT\n\n"])

                    for line in list(comFile.split("\n"))[5:]:
                        file.write(f"{line}\n")

                    file.writelines([
                        "--Link1--\n",
                        "%Chk=checkpoint.chk\n",
                        "#P HF/6-31G* SCF=Tight Geom=AllCheck Guess=Read\n",
                        "Pop=MK IOp(6/33=2, 6/41=10, 6/42=17)\n"
                    ])

    def generateJobScript(self):
        JobScripts().writeGausianJob(name = self.job.id, location=self.newPath)

    def submit(self):
        self.job.updateJob(runningtasks = self.job.getNextTask()[0])
        print(f"Submitted Gaussian job {self.job.name}")
        return super().submit(self.newPath)
        
    def isFinished(self):
        hasFinished, succesfull = False, True
        tail = self._readTail(self.newPath, "*.log")
        hasFinished = "Normal termination of Gaussian" in tail
        #succesfull = "****ORCA TERMINATED NORMALLY****" in str(tail)
            
        if hasFinished:
            if succesfull:
                self.job.updateJob(finnishedtasks = self.job.getRunningTask()[0])
                print(f"Gaussian Job {self.job.name} has finished succesfull")
            else:
                self.job.updateJob(failedtasks = self.job.getRunningTask()[0])
                print(f"Gaussian Job {self.job.name} run into a problem")
        else:
            print(f"Gaussian Job {self.job.name} is still running")


    def _setupMolecule(self):

        inputVars = Reader(f"{self.job.location}Input")

        molecule = Molecule(f"{self.newPath}/orca_opt.xyz")
        molecule.removeAtom(inputVars.getKeyword("RemoveAtomNr")-1)

        fragment = Molecule(f"{self.job.location}fragment.xyz")
        fragment.removeAtom(inputVars.getKeyword("RemoveAtomFragmentNr")-1)  #Everywhere -1 because index starts at 0
        
        combinedMolecules = Molecule.combineMolecules(fragment.getMol(), molecule.getMol(), (inputVars.getKeyword("CombineFragmentAt")-1,inputVars.getKeyword("CombineAtomAt")-1))

        di = self._findSubstring(smilesString="*N=N*", inStructure= combinedMolecules.getXYZBlock())[0]
        
        combinedMolecules.setDihedralAgle(180, np.array(di)-1)#Rotate the molecule to the right dihedral
        comFile = combinedMolecules.getCOMBlock()
        return comFile

if __name__ == "__main__":
    sys.path.append("Modules/Misc")
    from job import Job
    job = Job(name = "Test", id = 666, location="Calculations/TESTING/", tasks={"Gaussian":1})
    
    task = Gaussian_opt(job)
    task.writeInputFile()
    #task._setupMolecule()
    #task.moveFiles()
    #task.generateJobScript()
    #task.submit()
