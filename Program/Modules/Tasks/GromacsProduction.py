import os, shutil, sys, glob
sys.path.append("Modules/Misc")
from Sbatch import JobScripts
from task import Task
from writePlumed import writePlumed
from InputFileReader import Reader
import numpy as np
import subprocess
from writeMDP import writeMdpFile

class GromacsProd(Task):
    def __init__(self,job):
        Task.__init__(self, job)
        self.newPath = f"{self.job.location}Gromacs"
        self.executionOrder = [self.writeInputFile,
                               self.generateJobScript,
                               self.submit]

    def writeInputFile(self):
        self._getEnergys()
        self._writeTableFourier()

        temp = Reader(f"{self.job.location}Input").getKeyword("temp")
        with open(f"{self.newPath}/prod.mdp", "w") as file:
            file.write(writeMdpFile(job_type="equilibration", temp=temp))
        
        with open(f"{self.job.location}Gromacs/System.gro","r") as file:
            structure = file.read()
        dihedral = self._findSubstring(smilesString="*N=N*" ,inStructure=structure, inFormat="gro")[6]
        with open(f"{self.newPath}/plumed.dat", "w") as file:
            if Reader(f"{self.job.location}Input").getKeyword("tasks")[0] == "rates":
                file.write(writePlumed(dihedral, METAD=True, RATES=True))
            else:
                file.write(writePlumed(dihedral, METAD=True))

    def generateJobScript(self):
        with open(f"{self.newPath}/prod.sh","w") as file:
            file.writelines([
                "#!/usr/local_rwth/bin/zsh\n",
                "module load GCC/11.2.0 OpenMPI/4.1.1 GROMACS/2021.5-PLUMED-2.8.0\n",
                "gmx grompp -f prod.mdp -c nvt.gro -r nvt.gro -p System.top -o prod.tpr\n",
            ])
        os.chmod(f"{self.newPath}/prod.sh", 0o755)
        JobScripts().writeGromacsJob(name = self.job.id, location=self.newPath)
         
    def isFinished(self):
        tail = self._readTail(self.newPath, file = "prod.log")
        hasFinished = "Finished mdrun" in str(tail)
        succesfull = "Constraint error in algorithm" not in str(tail)

        if hasFinished:
            if succesfull:
                self.job.updateJob(finnishedtasks = self.job.getRunningTask()[0])
                print(f"Gromacs Job {self.job.name} has finished succesfull")
            else:
                self.job.updateJob(failedtasks = self.job.getRunningTask()[0])
                print(f"Gromacs Job {self.job.name} run into a problem")
        else:
            print(f"Gromacs Job {self.job.name} is still running")

    def submit(self):
        if Reader(f"{self.job.location}Input").getKeyword("tasks")[0] == "rates":
            with open(f"{self.newPath}/getFrames.sh","w") as file:
                file.writelines([
                    "#!/usr/local_rwth/bin/zsh\n",
                    "module load GCC/11.2.0 OpenMPI/4.1.1 GROMACS/2021.5-PLUMED-2.8.0\n",
                    'echo "0 0" | gmx traj -f nvt.trr -s nvt.tpr -oxt allNVT.gro -dt 4 -b 400\n'
                ])
            os.chmod(f"{self.newPath}/getFrames.sh", 0o755)
            ret = subprocess.run(f"./getFrames.sh",
                capture_output = True,
                text = True,
                cwd=f"{self.job.location}Gromacs")
            
            groFiles = []
            with open(f"{self.job.location}Gromacs/allNVT.gro", "r") as file:
                out = ""
                for line in file.readlines():
                    if "frame" in line and out != "":
                        groFiles.append(out)
                        out = ""
                    out += line
                groFiles.append(out)
            
            groFiles =np.random.choice(groFiles, 30, replace=False)

            os.makedirs(f"{self.job.location}Gromacs_Rates/Base")
            files = ["run_job.sh", "table_fourier.itp", "table_d0.xvg", "posre.itp", "plumed.dat", "System.top", "prod.mdp", "prod.sh"]
            for file in files:
                shutil.copy(f"{self.newPath}/{file}", f"{self.job.location}Gromacs_Rates/Base")
            for i in range(len(groFiles)):
                shutil.copytree(f"{self.job.location}Gromacs_Rates/Base", f"{self.job.location}Gromacs_Rates/{i+1}")
                with open(f"{self.job.location}Gromacs_Rates/{i+1}/nvt.gro", "w") as file:
                    file.write(groFiles[i])

                ret = subprocess.run(f"./prod.sh",
                        capture_output = True,
                        text = True,
                        cwd=f"{self.job.location}Gromacs_Rates/{i+1}")
                super().submit(f"{self.job.location}Gromacs_Rates/{i+1}")
        else:
            ret = subprocess.run(f"./prod.sh",
                        capture_output = True,
                        text = True,
                        cwd=self.newPath)
            print(f"Setup for Gromacs production job {self.job.name} finished with code {ret.returncode}")
            self.job.updateJob(runningtasks = self.job.getNextTask()[0])
            print(f"Submitted Gromacs Production job {self.job.name}")
            return super().submit(self.newPath)
        

    def _writeTableFourier(self):
        with open(f"{self.job.location}Amber/System.gro","r") as file:
            structure = file.read()
        dihedral = self._findSubstring(smilesString="*N=N*" ,inStructure=structure, inFormat="gro")[0]
        
        with open(f"{self.newPath}/table_fourier.itp", "w") as file:
            file.writelines([
                "; ai    aj    ak    al  funct   n   k\n",
                f"{dihedral[0]}   {dihedral[1]}   {dihedral[2]}   {dihedral[3]}       8       0   1   \n"  
            ])

    def _getEnergys(self):
            from scipy.interpolate import CubicSpline
            energyFiles = []
            subfolders = ["singlet_left", "singlet_right", "triplet_left", "triplet_right"]
            for subFolder in subfolders:
                energyFiles.append(glob.glob(f"{self.job.location}Orca_Dihedral/{subFolder}/*relaxscanscf.dat")[0])

            sing_left = np.genfromtxt(energyFiles[0])
            sing_right = np.genfromtxt(energyFiles[1])
            trip_left = np.genfromtxt(energyFiles[2])
            trip_right = np.genfromtxt(energyFiles[3])

            sing = np.append(sing_left[::-1].T, sing_right[1:].T, axis=1)
            trip = np.append(trip_left[::-1].T, trip_right[1:].T, axis=1)

            energy_combined = np.where(sing[1]-trip[1] < 0, sing[1], trip[1])
            phi = sing[0]


            phi_renormalized = phi - phi[0]
            phi_ges = np.append(-phi_renormalized[::-1], phi_renormalized[1:])

            E_normal_in_kJ = (energy_combined - energy_combined.min()) *2625.5  ##Eh in kJ/mol umrechenen
            E_ges = np.append(E_normal_in_kJ, E_normal_in_kJ[-2::-1])

            cs = CubicSpline(phi_ges, E_ges, bc_type='periodic')

            y = CubicSpline.__call__(cs, x = phi_ges, nu=1)
            y_minus = y[::-1]

            np.savetxt(f"{self.newPath}/table_d0.xvg", np.column_stack((phi_ges, E_ges, y_minus)), fmt="%12.8f\t %12.8f\t %12.8f")


if __name__ == "__main__":
    import sys
    sys.path.append("Modules/Misc")
    from job import Job
    job = Job(name = "Test", id = 666, location="Calculations/TESTING/", tasks={"Amber":1})

    task = GromacsProd(job)
    #task.moveFiles()
    #task.writeInputFile()
    #task.generateJobScript()
    # task.submit()
    # with open(f"{task.job.location}Gromacs/System.gro","r") as file:
    #     structure = file.read()
    # dihedral = task._findSubstring(smilesString="*N=N*" ,inStructure=structure, inFormat="gro")
    # print(dihedral)
    #task._changeTOPFile()

    
