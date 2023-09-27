import glob, os, sys, shutil
sys.path.append("Modules/")
from excel import Excel
from job import Job
from jobFinished import HasFinished
import writeJobFile

def submitJob(folder):
    import subprocess
    #shell = subprocess.run(["sbatch", "run_job.sh"], cwd=folder)
    shell = subprocess.run(["echo", "TEST"], cwd=folder)
    print(shell.returncode)

def setupNewCalculation(NewFolder):
    firstCalculationfolder = "Orca_opt"
    with open(f"{NewFolder}Input", "r") as InputFile:
        Name = InputFile.readline().split("=")[1].strip()
    with Excel() as schedule:
        schedule.createJob(Name, location = f"Calculations/{Name}/") #Setup New Job in the Excel Spreadsheet

    os.makedirs(f"Calculations/{Name}/{firstCalculationfolder}") # Create new SubFolders
    xyzFilePath = glob.glob(f"{NewFolder}*.xyz")[0]  #Get Path to the xyz File
    shutil.copy(xyzFilePath, f"Calculations/{Name}/{firstCalculationfolder}/start_molecule.xyz") #Copy xyz File to new directory
    shutil.copy(f"{NewFolder}Input", f"Calculations/{Name}/Input")
    shutil.rmtree(NewFolder, ignore_errors = False) # Delete Input Folder

#Handle new Jobs
for NewFolder in glob.glob("Input/*/"):
    setupNewCalculation(NewFolder)

#Get new Jobs to set up
with Excel() as scheduler:
    jobs = scheduler.readJobs()

for job in jobs:
    for runningJob in job.getRunningJobs():
        match runningJob:
            case "Orca_Opt":
                if HasFinished.orca(f"{job.location}{runningJob}/"):
                    job.updateJob(Orca_Opt = 1, Orca_Dihedral = 3, Gaussian = 3)
                
            case "Orca_Dihedral":
                if HasFinished.orcaDihedral(f"{job.location}{runningJob}/"):
                    job.updateJob(Orca_Dihedral = 1)
            case "Gaussian":
                if HasFinished.gaussian(f"{job.location}{runningJob}/"):
                    job.updateJob(Gaussian = 1)

    with Excel() as scheduler:
        scheduler.updateRow(job)
    