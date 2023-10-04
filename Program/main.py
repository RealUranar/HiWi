import glob, os, sys, shutil
sys.path.append("Modules/Jobs")
sys.path.append("Modules/Tasks")
from excel import Excel
from job import Job
from OldProgram.jobFinished import HasFinished
from Orca_Opt import Orca_opt
from Orca_Dihedral import Orca_Dihedral
from Gaussian import Gaussian_opt
from GromacsEnergy import GromacsEnergy
from GromacsEquillibration import GromacsEquill
from GromacsProduction import GromacsProd
from Amber import Amber

# Setting up the new Calculation in Excel spread sheet
def setupNewCalculation(NewFolder):
    with open(f"{NewFolder}Input", "r") as InputFile:
        Name = InputFile.readline().split("=")[1].strip()
    with Excel() as schedule:
        schedule.createJob(Name, location = f"Calculations/{Name}/") #Setup New Job in the Excel Spreadsheet

    os.makedirs(f"Calculations/{Name}/") #Create new Folder
    files = glob.glob(f"{NewFolder}/*")
    for file in files:
        shutil.copy(file, f"Calculations/{Name}/")  #Copy files temporarily
    shutil.rmtree(NewFolder, ignore_errors = False) #Remove input folder

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
                if HasFinished.orca(f"{job.location}{runningJob}/")[0]: #If it is done
                    if HasFinished.orca(f"{job.location}{runningJob}/")[1]: #If it is succesfull
                        job.updateJob(Orca_Opt = 1, Orca_Dihedral = 3)
                    else:
                        job.updateJob(Orca_Opt = -1)

            case "Orca_Dihedral":
                if HasFinished.orcaDihedral(f"{job.location}{runningJob}/")[0]:
                    if HasFinished.orcaDihedral(f"{job.location}{runningJob}/")[1]:
                        job.updateJob(Orca_Dihedral = 1, Gaussian = 3)
                    else:
                        job.updateJob(Orca_Dihedral = -1)

            case "Gaussian":
                if HasFinished.gaussian(f"{job.location}{runningJob}/")[0]:
                    if HasFinished.gaussian(f"{job.location}{runningJob}/")[1]:
                        job.updateJob(Gaussian = 1, Gromacs = 3)
                    else:
                        job.updateJob(Gaussian = -1)

            case "Gromacs":
                if HasFinished.gromacs(f"{job.location}{runningJob}/")[0]:
                    if HasFinished.gromacs(f"{job.location}{runningJob}/")[1]:
                        job.updateJob(Gromacs = 1)
                    else:
                        job.updateJob(Gromacs = -1)

    for nextJob in job.getNextJob():
        match nextJob:
            case "Orca_Opt":
                print(f"Starting Orca_Opt for Job id {job.id}")
                task = Orca_opt(job)
                task.moveFiles()
                task.writeInputFile()
                task.generateJobScript()
                task.submit()
                job.updateJob(Orca_Opt = 2)

            case "Orca_Dihedral":
                task = Orca_Dihedral(job)
                task.moveFiles()
                task.writeInputFile()
                task.generateJobScript()
                task.submit()
                job.updateJob(Orca_Dihedral = 2)
                print(f"Starting Orca_Dihedral for Job id {job.id}")

            case "Gaussian":
                task = Gaussian_opt(job)
                task.moveFiles()
                task.writeInputFile()
                task.generateJobScript()
                task.submit()
                job.updateJob(Gaussian = 2)
                print(f"starting Gaussian calc for Job id {job.id}")

            case "Gromacs":
                task = Amber(job)
                task.moveFiles()
                task.writeInputFile()
                task.submit()

                task = GromacsEnergy(job)
                task.moveFiles()
                task.writeInputFile()
                task.generateJobScript()
                task.submit()
                job.updateJob(Gromacs = 2)

                task = GromacsEquill(job)
                task.writeInputFile()
                task.generateJobScript()
                task.submit()

                task = GromacsProd(job)
                task.writeInputFile()
                task.generateJobScript()
                #task.submit()
                print(f"starting Gromacs calc for Job id {job.id}")


    with Excel() as scheduler:
        scheduler.updateRow(job)
    
