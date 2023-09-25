import glob, os, sys, shutil
sys.path.append("Modules/")
import excel

def setupNewCalculation(NewFolder):
    with open(f"{NewFolder}Input", "r") as InputFile:
        Name = InputFile.readline().split("=")[1].strip()

    with excel.Excel() as schedule:
        schedule.createJob(Name) #Setup New Job in the Excel Spreadsheet

    os.makedirs(f"Calculations/{Name}/Orca") # Create ne SubFolders
    xyzFilePath = glob.glob(f"{NewFolder}*.xyz")[0]  #Get Path to the xyz File
    shutil.copy(xyzFilePath, f"Calculations/{Name}/Orca/start_molecule.xyz") #Copy xyz File to new directory
    shutil.copy(f"{NewFolder}Input", f"Calculations/{Name}/Input")
    shutil.rmtree(NewFolder, ignore_errors = True) # Delete Input Folder

#Handle new Jobs
for NewFolder in glob.glob("Input/*/"):
    print(f"TEST1 {NewFolder}")
    setupNewCalculation(NewFolder)

#Get new Jobs to set up
with excel.Excel() as scheduler:
    jobs = scheduler.getNextJobs()

#print(jobs)
for id in jobs.keys():
    for job in jobs[id]["Next_Job"]:
        match job:
            case "ERROR!":
                print(f"Job {id} has crashed")
            case "RUNNING":
                print(f"Job {id} is running")
            case "Orca_Opt":
                print(f"Starting Orca opt for Job {id}")
            case "Orca_dihedral":
                print(f"Starting Orca dihedral scan for Job {id}")
            case "Gaussian":
                print(f"Starting gaussian opt for Job {id}")
            case "Gromax":
                print(f"Starting Gromax for Job {id}")
            case _:
                print("Something went really wrong")
    print()


