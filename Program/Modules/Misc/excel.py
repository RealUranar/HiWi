import glob, sys
sys.path.append("Modules/")
from dataclasses import dataclass
import pandas as pd
from job import Job


class Excel():
    def __init__(self, fileName = "schedule.xlsx") -> None:
        self.fileName = fileName
        pass

    def __enter__(self):
        excelFile = glob.glob(self.fileName)
        if excelFile == []:
            self.table = pd.DataFrame(columns = ["location","ID", "Name", "Orca_Opt", "Orca_Dihedral", "Gaussian", "Amber", "GromacsEnergy", "GromacsEquil", "GromacsProduction"])
        else:
            self.table = pd.read_excel(excelFile[0], index_col="pandas_Index")
        return self
    
    def __exit__(self, *args):
        self._ConditionalFormatting()
        self.table.to_excel(self.fileName ,index_label="pandas_Index", index=1)
        
    def createJob(self, name, location):
        ID = 1 if len(self.table["ID"]) == 0 else self.table.index.max() + 2
        new = pd.DataFrame(data = {"location": location,
                                "ID": pd.Series([ID], index= [ID-1]),
                                "Name": name,
                                "Orca_Opt": 3,
                                "Orca_Dihedral": 0,
                                "Gaussian" : 0,
                                "Amber" : 0,
                                "GromacsEnergy":0,
                                "GromacsEquil" : 0,
                                "GromacsProduction": 0})
        self.table = pd.concat([self.table, new], ignore_index=False)


    def updateRow(self, job:Job):
        ID = job.id
        self.table.iloc[ID-1,3:] = job.getPanda()

        
    def deleteJob(self, index):
        try:
            self.table = self.table.drop([index-1])
        except:
            print("Column does not exist!")


    def readJobs(self) ->list[Job]:
        jobs = []
        for index in self.table.index:
            data = self.table.iloc[index]
            jobs.append(Job(name = data["Name"], id=data["ID"], location=data["location"] ,tasks= dict(data[3:])))
        return jobs
        
            
    def _ConditionalFormatting(self):
        def colorCode(row):
            return ["background-color: green" if val == 1   # Done
                    else "background-color: yellow" if val == 2   #Working
                    else "background-color: blue" if val == 3  #Ready
                    else "background-color: red" if val == -1  #Error
                    else '' for val in row]
        self.table.style.apply(colorCode, axis=1, subset=pd.IndexSlice[:, ["Orca_Opt", "Orca_Dihedral", "Gaussian", "Amber", "GromacsEnergy", "GromacsEquil", "GromacsProduction"]])
        #self.table = self.table.style.apply(colorCode, axis=1, subset=pd.IndexSlice[:, ["Orca_Opt", "Orca_Dihedral", "Gaussian", "Gromax"]])
        return

if __name__ == "__main__":
    with Excel("test.xlsx") as scheduler:
        scheduler.createJob("Test", "TEST2")
        job = scheduler.readJobs()[0]
        job.updateJob(Orca_Opt = 120)
        scheduler.updateRow(job)

        





