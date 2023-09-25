import glob, sys
sys.path.append("~/miniconda3/lib/python3.11/site-packages/openpyxl/")
sys.path.append("~/miniconda3/lib/python3.11/site-packages/")
import pandas as pd

class Excel():
    def __init__(self, fileName = "schedule.xlsx") -> None:
        self.fileName = fileName
        pass

    def __enter__(self):
        excelFile = glob.glob(self.fileName)
        if excelFile == []:
            self.table = pd.DataFrame(columns = ["ID", "Name", "Orca_Opt", "Orca_dihedral", "Gaussian", "Gromax"])
        else:
            self.table = pd.read_excel(excelFile[0], index_col="pandas_Index")
        return self
    
    def __exit__(self, *args):
        self._ConditionalFormatting()
        self.table.to_excel(self.fileName ,index_label="pandas_Index", index=1)
        
    def createJob(self, name):
        ID = 1 if len(self.table["ID"]) == 0 else self.table.index.max() + 2
        new = pd.DataFrame(data = {"ID": pd.Series([ID], index= [ID-1]),
                                "Name": name,
                                "Orca_Opt": 3,
                                "Orca_dihedral": 0,
                                "Gaussian" : 3,
                                "Gromax" : 0})
        self.table = pd.concat([self.table, new], ignore_index=False)

    def UpdateJob(self,ID, location, value):
        self.table.at[ID-1, location] = value

    def deleteJob(self, index):
        try:
            self.table = self.table.drop([index-1])
        except:
            print("Column does not exist!")
            
    def getNextJobs(self):
        Jobs = {}
        for row in self.table.iterrows():
            jobs = []
            for job in row[1].keys()[2:]:
                if row[1][job] == 2:
                    jobs.append("RUNNING")
                    break
                elif row[1][job] == -1:
                    jobs.append("ERROR!")
                    break
                elif row[1][job] == 3:
                    jobs.append(job)

            Jobs[row[1]["ID"]] = {"Name": row[1]["Name"],
                                  "Next_Job": jobs}
        return Jobs
    
    def _ConditionalFormatting(self):
        def colorCode(row):
            return ["background-color: green" if val == 1
                    else "background-color: yellow" if val == 2
                    else "background-color: blue" if val == 3
                    else "background-color: red" if val == -1
                    else '' for val in row]
        self.table = self.table.style.apply(colorCode, axis=1, subset=pd.IndexSlice[:, ["Orca_Opt", "Orca_dihedral", "Gaussian", "Gromax"]])
        return

if __name__ == "__main__":
    with Excel("test.xlsx") as scheduler:
        scheduler.createJob("Test")
        





