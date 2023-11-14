import numpy as np
class Reader():
    def __init__(self,fileName):
        with open(fileName) as file:
            content = file.read()

        keyWordsDict = {}
        for line in content.split("\n"):
            if line.find("#") != -1:
                line = line[:line.find("#")]
            if line == "":
                continue
            
            key, *arg = line.replace(" ", "").split("=")
            if len(arg) == 0:
                arg = key
            else:
                arg = arg[0]

            try:
                arg = int(arg)
            except:
                pass

            keyWordsDict[key.lower()] = arg
        
        self.keyWordsDict = keyWordsDict

    def getKeyword(self,keyword:str):
        basicTasks = [
                     "Orca_opt",
                    "Orca_Dihedral",
                    "Gaussian_opt",
                    "Amber",
                    "GromacsEnergy",
                    "GromacsEquill",
                    "GromacsProd"
                ]
        if keyword.lower() == "tasks":
            tasks = self.keyWordsDict[keyword.lower()].split(",")
            if tasks[0] == "Rates":
                return "Rates", basicTasks
            elif tasks[0] == "Barrier":
                return "Barrier", basicTasks
            else:
                return None, tasks
        try:
            return  self.keyWordsDict[keyword.lower()]
        except KeyError:
            return None
        


if __name__ == "__main__":
    inputFile ="Calculations/TESTING/Input"
    reader = Reader(inputFile)
    ret = reader.getKeyword("Tasks")
    print(ret)