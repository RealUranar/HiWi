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
        try:
            return  self.keyWordsDict[keyword.lower()]
        except KeyError:
            return False
        




if __name__ == "__main__":
    reader = Reader("Input")
    ret = reader.getKeyword("Name")
    print(ret)