class Reader():
    def readInputFile(self, fileName):
        with open(fileName) as file:
            content = file.read()

        keyWords = {}

        for line in content.split("\n"):
            if line.find("#") != -1:
                line = line[:line.find("#")]
            if line == "":
                continue
            
            line = line.replace(" ", "").split("=")

            keyWords[line[0]] = line[1]
        
        return keyWords




if __name__ == "__main__":
    reader = Reader()
    ret = reader.readInputFile("Input")
    print(ret)