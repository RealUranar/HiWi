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
            
            key, arg = line.replace(" ", "").split("=")
            try:
                arg = int(arg)
            except:
                pass

            keyWords[key] = arg
        
        return keyWords




if __name__ == "__main__":
    reader = Reader()
    ret = reader.readInputFile("Input")
    print(ret)