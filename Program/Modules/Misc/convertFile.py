def convertFile(filePath, inFormat = "gro", outFormat = "xyz"):
    with open(filePath, "r") as file:
        inFile = file.read()

    from openbabel import openbabel
    obConversion = openbabel.OBConversion()
    obConversion.SetInAndOutFormats(inFormat, outFormat)
    outFile = openbabel.OBMol()
    obConversion.ReadString(outFile, inFile)
    return obConversion.WriteString(outFile)