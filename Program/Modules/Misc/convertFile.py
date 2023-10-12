def convertFile(inFile:str, inFormat = "gro", outFormat = "xyz"):
    """_summary_

    Args:
        inFile (_type_): A string representing the whole input File
        inFormat (str, optional): _description_. Defaults to "gro".
        outFormat (str, optional): _description_. Defaults to "xyz".

    Returns:
        _type_: _description_
    """    
    from openbabel import openbabel
    obConversion = openbabel.OBConversion()
    obConversion.SetInAndOutFormats(inFormat, outFormat)
    outFile = openbabel.OBMol()
    obConversion.ReadString(outFile, inFile)
    return obConversion.WriteString(outFile)