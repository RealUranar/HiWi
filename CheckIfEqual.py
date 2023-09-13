import numpy as np


file1 = "lobsteroutQE"
file2 = "lobsteroutIntern"

def FindKeyword(key, lines):
    for i,line in enumerate(list(lines)):
        if line.find(key) != -1:
            return i

def ReadFile(FILENAME):
    with open(FILENAME, "r") as file:
        content =  file.readlines()
        START = FindKeyword("HALLO ICH BIN EIN TEST!!!", content)
        content = content[START:]
        return content
        # content = file.read()
        # start = content.find("HALLO ICH BIN EIN TEST!!!")
        # fileQE = content[start+27:].split("\n")
        # fileQE = fileQE[:len(fileQE)-4]


QE_File = ReadFile(file1)[:200]
INTERN_File = ReadFile(file2)[:200]

def StringToComplex(line):
    try:
        line_new = line.replace("(", "").replace(")", "").split(",")
        i, j = float(line_new[0]), float(line_new[1])
        return complex(i,j)
    except:
        return line



with open("Unequal.txt", "w") as file:
    for QE, INTERN in zip(QE_File, INTERN_File):
        if QE != INTERN:
            QE_neu = [StringToComplex(i) for i in QE.split()]
            INTERN_neu = [StringToComplex(i) for i in INTERN.split()]
            if type(QE_neu[0]) == type(complex(1,2)):
                for Q, I in zip(QE_neu, INTERN_neu):
                    if Q != I:
                        if np.abs(Q-I) > 0.0000001:
                            file.write(f'Fehler: {Q-I}\n Bei: {Q} != {I}\n\n')
            else:
                if QE.split() != INTERN.split():
                    file.write(f'{QE} != {INTERN}\n\n')
