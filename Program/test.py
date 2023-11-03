import glob, os, sys, shutil

sys.path.append("Modules/Misc")
from InputFileReader import Reader
from excel import Excel
from job import Job

print(Reader("Calculations/TESTING/Input").getKeyword("calcrates"))

if val := Reader("Calculations/TESTING/Input").getKeyword("calcrates"):
    print(val)