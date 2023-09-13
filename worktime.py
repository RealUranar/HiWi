import time
import pandas as pd


def timeToString(time):
    return f"{time[2]}.{time[1]}.{time[0]} {time[3]}:{time[4]}"

fileName = "test.ods"

print(timeToString(time.localtime()))
print(pd.to_datetime(timeToString(time.localtime())))
df = pd.read_excel(fileName)
print(df["Start"])