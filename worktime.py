import time
import pandas as pd

def timeToString(time):
    return f"{time[2]}.{f'{time[1]}'.rjust(2,'0')}.{time[0]} {time[3]}:{time[4]}:{time[5]}"


fileName = "work.xlsx"

try:
    df = pd.read_excel(fileName)
    index = df.shape[0]
    df.loc[index,"Start"] = pd.to_datetime(timeToString(time.localtime()))
    while True:
        time.sleep(10)
finally:
    df.loc[index,"End"] = pd.to_datetime(timeToString(time.localtime()))
    with pd.ExcelWriter(
        path = fileName,
        datetime_format="dd.mm.yy hh:mm:ss",
        date_format="dd.mm.yy",
        mode = "a",
        if_sheet_exists="overlay"
                        ) as writer:
        df.to_excel(writer, sheet_name="Times"  ,header = True,index =False)



