import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import datetime
import os

def central_difference(x, t):
    n = len(x)
    dx_dt = np.zeros(n)
    for i in range(1, n-1):
        dx = x[i+1] - x[i-1]
        dt = (t[i+1] - t[i-1])
        dx_dt[i] = dx / dt
    dt0 = (t[1] - t[0])
    dx_dt[0] = (x[1] - x[0]) / dt0
    dt_n = (t[n-1] - t[n-2])
    dx_dt[n-1] = (x[n-1] - x[n-2]) / dt_n
    return dx_dt


def parse(allData):
    for i in range(0,int(len(allData))):

        if(i%2==0):
            m = pd.read_csv(path+allData[i])
            f = pd.read_csv(path+allData[i+1])
            
            t = ((pd.to_datetime(m['Timestamp'].values)).asi8 - (pd.to_datetime(m['Timestamp'].values)).asi8[0])/10**9
            x = (m['Position'].values.astype(float)-3500)*-0.074/(18000-3500)
            v = np.ones_like(t)*0.3
            F = f['Force'].values.astype(float)
            dxdt = central_difference(x, t)

            minimum = min(len(t), len(F))

            F = F[:minimum]
            t = t[:minimum]
            x = x[:minimum]
            dxdt = dxdt[:minimum]
            v = v[:minimum]

            data = {
                'Time': t,
                'Position': x,
                'Velocity': dxdt,
                'Voltage': v,
                'Force' : F
            }
            df = pd.DataFrame(data)


            df.to_csv(path_merged+allData[i][0:len(allData[i])-4]+"_merged.csv", index=False)

# def clean(allData,num):
#     dataPath = os.getcwd()
#     MR_data = pd.read_csv(dataPath+'/v325000_a150_121mv_merged.csv')

#     dup_mask = MR_data['Time'].diff().shift(-1) < 0.01
#     MR_data = MR_data[~dup_mask].head(num)

#     MR_data.to_csv(dataPath+"/v325000_a150_121mv_clean.csv")


if __name__ == "__main__":
    path = os.getcwd()+"/data_collection_031124/"
    allData = os.listdir(path)
    path_merged = os.getcwd()+"/data_merged_031124/"

    if(not os.path.exists(path_merged)):
       os.mkdir(path_merged)
       
    parse(allData)
    