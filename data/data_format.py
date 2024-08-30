import numpy as np
import pandas as pd
import os

num = 30000
dataPath = os.getcwd()
MR_data = pd.read_csv(dataPath+'/data_merged_030824/v250000_a50_185mv_merged.csv')

dup_mask = MR_data['Time'].diff().shift(-1) < 0.01
MR_data = MR_data[~dup_mask].head(num)

MR_data.to_csv(dataPath+"/v250000_a50_185mv_clean.csv")
