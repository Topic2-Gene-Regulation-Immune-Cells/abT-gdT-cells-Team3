import pandas as pd
import numpy as np
import matplotlib.pyplot as p
import os 

p.ioff()

x = np.arange (0, 2*np.pi, 0.1)
y = np.sin(x)
# p.figure (1)
# p.clf ()
# p.plot (x,y)
# p.show ()

pd_data = pd.DataFrame ({'X':x, 'Y':y})
print(pd_data.head())
data_folder = 'D:\\Uni\\4. Semester\\Bioinfo\\datasets'
file_name = os.path.join (data_folder, 'processed atac seq data and called peaks.csv')
atac_data = pd.read_csv(file_name, sep=',')
print(atac_data.head())

from test2 import myprint
myprint ()