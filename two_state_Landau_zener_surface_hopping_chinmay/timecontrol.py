import os
import shutil
import time
import subprocess
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import fileinput
##################################################################
tar="LZ_var_U"
no_fold=30
inputfil="AFSSH.inp"
script="./run"
folds_in_tar=3
cwd=os.getcwd()


for i in range(1,folds_in_tar+1):
    dest=os.path.join(".",tar,str(i))
    shutil.copy(os.path.join(dest,inputfil),cwd)

    subprocess.call(script)

    for j in range(1,no_fold+1):
        end_path=os.path.join(".",str(j),"ended")
        print(end_path)
        while not os.path.exists(end_path):
            time.sleep(2)
    
    print(os.getcwd(),'a')
    with fileinput.FileInput(inputfil, inplace=True) as file:
        for line in file:
            print(line.replace("2 !! iflow", "3 !! iflow"), end='')
    
    os.chdir(cwd)
    subprocess.call(script)
    

    print(os.getcwd(),'b')

    avr=np.zeros(2)

    for i in range(1,no_fold+1):
        avr=avr+np.loadtxt("./"+str(i)+"/fort.1400")
    
    with open('plz_table.txt', 'a') as f:
            f.write("\n"+str(tar)+" how often is pLZ>0.2:"+str(avr[0]*100/avr[1])+"%")
    shutil.copy("pop.out",dest)
    shutil.copy("current.out",dest)
    
    



######################################################################



