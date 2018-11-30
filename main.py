# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import os
import glob
import math
import matplotlib.pyplot as plt
from tkinter import filedialog
from tkinter import *
import subprocess

# DESIGNER folder and path
designerDir = '/Volumes/vdrive/helpern_users/dhiman_s/Projects/DESIGNER_Evaluation/03_Analysis/Scripts/Diffusion-Kurtosis-Imaging/designer'

# Locate DWI root folder
root = Tk()
root.withdraw()
dwi_root = filedialog.askdirectory(parent=root,initialdir="/",title='Please select a directory')
print(dwi_root)

# Obtain list of subjects
subj_list = [d for d in os.listdir(dwi_root) if os.path.isdir(os.path.join(dwi_root, d))]

numSubj = len(subj_list)

print("Preprocessing " + str(numSubj) + " subjects")

os.chdir(designerDir)

subprocess.run('./DESIGNER.py',[dwi_root + '/' + subj_list[0] + '/' + subj_list[0] + '_dwi.nii.gz'])


# Call DESIGNER with input arguments
for x in subj_list
subprocess.call(["python", "DESIGNER.py", [dwi_root + '/' + subj_list[0] + '/' + subj_list[0] + '_dwi.nii.gz'], cwd=designerDir)