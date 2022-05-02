#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
heat_data_saving.py   Apr 2022

This file takes the processed data in through google colab and joins it
to the files processed in the computer

@author: naiacasina
"""
import pandas as pd
import numpy as np
import os
import pickle

os.chdir('/Users/naiacasina/Documents/IDEA SECOND/Sem 3/ENVS/Codes and Data/Python/colab/') 
     
heat_June = np.load("heat_June", allow_pickle=True)
heat_May = np.load("heat_May", allow_pickle=True)
heat_Apr = np.load("heat_Apr", allow_pickle=True)
heat_Nov = np.load("heat_Nov", allow_pickle=True)
heat_Dec = np.load("heat_Dec", allow_pickle=True)
heat_JuneMay = np.add(heat_June,heat_May)
heat_AprJune = np.add(heat_JuneMay,heat_Apr)
heat_data = np.add(heat_AprJune,heat_Nov)
heat_data_final = np.add(heat_data,heat_Dec)

os.chdir('/Users/naiacasina/Documents/IDEA SECOND/Sem 3/ENVS/Codes and Data/Heat/Results/')       
heat_data_final.dump('heat_colab')



