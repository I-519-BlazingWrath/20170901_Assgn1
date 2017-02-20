from __future__ import division
import sys;import math; import re; import subprocess; import os
from random import choice, randrange
import matplotlib.pyplot as plt; import numpy as np; import scipy as sc; from scipy import stats; import scipy.stats.distributions
import csv; import pandas as pd; import itertools


'''def chooseseason(month):
    row = ''
    if 1 <= month <= 3:
        row = 'q1'
    elif 4 <= month <= 6:
        row = 'q2'
    elif 7 <= month <= 9:
        row = 'q3'
    elif 10 <= month <= 12:
        row = 'q4'
    else:
        row = "??"
    return row'''
def chooseseason(month):
    row = ''
    if 4 <= month <= 9:
        row = 'q1'
    elif 1 <= month <= 3 or 10 <= month <= 12:
        row = 'q2'
    else:
        row = "??"
    return row

def makets(thedf, OP):
    thedf['season'] = pd.Series([chooseseason(x) for x in thedf['month']], index=thedf.index)#add a column to the df with
    #month translated into season
    groups = thedf.groupby(['year', 'season'])#group the df by year and season
    temp = list(range(0, len(groups), 1))#make a list from 0 to number.of.groups
    qdf = pd.DataFrame(columns=['year', 'season'], index=temp)#intialize an empty df with the right column names
    for index,guy in enumerate(groups):#for each group in the grouped data
        qdf.loc[index]=[guy[0][0],guy[0][1]]#populate a row of the empty df with the year and season
    #for index, row in enumerate(qdf.iterrows()):
     #   row = [groups[index][0][0], groups[index][0][1]]
    qdf['co2avg']=pd.Series([groups['co2flux'].mean()[x] for x in temp], index=qdf.index)#now add a column with avg
    #soil respiration for each year&season group

    OP.write(pd.DataFrame.to_csv(qdf))
    return qdf

mydf = pd.read_csv("C:\\Users\\rmoge\\GitHub\\QB2017_DivPro\\Data\\HF_soilresp.csv", infer_datetime_format=True)

cdf = mydf.loc[mydf['trt'] == "C"]
# for row in
hdf = mydf.loc[mydf['trt'] == "H"]
hndf = mydf.loc[mydf['trt'] == "HN"]
ndf = mydf.loc[mydf['trt'] == "N"]


try:
    cco2=open("C:\\Users\\rmoge\\GitHub\\QB2017_DivPro\\Data\\cflux.csv",'w')
    hco2=open("C:\\Users\\rmoge\\GitHub\\QB2017_DivPro\\Data\\hflux.csv",'w')
    hnco2=open("C:\\Users\\rmoge\\GitHub\\QB2017_DivPro\\Data\\hnflux.csv",'w')
    nco2=open("C:\\Users\\rmoge\\GitHub\\QB2017_DivPro\\Data\\nflux.csv",'w')
    c=makets(cdf,cco2); h=makets(hdf,hco2); hn= makets(hndf,hnco2); n= makets(ndf,nco2)
    
finally:
    cco2.close(); hco2.close(), hnco2.close(), nco2.close()

grouped = ndf.groupby(['year', 'month'])
meanspermonth = grouped['co2flux'].mean()
tester = list(range(0, len(grouped), 1))#make a list from 0 to number.of.groups


ndfg = ndf.groupby(['year', 'month']).groups
type(ndf.groupby(['year', 'month']).groups)
sl = (ndf['month'])
ndff = ndf
ndff['e'] = pd.Series([chooseseason(x) for x in ndff['month']], index=ndff.index)

print(chooseseason(6))