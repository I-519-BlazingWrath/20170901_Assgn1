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
    qdf['Nminavg']=pd.Series([groups['netnh4'].mean()[x] for x in temp], index=qdf.index)#now add a column with avg
    #soil respiration for each year&season group

    OP.write(pd.DataFrame.to_csv(qdf))
    return qdf

#mydf = pd.read_csv("C:\\Users\\rmoge\\GitHub\\QB2017_DivPro\\Data\\HF_nmin.csv", infer_datetime_format=True)
mydf = pd.read_csv("C:\\Users\\rmoge\\GitHub\\QB2017_DivPro\\Data\\HF_nmin_NA.omit_RZM.csv", infer_datetime_format=True)
'''cdfna = mydf.loc[mydf['trt'] == "C"]
cdf=cdfna.loc[cdfna['netnh4']!="nan"]
hdf = mydf.loc[mydf['trt'] == "H" and mydf['netnh4'] != "NA"]
hndf = mydf.loc[mydf['trt'] == "HN" and mydf['netnh4'] != "NA"]
ndf = mydf.loc[mydf['trt'] == "N" and mydf['netnh4'] != "NA"]'''

cdf = mydf.loc[mydf['trt'] == "C"]
hdf = mydf.loc[mydf['trt'] == "H"]
hndf = mydf.loc[mydf['trt'] == "HN"]
ndf = mydf.loc[mydf['trt'] == "N"]
mincdf=cdf.loc[cdf['hor']=='mineral']
minhdf=hdf.loc[hdf['hor']=='mineral']
minhndf=hndf.loc[hndf['hor']=='mineral']
minndf=ndf.loc[ndf['hor']=='mineral']
orgcdf=cdf.loc[cdf['hor']=='organic']
orghdf=hdf.loc[hdf['hor']=='organic']
orghndf=hndf.loc[hndf['hor']=='organic']
orgndf=ndf.loc[ndf['hor']=='organic']


'''for i in (cdfna['netnh4']):
    print(i)
    print(type(i))
'''

try:
    orgNminc=open("C:\\Users\\rmoge\\GitHub\\QB2017_DivPro\\Data\\orgNminc.csv",'w')
    orgNminh=open("C:\\Users\\rmoge\\GitHub\\QB2017_DivPro\\Data\\orgNminh.csv",'w')
    orgNminhn=open("C:\\Users\\rmoge\\GitHub\\QB2017_DivPro\\Data\\orgNminhn.csv",'w')
    orgNminn=open("C:\\Users\\rmoge\\GitHub\\QB2017_DivPro\\Data\\orgNminn.csv",'w')
    minNminc=open("C:\\Users\\rmoge\\GitHub\\QB2017_DivPro\\Data\\minNminc.csv",'w')
    minNminh=open("C:\\Users\\rmoge\\GitHub\\QB2017_DivPro\\Data\\minNminh.csv",'w')
    minNminhn=open("C:\\Users\\rmoge\\GitHub\\QB2017_DivPro\\Data\\minNminhn.csv",'w')
    minNminn=open("C:\\Users\\rmoge\\GitHub\\QB2017_DivPro\\Data\\minNminn.csv",'w')
    
    orgc=makets(orgcdf,orgNminc); orgh=makets(orghdf,orgNminh); orghn=makets(orghndf,orgNminhn); orgn=makets(orgndf,orgNminn)
    minc=makets(mincdf,minNminc); minh=makets(minhdf,minNminh); minhn=makets(minhndf,minNminhn); minn=makets(minndf,minNminn)
    
finally:
    orgNminc.close(); orgNminh.close(), orgNminhn.close(), orgNminn.close()
    minNminc.close(); minNminh.close(), minNminhn.close(), minNminn.close()

'''grouped = ndf.groupby(['year', 'month'])
meanspermonth = grouped['co2flux'].mean()
tester = list(range(0, len(grouped), 1))#make a list from 0 to number.of.groups


ndfg = ndf.groupby(['year', 'month']).groups
type(ndf.groupby(['year', 'month']).groups)
sl = (ndf['month'])
ndff = ndf
ndff['e'] = pd.Series([chooseseason(x) for x in ndff['month']], index=ndff.index)

print(chooseseason(6))'''




