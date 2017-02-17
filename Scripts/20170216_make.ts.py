from __future__ import division
import sys; import math; import re; import subprocess; import os
from random import choice, randrange
import matplotlib.pyplot as plt; import numpy as np; import scipy as sc; from scipy import stats; import scipy.stats.distributions

def go(IP):

    geneslist=[]
    headerslist=()
    for index,line in enumerate(IP):
        if index==0:#skip the first line of headers
            headerslist+=tuple(line.strip("\n").split("\t"))
            
        
        elif index!=0:
             geneslist+=[line.strip("\n").split(",")]
    return geneslist,headerslist


def chooseseason(month):
    row=''
    if 4<=month<=9:
        row='q1'
    elif 1<=month<=3 or 10<=month<=12:
        row='q2'
    else:
        row="??"
    return row
	
def makets(thedf):
	thedf['season'] = pd.Series([chooseseason(x) for x in thedf['month']], index=thedf.index)
	groups=thedf.groupby(['year', 'season'])
	temp=[range(0,len(groups),1)]
	qdf=pd.DataFrame(columns=['year','season'],index=temp)
	for index,row in enumerate(qdf):
		row=[groups[index][0][0], groups[index][0][1]]
	