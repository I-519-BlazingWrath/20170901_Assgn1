# -*- coding: utf-8 -*-
"""
Created on Fri Aug 25 08:19:20 2017

@author: rmoge
"""
from __future__ import division
import re, sys, math, operator,random, copy,collections; import numpy as np 
from itertools import groupby; import pprint as pp
import matplotlib.pyplot as plt
from Bio import SeqIO; from Bio import Seq;#from Bio.Alphabet import generic_dna, generic_protein
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC

def chooseRandomComposition():
    empty={"A":0, "C":0, "G":0,"T":0};compD=copy.deepcopy(empty)
    for b in empty.keys():
        empty[b]=random.random()
    myL=[empty[i] for i in ("A","C","G","T")]
    correct=(1/sum(myL))
    for b in compD.keys():
        compD[b]=(empty[b]*correct)
    return compD
composition=chooseRandomComposition()
    

def makeRandomSeq(composition={},length=0):
    
    # composition = {"A":.15, "C":.25, "G":.1,"T":.5}
    seq=''
    ab=composition.keys()

    bound=0
    bounds=[]
    for char in ab:
        bound+=composition[char]
        bounds.append(bound)
    print(bounds)

    assert 0.99<bounds[-1]<1.01#sanity check. add this to your 567 code...This saves you from spending months try to interpp
    #a garbage output    If this fails, it will spit out an assertion error and refuse to run---thereby avoiding spewing bogus
    #input

    for i in range(length):
        roll=random.random()
        for char,bound in zip(ab,bounds):#<--LEARN HOW TO USE ZIP, DUDE
            if roll<=bound:
                seq+=char
                break#once you find that it is a C, don't ALSO add a G and a T
    #print(seq)
    return seq
#composition = {"A":.15, "C":.25, "G":.1,"T":.5}
#seq=Seq(makeRandomSeq(composition,1590390229),Bio.Alphabet.generic_dna)
seq=Seq(makeRandomSeq(composition,1000),Bio.Alphabet.generic_dna)
record=SeqIO.SeqRecord(seq, 'random_sequence')
SeqIO.write(record,'random.fasta','fasta')

def getRandomReads(genomeIN='',readnumber=0,readlen=0):
    readL=[];genome=str(genomeIN)
    for i in range(readnumber):
        start=random.randint(0,(len(genome)-readlen+1))
        tempread=genome[start:(start+readlen+1)]
        print(tempread)
        readL+=[tempread]
    pp.pprint(readL)
    my_seqs=[]
    for index,s in enumerate(readL):
        my_seqs.append(SeqRecord(Seq(s,Bio.Alphabet.generic_dna),id=('randomread'+str(index))))
    return my_seqs

#ourreads=getRandomReads(seq,220887531,216)
ourreads=getRandomReads(seq,100,150)
pp.pprint(ourreads)

record=SeqIO.SeqRecord(seq, 'random_sequence')
SeqIO.write(record,'random.fasta','fasta')
for index,r in enumerate(ourreads):
    pass
with open('randomreads.fasta','w') as OPfile:
    SeqIO.write(ourreads,OPfile,'fasta')
    OPfile.write("\n")
    '''for index,r in enumerate(ourreads):
        SeqIO.write(SeqIO.SeqRecord(r,"random_read"+str(index)),OPfile,'fasta')
'''