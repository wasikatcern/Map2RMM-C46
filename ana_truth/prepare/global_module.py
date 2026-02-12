import random
import json,sys
sys.path.append("modules/")

import ROOT;
from ROOT import TTree;
from array import *
from math import *
import math,sys,os
from array import array
from decimal import Decimal
import numpy
import random
import sys,zipfile,json,math


### This is common for all
with open('../config.json') as json_file:
    data = json.load(json_file)
    maxNumber=int(data['maxNumber'])
    maxTypes=int(data['maxTypes'])
    mSize=int(data['mSize'])
print ("maxNumber=",maxNumber," maxTypes=",maxTypes," mSize=",mSize)
mSize=maxTypes*maxNumber+1;

# dijet invariant mass
x=1+0*maxNumber+1  # X position  
y=1+0*maxNumber    # Y position 
mjj=(x,y) #  index of Mjj  matrix ellement 

# PT of first jet
x=1+0*maxNumber  # X position  
y=1+0*maxNumber  # Y position 
pT=(x,y) #  index of Mjj  matrix ellement 


#  bb mass 
x=1+1*maxNumber+1
y=1+1*maxNumber
mbb=(x,y)

#  bj mass 
x=1+1*maxNumber+1
y=1+0*maxNumber
mbj=(x,y) 

# e+e- 
x=1+2*maxNumber+1
y=1+2*maxNumber
mee=(x,y)

# mu+mu 
x=1+3*maxNumber+1
y=1+3*maxNumber
mmumu=(x,y)


### This list contains excluded values for Z-score calculation
### We excluding pT of leading jet, Mjj and mbb
# excluded_val= ( pT, mjj, mbb)
excluded_val= (mjj, mbb, mbj)
print ("Excluded cells=",excluded_val) 

#### Exclusion values for RMM matrix #############
###################################


# dijet invariant mass
x=2 # X position  
y=1 # Y position 
inx1=x*mSize+y; #  index of hxw matrix ellement 

# pT1 
x=1 # X position  
y=1 # Y position 
inx2=x*mSize+y; #  index of hxw matrix ellement 

# Mjj for for light-jet + b-jets
x=1+maxNumber # X position  
y=1 # Y position 
inx3=x*mSize+y; #  index of hxw matrix ellement 

# pT1 for for b-jets
x=1+maxNumber # X position  
y=1+maxNumber # Y position 
inx4=x*mSize+y; #  index # pT for for b-jets

# Mjj for for 2-b jets 
x=2+maxNumber # X position  
y=1+maxNumber # Y position 
inx5=x*mSize+y; #  index of hxw matrix ellement 


# exlusion matrix for RMM in terms of indexes (how it is packed) 
excluded=(inx1,inx2,inx3,inx4,inx5)




