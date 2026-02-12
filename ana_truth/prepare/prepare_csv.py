# S.Chekanov
# Convert data or MC to zip

import os,sys, math, csv
from global_module import *

print ('Number of arguments:', len(sys.argv), 'arguments.')
print ('Argument List:', str(sys.argv))
n = len(sys.argv)
if (n != 4):
      print ("No argumentts!")
      sys.exit()

inputData=sys.argv[1]
outputData=sys.argv[2]
label=str(sys.argv[3]) 


# if data and only 10%
onlyFraction=1.0
if (outputData.find("10percent")>-1):
                     onlyFraction=0.1
if (outputData.find("1percent")>-1):
                     onlyFraction=0.01
if (outputData.find("01percent")>-1):
                     onlyFraction=0.001
if (outputData.find("001percent")>-1):
                     onlyFraction=0.0001

# shift by 1 event when we wee _1
shift=0
if (outputData.find("_1")>-1):
                     shift=1
# shift by 1 event when we wee _2
if (outputData.find("_2")>-1):
                     shift=2
# shift by 1 event when we wee _3
if (outputData.find("_3")>-1):
                     shift=3


print("Process fraction =", onlyFraction*100,"%")

proc=[ inputData ]
outRMM=[ outputData ]

# max number of events 2M for training and testing 
MaxTrain=-1
print("Max number of events=",MaxTrain)

figdir="figs/"
epsfig=figdir+__file__.replace(".py",".eps")

# output
csv_file=open(str(outRMM[0]), "w")
#writer = csv.writer(csv_file, delimiter=',')

columnsX="Run,Event,Weight"
for i in range(1,mSize*mSize+1):
       columnsX=columnsX+",V_"+str(i) 
csv_file.write(columnsX+",Label\n")

rfile=[]
for i in proc:
     rfile.append(ROOT.TFile.Open(i))
     print(i) 

test=[]
pat=[]
validation=[]
ev=0

# Function to convert  
def listToString(s): 
    # initialize an empty string
    nn=0
    str1 = "" 
    # traverse in the string  
    for ele in s: 
        ST='%.7E' % Decimal(ele)
        if (ele == 0): ST="0" # compact 0 representation 
        if nn==0: str1 = ST
        else: str1 = str1+","+ST
        nn=nn+1
    # return string  
    return str1 


# we also implement a shift by +1 to have statistically differents sample
me=shift 
if (onlyFraction<1.0):
   print("This is another 1% with shift =",me)


print("Start loop")
for i in range(len(proc)):
   ev=0
   for event in rfile[i].inputNN:
      Trun=1; #  event.run
      me=me+1;
      Tevent=me; # event.event
      if (onlyFraction<1.0):
          takeevent=False;
          #if (me%int(1/onlyFraction) == 0): takeevent=True 
          if (Tevent%int(1/onlyFraction) == 0): takeevent=True 
          if (takeevent == False): continue 

          #ran=random.uniform(0, 1)
          #if (ran>onlyFraction): continue


      NN=(event.proj).size()
      a=event.proj
      inx1=event.proj_index1
      inx2=event.proj_index2
      Tweight=event.weight # for MC with weigths
      emptyMatrix = numpy.zeros(shape=(mSize,mSize))

      ST='%.6E' % Decimal(Tweight)
      if (Tweight==1):  ST="1"
      pos=str(Trun)+","+str(Tevent)+","+str(ST)

      for i3 in range(NN):
                w=inx1[i3];
                h=inx2[i3];
                val=float(a[i3])
                # do some exclusion if needed
                #if (h,w) in excluded_val: 
                #                     # print "Position",[h,w]," is excluded"
                #                     continue; 
                emptyMatrix[w][h] = val

      data=(emptyMatrix.flatten()).tolist()

      sdata=listToString(data)
      csv_file.write(pos+","+sdata+","+label+"\n") # last number is the label (MC) 

      ev=ev+1;
      if ev<10000 and ev%100==0:  print ("Event=", ev, " records=",NN," total processed=",me)
      if (ev%1000==0): print ("Event=", ev, " records=",NN," total processed=",me)
      if (MaxTrain>0):
                  if (ev>MaxTrain): break;

print ("Close csv_file") 
csv_file.close()
