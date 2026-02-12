import random
import sys
sys.path.append("modules/")

# import atlas styles
from array import *
from math import *
from ROOT import TH1D,TF1,TCanvas,TColor,TPostScript,TProfile2D,THStack,TRandom3,TFile,TLatex,TLegend,TPaveText,TGraphErrors,kRed,kBlue,kGreen,kCyan,kAzure,kYellow,kTRUE
import math,sys,os
import shapiro
import ROOT
from array import array
from decimal import Decimal
import sys,zipfile,json,math
from ROOT import gROOT, gPad, gStyle, gRandom
import math,sys,os 
import numpy
import random
import sys,zipfile,json,math
from array import *
from math import *
import json


# CM energy of collisions
CMS=13600.0

## masses labels 
masseslab=["jj","jb","bb","je","j\;\mu","j\;\gamma","b\;\ell","b\;\mu","b\;\gamma"]
## masses data in histograms
massesdata=["jj","jb","bb","je","jm","jg","be","bm","bg"]


# triggers
labels=["MET", "\ell", "2 \ell ", "1 #gamma ", "2 #gamma ", "1j ", "4j "]
labels_map={1:"MET",2:"l",3:"2 l",4:"1 #gamma ",5:"2 #gamma ",6:"1j ",7:"4j "}
labels_color=[1, 2, 4, 6, 8, 44, 1]
labels_style=[24, 20, 21, 22, 23, 20, 34]

xmapcolor={}
xmapcolor[500]=31
xmapcolor[700]=32
xmapcolor[1000]=33
xmapcolor[1500]=34
xmapcolor[2000]=36


# MG5 cross sections in pb 
# X->HH from MG5 in fb (not pb)  
# note that these cross section do not take into account 
# branching ratio of one H -> ZZ and second to H->bb.
xfactor=1000*(1.0/3.687)   # fb x branching ratio  
mg5xcross={}
mg5xcross[500]=12.4*xfactor 
mg5xcross[700]=2.0*xfactor 
mg5xcross[1000]=0.2*xfactor 
mg5xcross[1500]=0.01*xfactor 
mg5xcross[2000]=0.001*xfactor 


# factors for "SH" to get the same significance
xfactor=1000*(1.0/3.687)   # fb x branching ratio
mg5Sxcross={}
#mg5Sxcross[500]=1.4*xfactor
#mg5Sxcross[700]=0.187*xfactor
#mg5Sxcross[1000]=0.015*xfactor
#mg5Sxcross[1500]=0.00065*xfactor
#mg5Sxcross[2000]=0.000065*xfactor

mg5Sxcross[500]=12.4*xfactor
mg5Sxcross[700]=2.0*xfactor
mg5Sxcross[1000]=0.2*xfactor
mg5Sxcross[1500]=0.01*xfactor
mg5Sxcross[2000]=0.001*xfactor



#########################################################
# cut to select outlier events

# expcted lumin in pb
ExpectedLumiFB=140


# MC region for 10 pb working point ("data limit") 
# It is taken from PRL paper 2024
CutOutlier_10PB=-9.10 

# WP region for mean-3*RMS pb working point 
CutOutlier_3RMS=-7.85

# main cuts
# CutOutlier=CutOutlier_10PB
CutOutlier=CutOutlier_3RMS



mcTT="t#bar{t}+single t"
mcWZ="W/Z+jet"


mjjBinsL = [0, 14,28,43,58,72,86,99,112,125,138,151,164,177,190, 203, 216, 229, 243, 257, 272, 287, 303, 319, 335, 352, 369, 387, 405, 424, 443, 462, 482, 502, 523, 544, 566, 588, 611, 634, 657, 681, 705, 730, 755, 781, 807, 834, 861, 889, 917, 946, 976, 1006, 1037, 1068, 1100, 1133, 1166, 1200, 1234, 1269, 1305, 1341, 1378, 1416, 1454, 1493, 1533, 1573, 1614, 1656, 1698, 1741, 1785, 1830, 1875, 1921, 1968, 2016, 2065, 2114, 2164, 2215, 2267, 2320, 2374, 2429, 2485, 2542, 2600, 2659, 2719, 2780, 2842, 2905, 2969, 3034, 3100, 3167, 3235, 3305, 3376, 3448, 3521, 3596, 3672, 3749, 3827, 3907, 3988, 4070, 4154, 4239, 4326, 4414, 4504, 4595, 4688, 4782, 4878, 4975, 5074, 5175, 5277, 5381, 5487, 5595, 5705, 5817, 5931, 6047, 6165, 6285, 6407, 6531, 6658, 6787, 6918, 7052, 7188, 7326, 7467, 7610, 7756, 7904, 8055, 8208, 8364, 8523, 8685, 8850, 9019, 9191, 9366, 9544, 9726, 9911, 10100, 10292, 10488, 10688, 10892, 11100, 11312, 11528, 11748, 11972, 12200, 12432, 12669, 12910, 13156];

mjjBins = array("d", mjjBinsL)


from os.path import exists

confile='data/config.json'
file_exists = exists( confile )

if (file_exists):
 with open( confile ) as json_file:
    data = json.load(json_file)
    maxNumber=int(data['maxNumber'])
    maxTypes=int(data['maxTypes'])
    mSize=int(data['mSize'])
    print("Read from file ",confile)
else:
    maxNumber=10
    maxTypes=5
    mSize=5


print ("maxNumber=",maxNumber," maxTypes=",maxTypes," mSize=",mSize) 
mSize=maxTypes*maxNumber+1;

######################### define position of invarinat masses from RMM #################

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
x=1+1*maxNumber
y=1+0*maxNumber
mbj=(x,y) 

# mu+mu 
x=1+2*maxNumber+1
y=1+2*maxNumber
mmumu=(x,y)

# e+e 
x=1+3*maxNumber+1
y=1+3*maxNumber
mee=(x,y)

# j+mu 
x=1+2*maxNumber
y=1+0*maxNumber
mjmu=(x,y)

# j+e 
x=1+3*maxNumber
y=1+0*maxNumber
mje=(x,y)

# j+gamma 
x=1+4*maxNumber
y=1+0*maxNumber
mjg=(x,y)

# b+mu 
x=1+2*maxNumber
y=1+1*maxNumber
mbmu=(x,y)

# b+e  
x=1+3*maxNumber
y=1+1*maxNumber
mbe=(x,y)

# b+gamma 
x=1+4*maxNumber
y=1+1*maxNumber
mbg=(x,y)

############# end invariant mass definitions using RMM ############

### This list contains excluded values for Z-score calculation
### We excluding pT of leading jet, Mjj and mbb
# excluded_val= ( pT, mjj, mbb)
# excluded_val= (mjj, mbb)
# print ("Excluded cells=",excluded_val ) 

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


# get a labale for trigger 
def getTriggerLabel(trig_type):
  Tlab="1 lep"
  # the algoritm  the algorithm: HighsPt cut x 3 
  if (int(trig_type)==1):
             Tlab="T1:\; MET"
  if (int(trig_type)==2): # 1 lepton 
             Tlab="T2:\; 1 \ell"
  if (int(trig_type)==3): # 2 lepton pT>25 GeV  
             Tlab="T3:\; 2 \ell"
  if (int(trig_type)==4): #  single photon 
             Tlab="T4:\; 1 \gamma"
  if (int(trig_type)==5): #  2 photon    
             Tlab="T5:\; 2 \gamma"
  if (int(trig_type)==6): #  single jet (pT>500 GeV)  
             Tlab="T6:\; 1 jet"
  if (int(trig_type)==7): #  4 jet  (lead 200 GeV) 
             Tlab="T7:\; 4 jets"
  return Tlab



# Save mathplot in CSV file
import csv
def SavePlotXY(xfile,lines, Xlab="X", Ylab="Y"):
   #NrLines=len(lines)
   #print("Nr of lines",NrLines)
   print("Save plot in CSV ",xfile);
   with open(xfile, 'w') as myfile:
            data=lines[0].get_data()
            writer = csv.writer(myfile)
            writer.writerow([Xlab, Ylab])
            for i in range(len(data[0])):
                writer.writerow([data[0][i], data[1][i]])

from numpy import savetxt
def SaveNumpyData(xfile,lines):
    savetxt(xfile,lines,delimiter=',')

# do counting statistics for >3 events
# for less, increase the error 
# Set to 0 if Nr of entries less than 1 
def countingErrors(hhh):
  for i in range(1, hhh.GetNbinsX()):
     D = hhh.GetBinContent(i);
     if (D>3.0): 
              hhh.SetBinError(i, sqrt(D) )
     if (D>0.999 and D<=3): 
              hhh.SetBinError(i, 2*sqrt(D) )
     if (D<1.0):
              hhh.SetBinError(i,0)
              hhh.SetBinContent(i,0)


# calculate loss cut for 10% and 1 % of data
def findCutvalues(hhh):
  for i in range(40, 150, 1):
     cutval = i/10
     intigral = hhh.Integral(hhh.FindBin(cutval), hhh.FindBin(0))
     ratio = intigral/hhh.Integral()
     if (ratio == 0.1 or (math.isclose(ratio, 0.1, abs_tol = 0.01)==1)):
        Xcut1=cutval
     elif (ratio == 0.01 or (math.isclose(ratio, 0.01, abs_tol = 0.01)==1)):
        Xcut2=cutval
  return [Xcut1,Xcut2]
    

# do counting statistics only if errors are smaller than counting 
def countingErrorsCorrect(hhh):
  for i in range(1, hhh.GetNbinsX()):
     D = hhh.GetBinContent(i);
     E=  hhh.GetBinError(i);
     Expected= 0;
     if (D>1):  Expected= sqrt(D)
     if (E< Expected and D>1):
              hhh.SetBinError(i, sqrt(D) )


def SavePlotHisto(xfile,ax):
   print("Save histogram in CSV ",xfile);
   p = ax.patches
   with open(xfile, 'w') as myfile:
            writer = csv.writer(myfile)
            writer.writerow(["Xlow", "Height"])
            for i in range(len(p) ):
                lower_left_corner=p[i].get_xy() 
                #writer.writerow([ lower_left_corner[0], p[i].get_width(), p[i].get_height()  ])
                #writer.writerow([ lower_left_corner[0], p[i].get_height()  ])
                writer.writerow([ lower_left_corner[0], lower_left_corner[1]  ])
## draw axis
def drawXAxis(sf,gPad,XMIN,YMIN,XMAX,YMAX,nameX,nameY,showXAxis=True, showYAxis=True):
 h=gPad.DrawFrame(XMIN,YMIN,XMAX,YMAX);
 ay=h.GetYaxis();
 ay.SetLabelFont(42)

 if (sf==1):
             ay.SetLabelSize(0.05)
             ay.SetTitleSize(0.06)

 if (sf==2 or sf==3):
             ay.SetLabelSize(0.10)
             ay.SetTitleSize(0.3)
 if (sf==20):
             ay.SetLabelSize(0.18)
 if (sf==30):
             ay.SetLabelSize(0.12)
# ay.SetTitleSize(0.1)
 ay.SetNdivisions(505);
 if (sf==1): ay.SetTitle( nameY )
 # ay.Draw("same")
 ax=h.GetXaxis();
 if (sf==1 or sf==2): ax.SetTitle( nameX );
 if (sf==30): ax.SetTitle( nameX );
 ax.SetTitleOffset(1.18)
 ay.SetTitleOffset(0.8)

 ax.SetLabelFont(42)
 # ax.SetTitleFont(42)
 ay.SetLabelFont(42)
 # ay.SetTitleFont(42)
 ax.SetLabelSize(0.12)
 ax.SetTitleSize(0.14)

 if (showXAxis==False):
         ax.SetLabelSize(0)
         ax.SetTitleSize(0)
 if (showYAxis):
          ay.SetLabelSize(0)
          ay.SetTitleSize(0)

 #ay.SetTitleSize(0.14)
 if (sf==30):
          ax.SetLabelSize(0.12)
          ax.SetTitleSize(0.12)
 if (sf==2 or sf==3):
             ay.SetLabelSize(0.12)
             ay.SetTitleSize(0.2)

 ax.Draw("same");
 ay.Draw("same");
 return h

def style3par(back):
     back.SetNpx(100); back.SetLineColor(4); back.SetLineStyle(1)
     back.SetLineWidth(2)
     back.SetParameter(0,4.61489e-02)
     back.SetParameter(1,1.23190e+01)
     back.SetParameter(2,3.65204e+00)

     #back.SetParameter(3,-6.81801e-01)
     #back.SetParLimits(0,0,100)
     # back.SetParLimits(1,0,12)
     #back.SetParLimits(2,-100,100)
     return back

def style5par(back):
     back.SetNpx(200); back.SetLineColor(4); back.SetLineStyle(1)
     back.SetLineWidth(2)
     back.SetParameter(0,6.0e+10)
     back.SetParameter(1,80)
     back.SetParameter(2,40)
     back.SetParameter(3,11)
     back.SetParameter(4,1.0)
     #back.SetParLimits(0,0,10000)
     #back.SetParLimits(1,0,100000000)
     #back.SetParLimits(2,-10000,10000)
     #back.SetParLimits(3,-400,400)
     return back

# get width of the bin near the mass
def getBinWidth(bins,peak):
    imean = bins.FindBin(peak)
    return bins.GetBinCenter(imean+1) - bins.GetBinCenter(imean);


# get run from MC sample
def getRun(sample):
    parts=sample.split(".")
    run=parts[1]
    return run 

# residual plots: input histoogram, function, file name 
# http://sdittami.altervista.org/shapirotest/ShapiroTest.html
from module_functions import  Gauss
def showResiduals(hh,func,fname, MyMinX=-12,  MyMaxX=12, isKS=True):
   print ("showResiduals: Calculate residuals.")
   MyBins=100
   res=TH1D("Residuals","Residuals",MyBins,MyMinX,MyMaxX);
   res.SetTitle("")
   res.SetStats(1)
   res.SetLineWidth(2)
   res.SetMarkerColor( 1 )
   res.SetMarkerStyle( 20 )
   res.SetMarkerSize( 0.8 )
   res.SetFillColor(42)
   nameX="D_{i} - F_{i} / #Delta D_{i}"
   nameY="Entries"
   FitMin=func.GetXmin()
   FitMax=func.GetXmax()
   print ("Fit min=",FitMin,"  max=",FitMax)
   nres=0.0
   residuals=[]
   for i in range(1,hh.GetNbinsX()):
     center=hh.GetBinCenter(i)
     if (hh.GetBinContent(i)>0 and center>FitMin and center<FitMax):
       center=hh.GetBinCenter(i)
       x=hh.GetBinCenter(i)
       D = hh.GetBinContent(i);
       Derr = hh.GetBinError(i);
       B = func.Eval(center);
       frac=0
       if Derr>0:
          frac = (D-B)/Derr
       residuals.append(frac)
       res.Fill(frac)
       nres=nres+1.0
   res.SetStats(1)
   bcallable=Gauss()
   back=TF1("back",bcallable,MyMinX,MyMaxX,3);
   back.SetNpx(200); back.SetLineColor(4); back.SetLineStyle(1)
   back.SetLineWidth(2)
   back.SetParameter(0,10)
   back.SetParameter(1,0)
   back.SetParameter(2,1.0)
   back.SetParLimits(2,0.1,1000)
   #back.SetParLimits(0,0.01,10000000)
   #back.SetParLimits(1,-5.0,5.0)
   #back.SetParLimits(2,0.0,5.0)

   #back.FixParameter(1,0)
   #back.FixParameter(2,1.0)
   nn=0
   chi2min=10000
   parbest=[]
   for i in range(10):
     fitr=res.Fit(back,"SMR0")
     print ("Status=",int(fitr), " is valid=",fitr.IsValid())
     if (fitr.IsValid()==True):
             chi2=back.GetChisquare()/back.GetNDF()
             if chi2<chi2min:
                    nn=nn+1
                    if nn>3:
                           break;
                    back.SetParameter(0,random.randint(0,10))
                    back.SetParameter(1,random.randint(-1,1))
                    back.SetParameter(2,random.randint(0,2.0))
                    par = back.GetParameters()

   #fitr=res.Fit(back,"SMR0")
   fitr.Print()
   print ("Is valid=",fitr.IsValid())

   par = back.GetParameters()
   err=back.GetParErrors()
   chi2= back.GetChisquare()
   ndf=back.GetNDF()
   print ("Chi2=", chi2," ndf=",ndf, " chi2/ndf=",chi2/ndf)
   prob=fitr.Prob();
   print ("Chi2 Probability=",fitr.Prob());
   # make reference for normal
   norm_mean=0
   norm_width=1
   normal=TH1D("Normal with sig=1","Reference normal",MyBins,MyMinX,MyMaxX);
   normal.SetLineWidth(3)
   normal.SetLineColor(2)
   # normal.SetFillColor( 5 )

   maxe=5000
   r = TRandom3()
   for i in range (maxe) :
          xA = r.Gaus(norm_mean, norm_width)
          normal.Fill(xA)
   norM=nres/maxe
   normal.Scale(norM)
   pKSbinned = res.KolmogorovTest(normal)
   KSprob="KS prob ="+"{0:.2f}".format(pKSbinned)
   if (isKS): print (KSprob)

   shapiro_prob=shapiro.ShapiroWilkW(residuals)
   Shapiro="ShapiroWilk ={0:.2f}".format(shapiro_prob)
   print (Shapiro)

   gROOT.SetStyle("ATLAS");
   gStyle.SetOptStat(220002210);
   gStyle.SetStatW(0.32)
   c2=TCanvas("c","BPRE",10,10,600,540);
   c2.Divide(1,1,0.008,0.007);
   c2.SetBottomMargin(0.1)
   c2.SetTopMargin(0.05)
   c2.SetRightMargin(0.02)
   c2.SetLeftMargin(0.10)

   binmax = normal.GetMaximumBin();
   Ymax=normal.GetBinContent(normal.FindBin(0));
   for i in range(1,res.GetNbinsX()):
      if res.GetBinContent(i)>Ymax: Ymax=res.GetBinContent(i);
   Ymax=1.15*Ymax;

   #h=gPad.DrawFrame(MyMinX,0,MyMaxX,Ymax)

   ps2 = TPostScript( fname,113)

   res.SetStats(1)
   gStyle.SetOptStat(220002200);
   gStyle.SetStatW(0.32)

   res.SetAxisRange(0, Ymax,"y");
   res.Draw("histo")
   back.Draw("same")
   if (isKS): normal.Draw("histo same")
   leg2=TLegend(0.11, 0.6, 0.39, 0.90);
   leg2.SetBorderSize(0);
   leg2.SetTextFont(62);
   leg2.SetFillColor(10);
   leg2.SetTextSize(0.04);
   leg2.AddEntry(res,"Residuals","f")
   leg2.AddEntry(back,"Gauss fit","l")
   mean= "mean="+"{0:.2f}".format(par[1])
   mean_err= "#pm "+"{0:.2f}".format(err[1])
   sig= "#sigma="+"{0:.2f}".format(par[2])
   sig_err= "#pm "+"{0:.2f}".format(err[2])
   leg2.AddEntry(back,mean+mean_err,"")
   leg2.AddEntry(back,sig+sig_err,"")
   leg2.AddEntry(back,"#chi^{2}/ndf="+"{0:.2f}".format(chi2/ndf)+"(p="+"{0:.2f})".format(prob),"")
   leg2.AddEntry(back,Shapiro,"")
   if (isKS): leg2.AddEntry(normal,"Normal (#sigma=1)","l")
   if (isKS): leg2.AddEntry(back,KSprob,"")

   leg2.Draw("same");
   ax=res.GetXaxis();
   ax.SetTitle( nameX );
   ay=res.GetYaxis();
   ay.SetTitle( nameY );
   ax.SetTitleOffset(1.0); ay.SetTitleOffset(1.0)
   ax.Draw("same")
   ay.Draw("same")
   gPad.RedrawAxis()
   c2.Update()
   ps2.Close()
   c2.Close();
   print (fname, " done")

# get a histogram from input file (simple), but use lumi weigths if needed 
def getHistogram(inputroot,histo_name,LumiWeight=1):
      global lumi
      print(inputroot,histo_name)
      rf=TFile(inputroot)
      #rf.ls()
      tmp=rf.Get(histo_name)
      tmp=tmp.Clone()
      tmp.SetDirectory(0)
      tmp.Scale(LumiWeight)
      rf.Close()
      #tmp.Scale(LumiWeight)
      #countingErrors(tmp)
      return tmp


def TH1Error2Zero(h1):
    for i in range(h1.GetNbinsX()):
        y = h1.GetBinContent(i+1)
        ey = h1.GetBinError(i+1)
        x = h1.GetBinCenter(i+1)
        ex = h1.GetBinWidth(i+1)
        h1.SetBinContent(i+1, y)
        h1.SetBinError(i+1, ey)
        if (y > 0):
            h1.SetBinError(i+1, 0)
        else:
            h1.SetBinError(i+1, 0)

# macro to divide by bin width taking into account Nr of cores
def getBinSize(fdata):
   xbins=fdata.Get("bins_m")
   nnn=fdata.Get("cpucores"); 
   nCPU=nnn.GetBinContent(2)
   print ("Nr of cores=",nCPU)
   xbins.Scale(1.0/nCPU)
   bins=xbins.Clone();
   bins.SetDirectory(0)
   TH1Error2Zero(bins)
   #bins.Print("all")
   return bins


def GetZVal (p, excess) :
  #the function normal_quantile converts a p-value into a significance,
  #i.e. the number of standard deviations corresponding to the right-tail of 
  #a Gaussian
  if excess :
    zval = ROOT.Math.normal_quantile(1-p,1);
  else :
    zval = ROOT.Math.normal_quantile(p,1);

  return zval


# calculate significance
def signif(S,B):
     # print("Calculate significance for B=",B," and S=",S)
     if (S<10e-6): return S/math.sqrt(B) # avoid domain problem 
     else:  return math.sqrt(2* ( (S+B)*math.log(1+ float(S)/B)-S ) )

# get significance near the mass "peak" (gev)  using bkg and sig histograms, for a list of masses
def getSignificances(bkg, sig, peak=120, jsout="out/significance.js"):
     Idx=bkg.FindBin(peak)
     print("3 bins with maximum entries=",Idx)
     signif_dic={}
     for i in range(len(sig)):
         his=sig[i]
         N_signal = his.GetBinContent( Idx )+his.GetBinContent( Idx-1 )+his.GetBinContent( Idx+1 );
         N_bkg=bkg.GetBinContent( Idx )+bkg.GetBinContent( Idx-1 )+bkg.GetBinContent( Idx+1 );
         Sign=signif(N_signal, N_bkg)
         print("Mass=",his.GetTitle(), " Sign=",Sign)
         signif_dic[his.GetTitle()]=float(Sign)
     with open(jsout, "w") as f:
      f.write(json.dumps(signif_dic, indent=4)) 
     print("Write significances to=",jsout)

