# Plot loss from anomaly detection

import glob,sys
sys.path.append("modules/")
from AtlasStyle import *
from AtlasUtils import *
from global_module import *

from ROOT import TH1D,TF1,TProfile2D,TEllipse, THStack,TRandom3,TFile,TLatex,TLegend,TPaveText,TGraphErrors,kRed,kBlue,kGreen,kCyan,kAzure,kYellow,kTRUE
import ROOT


print ('Number of arguments:', len(sys.argv), 'arguments.') 
print ('Argument List:', str(sys.argv))
print ('Use as: script.py -b 0 (or 1,2)') 
myinput="interactive"
if (len(sys.argv) ==2):
   myinput = sys.argv[1]
print ("Mode=",myinput) 


gROOT.Reset()
figdir="figs/"
fname=os.path.basename(__file__)
#fname=fname.replace(".py","_"+trig_type+".py");
epsfig=figdir+(fname).replace(".py",".eps")

nameX="log (loss)"
nameY="Events"
Ymin=0.01 
Ymax=10000000000
Xmin=-12
Xmax=-4

######################################################
gROOT.SetStyle("Plain");
c1=TCanvas("c_massjj","BPRE",10,10,600,500);
c1.Divide(1,1,0.008,0.007);
ps1 = TPostScript( epsfig,113)

c1.cd(1);
gPad.SetLogy(1)
gPad.SetLogx(0)
gPad.SetTopMargin(0.05)
gPad.SetBottomMargin(0.12)
gPad.SetLeftMargin(0.14)
gPad.SetRightMargin(0.04)


h=gPad.DrawFrame(Xmin,Ymin,Xmax,Ymax);
h.Draw()


name="Loss"
xfile="out/tev13.6pp_pythia8_ttbar_2lep_ADFilter.root"
xfile=TFile( xfile )
ttbar=xfile.Get(name)
ttbar.SetTitle("")
ttbar.SetStats(0)
ttbar.SetLineWidth(2)
ttbar.SetLineColor( 1 )
ttbar.SetMarkerColor( 1 )
ttbar.SetMarkerSize( 1.1 )
ttbar.SetLineStyle(2)
#ttbar.SetFillColor( 6 )
ttbar.SetAxisRange(Ymin, Ymax,"y");
ttbar.SetAxisRange(Xmin, Xmax,"x");


# expected events
cross=xfile.Get("cross");
xsec=cross.GetBinContent(1)
lumi=float(cross.GetBinContent(5))
print("Cross=",xsec," lumi=",lumi)
CurrentLumuFB=lumi/1000.0
Scale=ExpectedLumiFB/CurrentLumuFB;
ttbar.Scale(Scale)
lumi=ExpectedLumiFB*1000;



ttbar.Draw("same histo")

dMean=ttbar.GetMean()
dData=dMean+3*ttbar.GetRMS()
print("Anomaly region defined at mean-3*RMS as =",dData)


# show hh BSM
bmass=[500,700,1000,1500,2000]
bhitos={} 
for b in range(len(bmass)):
     mass=bmass[b]
     xfile1="out/pythia8_X"+str(mass)+"GeV_SH2bbll_ADFilter.root"
     xfile1=TFile( xfile1 )
     bsm=xfile1.Get(name)
     bsm.SetDirectory(0)
     bsm.SetTitle("")
     bsm.SetStats(0)
     bsm.SetLineWidth(2)
     bsm.SetLineColor( 1 )
     bsm.SetMarkerColor( 1 )
     bsm.SetMarkerSize( 1.1 )

     crossHisto=xfile1.Get("cross");
     crossZ=mg5Sxcross[mass]
     nevents=float(crossHisto.GetBinContent(2))
     CurrentLumuZ=nevents/crossZ # lumi in fb-1
     print("Cross=",crossZ,"fb lumi=",CurrentLumuZ," name=",name, "gen events=",nevents)
     Scale=ExpectedLumiFB/CurrentLumuZ;
     bsm.Scale(Scale)


     if mass in xmapcolor:
                  colo=xmapcolor[mass]
                  #bsm.SetFillColor(colo)
                  bsm.SetLineColor( colo )
     bsm.SetAxisRange(Ymin, Ymax,"y");
     bsm.SetAxisRange(Xmin, Xmax,"x");
     bsm.Draw("same histo")
     bhitos[mass]=bsm
     xfile1.Close()



Xcut=CutOutlier
#xsum14000=data.Integral(data.FindBin(Xcut), data.FindBin(0));
#print("Summ=",xsum14000, " for cut=",Xcut)
x1=c1.XtoPad(Xcut)
x2=c1.XtoPad(Xcut)
ar6=TArrow(x1,Ymin,x2,c1.YtoPad(50000),0.05,">");
ar6.SetLineWidth(3)
ar6.SetLineStyle(1)
ar6.SetLineColor(2)
ar6.Draw("same")

ax=h.GetXaxis(); ax.SetTitleOffset(0.8)
ax.SetTitle( nameX );
ay=h.GetYaxis(); ay.SetTitleOffset(0.8)
ay.SetTitle( nameY );
ax.SetTitleOffset(1.1); ay.SetTitleOffset(1.6)
ax.Draw("same")
ay.Draw("same")

leg2=TLegend(0.16, 0.62, 0.82, 0.93);
leg2.SetBorderSize(0);
leg2.SetTextFont(62);
leg2.SetFillColor(10);
leg2.SetTextSize(0.04);
leg2.AddEntry(ttbar,"PYTHIA8 SM (t#bar{t})","lf")
#leg2.AddEntry(ttbar,"M(S)=M(X)/2","")
for b in range(len(bmass)):
     mass=bmass[b]
     leg2.AddEntry(bhitos[mass],"X("+str(mass)+ ")#rightarrow SH","lf")

#leg2.AddEntry(bsm,"SSM\; W\,{\\prime}(3 TeV) \\rightarrow Z\,{\\prime} (2\, TeV) W (all\, decays)","lf")
leg2.Draw("same")

leg3=TLegend(0.6, 0.62, 0.92, 0.80);
leg3.SetBorderSize(0);
leg3.SetTextFont(62);
leg3.SetFillColor(10);
leg3.SetTextSize(0.04);
leg3.AddEntry(ar6,"AD filter","l")
leg3.Draw("same")


# myText(0.65,0.52,1,0.04,UsedData0)
myText(0.6,0.65,1,0.04,"RMM input with:")
myText(0.6,0.6,1,0.04,"2 leptons pT>15 GeV")
myText(0.6,0.55,1,0.04,"2-bjets pT>20 GeV")
# ATLASLabel(0.19,0.89,0.14,0.03)

print (epsfig) 
gPad.RedrawAxis()
c1.Update()
ps1.Close()
if (myinput != "-b"):
              if (input("Press any key to exit") != "-9999"):
                         c1.Close(); sys.exit(1);




