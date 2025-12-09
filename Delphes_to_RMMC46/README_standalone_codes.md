## Standalone codes (Start from Delphes files or RMM files):

Inside "Delphes_to_RMMC46" directory, we have codes to :

1. Convert Delphes files into skimmed files with fewer inforamtion:
   
root -l -b 'skim_Delphes.C("Hplus_1800GeV_SLHA2_delphes.root", "skimmed_delphes.root", 13000, 100)'

2. Convert skimmed Delphes files into RMM matrix, like the folllowing :
root -l -b -q 'make_RMMs.C("skimmed_delphes.root","rmm_events_100.csv",100)'

root -l -b -q 'make_RMMs.C' --args --in skimmed_delphes.root --out rmm_events_100.csv --nevents 100

root -l -b -q 'make_RMMs.C' --args --in skimmed_delphes.root --out rmm.csv --selroot my_selected.root --iconfig 0

4. Convert RMM into RMM-C46 :
   
python make_RMM_C46.py --csv out/pythia8_X1500GeV_HH2bbll_data100percent.csv.gz --out hh_1500_c46_frob.csv

python make_RMM_C46.py --csv out/pythia8_X1500GeV_HH2bbll_data100percent.csv.gz --out hh_1500_c46_add.csv --style add 
