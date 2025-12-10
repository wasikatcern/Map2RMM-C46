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

# Exact commands used to convert all RMMs into RMM-C46:

python make_RMM_C46.py --csv out/pythia8_X500GeV_HH2bbll_data100percent.csv.gz --out out_C46/hh_500_c46_frob.csv

python make_RMM_C46.py --csv out/pythia8_X700GeV_HH2bbll_data100percent.csv.gz --out out_C46/hh_700_c46_frob.csv

python make_RMM_C46.py --csv out/pythia8_X1000GeV_HH2bbll_data100percent.csv.gz --out out_C46/hh_1000_c46_frob.csv

python make_RMM_C46.py --csv out/pythia8_X1500GeV_HH2bbll_data100percent.csv.gz --out out_C46/hh_1500_c46_frob.csv

python make_RMM_C46.py --csv out/pythia8_X2000GeV_HH2bbll_data100percent.csv.gz --out out_C46/hh_2000_c46_frob.csv


python make_RMM_C46.py --csv out/pythia8_X500GeV_SH2bbll_data100percent.csv.gz --out out_C46/sh_500_c46_frob.csv

python make_RMM_C46.py --csv out/pythia8_X700GeV_SH2bbll_data100percent.csv.gz --out out_C46/sh_700_c46_frob.csv

python make_RMM_C46.py --csv out/pythia8_X1000GeV_SH2bbll_data100percent.csv.gz --out out_C46/sh_1000_c46_frob.csv

python make_RMM_C46.py --csv out/pythia8_X1500GeV_SH2bbll_data100percent.csv.gz --out out_C46/sh_1500_c46_frob.csv

python make_RMM_C46.py --csv out/pythia8_X2000GeV_SH2bbll_data100percent.csv.gz --out out_C46/sh_2000_c46_frob.csv


python make_RMM_C46.py --csv out/tev13.6pp_pythia8_ttbar_2lep_data10percent.csv.gz --out out_C46/wzjet_c46_frob.csv

python make_RMM_C46.py --csv out/tev13.6pp_pythia8_wzjet_2lep_data10percent.csv.gz --out out_C46/tt_c46_frob_10.csv
