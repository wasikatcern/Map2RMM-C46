## Map2RMM-C46 (Compact 46)

Map2RMM-C46 is a C++ library and example program that converts Monte-Carlo collision events into a 46-dimensional compact representation of the Rapidity–Mass Matrix (RMM), called RMM-C46.

Unlike the original Map2RMM framework, this version:

does not write or store the full RMM,

does not use sparse RMM branches (proj, proj_index1, proj_index2),

builds the full RMM only in memory,

immediately compresses it into a 46-element feature vector, and

writes only the compact C46 vector to a ROOT TTree (inputNN.c46).

This makes the output highly compact and directly usable for machine-learning pipelines.

The example included here reads ProMC Monte-Carlo files from HepSim, reconstructs physics objects (jets, muons, electrons, photons, MET), builds an RMM internally, compresses it to C46, and writes the results to a ROOT output file.

## About the RMM-C46 Representation

The standard RMM is an N×N matrix encoding:

MET

Transverse energies (ET)

Transverse masses (mT)

Rapidity-based angular variables

Pairwise invariant masses

Map2RMM-C46 compresses this matrix into 46 physics-motivated features:

1 MET term

15 ET / mT / hL terms

15 “h-terms” (rapidity-based blocks)

15 “m-terms” (mass-based blocks)

Up to 5 object types are supported.
The example uses:

jets, muons, electrons, photons

Any unused type is padded with zeros.

You may configure multiplicity and number of types in example.cc:

const int maxNumber = 5;  // multiplicity per object type
const int maxTypes  = 4;  // jets, muons, electrons, photons

## Installation

These instructions apply to any Linux system with GCC ≥ 4.

# 1. Install dependencies

You must have:

ProMC

ROOT

FastJet

Verify:

echo $PROMC
echo $ROOTSYS
echo $FASTJET


Each should print a valid directory.

# 2. Build the Map2RMM-C46 library

Go to:

map2rmm46/


Compile:

make


This generates:

map2rmm46/lib/libmap2rmm.so
map2rmm46/lib/libmap2rmm_static.a

# 3. Build the example analysis program

Return to the project root directory and run:

make


This links against:

FastJet

ROOT

ProMC

map2rmm46

The resulting executable is:

./example

# 4. Download ProMC data

Place files into the data/ directory.
For example:

hs-get tev100_higgs_ttbar_mg5 data


See HepSim documentation for more datasets.

# 5. Run the Map2RMM-C46 analysis

Create a file listing input ProMC files (e.g. data.in), then run:

./example data.in output.root


This produces:

physics histograms

a ROOT TTree inputNN containing the compact 46-dimensional feature vectors

Output Format

The resulting output.root contains:

TTree: inputNN
   id    (Int_t)
   c46   (std::vector<Double32_t>)   # 46-dimensional RMM-C46 vector


There is no full RMM stored.

There are no sparse matrices (proj, proj_index1, proj_index2).

The entire RMM is compressed into the single branch c46.

Preparing ML-Ready CSV Files

A new prepare.py script is included.

It converts the ROOT file into a CSV file containing:

id, C46_1, C46_2, ..., C46_46


## Usage:

python prepare.py output.root c46_output.csv


Example output:

id,C46_1,C46_2,...,C46_46
0,0.0041,0.8123,...,0.1231
0,0.0009,1.2345,...,0.9912
...


This CSV is directly compatible with PyTorch, TensorFlow, scikit-learn, XGBoost, etc.

## Reference

Original RMM paper:
S. Chekanov,
Imaging particle collision data for event classification using machine learning,
NIMA A931 (2019) 92
https://arxiv.org/abs/1805.11650

Further applications:
S. Chekanov,
Machine learning using rapidity–mass matrices for event classification problems in HEP,
ANL-HEP-147750
https://arxiv.org/abs/1810.06669

Author

Original RMM was developed by S. Chekanov (ANL) and updated for RMM-C46 by Wasikul Islam (Univ. of Wisconsin-Madison).
This version extends the framework by integrating a native RMM → C46 compression engine directly into the C++ backend, optimized for modern high-energy-physics ML workflows.

