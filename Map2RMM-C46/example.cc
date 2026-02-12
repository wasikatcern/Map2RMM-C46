 /***************************************************************************
 *  How to use ProMC files from HepSim, and how to build anti-KT jets 
 *  S.Chekanov (ANL) chekanov@anl.gov
 *  A library for HEP events storage and processing based on Google's PB   
 *  The project web site: http://atlaswww.hep.anl.gov/hepsim/
****************************************************************************/

#include<iostream>
#include<fstream>
#include<stdlib.h>
#include <sys/stat.h>
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TProfile2D.h>
#include <TProfile.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLine.h>
#include <TMath.h>
#include <TGraphErrors.h>
#include <TVector3.h>
#include <TLorentzVector.h>

#include "Pythia8/Pythia.h"

#include "ProMCBook.h"
#include "ProMCStat.h"
#include "ProMC.pb.h"

#include "map2rmm46/inc/SystemOfUnits.h"
#include "map2rmm46/inc/LParticle.h"
#include "map2rmm46/inc/CParticle.h"

using namespace fastjet;

using namespace std;
using namespace promc;

// project event
float**  map2rmm(const float CMS, const int maxN, const int maxNumberTypes,
                      const vector<LParticle> missing,
                      const vector<LParticle> jets,
                      const vector<LParticle> muons,
                      const vector<LParticle> electrons,
                      const vector<LParticle> photons);

std::vector<float> map2rmm46(const float CMS,
                             const int   maxN,
                             const int   maxNumberTypes,
                             const std::vector<LParticle>& missing,
                             const std::vector<LParticle>& jets,
                             const std::vector<LParticle>& muons,
                             const std::vector<LParticle>& electrons,
                             const std::vector<LParticle>& photons,
                             bool useFrob = true);


// find all files inside a directory
std::vector<std::string> open(std::string name = "data.in") {

        vector<std::string> ntup;
        ifstream myfile;
        myfile.open(name.c_str(), ios::in);

        if (!myfile) {
                cerr << " -> Can't open input file:  " << name << endl;
                exit(1);
        } else {
                cout << "-> Read data file=" << name << endl;
        }

        string temp;
        while (myfile >> temp) {
                //the following line trims white space from the beginning of the string
                temp.erase(temp.begin(), std::find_if(temp.begin(), temp.end(), not1(ptr_fun<int, int>(isspace))));
                if (temp.find("#") == 0) continue;
                ntup.push_back(temp);
        }
        cout << "-> Number of files=" << ntup.size()  << endl;
        myfile.close();

        for (unsigned int i=0; i<ntup.size(); i++) {
                cout << ".. file to analyse="+ntup[i] << endl;
        }

        return ntup;
};




int main(int argc, char **argv) {

        if (argc != 3) {
                cout << " Usage: ./example  input.list output.root " << endl;
                exit(0);
        }

        vector<string> names_debug = {"met","j","m","e","g"};

        // input file with paths to ProMC files
        string inputfile = argv[1];

        // output ROOT file
        string outputfile = argv[2];

        //int nevhisto=100; // draw this many events for the RMM
        int nevhisto=0; // disable individual RMM plots in this example

        // center-of-mass energy
        float CMS=13000.0;

        // maximum number of objects for each type
        int maxNumber=5;

        // number of object types
        int maxTypes=4; // jets, muons, electrons, photons

        cout << "Project using max number of each object=" << maxNumber << endl;
        cout << "Number of particle types==" << maxTypes << endl;
        const int mSize=maxTypes*maxNumber+1;
        // non-zero in triangular matrix
        int NonZero=((1+mSize)*mSize)/2;


        std::vector<string> Names1;
        Names1.push_back(names_debug[0]);

        for (int h = 1; h < maxTypes+1; h++) {
                for (int i = 1; i <  maxNumber+1; i++) {
                        ostringstream ss;
                        ss << i;
                        Names1.push_back(names_debug[h]+"_{"+ss.str()+"}");
                }
        }


        std::vector<string> Names2;
        for (unsigned int i=0; i<Names1.size(); i++) {
                cout << "Name=" << i << " " << Names1.at(i) << endl;
                Names2.push_back(Names1.at(i));
        }


        struct stat st;
        if(stat(outputfile.c_str(), &st) == 0) {
                cout << "File already exist. Exit." << endl;
                return 0;
        }


        TFile *RootFile = new TFile(outputfile.c_str(), "recreate", "Histograms"); // create output file

        //jet histograms
        TH1D * h_jet_pt = new TH1D("h_jet_pt", "jet pT",100,0,1000);
        TH1D * h_muon_pt = new TH1D("h_muon_pt", "muon pT",100,0,1000);
        TH1D * h_elec_pt = new TH1D("h_elec_pt", "electron pT",100,0,1000);
        TH1D * h_phot_pt = new TH1D("h_phot_pt", "photon pT",100,0,1000);

        TH1D * h_met = new TH1D("met_pt", "MET",100,0,1000);
        TH1D * h_metphi = new TH1D("met_phi", "MET phi",32,0,3.2);

        TH1D * h_mjj = new TH1D("m_jj", "dijet mass",50,0,3000);
        TH1D * h_mmu = new TH1D("m_mumu", "dimuon mass",50,0,3000);
        TH1D * h_mee = new TH1D("m_ee", "dielectron mass",50,0,3000);
        TH1D * h_mgamgam = new TH1D("m_gamgam", "diphoton mass",50,0,3000);

        TH1D * h_ajetjet = new TH1D("a_jetjet", "angle jet-jet",50,0,3.2);
        TH1D * h_amuonmuon = new TH1D("a_muonmuon", "angle muon-muon",50,0,3.2);
        TH1D * h_aelectron = new TH1D("a_ee", "angle e-e",50,0,3.2);
        TH1D * h_agamgam = new TH1D("a_gamgam", "angle gam-gam",50,0,3.2);

        TH1D * h_ajetjet_cos = new TH1D("a_jetjet_cos", "angle jet-jet 1-cos",50,0,1.0);
        TH1D * h_amuonmuon_cos = new TH1D("a_muonmuon_cos", "angle muon-muon 1-cos",50,0,1.0);
        TH1D * h_aelectron_cos = new TH1D("a_ee_cos", "angle e-e 1-cos",50,0,1.0);
        TH1D * h_agamgam_cos = new TH1D("a_gamgam_cos", "angle gam-gam 1-cos",50,0,1.0);

        TH1D * h_matAng = new TH1D("matAng", "matrix part with angles",500,0,1.0);
        TH1D * h_matPt = new TH1D("matPt", "matrix part (diag) with pT",500,0,1.0);
        TH1D * h_matMass = new TH1D("matMass", "matrix part with masses",500,0,1.0);
        TH1D * h_matAngNorm = new TH1D("matAngNorm", "matrix part with angles (norm)",500,0,1.0);

        TH2D * h_proj = new TH2D("proj", "values", mSize, 0, (double)(mSize), mSize, 0, (double)(mSize));
        TH2D * h_events = new TH2D("events", "events", mSize, 0, (double)(mSize), mSize, 0, (double)mSize);

        TProfile2D * h_prof = new TProfile2D("profile", "profile", mSize, 0, (double)(mSize), mSize, 0, (double)mSize);
        h_prof->SetOption("COLZ TEXT");

        TProfile * h_info = new TProfile("info","info histogram", 100,0,100);
        h_info->Fill("Events weighted",0);
        h_info->Fill("Nr of dijets",0);
        h_info->Fill("Dijet section",0);
        h_info->Fill("Cross section",0);

        TH1D * h_debug = new TH1D("debug", "debug",2,0,2);

        TH2D * h_proje[ nevhisto];
        for(int i=0; i<nevhisto; i++) h_proje[i] = new TH2D(Form("event_%02d",i),Form("event_%02d",i), mSize, 0, (double)(mSize), mSize, 0, (double)mSize);

        // debug histogram
        h_debug->Fill("Start",1.0);

        // initialize RMM binning to zeroes
        for (int h = 0; h < mSize; h++) {
                for (int w = 0; w < mSize; w++)
                {
                        int i1=h;
                        int i2=w;
                        h_proj->Fill((Names2.at(i1)).c_str(),  (Names1.at(i2)).c_str(),0);
                        h_prof->Fill((Names2.at(i1)).c_str(),  (Names1.at(i2)).c_str(),0);
                        h_events->Fill((Names2.at(i1)).c_str(), (Names1.at(i2)).c_str(),0);
                        for(int i=0; i<nevhisto; i++) h_proje[i]->Fill((Names2.at(i1)).c_str(), (Names1.at(i2)).c_str(),0);
                }}

        TTree *  m_tree  = new TTree("inputNN","inputNN");
        m_tree->SetAutoFlush(100000);
        Int_t m_id;
        std::vector<Double32_t> m_c46;

//  std::vector<Double32_t> m_proj; // removed for Map2RMM46
//         std::vector<UInt_t> m_proj_index1; // removed for Map2RMM46
//         std::vector<UInt_t> m_proj_index2; // removed for Map2RMM46

//  std::vector<Double32_t> m_multi; // removed for Map2RMM46
        m_tree->Branch("id",  &m_id);
        m_tree->Branch("c46",   &m_c46);
//         m_tree->Branch("proj",   &m_proj); // removed for Map2RMM46
//         m_tree->Branch("proj_index1",   &m_proj_index1); // removed for Map2RMM46
//         m_tree->Branch("proj_index2",   &m_proj_index2); // removed for Map2RMM46

//  m_tree->Branch("multiplicity",   &m_multi); // removed for Map2RMM46

        /*
                Int_t id; 
                Double32_t g,e,m,j,gg,eg,mg,jg,ee,me,je,mm,jm,jj;
                const char* tupleNames = "id:g:e:m:j:gg:eg:mg:jg:ee:me:je:mm:jm:jj" ;
                NTuple *nt = new NTuple("RMM", "RMM",tupleNames,&id,&g,&e,&m,&j,&gg,&eg,&mg,&jg,&ee,&me,&je,&mm,&jm,&jj, 16);
                //nt->project(arr, 16 );
        */

        ProMCBook *epbook = NULL;
        ProMCEvent promc_evt;

        vector<string> flist = open(inputfile);

        double Lumi=1000; // 1/fb
        double Nevents=0;
        double xcross=0;
        int nfiles=0;
        double cross=0;
        int ntot=0;

        double R=0.4; //  anti-kt R parameter
        double pTMin = 20.0; // min pT
        double etaMax = 4.5; // max eta

        double lumi = 0;

        int nDiJets=0;

        // loop over files
        for(unsigned int ff=0;ff<flist.size();ff++) {

                string file=flist[ff];
                cout << "-> Open file=" << file << endl;
                epbook = new ProMCBook(file.c_str(),"r");

                // histogram with Cross section
                TH1F * h_cros = new TH1F("h_cros","Cross section",1,0,1);

                // total cross section * efficiency
                ProMCStat stat = epbook->getStatistics();
                float xs = stat.cross_section();
                float md = stat.cross_section_error();
                float sumw = stat.sum_of_weights();
                float sumw2 = stat.sum_of_weights2();
                int nevents=int(stat.ntried());
                cout << " Cross section      = " << xs << " +- " << md << " pb" << endl;
                cout << " Number of events  = " << nevents << endl;
                cout << " Sum of weights    = " << sumw << endl;
                cout << " Sum of weights^2  = " << sumw2 << endl;

                // we will use this
                // xs = sum_w * <sigma>/ntried
                //double cross_section = sumw * xs / double(ntried);

                double wgt=xs/double(nevents); // each event weight
                h_cros->Fill(0.5,xs);

                double wgt_lumi=wgt*Lumi; // event scaling for 1/fb

                lumi = lumi + h_cros->GetMean()*Lumi;
                cout << " -> File luminosity=" << h_cros->GetMean()*Lumi << " pb-1" << endl;
                cout << " -> Files luminosity (after this file) =" << lumi << " pb-1" << endl;

                delete h_cros;

                string desc=epbook->getDescription();
                cout << desc << endl;

                int nevent=epbook->sizeEvents();

                // debug info
                h_info->Fill("Cross section",xs);
                h_info->Fill("Events weighted",nevents);

                // event loop
                int event=0;
                while (epbook->next()) {

                        if (nevent>10) {
                                if (event%5000==0) cout << "Selected events=" << event << " (out of " << nevent << ")" << endl;
                        }

                        epbook->read(promc_evt);

                        const ProMCHeader &header = promc_evt.header();
                        const ProMCEvent_Lumi &lumi_hdr = header.lumi();

                        // get event with units
                        const ProMCEvent &eve = promc_evt.event();
                        const ProMCEvent_Units &m_unit = header.momentum_unit();
                        const ProMCEvent_Units &l_unit = header.length_unit();

                        double  mom_unit = pow(10.0, m_unit.mult());
                        double  len_unit = pow(10.0, l_unit.mult());

                        // Event weight
                        double weight = lumi_hdr.weight();

                        // particles
                        int m_boson1 = 0;
                        int m_boson2 = 0;
                        int m_lep1=0,m_lep2=0;
                        int m_nu1=0,m_nu2=0;

                        int nlist=0;
                        TLorentzVector boson1,boson2;
                        TLorentzVector lep1,lep2;
                        TLorentzVector nu1,nu2;

                        // read jets, muons, electrons, photons, MET
                        vector<LParticle> jets,muons,electrons,photons,missing;
                        jets.clear();
                        muons.clear();
                        electrons.clear();
                        photons.clear();
                        missing.clear();

                        //anti-kt jets
                        vector<fastjet::PseudoJet> particlesForJets;
                        particlesForJets.clear();

                        // loop over particles
                        for (int i = 0; i < eve.particles_size(); i++) {

                                const ProMCEvent_Particle &pa = eve.particles(i);
                                int id=pa.id();
                                int status=pa.status();

                                double px=pa.px()/mom_unit;
                                double py=pa.py()/mom_unit;
                                double pz=pa.pz()/mom_unit;
                                double ee=pa.energy()/mom_unit;

                                TLorentzVector l(px,py,pz,ee);

                                if (status != 1) continue; // final states only

                                if (fabs(l.Eta())>etaMax) continue;

                                int PDG = pa.pdg_id();

                                // neutrinos
                                if (PDG==12 || PDG==14 || PDG==16 ) {
                                        // missing energy
                                        LParticle p;
                                        p.SetP(l);
                                        p.SetType(1);
                                        missing.push_back(p);
                                        continue;
                                }

                                // muons
                                if (fabs(PDG)==13) {
                                        if (l.Pt()>pTMin) {
                                                LParticle p;
                                                p.SetP(l);
                                                p.SetType(1);
                                                muons.push_back(p);
                                                h_muon_pt->Fill(l.Pt());
                                        }
                                        continue;
                                }

                                // electrons
                                if (fabs(PDG)==11) {
                                        if (l.Pt()>pTMin) {
                                                LParticle p;
                                                p.SetP(l);
                                                p.SetType(1);
                                                electrons.push_back(p);
                                                h_elec_pt->Fill(l.Pt());
                                        }
                                        continue;
                                }

                                // photons
                                if (fabs(PDG)==22) {
                                        if (l.Pt()>pTMin) {
                                                LParticle p;
                                                p.SetP(l);
                                                p.SetType(1);
                                                photons.push_back(p);
                                                h_phot_pt->Fill(l.Pt());
                                        }
                                        continue;
                                }

                                // everything else goes into jets
                                if (l.Pt()>0.0) {
                                        fastjet::PseudoJet particle(px,py,pz,ee);
                                        particle.set_user_index(PDG);
                                        particlesForJets.push_back(particle);
                                }
                        } // end loop over particles


                        // build jets with anti-kt
                        fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, R);
                        fastjet::ClusterSequence clustSeq(particlesForJets, jet_def);
                        vector<fastjet::PseudoJet> jets_akt = sorted_by_pt(clustSeq.inclusive_jets(pTMin));

                        // fill LParticle jets from PseudoJets
                        for (unsigned int i=0; i<jets_akt.size(); i++) {
                                TLorentzVector l(jets_akt[i].px(),jets_akt[i].py(),jets_akt[i].pz(),jets_akt[i].E());
                                if (fabs(l.Eta())>etaMax) continue;
                                LParticle p;
                                p.SetP(l);
                                p.SetType(1);
                                jets.push_back(p);
                                h_jet_pt->Fill(l.Pt());
                        }


                        // check dijets
                        if (jets.size()>1) {
                                for (unsigned int i1=0; i1<jets.size()-1; i1++){
                                        LParticle LPP1=jets.at(i1);
                                        TLorentzVector LP1=LPP1.GetP();
                                        for (unsigned int i2=i1+1; i2<jets.size(); i2++){
                                                LParticle LPP2=jets.at(i2);
                                                TLorentzVector LP2=LPP2.GetP();
                                                TLorentzVector LPP=LP1+LP2;
                                                h_mjj->Fill(LPP.M());

                                                double ang=LP1.Angle(LP2.Vect());
                                                double cang=1-cos(ang);
                                                h_ajetjet->Fill(ang);
                                                h_ajetjet_cos->Fill(cang);
                                        }
                                }
                        }

                        // check dimuons
                        if (muons.size()>1) {
                                for (unsigned int i1=0; i1<muons.size()-1; i1++){
                                        LParticle LPP1=muons.at(i1);
                                        TLorentzVector LP1=LPP1.GetP();
                                        for (unsigned int i2=i1+1; i2<muons.size(); i2++){
                                                LParticle LPP2=muons.at(i2);
                                                TLorentzVector LP2=LPP2.GetP();
                                                TLorentzVector LPP=LP1+LP2;
                                                h_mmu->Fill(LPP.M());

                                                double ang=LP1.Angle(LP2.Vect());
                                                double cang=1-cos(ang);
                                                h_amuonmuon->Fill(ang);
                                                h_amuonmuon_cos->Fill(cang);
                                        }
                                }
                        }

                        // check dielectrons
                        if (electrons.size()>1) {
                                for (unsigned int i1=0; i1<electrons.size()-1; i1++){
                                        LParticle LPP1=electrons.at(i1);
                                        TLorentzVector LP1=LPP1.GetP();
                                        for (unsigned int i2=i1+1; i2<electrons.size(); i2++){
                                                LParticle LPP2=electrons.at(i2);
                                                TLorentzVector LP2=LPP2.GetP();
                                                TLorentzVector LPP=LP1+LP2;
                                                h_mee->Fill(LPP.M());

                                                double ang=LP1.Angle(LP2.Vect());
                                                double cang=1-cos(ang);
                                                h_aelectron->Fill(ang);
                                                h_aelectron_cos->Fill(cang);
                                        }
                                }
                        }

                        // check diphotons
                        if (photons.size()>1) {
                                for (unsigned int i1=0; i1<photons.size()-1; i1++){
                                        LParticle LPP1=photons.at(i1);
                                        TLorentzVector LP1=LPP1.GetP();
                                        for (unsigned int i2=i1+1; i2<photons.size(); i2++){
                                                LParticle LPP2=photons.at(i2);
                                                TLorentzVector LP2=LPP2.GetP();
                                                TLorentzVector LPP=LP1+LP2;
                                                h_mgamgam->Fill(LPP.M());

                                                double ang=LP1.Angle(LP2.Vect());
                                                double cang=1-cos(ang);
                                                h_agamgam->Fill(ang);
                                                h_agamgam_cos->Fill(cang);
                                        }
                                }
                        }

                        // build MET from neutrinos if missing vector is empty
                        TLorentzVector lmet(0,0,0,0);
                        double pxsum=0, pysum=0;
                        for (unsigned int i=0; i<missing.size(); i++) {
                                TLorentzVector lp = missing[i].GetP();
                                pxsum += lp.Px();
                                pysum += lp.Py();
                        }
                        double met_pt=sqrt(pxsum*pxsum+pysum*pysum);
                        double met_phi=0;
                        if (met_pt>0) met_phi=atan2(pysum,pxsum);

                        if (met_pt<20.0) {
                                met_pt=0;
                                met_phi=0;
                        }

                        lmet.SetPtEtaPhiM(met_pt, 0, met_phi, 0);
                        h_met->Fill(met_pt);
                        lmet.SetPz(0);
                        LParticle pmet;
                        pmet.SetP(lmet);
                        pmet.SetType(1);
                        missing.push_back(pmet);


                        m_c46.clear();
//                         m_proj_index1.clear(); // removed for Map2RMM46
//                         m_proj_index2.clear(); // removed for Map2RMM46

//                      m_multi.clear(); // removed for Map2RMM46
                        //m_m=muons.size();
                        //m_e=electrons.size();
                        //m_g=photons.size();
                        //m_j=jets.size();
                        m_id=0;


                        // Require at least something non-trivial in the event
                        if (met_pt==0 && jets.size()==0 && muons.size()==0 && electrons.size() == 0 && photons.size()==0) continue;


                        // === Build 46-dimensional RMM-C46 directly (no intermediate RMM file stored) ===
                        std::vector<float> c46 = map2rmm46(CMS,
                                                          maxNumber,
                                                          maxTypes,
                                                          missing,
                                                          jets,
                                                          muons,
                                                          electrons,
                                                          photons,
                                                          true);  // true = use Frobenius aggregation

                        // copy to Double32_t vector for ROOT I/O
                        m_c46.clear();
                        for (std::size_t ic = 0; ic < c46.size(); ++ic) {
                                m_c46.push_back(static_cast<Double32_t>(c46[ic]));
                        }

                        // optional: you can set a meaningful event id here if desired
                        m_id = 0;


                        event++;
                        m_tree->Fill();

                } // end event loop


                ProMCStat stat2 = epbook->getStatistics();
                cross=stat2.cross_section_accumulated();
                epbook->close(); // close
                nfiles++;
                xcross=xcross+cross;


        } // end loop over all files



        xcross=xcross/(double)nfiles; // average cross for all files
        cout << "Total events=" << ntot << endl;
        cout << "Total files=" << nfiles << endl;
        cout << "Total cross section=" << xcross << " pb" << endl;

        h_info->Fill("Cross section",xcross);

        h_info->Fill("Nr of dijets",nDiJets); // calibration check
        h_info->Fill("Dijet section",nDiJets/lumi); // calibration check
        cout << " Nr of dijets=" << nDiJets << endl;
        cout <<" Observed cross section=" << nDiJets/lumi << " pb " << endl;


        RootFile->Write();
        RootFile->Print();
        RootFile->Close();

        cout << "Writing ROOT file "+ outputfile << endl;

        return 0;
}

