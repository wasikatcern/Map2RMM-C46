/**************************************************
*       Project event into a matrix	
*	For details see http://arxiv.org/abs/1002.3982
**************************************************/

using namespace std; 
#include<iostream>
#include "TMath.h"
#include "TLorentzVector.h"
#include "SystemOfUnits.h"
#include "LParticle.h"
#include "CParticle.h"
#include <string>
#include <vector>

static const int pairs[15][2] = {
    {0,0},
    {1,0},{1,1},
    {2,0},{2,1},{2,2},
    {3,0},{3,1},{3,2},{3,3},
    {4,0},{4,1},{4,2},{4,3},{4,4}
};



/*

 We build matrix like this assuming maxN=2:
 j -jet
 m- muons
 e- electron
 g -photon

 All variables are ratio to CMS 

 
     |   MET    j1       j2       mu1     mu2     e1      e2     g1      g2
 ------------------------------------------------------------------------- 
 MET |  MT^2   MT(j1)  MT(j2)   MT(m1)   MT(m2)  MT(e1)  MT(e2) MT(g1) MT(g2)  
     |      
 j1  |  h_j1     et(j1) m(j1,j2) 
     |  
 j2  |  h_j2               et(j1)
     |
 m1  |                         et(m1) 
     | 
 m2  |                                 et(m2) 
     |
 e1  |   
     |
 e2  |   
     |
 g1  |   
     |
 g2  |  


*/


// aggregate the values..
double getFrobeniusAggregationET(const unsigned int INCR, const unsigned int Nparticles,  const unsigned int maxNumber, float** outMatrix){
  double vv = 0;
  for (unsigned int k1 = 0; k1<Nparticles; k1++) {
     double v=outMatrix[k1+INCR*maxNumber+1][k1+INCR*maxNumber+1];
     vv = vv + v*v;
  }
  vv=TMath::Sqrt(   vv );
  return vv;

};	

// aggregate the values..
double getFrobeniusAggregationMET(const unsigned int INCR, const unsigned int Nparticles, const unsigned int maxNumber,  float** outMatrix){
  double vv = 0;
  for (unsigned int k1 = 0; k1<Nparticles; k1++) {
     double v=outMatrix[0][k1+INCR*maxNumber+1];
     vv = vv + v*v;
  }
  vv=TMath::Sqrt(   vv );
  return vv;

};


// aggregate the values..
double getFrobeniusAggregationLong(const unsigned int INCR, const unsigned int Nparticles, const unsigned int maxNumber,  float** outMatrix){
  double vv = 0;
  const float angNorm =  0.15; // reduce (ad-hoc to keep it under 1)
  for (unsigned int k1 = 0; k1<Nparticles; k1++) {
     double v=angNorm*outMatrix[k1+INCR*maxNumber+1][0];
     vv = vv + v*v;
  }
  vv=TMath::Sqrt(   vv );
  return vv;

};

/*
// aggregate the rapidity values..
double getFrobeniusAggregationRapidityDiff(const unsigned int INCR1, const unsigned int INCR2,  
		const unsigned int Nparticles1, const unsigned int Nparticles2, const unsigned int maxNumber,   float** outMatrix){
  double vv = 0;
  const float angNorm =  0.15; // reduce (ad-hoc to keep it under 1)

  for (unsigned int k1 = 0; k1<Nparticles1; k1++) {
      for (unsigned int k2 = 0; k2<Nparticles2; k2++) { 
       if (k1>k2) {
       double v= angNorm*outMatrix[k1+INCR1*maxNumber+1][k2+INCR2*maxNumber+1]; 
       vv += v*v;
       }
    }
  }
  vv=TMath::Sqrt(   vv );
  return vv;

};

// aggregate the inv mass  values..
double getFrobeniusAggregationMass(const unsigned int INCR1, const unsigned int INCR2, 
                const unsigned int Nparticles1, const unsigned int Nparticles2, const unsigned int maxNumber,   float** outMatrix){
  double vv = 0;
  cout << "Nparticles1=" << Nparticles1 << " Nparticles2=" << Nparticles2 << endl; 
  
  for (unsigned int k1 = 0; k1<Nparticles1; k1++) {
      for (unsigned int k2 = 0; k2<Nparticles2; k2++) {
       if (k1<k2) {
       double v=outMatrix[k1+INCR1*maxNumber+1][k2+INCR2*maxNumber+1];
       cout << "C# index k1=" << k1 << " k2=" << k2 << " MATRIX=" << "("<< k1+INCR1*maxNumber+1 << ","<< k2+INCR2*maxNumber+1 << ") value="<< v  << endl;
       vv += v*v;
       }
    }
  }
  cout << "Sum="<< vv << "  sqrt="<< TMath::Sqrt(   vv ) << endl; 
  vv=TMath::Sqrt(   vv );
  return vv;
};
*/

// ============================================================
// Aggregate invariant-mass values for RMM-C46 (Updated Version)
// ============================================================

// aggregate invariant-mass values (Updated C46 definition)
double getFrobeniusAggregationMass(
    const unsigned int INCR1,
    const unsigned int INCR2,
    const unsigned int /*Nparticles1*/,   // intentionally ignored
    const unsigned int /*Nparticles2*/,   // intentionally ignored
    const unsigned int maxNumber,
    float** outMatrix)
{
    double sum = 0.0;
    const bool sameType = (INCR1 == INCR2);

    for (unsigned int k1 = 0; k1 < maxNumber; ++k1) {
        for (unsigned int k2 = 0; k2 < maxNumber; ++k2) {

            // Same-type blocks: strict upper triangle only
            if (sameType && k1 >= k2) continue;

            // Cross-type blocks: full block (NO triangle cut)
            const double v =
                outMatrix[k1 + INCR1 * maxNumber + 1]
                         [k2 + INCR2 * maxNumber + 1];

            sum += v * v;
        }
    }

    return TMath::Sqrt(sum);
}

// ============================================================
// Aggregate rapidity-difference values for RMM-C46 (Updated Version)
// ============================================================
double getFrobeniusAggregationRapidityDiff(
    const unsigned int INCR1,
    const unsigned int INCR2,
    const unsigned int /*Nparticles1*/,   // intentionally ignored
    const unsigned int /*Nparticles2*/,   // intentionally ignored
    const unsigned int maxNumber,
    float** outMatrix)
{
    double sum = 0.0;
    const bool sameType = (INCR1 == INCR2);

    for (unsigned int k1 = 0; k1 < maxNumber; ++k1) {
        for (unsigned int k2 = 0; k2 < maxNumber; ++k2) {

            // Same-type blocks: strict lower triangle only
            if (sameType && k1 <= k2) continue;

            // Cross-type blocks: full block (NO triangle cut)
            const double v =
                outMatrix[k1 + INCR1 * maxNumber + 1]
                         [k2 + INCR2 * maxNumber + 1];

            sum += v * v;
        }
    }

    return TMath::Sqrt(sum);
}




// return angle product  
float getAngle(const float CMS, const TLorentzVector p1, const TLorentzVector p2){
      //float ang=p1.Angle(p2.Vect());
      //ang=(1-TMath::Cos(ang));
      //double  ene=sqrt(p1.E()*p2.E())/CMS; 
      double y1=p1.Rapidity();
      double y2=p2.Rapidity();
      double HL=TMath::CosH( 0.5*(y2-y1) )-1;
      return float(HL);
      }

// get masses
float getMass(const float CMS, const TLorentzVector p1, const TLorentzVector p2){
      TLorentzVector pp=p1+p2;
      float xmass=pp.M()/CMS; 
      //if (xmass>1) cout << "xmass=" << xmass << " CMS=" << CMS << endl;
      return xmass;
      } 

// get cosh(y)-1
// y=0 correspongs 0 
// y=4 correspondsi 2.7/100
float getHL(const TLorentzVector p1){
      double y=p1.Rapidity();
      double HL=TMath::CosH(y)-1; 
      //cout << "y=" << y << " HL=" << HL << endl;
      return (float)(HL);
      }

// return transverse mass using experimental method 
float getMT(const TLorentzVector met, const TLorentzVector jet){

     /* massless approximation 
     double DeltaPhi = jet.DeltaPhi(met);
     double s=2*jet.Et()*met.Et()*(1-TMath::Cos( DeltaPhi ));
     // Mt for massless
     double Mt_massless=0;
     if (s>0) Mt_massless=TMath::Sqrt(s);
     */

     // Mt exact 
     double ss= (jet.Et()+met.Et()) * (jet.Et() +met.Et()) -
                ( (jet.Px()+met.Px())*(jet.Px() +met.Px()) ) - ( (jet.Py()+met.Py())*(jet.Py() +met.Py()) );
     double Mt_exact=0;
     if (ss>0) Mt_exact= TMath::Sqrt(ss);

     // Mt corrected
     //double mW_t = TMath::Sqrt(jet.Et()*jet.Et() + met.Et()*met.Et() - 2*jet.Et()*met.Et()*TMath::Cos(jet.DeltaPhi(met)));
     //cout << "Mt_massless=" << Mt_massless << " Mt from sum=" << (met+jet).Mt() << " Mt_exact=" << Mt_exact << endl; 
     return (float)Mt_exact;
     }


// create matrix with event projection
// maxN max number of particles of a given type
// maxNumberTypes total number of types
float*  projectevent_c46(const float CMS, const int maxN, const int maxNumberTypes, 
                      const vector<LParticle> missing,   
                      const vector<LParticle> jets, 
                      const vector<LParticle> bjets,
                      const vector<LParticle> muons, 
                      const vector<LParticle> electrons, 
                      const vector<LParticle> photons) {


            //just to make sure, sort 
            /*
            if (jets.size()>1)  std::sort(jets.begin(), jets.end(), greater<LParticle>() ) ;
            if (muons.size()>1) std::sort(muons.begin(), muons.end(), greater<LParticle>() ) ;
            if (electrons.size()>1) std::sort(electrons.begin(), electrons.end(), greater<LParticle>() ) ;
            if (photons.size()>1) std::sort(photons.begin(), photons.end(), greater<LParticle>() ) ;
            */

            const int maxNumber=maxN; // max number for each object
            const int maxTypes=maxNumberTypes;  // max numbers of types
            const int maxSize=maxNumber*maxTypes+1; // + MET 

            const int height=maxSize;
            const int width=maxSize;

	    // initialize the standard RMM
            float** outMatrix = 0;
            outMatrix = new float*[height];
            for (int h=0; h<height; h++) outMatrix[h] = new float[width]; 

            unsigned  int INCR=0;         
            for(int i=0; i<maxSize; i++)
                 for(int j=0; j<maxSize; j++)  outMatrix[i][j]=0;


                  // invariant masses of same objects 
                  INCR=0; 
                  int ii=0;
                  LParticle MET=(LParticle)missing.at(ii);
                  TLorentzVector LMET=MET.GetP();
                  outMatrix[ii+INCR*maxNumber][ii+INCR*maxNumber]=LMET.Et()/CMS;


                  INCR=0; 
                  unsigned  int mJets=maxNumber;
                  if (jets.size()<mJets) mJets=jets.size();
                  if (mJets>0) {
                  for (unsigned int k1 = 0; k1<mJets; k1++) {
                    LParticle LPP1=(LParticle)jets.at(k1);
                    TLorentzVector LP1=LPP1.GetP();
                    if (LMET.Et()>0) outMatrix[0][k1+INCR*maxNumber+1]=getMT(LMET, LP1)/CMS;
                    outMatrix[k1+INCR*maxNumber+1][0]=getHL(LP1);

                    if (k1 == 0) {
                        outMatrix[k1+INCR*maxNumber+1][k1+INCR*maxNumber+1]=LP1.Et()/CMS;
                    } else if (k1>0)  {
                          LParticle LPP3=(LParticle)jets.at(k1-1);
                          TLorentzVector LP3=LPP3.GetP();
                          float imbalance=(LP3.Et() - LP1.Et())/(LP3.Et() + LP1.Et()); 
                          outMatrix[k1+INCR*maxNumber+1][k1+INCR*maxNumber+1]=imbalance;
                    }
                   
                    //  non-diagonal 
                    for (unsigned int k2 = 0; k2<mJets; k2++) {
                     LParticle LPP2=(LParticle)jets.at(k2);
                     TLorentzVector LP2=LPP2.GetP();
                     if (k1<k2)    outMatrix[k1+INCR*maxNumber+1][k2+INCR*maxNumber+1]=getMass(CMS,LP1,LP2);
                     if (k1>k2)    outMatrix[k1+INCR*maxNumber+1][k2+INCR*maxNumber+1]=getAngle(CMS,LP1,LP2);
                    }
                  }
                 }




                  INCR=1;
                  unsigned  int mJetsB=maxNumber;
                  if (bjets.size()<mJetsB) mJetsB=bjets.size();
                  if (mJetsB>0) {
                  for (unsigned int k1 = 0; k1<mJetsB; k1++) {
                    LParticle LPP1=(LParticle)bjets.at(k1);
                    TLorentzVector LP1=LPP1.GetP();
                    if (LMET.Et()>0) outMatrix[0][k1+INCR*maxNumber+1]=getMT(LMET, LP1)/CMS;
                    outMatrix[k1+INCR*maxNumber+1][0]=getHL(LP1);

                    if (k1 == 0) {
                        outMatrix[k1+INCR*maxNumber+1][k1+INCR*maxNumber+1]=LP1.Et()/CMS;
                    } else if (k1>0)  {
                          LParticle LPP3=(LParticle)bjets.at(k1-1);
                          TLorentzVector LP3=LPP3.GetP();
                          float imbalance=(LP3.Et() - LP1.Et())/(LP3.Et() + LP1.Et());
                          outMatrix[k1+INCR*maxNumber+1][k1+INCR*maxNumber+1]=imbalance;
                    }

                    //  non-diagonal 
                    for (unsigned int k2 = 0; k2<mJetsB; k2++) {
                     LParticle LPP2=(LParticle)bjets.at(k2);
                     TLorentzVector LP2=LPP2.GetP();
                     if (k1<k2)    outMatrix[k1+INCR*maxNumber+1][k2+INCR*maxNumber+1]=getMass(CMS,LP1,LP2);
                     if (k1>k2)    outMatrix[k1+INCR*maxNumber+1][k2+INCR*maxNumber+1]=getAngle(CMS,LP1,LP2);
                    }
                  }
                 }



                  INCR=2;
                  unsigned  int mMuons=maxNumber;
                  if (muons.size()<mMuons) mMuons=muons.size();
                  if (mMuons>0) {
                  for (unsigned int k1 = 0; k1<mMuons; k1++) {
                    LParticle LPP1=(LParticle)muons.at(k1);
                    TLorentzVector LP1=LPP1.GetP();
                    if (LMET.Et()>0) outMatrix[0][k1+INCR*maxNumber+1]=getMT(LMET, LP1)/CMS; 
                    outMatrix[k1+INCR*maxNumber+1][0]=getHL(LP1);

                    if (k1 == 0) {
                        outMatrix[k1+INCR*maxNumber+1][k1+INCR*maxNumber+1]=LP1.Et()/CMS;
                    } else if (k1>0)  {
                          LParticle LPP3=(LParticle)muons.at(k1-1);
                          TLorentzVector LP3=LPP3.GetP();
                          float imbalance=(LP3.Et() - LP1.Et())/(LP3.Et() + LP1.Et());
                          outMatrix[k1+INCR*maxNumber+1][k1+INCR*maxNumber+1]=imbalance;
                    }


                    for (unsigned int k2 = 0; k2<mMuons; k2++) {
                     LParticle LPP2=(LParticle)muons.at(k2);
                     TLorentzVector LP2=LPP2.GetP();
                     if (k1<k2)    outMatrix[k1+INCR*maxNumber+1][k2+INCR*maxNumber+1]=getMass(CMS,LP1,LP2);
                     if (k1>k2)    outMatrix[k1+INCR*maxNumber+1][k2+INCR*maxNumber+1]=getAngle(CMS,LP1,LP2);
                  }}
                 }


                 INCR=3;
                 unsigned  int mEle=maxNumber;
                  if (electrons.size()<mEle) mEle=electrons.size();
                  if (mEle>0) {
                  for (unsigned int k1 = 0; k1<mEle; k1++) {
                    LParticle LPP1=(LParticle)electrons.at(k1);
                    TLorentzVector LP1=LPP1.GetP();
                    if (LMET.Et()>0) outMatrix[0][k1+INCR*maxNumber+1]=getMT(LMET, LP1)/CMS; 
                    outMatrix[k1+INCR*maxNumber+1][0]=getHL(LP1);

                    if (k1 == 0) {
                        outMatrix[k1+INCR*maxNumber+1][k1+INCR*maxNumber+1]=LP1.Et()/CMS;
                    } else if (k1>0)  {
                          LParticle LPP3=(LParticle)electrons.at(k1-1);
                          TLorentzVector LP3=LPP3.GetP();
                          float imbalance=(LP3.Et() - LP1.Et())/(LP3.Et() + LP1.Et());
                          outMatrix[k1+INCR*maxNumber+1][k1+INCR*maxNumber+1]=imbalance;
                    }

                    // non-diagonal
                    for (unsigned int k2 = 0; k2<mEle; k2++) {
                     LParticle LPP2=(LParticle)electrons.at(k2);
                     TLorentzVector LP2=LPP2.GetP();
                     if (k1<k2)    outMatrix[k1+INCR*maxNumber+1][k2+INCR*maxNumber+1]=getMass(CMS,LP1,LP2);
                     if (k1>k2)    outMatrix[k1+INCR*maxNumber+1][k2+INCR*maxNumber+1]=getAngle(CMS,LP1,LP2);
                  }}
                 }


                 INCR=4;
                 unsigned  int mPho=maxNumber;
                  if (photons.size()<mPho) mPho=photons.size();
                  if (mPho>0) {
                  for (unsigned int k1 = 0; k1<mPho; k1++) {
                    LParticle LPP1=(LParticle)photons.at(k1);
                    TLorentzVector LP1=LPP1.GetP();
                    if (LMET.Et()>0) outMatrix[0][k1+INCR*maxNumber+1]=getMT(LMET, LP1)/CMS; 
                    outMatrix[k1+INCR*maxNumber+1][0]=getHL(LP1);

                    if (k1 == 0) {
                        outMatrix[k1+INCR*maxNumber+1][k1+INCR*maxNumber+1]=LP1.Et()/CMS;
                    } else if (k1>0)  {
                        LParticle LPP3=(LParticle)photons.at(k1-1);
                        TLorentzVector LP3=LPP3.GetP();
                        float imbalance=(LP3.Et() - LP1.Et())/(LP3.Et() + LP1.Et());
                        outMatrix[k1+INCR*maxNumber+1][k1+INCR*maxNumber+1]=imbalance;
                    }


                    // non-diagonal
                    for (unsigned int k2 = 0; k2<mPho; k2++) {
                     LParticle LPP2=(LParticle)photons.at(k2);
                     TLorentzVector LP2=LPP2.GetP();
                     if (k1<k2) outMatrix[k1+INCR*maxNumber+1][k2+INCR*maxNumber+1]=getMass(CMS,LP1,LP2); 
                     if (k1>k2) outMatrix[k1+INCR*maxNumber+1][k2+INCR*maxNumber+1]=getAngle(CMS,LP1,LP2); 
                  }}
                 }




                // diagonal ellements, such as jet-lep, etc. First, do jets:
                  unsigned  int INCR_SHIFT=0; // for jets =1 
                  for (unsigned int k1 = 0; k1<mJets; k1++) {
                   LParticle LPP1=(LParticle)jets.at(k1);
                   TLorentzVector LP1=LPP1.GetP();

                   INCR=1;
                   for (unsigned int k2 = 0; k2<mJetsB; k2++) {
                     LParticle LPP2=(LParticle)bjets.at(k2);
                     TLorentzVector LP2=LPP2.GetP();
                     TLorentzVector LPP=LP1+LP2;
                     outMatrix[k1+INCR_SHIFT*maxNumber+1][k2+INCR*maxNumber+1]=getMass(CMS,LP1,LP2);
                     outMatrix[k2+INCR*maxNumber+1][k1+INCR_SHIFT*maxNumber+1]=getAngle(CMS,LP1,LP2);
                  }

                   INCR=2;
                   for (unsigned int k2 = 0; k2<mMuons; k2++) {
                     LParticle LPP2=(LParticle)muons.at(k2);
                     TLorentzVector LP2=LPP2.GetP();
                     TLorentzVector LPP=LP1+LP2;
                     outMatrix[k1+INCR_SHIFT*maxNumber+1][k2+INCR*maxNumber+1]=getMass(CMS,LP1,LP2); 
                     outMatrix[k2+INCR*maxNumber+1][k1+INCR_SHIFT*maxNumber+1]=getAngle(CMS,LP1,LP2);  
                  }

                  INCR=3;  
                  for (unsigned int k2 = 0; k2<mEle; k2++) {
                     LParticle LPP2=(LParticle)electrons.at(k2);
                     TLorentzVector LP2=LPP2.GetP();
                     TLorentzVector LPP=LP1+LP2;
                     outMatrix[k1+INCR_SHIFT*maxNumber+1][k2+INCR*maxNumber+1]=getMass(CMS,LP1,LP2);  
                     outMatrix[k2+INCR*maxNumber+1][k1+INCR_SHIFT*maxNumber+1]=getAngle(CMS,LP1,LP2); 
                  }

                  INCR=4;
                  for (unsigned int k2 = 0; k2<mPho; k2++) {
                     LParticle LPP2=(LParticle)photons.at(k2);
                     TLorentzVector LP2=LPP2.GetP();
                     TLorentzVector LPP=LP1+LP2;
                     outMatrix[k1+INCR_SHIFT*maxNumber+1][k2+INCR*maxNumber+1]=getMass(CMS,LP1,LP2); 
                     outMatrix[k2+INCR*maxNumber+1][k1+INCR_SHIFT*maxNumber+1]=getAngle(CMS,LP1,LP2); 
                  }

                }


                // b-jets
                // diagonal ellements, such as jet-lep, etc. First, do jets:
                  INCR_SHIFT=1; // for bjets =2 
                  for (unsigned int k1 = 0; k1<mJetsB; k1++) {
                   LParticle LPP1=(LParticle)bjets.at(k1);
                   TLorentzVector LP1=LPP1.GetP();

                   INCR=2;
                   for (unsigned int k2 = 0; k2<mMuons; k2++) {
                     LParticle LPP2=(LParticle)muons.at(k2);
                     TLorentzVector LP2=LPP2.GetP();
                     TLorentzVector LPP=LP1+LP2;
                     outMatrix[k1+INCR_SHIFT*maxNumber+1][k2+INCR*maxNumber+1]=getMass(CMS,LP1,LP2);
                     outMatrix[k2+INCR*maxNumber+1][k1+INCR_SHIFT*maxNumber+1]=getAngle(CMS,LP1,LP2);
                  }

                  INCR=3;
                  for (unsigned int k2 = 0; k2<mEle; k2++) {
                     LParticle LPP2=(LParticle)electrons.at(k2);
                     TLorentzVector LP2=LPP2.GetP();
                     TLorentzVector LPP=LP1+LP2;
                     outMatrix[k1+INCR_SHIFT*maxNumber+1][k2+INCR*maxNumber+1]=getMass(CMS,LP1,LP2);
                     outMatrix[k2+INCR*maxNumber+1][k1+INCR_SHIFT*maxNumber+1]=getAngle(CMS,LP1,LP2);
                  }

                  INCR=4;
                  for (unsigned int k2 = 0; k2<mPho; k2++) {
                     LParticle LPP2=(LParticle)photons.at(k2);
                     TLorentzVector LP2=LPP2.GetP();
                     TLorentzVector LPP=LP1+LP2;
                     outMatrix[k1+INCR_SHIFT*maxNumber+1][k2+INCR*maxNumber+1]=getMass(CMS,LP1,LP2);
                     outMatrix[k2+INCR*maxNumber+1][k1+INCR_SHIFT*maxNumber+1]=getAngle(CMS,LP1,LP2);
                  }

                }





                  // muons
                  INCR_SHIFT=2; // for muons =3 
                  for (unsigned int k1 = 0; k1<mMuons; k1++) {
                   LParticle LPP1=(LParticle)muons.at(k1);
                   TLorentzVector LP1=LPP1.GetP();

                  INCR=3;
                  for (unsigned int k2 = 0; k2<mEle; k2++) {
                     LParticle LPP2=(LParticle)electrons.at(k2);
                     TLorentzVector LP2=LPP2.GetP();
                     TLorentzVector LPP=LP1+LP2;
                     outMatrix[k1+INCR_SHIFT*maxNumber+1][k2+INCR*maxNumber+1]=getMass(CMS,LP1,LP2); 
                     outMatrix[k2+INCR*maxNumber+1][k1+INCR_SHIFT*maxNumber+1]=getAngle(CMS,LP1,LP2); 
                  }

                  INCR=4;
                  for (unsigned int k2 = 0; k2<mPho; k2++) {
                     LParticle LPP2=(LParticle)photons.at(k2);
                     TLorentzVector LP2=LPP2.GetP();
                     TLorentzVector LPP=LP1+LP2;
                     outMatrix[k1+INCR_SHIFT*maxNumber+1][k2+INCR*maxNumber+1]=getMass(CMS,LP1,LP2); 
                     outMatrix[k2+INCR*maxNumber+1][k1+INCR_SHIFT*maxNumber+1]=getAngle(CMS,LP1,LP2); 
                  }

                }


                   // electrons 
                 INCR_SHIFT=3; // for electrons =3 
                 for (unsigned int k1 = 0; k1<mEle; k1++) {
                  LParticle LPP1=(LParticle)electrons.at(k1);
                  TLorentzVector LP1=LPP1.GetP();
                  INCR=4;
                  for (unsigned int k2 = 0; k2<mPho; k2++) {
                     LParticle LPP2=(LParticle)photons.at(k2);
                     TLorentzVector LP2=LPP2.GetP();
                     TLorentzVector LPP=LP1+LP2;
                     outMatrix[k1+INCR_SHIFT*maxNumber+1][k2+INCR*maxNumber+1]=getMass(CMS,LP1,LP2); 
                     outMatrix[k2+INCR*maxNumber+1][k1+INCR_SHIFT*maxNumber+1]=getAngle(CMS,LP1,LP2); 
                  }

                }


  // here we convert RMM to RMM-C46
  // The standard RMM is given by th 2D matrix: outMatrix
  float* outMatrix_c46 = new float[46];
  outMatrix_c46[0] = outMatrix[0][0]; // MET 
  // ET components 
  outMatrix_c46[1]= getFrobeniusAggregationET(0, mJets,  maxNumber,  outMatrix);
  outMatrix_c46[2]= getFrobeniusAggregationET(1, mJetsB, maxNumber,  outMatrix);
  outMatrix_c46[3]= getFrobeniusAggregationET(2, mMuons, maxNumber,  outMatrix);
  outMatrix_c46[4]= getFrobeniusAggregationET(3, mEle,   maxNumber,  outMatrix);
  outMatrix_c46[5]= getFrobeniusAggregationET(4, mPho,   maxNumber,  outMatrix);
  // transverse masses
  outMatrix_c46[6]= getFrobeniusAggregationMET(0, mJets,  maxNumber,  outMatrix);
  outMatrix_c46[7]= getFrobeniusAggregationMET(1, mJetsB, maxNumber,  outMatrix);
  outMatrix_c46[8]= getFrobeniusAggregationMET(2, mMuons, maxNumber,  outMatrix);
  outMatrix_c46[9]= getFrobeniusAggregationMET(3, mEle,   maxNumber,  outMatrix);
  outMatrix_c46[10]= getFrobeniusAggregationMET(4, mPho,   maxNumber,  outMatrix);
  // longit. 
  outMatrix_c46[11]= getFrobeniusAggregationLong(0, mJets,  maxNumber,  outMatrix);
  outMatrix_c46[12]= getFrobeniusAggregationLong(1, mJetsB, maxNumber,  outMatrix);
  outMatrix_c46[13]= getFrobeniusAggregationLong(2, mMuons, maxNumber,  outMatrix);
  outMatrix_c46[14]= getFrobeniusAggregationLong(3, mEle,   maxNumber,  outMatrix);
  outMatrix_c46[15]= getFrobeniusAggregationLong(4, mPho,   maxNumber,  outMatrix);


  // rapidity-difference zones (15 total)

 for (int p = 0; p < 15; ++p) {
    int ti = pairs[p][0];
    int tj = pairs[p][1];
    unsigned int Ni =
        (ti==0 ? mJets  :
         ti==1 ? mJetsB :
         ti==2 ? mMuons :
         ti==3 ? mEle   :
                 mPho);

    unsigned int Nj =

        (tj==0 ? mJets  :
         tj==1 ? mJetsB :
         tj==2 ? mMuons :
         tj==3 ? mEle   :
                 mPho);

    outMatrix_c46[16 + p] =
        getFrobeniusAggregationRapidityDiff(
            tj, ti,   // NOTE transpose for h-terms
            Nj, Ni,
            maxNumber,
            outMatrix
        );
   }




// invariant-mass zones (15 total)
  for (int p = 0; p < 15; ++p) {
    int ti = pairs[p][0];
    int tj = pairs[p][1];
    unsigned int Ni =
        (ti==0 ? mJets  :
         ti==1 ? mJetsB :
         ti==2 ? mMuons :
         ti==3 ? mEle   :
                 mPho);
    unsigned int Nj =
        (tj==0 ? mJets  :
         tj==1 ? mJetsB :
         tj==2 ? mMuons :
         tj==3 ? mEle   :
                 mPho);


    outMatrix_c46[31 + p] =
        getFrobeniusAggregationMass(
            ti, tj,
            Ni, Nj,
            maxNumber,
            outMatrix
        );

    }



/*
  // rapidity-difference zones (15 in total) 
  // jets
  outMatrix_c46[16]= getFrobeniusAggregationRapidityDiff(0,1, mJets,  mJets,  maxNumber,  outMatrix);
  outMatrix_c46[17]= getFrobeniusAggregationRapidityDiff(0,2, mJets,  mJetsB,  maxNumber,  outMatrix);
  outMatrix_c46[18]= getFrobeniusAggregationRapidityDiff(0,3, mJets,  mMuons,  maxNumber,  outMatrix); 
  outMatrix_c46[19]= getFrobeniusAggregationRapidityDiff(0,4, mJets,  mEle,  maxNumber,  outMatrix);
  outMatrix_c46[20]= getFrobeniusAggregationRapidityDiff(0,5, mJets,  mPho,  maxNumber,  outMatrix);
  // b-jets 
  outMatrix_c46[21]= getFrobeniusAggregationRapidityDiff(1,1, mJetsB, mJetsB,  maxNumber,  outMatrix); 
  outMatrix_c46[22]= getFrobeniusAggregationRapidityDiff(1,2, mJetsB, mMuons,  maxNumber,  outMatrix);
  outMatrix_c46[23]= getFrobeniusAggregationRapidityDiff(1,3, mJetsB, mEle,  maxNumber,  outMatrix);
  outMatrix_c46[24]= getFrobeniusAggregationRapidityDiff(1,4, mJetsB, mPho,  maxNumber,  outMatrix);
  // muons
  outMatrix_c46[25]= getFrobeniusAggregationRapidityDiff(2,2, mMuons, mMuons,  maxNumber,  outMatrix);
  outMatrix_c46[26]= getFrobeniusAggregationRapidityDiff(2,3, mMuons, mEle,  maxNumber,  outMatrix);
  outMatrix_c46[27]= getFrobeniusAggregationRapidityDiff(2,4, mMuons, mPho,  maxNumber,  outMatrix);
  // electrons 
  outMatrix_c46[28]= getFrobeniusAggregationRapidityDiff(3,3, mEle, mEle,  maxNumber,  outMatrix);
  outMatrix_c46[29]= getFrobeniusAggregationRapidityDiff(3,4, mEle, mPho,  maxNumber,  outMatrix);
  // photons
  outMatrix_c46[30]= getFrobeniusAggregationRapidityDiff(4,4, mPho, mPho,  maxNumber,  outMatrix);

  // invariant masses
  // jets
  cout << " " << endl;
  cout << "Matrix=31" << endl;
  outMatrix_c46[31]= getFrobeniusAggregationMass(0,0, mJets,  mJets,  maxNumber,  outMatrix);
  cout << "Matrix=32" << endl;
  outMatrix_c46[32]= getFrobeniusAggregationMass(0,1, mJets,  mJetsB,  maxNumber,  outMatrix);

  outMatrix_c46[33]= getFrobeniusAggregationMass(0,2, mJets,  mMuons,  maxNumber,  outMatrix);
  outMatrix_c46[34]= getFrobeniusAggregationMass(0,3, mJets,  mEle,  maxNumber,  outMatrix);
  outMatrix_c46[35]= getFrobeniusAggregationMass(0,4, mJets,  mPho,  maxNumber,  outMatrix);
  // b-jets 
  outMatrix_c46[36]= getFrobeniusAggregationMass(1,1, mJetsB, mJetsB,  maxNumber,  outMatrix);
  outMatrix_c46[37]= getFrobeniusAggregationMass(1,2, mJetsB, mMuons,  maxNumber,  outMatrix);
  outMatrix_c46[38]= getFrobeniusAggregationMass(1,3, mJetsB, mEle,  maxNumber,  outMatrix);
  outMatrix_c46[39]= getFrobeniusAggregationMass(1,4, mJetsB, mPho,  maxNumber,  outMatrix);
  // muons
  outMatrix_c46[40]= getFrobeniusAggregationMass(2,2, mMuons, mMuons,  maxNumber,  outMatrix);
  outMatrix_c46[41]= getFrobeniusAggregationMass(2,3, mMuons, mEle,  maxNumber,  outMatrix);
  outMatrix_c46[42]= getFrobeniusAggregationMass(2,4, mMuons, mPho,  maxNumber,  outMatrix);
  // electrons 
  outMatrix_c46[43]= getFrobeniusAggregationMass(3,3, mEle, mEle,  maxNumber,  outMatrix);
  outMatrix_c46[44]= getFrobeniusAggregationMass(3,4, mEle, mPho,  maxNumber,  outMatrix);
  // photons
  outMatrix_c46[45]= getFrobeniusAggregationMass(4,4, mPho, mPho,  maxNumber,  outMatrix);
  */


   cout << "\n Event "<< endl;
   for (int p = 0; p < 15; ++p) {
        int id=31+p;
        cout << "outMatrix_c46["<< id << "]=" << outMatrix_c46[id] << endl;
   }



  return outMatrix_c46; 


} 
