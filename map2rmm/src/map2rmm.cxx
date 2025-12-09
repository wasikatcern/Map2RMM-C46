/**************************************************
 * Project event into RMM and compress to C46
 **************************************************/

#include <vector>
#include <cmath>
#include <algorithm>
#include <iostream>

#include "TMath.h"
#include "TLorentzVector.h"
#include "SystemOfUnits.h"
#include "LParticle.h"
#include "CParticle.h"

using namespace std;

/*
  RMM layout (conceptual):

  Index 0: MET row/column.
  Then maxNumberTypes blocks of size maxN for each object type.
  For this implementation, the types are in the order:

    t = 0 : jets
    t = 1 : muons
    t = 2 : electrons
    t = 3 : photons
    (if maxNumberTypes < 4, we use fewer)

  The RMM is filled as follows:

  [0,0]      : MET_ET / CMS
  [0, obj]   : mT(met,obj) / CMS
  [obj, 0]   : hL(obj) = cosh(rapidity) - 1
  diagonal   : first object of each type => ET/CMS;
               others => ET-imbalance with respect to previous object
  upper-right blocks (including upper triangle of same type) : invariant mass / CMS
  lower-left blocks (including lower triangle of same type)  : rapidity-based angle (HL)

  This reproduces the standard RMM pattern: masses in upper-right,
  rapidity differences (or cosh(y)-1 variants) in lower-left.
*/

// ----------------------------------------------------------------------
// Scalar helpers
// ----------------------------------------------------------------------

// angle-like quantity between two 4-vectors (rapidity-based)
float getAngle(const float CMS, const TLorentzVector& p1, const TLorentzVector& p2) {
    (void)CMS; // currently unused, but kept for compatibility
    double y1 = p1.Rapidity();
    double y2 = p2.Rapidity();
    double HL = TMath::CosH(0.5 * (y2 - y1)) - 1.0;
    return static_cast<float>(HL);
}

// invariant mass / CMS
float getMass(const float CMS, const TLorentzVector& p1, const TLorentzVector& p2) {
    TLorentzVector pp = p1 + p2;
    float xmass = pp.M() / CMS;
    return xmass;
}

// cosh(y) - 1 (single-particle rapidity measure)
float getHL(const TLorentzVector& p1) {
    double y  = p1.Rapidity();
    double HL = TMath::CosH(y) - 1.0;
    return static_cast<float>(HL);
}

// transverse mass (exact, standard experimental definition)
float getMT(const TLorentzVector& met, const TLorentzVector& jet) {
    double et_sum = jet.Et() + met.Et();
    double px_sum = jet.Px() + met.Px();
    double py_sum = jet.Py() + met.Py();

    double ss = et_sum * et_sum - px_sum * px_sum - py_sum * py_sum;
    double Mt_exact = (ss > 0.0) ? TMath::Sqrt(ss) : 0.0;
    return static_cast<float>(Mt_exact);
}

// ----------------------------------------------------------------------
// Build full RMM matrix (float**), same interface as original Map2RMM
// ----------------------------------------------------------------------
//
//  CMS            : center-of-mass energy (GeV)
//  maxN           : max multiplicity per object type
//  maxNumberTypes : number of types (<=4 for this implementation)
//  missing        : vector with MET-like objects (we use missing[0] if present)
//  jets, muons,
//  electrons,
//  photons        : physics objects for each type
//
//  Returns: float** of size (maxNumberTypes*maxN+1)^2, all allocated with new[]
//           Caller must delete[] each row + delete[] the top pointer.
//

float** map2rmm(const float  CMS,
                const int    maxN,
                const int    maxNumberTypes,
                const vector<LParticle> missing,
                const vector<LParticle> jets,
                const vector<LParticle> muons,
                const vector<LParticle> electrons,
                const vector<LParticle> photons)
{
    const int maxNumber = maxN;
    const int maxTypes  = maxNumberTypes;
    const int maxSize   = maxNumber * maxTypes + 1;  // + MET row/col

    const int height = maxSize;
    const int width  = maxSize;

    // allocate and zero
    float** outMatrix = new float*[height];
    for (int i = 0; i < height; ++i) {
        outMatrix[i] = new float[width];
        for (int j = 0; j < width; ++j) {
            outMatrix[i][j] = 0.0f;
        }
    }

    // Prepare MET four-vector (if available)
    TLorentzVector LMET;
    if (!missing.empty()) {
        LMET = missing[0].GetP();
    } else {
        LMET.SetPxPyPzE(0, 0, 0, 0);
    }

    // MET cell
    outMatrix[0][0] = static_cast<float>(LMET.Et() / CMS);

    // Collect objects per type, in the order:
    // t = 0: jets
    // t = 1: muons
    // t = 2: electrons
    // t = 3: photons
    vector< vector<TLorentzVector> > typeObjs(4);

    auto fillVec = [&](const vector<LParticle>& src, vector<TLorentzVector>& dst) {
        dst.clear();
        int n = std::min<int>(maxNumber, src.size());
        dst.reserve(n);
        for (int i = 0; i < n; ++i) {
            TLorentzVector v = src[i].GetP();
            dst.push_back(v);
        }
    };

    fillVec(jets,      typeObjs[0]);
    fillVec(muons,     typeObjs[1]);
    fillVec(electrons, typeObjs[2]);
    fillVec(photons,   typeObjs[3]);

    // How many types do we actually use?
    const int nTypesUsed = std::min<int>(maxTypes, 4);

    // Helper to compute the starting index for a given type
    auto typeOffset = [&](int t) {
        return 1 + t * maxNumber;
    };

    // 1) Fill diagonal cells, row/col with MET (mT, hL)
    for (int t = 0; t < nTypesUsed; ++t) {
        const auto& objs = typeObjs[t];
        const int nObj   = static_cast<int>(objs.size());
        const int base   = typeOffset(t);

        if (nObj == 0) continue;

        for (int k = 0; k < nObj; ++k) {
            const TLorentzVector& obj = objs[k];
            const int row = base + k;
            const int col = row;

            // [0,row] : mT(met,obj)/CMS
            outMatrix[0][row] = getMT(LMET, obj) / CMS;

            // [row,0] : hL(obj)
            outMatrix[row][0] = getHL(obj);

            // diagonal: leading object => ET/CMS, others => ET imbalance
            if (k == 0) {
                outMatrix[row][col] = static_cast<float>(obj.Et() / CMS);
            } else {
                const TLorentzVector& prev = objs[k - 1];
                double num = prev.Et() - obj.Et();
                double den = prev.Et() + obj.Et();
                float imbalance = (den != 0.0) ? static_cast<float>(num / den) : 0.0f;
                outMatrix[row][col] = imbalance;
            }
        }
    }

    // 2) Fill pairwise blocks:
    //    - same type: upper triangle = mass, lower triangle = angle
    //    - different types:
    //        t1 < t2 -> upper-right block = mass
    //        t1 > t2 -> lower-left block = angle
    for (int t1 = 0; t1 < nTypesUsed; ++t1) {
        const auto& objs1 = typeObjs[t1];
        const int n1 = static_cast<int>(objs1.size());
        if (n1 == 0) continue;

        const int base1 = typeOffset(t1);

        for (int t2 = 0; t2 < nTypesUsed; ++t2) {
            const auto& objs2 = typeObjs[t2];
            const int n2 = static_cast<int>(objs2.size());
            if (n2 == 0) continue;

            const int base2 = typeOffset(t2);

            for (int i = 0; i < n1; ++i) {
                const TLorentzVector& p1 = objs1[i];
                const int r1 = base1 + i;

                for (int j = 0; j < n2; ++j) {
                    const TLorentzVector& p2 = objs2[j];
                    const int r2 = base2 + j;

                    if (t1 == t2) {
                        if (i == j) continue; // diagonal already set
                        if (i < j) {
                            // upper triangle: mass / CMS
                            outMatrix[r1][r2] = getMass(CMS, p1, p2);
                        } else {
                            // lower triangle: angle
                            outMatrix[r1][r2] = getAngle(CMS, p1, p2);
                        }
                    } else if (t1 < t2) {
                        // upper-right block: mass
                        outMatrix[r1][r2] = getMass(CMS, p1, p2);
                    } else {
                        // lower-left block: angle
                        outMatrix[r1][r2] = getAngle(CMS, p1, p2);
                    }
                }
            }
        }
    }

    return outMatrix;
}

// =====================================================================
//  RMM → RMM-C46 compression (C46 features)
// =====================================================================
//
// This implements a 46-dimensional summary of the RMM matrix:
//   1  MET term
//   5  ET terms     (one per type, padded up to 5)
//   5  mT (T) terms (first row segments, per type, padded up to 5)
//   5  hL (L) terms (first column segments, per type, padded up to 5)
//   15 h-terms      (lower blocks & triangles)
//   15 m-terms      (upper blocks & triangles)
//
// We support up to 5 logical types. For Map2RMM (4 types), any
// "5-th type" features are set to zero automatically.
//

static inline float agg_add_cpp(const std::vector<float>& vals) {
    double s = 0.0;
    for (double v : vals) s += v;
    return static_cast<float>(s);
}

static inline float agg_frob_cpp(const std::vector<float>& vals) {
    double s2 = 0.0;
    for (double v : vals) s2 += v * v;
    return static_cast<float>(std::sqrt(s2));
}

/// Build 46 C46 features from a full RMM matrix A.
/// A has size mSize x mSize, with:
///   mSize = maxNumberTypes * maxNumber + 1.
/// The result is written into `out` (size 46).
void buildC46FromRMM(float** A,
                     int maxNumber,
                     int maxNumberTypes,
                     std::vector<float>& out,
                     bool useFrob = true)
{
    out.assign(46, 0.0f);

    const int mSize = maxNumberTypes * maxNumber + 1;
    if (mSize <= 1) {
        std::cerr << "[buildC46FromRMM] matrix too small, mSize=" << mSize << std::endl;
        return;
    }

    auto agg = [&](const std::vector<float>& vals) {
        return useFrob ? agg_frob_cpp(vals) : agg_add_cpp(vals);
    };

    // We support up to 5 logical types.
    const int nTypes = std::min(maxNumberTypes, 5);

    struct Slice { int begin; int end; };
    std::vector<Slice> slices;
    slices.reserve(nTypes);
    for (int t = 0; t < nTypes; ++t) {
        Slice s;
        s.begin = 1 + t * maxNumber;
        s.end   = s.begin + maxNumber; // exclusive
        slices.push_back(s);
    }

    auto getSlice = [&](int tIndex) -> Slice {
        if (tIndex < nTypes) return slices[tIndex];
        // dummy empty slice
        return Slice{1, 1};
    };

    int idx = 0;

    // -----------------------------------------------------
    // 1) MET term
    // -----------------------------------------------------
    out[idx++] = static_cast<float>(A[0][0]);

    // -----------------------------------------------------
    // 2) 5 ET terms (diagonal of same-type blocks)
    // -----------------------------------------------------
    for (int t = 0; t < nTypes; ++t) {
        const Slice& s = slices[t];
        std::vector<float> vals;
        vals.reserve(maxNumber);
        for (int i = s.begin; i < s.end; ++i) {
            vals.push_back(static_cast<float>(A[i][i]));
        }
        out[idx++] = agg(vals);
    }
    // pad up to 5
    for (int t = nTypes; t < 5; ++t) {
        out[idx++] = 0.0f;
    }

    // -----------------------------------------------------
    // 3) 5 T terms (first row segments: mT blocks)
    // -----------------------------------------------------
    for (int t = 0; t < nTypes; ++t) {
        const Slice& s = slices[t];
        std::vector<float> vals;
        vals.reserve(maxNumber);
        for (int j = s.begin; j < s.end; ++j) {
            vals.push_back(static_cast<float>(A[0][j]));
        }
        out[idx++] = agg(vals);
    }
    for (int t = nTypes; t < 5; ++t) {
        out[idx++] = 0.0f;
    }

    // -----------------------------------------------------
    // 4) 5 L terms (first column segments: hL blocks)
    // -----------------------------------------------------
    for (int t = 0; t < nTypes; ++t) {
        const Slice& s = slices[t];
        std::vector<float> vals;
        vals.reserve(maxNumber);
        for (int i = s.begin; i < s.end; ++i) {
            vals.push_back(static_cast<float>(A[i][0]));
        }
        out[idx++] = agg(vals);
    }
    for (int t = nTypes; t < 5; ++t) {
        out[idx++] = 0.0f;
    }

    // -----------------------------------------------------
    // 5) 15 h-terms & 15 m-terms
    //
    // Order of (ti, tj) pairs (logical type indices 0..4):
    //  (0,0),
    //  (1,0), (1,1),
    //  (2,0), (2,1), (2,2),
    //  (3,0), (3,1), (3,2), (3,3),
    //  (4,0), (4,1), (4,2), (4,3), (4,4)
    //
    // h-terms:
    //  - same-type: strictly lower triangle of block
    //  - cross-type: lower-left block (tj rows, ti cols)
    //
    // m-terms:
    //  - same-type: strictly upper triangle of block
    //  - cross-type: upper-right block (ti rows, tj cols)
    // -----------------------------------------------------
    const int orderPairs[15][2] = {
        {0,0},
        {1,0},{1,1},
        {2,0},{2,1},{2,2},
        {3,0},{3,1},{3,2},{3,3},
        {4,0},{4,1},{4,2},{4,3},{4,4}
    };

    // --- h-terms ---
    for (int k = 0; k < 15; ++k) {
        int ti = orderPairs[k][0];
        int tj = orderPairs[k][1];

        if (ti >= nTypes || tj >= nTypes) {
            out[idx++] = 0.0f;
            continue;
        }

        Slice si = getSlice(ti);
        Slice sj = getSlice(tj);
        std::vector<float> vals;

        if (ti == tj) {
            // same-type: strictly lower triangle
            for (int i = si.begin; i < si.end; ++i) {
                for (int j = sj.begin; j < sj.end; ++j) {
                    if (i > j) vals.push_back(static_cast<float>(A[i][j]));
                }
            }
        } else {
            // cross-type: lower-left block (tj rows, ti cols)
            for (int i = sj.begin; i < sj.end; ++i) {
                for (int j = si.begin; j < si.end; ++j) {
                    vals.push_back(static_cast<float>(A[i][j]));
                }
            }
        }
        out[idx++] = agg(vals);
    }

    // --- m-terms ---
    for (int k = 0; k < 15; ++k) {
        int ti = orderPairs[k][0];
        int tj = orderPairs[k][1];

        if (ti >= nTypes || tj >= nTypes) {
            out[idx++] = 0.0f;
            continue;
        }

        Slice si = getSlice(ti);
        Slice sj = getSlice(tj);
        std::vector<float> vals;

        if (ti == tj) {
            // same-type: strictly upper triangle
            for (int i = si.begin; i < si.end; ++i) {
                for (int j = sj.begin; j < sj.end; ++j) {
                    if (i < j) vals.push_back(static_cast<float>(A[i][j]));
                }
            }
        } else {
            // cross-type: upper-right block (ti rows, tj cols)
            for (int i = si.begin; i < si.end; ++i) {
                for (int j = sj.begin; j < sj.end; ++j) {
                    vals.push_back(static_cast<float>(A[i][j]));
                }
            }
        }
        out[idx++] = agg(vals);
    }

    if (idx != 46) {
        std::cerr << "[buildC46FromRMM] internal error: filled "
                  << idx << " features instead of 46\n";
    }
}

// ----------------------------------------------------------------------
// Convenience wrapper: Map2RMM46
//
// This is what you should call from example.cc in your new
// Map2RMM46 code. It:
//
//   1) builds the full RMM with map2rmm(...),
//   2) compresses it to 46 features with buildC46FromRMM(...),
//   3) deletes the RMM matrix,
//   4) returns std::vector<float> of size 46.
// ----------------------------------------------------------------------

std::vector<float> map2rmm46(const float CMS,
                             const int   maxN,
                             const int   maxNumberTypes,
                             const std::vector<LParticle>& missing,
                             const std::vector<LParticle>& jets,
                             const std::vector<LParticle>& muons,
                             const std::vector<LParticle>& electrons,
                             const std::vector<LParticle>& photons,
                             bool useFrob /* = true */)
{
    // 1) Build full RMM
    float** projArray = map2rmm(CMS,
                                maxN,
                                maxNumberTypes,
                                missing,
                                jets,
                                muons,
                                electrons,
                                photons);

    // 2) Compress to C46
    std::vector<float> c46;
    buildC46FromRMM(projArray, maxN, maxNumberTypes, c46, useFrob);

    // 3) Free RMM
    const int mSize = maxNumberTypes * maxN + 1;
    for (int i = 0; i < mSize; ++i) {
        delete [] projArray[i];
    }
    delete [] projArray;

    // 4) Return C46
    return c46;
}

