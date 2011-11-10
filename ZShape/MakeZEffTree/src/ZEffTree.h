#include "TFile.h"
#include "TBranch.h"
#include "TTree.h"
#include <map>
#include <string>

#ifndef ZEffTree_h_included
#define ZEffTree_h_included

class ZEffTree {
    public:

        ZEffTree(TFile& f, bool writable);

        static int cutToBit(const std::string& cut);

        struct ZInfo {
            float eta[2];
            float phi[2];
            float pt[2];
            float mz, yz, qtz;
            int bits[2];
            int nverts;
            bool isSelected(int ielec, const std::string& bitname) const;
            void setBit(int ielec, const std::string& bitname, bool val);
        } gen, reco;

        void Fill() { 
            m_tree->Fill();
        }

        void Write() { 
            m_file.Write();
        }

        int Entries() { 
            return m_tree->GetEntries(); 
        }

        bool GetNextEvent() {
            if (nevt==Entries()) {
                return false; 
            }
            m_tree->GetEntry(nevt); 
            nevt++; 
            return true; 
        }

        void Clear() { 
            gen.eta[0] = 0; gen.eta[1] = 0; gen.phi[0] = 0; gen.phi[1] = 0; gen.pt[0] = 0; gen.pt[1] = 0; gen.mz = 0; gen.yz = 0; gen.qtz = 0; gen.bits[0] = 0; gen.bits[1] = 0; gen.nverts = 0;
            reco.eta[0] = 0; reco.eta[1] = 0; reco.phi[0] = 0; reco.phi[1] = 0; reco.pt[0] = 0; reco.pt[1] = 0; reco.mz = 0; reco.yz = 0; reco.qtz = 0; reco.bits[0] = 0; reco.bits[1] = 0; reco.nverts = 0;
        }

    private:
        TFile& m_file;
        TTree* m_tree;
        int nevt;
        TBranch* br_gen, *br_reco;

        void makeBranches(bool writable) {
            nevt=0;
            if (writable) {
                m_file.cd();
                m_tree=new TTree("ZEffs","Minnesota ZEffs");
                br_gen=m_tree->Branch("gen",&gen,"eta0/f:eta1:phi0:phi1:pt0:pt1:mz:yz:qtz:bits0/I:bits1:nverts");
                br_reco=m_tree->Branch("reco",&reco,"eta0/f:eta1:phi0:phi1:pt0:pt1:mz:yz:qtz:bits0/I:bits1:nverts");
            } else {
                m_tree=(TTree*)m_file.Get("ZEffs");
                m_tree->SetBranchAddress("gen",&gen);
                m_tree->SetBranchAddress("reco",&reco);
            }
        }

        /* Bit packing code*/
        /* 
        int m_bitbox;

        void clearBits() {
            m_bitbox = 0;
        }

        void setBit(int bitnum, bool val) {
            if (bitnum>=0 && bitnum<31) {
                int mask = 1 << bitnum;
                m_bitbox = val ? (m_bitbox | mask) : (m_bitbox & ~mask);    
            }
        };

        bool readBit(int bitnum) const {
            if (bitnum<0 || bitnum>=31) return false;
            int mask = 1 << bitnum;
            return m_bitbox & mask;
        };

        void flipBit(int bitnum) {
            if (bitnum>=0 && bitnum<31) {
                int mask = 1 << bitnum;
                m_bitbox ^= mask;
            }
        };

        void clearBit(int bitnum) {
            setBit(bitnum, false);
        };
        */

        /* Map from names to bitnums  */
        static std::map<std::string, int> cutToBits_;

        void prepBitmap();

};

#endif
