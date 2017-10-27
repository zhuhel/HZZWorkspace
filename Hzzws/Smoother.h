//    Description:  Smooth the distributions from the mini-tree
// 
#ifndef __Smoother_H__
#define __Smoother_H__
#include <string>

#include <RooArgSet.h>
#include <TH1.h>
#include <TTree.h>
#include <TChain.h>
#include "TVectorD.h"
#include "RooKeysPdf.h"
#include "RooNDKeysPdf.h"

using namespace std;

class Smoother{
   public:
        Smoother(const string& outname, std::vector<float> rho, std::string mirror);
        ~Smoother();

        void setMakeValidPlots(bool b){m_makeValidPlots=b;}

        void smooth(const string& input_name, 
                const string& oname, 
                const string& treename, 
                const RooArgSet& treeobs, 
                const string& cuts) const ;

        void setNickname(const string& nick) { m_nick = nick; }

   private:
        std::vector<float> m_rho;
        std::string m_mirror;
        std::string m_nick;
        TFile* outfile;
        bool m_makeValidPlots;

        void ReadMirror(RooKeysPdf::Mirror&) const;
};

#endif
