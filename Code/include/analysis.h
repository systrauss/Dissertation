//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Jan 11 12:37:54 2018 by ROOT version 5.34/19
// from TTree evtTree/133Ba_90deg_facingdown_trig(coin1(siliall+ge1+ge2) + pulser)
// found on file: rootFiles/Data/run1.root
//////////////////////////////////////////////////////////

#ifndef analysis_h
#define analysis_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include <TObject.h>

// Fixed size dimensions of array or collections stored in the TTree if any.

class analysis {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
 //eventData       *event;
   UInt_t          fUniqueID;
   UInt_t          fBits;
   UShort_t        sili_time[6];
   UShort_t        pulser_time;
   UShort_t        ge_time[2];
   UShort_t        n1_time;
   UShort_t        n2_time;
   UShort_t        bgo6_time;
   UShort_t        bgo1_time;
   UShort_t        bgo5_time;
   UShort_t        bgo4_time;
   UShort_t        buncher_time;
   Float_t         sili_en[6];
   UShort_t        sili_en_raw[6];
   UShort_t        ge_en_raw[2];
   Float_t         ge_en[2];
   UShort_t        n1_en[2];
   UShort_t        n2_en[2];
   UShort_t        pulser_en;
   UShort_t        veto_timing_en;
   UShort_t        bgo6_en;
   UShort_t        bgo5_en;
   UShort_t        bgo1_en;
   UShort_t        bgo4_en;

   // List of branches
   TBranch        *b_event_fUniqueID;   //!
   TBranch        *b_event_fBits;   //!
   TBranch        *b_event_sili_time;   //!
   TBranch        *b_event_pulser_time;   //!
   TBranch        *b_event_ge_time;   //!
   TBranch        *b_event_n1_time;   //!
   TBranch        *b_event_n2_time;   //!
   TBranch        *b_event_bgo6_time;   //!
   TBranch        *b_event_bgo1_time;   //!
   TBranch        *b_event_bgo5_time;   //!
   TBranch        *b_event_bgo4_time;   //!
   TBranch        *b_event_buncher_time;   //!
   TBranch        *b_event_sili_en;   //!
   TBranch        *b_event_sili_en_raw;   //!
   TBranch        *b_event_ge_en_raw;   //!
   TBranch        *b_event_ge_en;   //!
   TBranch        *b_event_n1_en;   //!
   TBranch        *b_event_n2_en;   //!
   TBranch        *b_event_pulser_en;   //!
   TBranch        *b_event_veto_timing_en;   //!
   TBranch        *b_event_bgo6_en;   //!
   TBranch        *b_event_bgo5_en;   //!
   TBranch        *b_event_bgo1_en;   //!
   TBranch        *b_event_bgo4_en;   //!

   analysis(TTree *tree=0);
   virtual ~analysis();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop(const char* fileName, int nRunNum, bool bTripleCoin);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef analysis_cxx
analysis::analysis(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("rootFiles/Data/run1.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("rootFiles/Data/run1.root");
      }
      f->GetObject("evtTree",tree);

   }
   Init(tree);
}

analysis::~analysis()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t analysis::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t analysis::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void analysis::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("fUniqueID", &fUniqueID, &b_event_fUniqueID);
   fChain->SetBranchAddress("fBits", &fBits, &b_event_fBits);
   fChain->SetBranchAddress("sili_time[6]", sili_time, &b_event_sili_time);
   fChain->SetBranchAddress("pulser_time", &pulser_time, &b_event_pulser_time);
   fChain->SetBranchAddress("ge_time[2]", ge_time, &b_event_ge_time);
   fChain->SetBranchAddress("n1_time", &n1_time, &b_event_n1_time);
   fChain->SetBranchAddress("n2_time", &n2_time, &b_event_n2_time);
   fChain->SetBranchAddress("bgo6_time", &bgo6_time, &b_event_bgo6_time);
   fChain->SetBranchAddress("bgo1_time", &bgo1_time, &b_event_bgo1_time);
   fChain->SetBranchAddress("bgo5_time", &bgo5_time, &b_event_bgo5_time);
   fChain->SetBranchAddress("bgo4_time", &bgo4_time, &b_event_bgo4_time);
   fChain->SetBranchAddress("buncher_time", &buncher_time, &b_event_buncher_time);
   fChain->SetBranchAddress("sili_en[6]", sili_en, &b_event_sili_en);
   fChain->SetBranchAddress("sili_en_raw[6]", sili_en_raw, &b_event_sili_en_raw);
   fChain->SetBranchAddress("ge_en_raw[2]", ge_en_raw, &b_event_ge_en_raw);
   fChain->SetBranchAddress("ge_en[2]", ge_en, &b_event_ge_en);
   fChain->SetBranchAddress("n1_en[2]", n1_en, &b_event_n1_en);
   fChain->SetBranchAddress("n2_en[2]", n2_en, &b_event_n2_en);
   fChain->SetBranchAddress("pulser_en", &pulser_en, &b_event_pulser_en);
   fChain->SetBranchAddress("veto_timing_en", &veto_timing_en, &b_event_veto_timing_en);
   fChain->SetBranchAddress("bgo6_en", &bgo6_en, &b_event_bgo6_en);
   fChain->SetBranchAddress("bgo5_en", &bgo5_en, &b_event_bgo5_en);
   fChain->SetBranchAddress("bgo1_en", &bgo1_en, &b_event_bgo1_en);
   fChain->SetBranchAddress("bgo4_en", &bgo4_en, &b_event_bgo4_en);
   Notify();
}

Bool_t analysis::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void analysis::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t analysis::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef analysis_cxx
