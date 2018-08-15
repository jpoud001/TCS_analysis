#include <TChain.h>


void make()
{
  TChain *ch = new TChain("tr1", "rga_analysis");
  ch->Add("./skim_TCS_JPsi_*_All.root/tr1"); 
  ch->MakeClass("rga_analysis");
}
