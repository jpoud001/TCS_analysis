#define rga_analysis_cxx
#include "rga_analysis.h"
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <sstream>
#include <string>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <math.h>
#include <TH1D.h>
#include <TH1F.h>
#include <TH2D.h>
#include <TF1.h>
#include <TH2F.h>
#include <TFile.h>
#include <iomanip>
#include <time.h>
#include <TVector3.h>
#include <TMath.h>
#include "TCutG.h"
#include <TLine.h>
#include <TLorentzVector.h>
#include <stdlib.h>
#include <vector>
#include <cmath>
#include <TGraphErrors.h>
#include "TLegend.h"
#include <map>

typedef std::unordered_map<int,vector<int>> vectorMap;

void rga_analysis::Loop()
{
    TCanvas *c1 = new TCanvas("c1","c1",2500,1500);
    TCanvas *c2 = new TCanvas("c2","c2",2500,2500);
    TCanvas *c3 = new TCanvas("c3","c3",2500,2500);
    TCanvas *c4 = new TCanvas("c4","c4",2500,2500);
    TCanvas *c5 = new TCanvas("c5","c5",2500,2500);
    
    //gStyle->SetOptTitle(0);
    gStyle->SetTitleAlign(23);
    gStyle->SetLineWidth(3);
    
    //TH1D *h_rec_charge = new TH1D("h_rec_charge","rec_charge",100,0,0);
    //h_rec_charge->SetXTitle("Rec.Charge");
    //h_rec_charge->SetLineWidth(2);
    //h_rec_charge->SetCanExtend(TH1::kXaxis);
    
    TH1D *h_rec_momen_elec = new TH1D("h_rec_momen_elec","Rec.Elec.Momentum",375,0.0,15.0);
    TH1D *h_rec_momen_posit = new TH1D("h_rec_momen_posit","Rec.Posit.Momentum",375,0.0,15.0);
    TH1D *h_rec_momen_proton = new TH1D("h_rec_momen_proton","Rec.Proton.Momentum",375,0.0,15.0);
    
    
    TH1D *h_rec_theta_elec = new TH1D("h_rec_theta_elec","Rec.Elec.theta",360,0.0,180.0);
    TH1D *h_rec_theta_posit = new TH1D("h_rec_theta_posit","Rec.Posit.theta",360,0.0,180.0);
    TH1D *h_rec_theta_proton = new TH1D("h_rec_theta_proton","Rec.Proton.theta",360,0.0,180.0);
    
    TH1D *h_rec_phi_elec = new TH1D("h_rec_phi_elec","Rec.Elec.phi",360,-180.0,180.0);
    TH1D *h_rec_phi_posit = new TH1D("h_rec_phi_posit","Rec.Posit.phi",360,-180.0,180.0);
    TH1D *h_rec_phi_proton = new TH1D("h_rec_phi_proton","Rec.Proton.phi",360,-180.0,180.0);
    
    TH1D *h_elec_sampl_frac= new TH1D("h_elec_sampl_frac","Elec.Sampling.fraction",100,0.0,0.4);
    TH1D *h_posit_sampl_frac= new TH1D("h_posit_sampl_frac","Posit.Sampling.fraction",100,0.0,0.4);
    
    TH2D *h_theta_phi_elec= new TH2D("h_theta_phi_elec","Rec.Theta.vs.Phi.elec",360,-180.0,180.0,360,0.0,180.0);
    TH2D *h_theta_phi_posit= new TH2D("h_theta_phi_posit","Rec.Theta.vs.Phi.posit",360,-180.0,180.0,360,0.0,180.0);
    TH2D *h_theta_mom_elec= new TH2D("h_theta_mom_elec","Rec.Theta.vs.Mom.elec",375,0.0,15.0,360,0.0,180.0);
    TH2D *h_theta_mom_posit= new TH2D("h_theta_mom_posit","Rec.Theta.vs.Mom.posit",375,0.0,15.0,360,0.0,180.0);
    
    TH2D *h_beta_vs_mom= new TH2D("h_beta_vs_mom","REC.beta.vs.Mom",1000,0,0,375,0,0);
    
    /*
    // add to get default plot
    h_rec_momen_elec->SetCanExtend(TH1::kXaxis);
    h_rec_momen_elec->SetBinsLength(-1);
    h_rec_momen_posit->SetCanExtend(TH1::kXaxis);
    h_rec_momen_posit->SetBinsLength(-1);
    h_rec_momen_proton->SetCanExtend(TH1::kXaxis);
    h_rec_momen_proton->SetBinsLength(-1);
    */
    
    
    //Function used
    vectorMap loadMapIndex(vector<int>* branch_pindex);
    double cal_energy(vector<double>* cal_ener_branch,vectorMap Cmap,int keyval);
    
    if (fChain == 0) return;
    
    Long64_t nentries = fChain->GetEntriesFast();
    
    Long64_t nbytes = 0, nb = 0;
    
    for (Long64_t jentry=0; jentry<nentries;jentry++) { //start loop over all the events/jentry
        
        //if(jentry>1500) break; //for debugging
        
        double elec_cal_energy=0.0, posit_cal_energy=0.0;
        
        Long64_t ientry = LoadTree(jentry);
        if (ientry < 0) break;
        nb = fChain->GetEntry(jentry);   nbytes += nb;
        
        //----------------------------------------
        //------load map for each entry/event for pindex-------
        vectorMap calomap = loadMapIndex(REC_Calorimeter_pindex);
        //----------------------------------------
        
        // loop over vectors elements of "particle" branch/bank
        for(int vecEntry=0; vecEntry< REC_Particle_pid->size();++vecEntry){
            //h_rec_charge->Fill(REC_Particle_charge->at(vecEntry));
           
            // Looking for the map between particle and calorimeter
            
            auto calorimeter_energy=cal_energy(REC_Calorimeter_energy,calomap,vecEntry);
            
            //#########################################
            // Choosing different particle ID and getting 4momentum
            auto PID = REC_Particle_pid->at(vecEntry);
            auto px = REC_Particle_px->at(vecEntry);
            auto py = REC_Particle_py->at(vecEntry);
            auto pz = REC_Particle_pz->at(vecEntry);
            TVector3 Mom3Vector(px,py,pz);
            
            auto beta= REC_Particle_beta->at(vecEntry);
            h_beta_vs_mom->Fill(Mom3Vector.Mag(),beta);
            
            if(PID == 11){
                TVector3 elecMom3Vector(Mom3Vector);
                h_rec_momen_elec->Fill(elecMom3Vector.Mag());
                h_rec_theta_elec->Fill(elecMom3Vector.Theta()*TMath::RadToDeg());
                h_rec_phi_elec->Fill(elecMom3Vector.Phi()*TMath::RadToDeg());
                //TLorentzVector elecMom4Vector(elecMom3Vector,sqrt(elecMom3Vector.Mag2()+0.0005*0.0005));
                h_theta_phi_elec->Fill(elecMom3Vector.Phi()*TMath::RadToDeg(),elecMom3Vector.Theta()*TMath::RadToDeg());
                h_theta_mom_elec->Fill(elecMom3Vector.Mag(),elecMom3Vector.Theta()*TMath::RadToDeg());
                
                if(calomap.find(vecEntry)!=calomap.end()){
                    h_elec_sampl_frac->Fill(calorimeter_energy/elecMom3Vector.Mag());
                }
            } else if(PID== -11){
                TVector3 posMom3Vector(Mom3Vector);
                h_rec_momen_posit->Fill(posMom3Vector.Mag());
                h_rec_theta_posit->Fill(posMom3Vector.Theta()*TMath::RadToDeg());
                h_rec_phi_posit->Fill(posMom3Vector.Phi()*TMath::RadToDeg());
                //TLorentzVector posMom4Vector(posMom3Vector,sqrt(posMom3Vector.Mag2()+0.0005*0.0005));
                h_theta_phi_posit->Fill(posMom3Vector.Phi()*TMath::RadToDeg(),posMom3Vector.Theta()*TMath::RadToDeg());
                h_theta_mom_posit->Fill(posMom3Vector.Mag(),posMom3Vector.Theta()*TMath::RadToDeg());
                
                if(calomap.find(vecEntry)!=calomap.end()){
                    h_posit_sampl_frac->Fill(calorimeter_energy/posMom3Vector.Mag());
                }
            }else if(PID== 2212){
                TVector3 protMom3Vector(Mom3Vector);
                h_rec_momen_proton->Fill(protMom3Vector.Mag());
                h_rec_theta_proton->Fill(protMom3Vector.Theta()*TMath::RadToDeg());
                h_rec_phi_proton->Fill(protMom3Vector.Phi()*TMath::RadToDeg());
                //TLorentzVector protMom4Vector(protMom3Vector,sqrt(protMom3Vector.Mag2()+0.938*0.938));
            }else if(PID== 22){
                TVector3 photMom3Vector(Mom3Vector);
                
                //TLorentzVector photMom4Vector(photMom3Vector,sqrt(photMom3Vector.Mag2()));
            }else{
            }
            //###########################################
            
        }//end of for vecEntry of particle branch/bank loop
        
    } // end of the for loop over all events/entry "jentry<=nentries"
    
    
    //c1->cd(); h_rec_charge->Draw();
    //c1->Print("figs/rec_charge.png"); c1->Clear();
    
    c1->Divide(3,1);
    c1->cd(1); h_rec_momen_elec->SetLineWidth(3); h_rec_momen_elec->Draw();
    c1->cd(2); h_rec_momen_posit->SetLineWidth(3); h_rec_momen_posit->Draw();
    c1->cd(3); h_rec_momen_proton->SetLineWidth(3); h_rec_momen_proton->Draw();
    c1->Print("figs/2rec_mom.png"); //c1->Clear();
    
    c2->Divide(3,2);
    c2->cd(1); h_rec_theta_elec->SetLineWidth(3); h_rec_theta_elec->Draw();
    c2->cd(2); h_rec_theta_posit->SetLineWidth(3); h_rec_theta_posit->Draw();
    c2->cd(3); h_rec_theta_proton->SetLineWidth(3); h_rec_theta_proton->Draw();
    c2->cd(4); h_rec_phi_elec->SetLineWidth(3); h_rec_phi_elec->Draw();
    c2->cd(5); h_rec_phi_posit->SetLineWidth(3); h_rec_phi_posit->Draw();
    c2->cd(6); h_rec_phi_proton->SetLineWidth(3); h_rec_phi_proton->Draw();
    c2->Print("figs/2rec_theta_phi.png");
    
    c3->Divide(2,1);
    c3->cd(1); h_elec_sampl_frac->SetLineWidth(3); h_elec_sampl_frac->Draw();
    c3->cd(2); h_posit_sampl_frac->SetLineWidth(3); h_posit_sampl_frac->Draw();
    c3->Print("figs/2sampling_frac.png");
    
    c4->Divide(2,2);
    c4->cd(1); h_theta_phi_elec->Draw("COLZ");
    c4->cd(2); h_theta_phi_posit->Draw("COLZ");
    c4->cd(3); h_theta_mom_elec->Draw("COLZ");
    c4->cd(4); h_theta_mom_posit->Draw("COLZ");
    c4->Print("figs/2theta_phi_mom.png");
    
    c5->cd();
    h_beta_vs_mom->Draw("COLZ");
    c5->Print("figs/2beta_vs_mom.png");
    
} // End of analysis code/script






//----------------------------------------
//   Declaration of functions
//----------------------------------------

//function definition for the mapping (link different REC_branches)
vectorMap loadMapIndex(vector<int>* branch_pindex){
    vectorMap idxmap;
    idxmap.clear();
    if(branch_pindex!=NULL){
        for(int vecEntryPos=0; vecEntryPos<branch_pindex->size(); ++vecEntryPos){
            int vecElement = branch_pindex->at(vecEntryPos);
                idxmap[vecElement].push_back(vecEntryPos);
        }
    }
    return idxmap;
}

//funtion for calculating total cal_energy(PCal+ECin+ECout)
double cal_energy(vector<double>* cal_ener_branch,vectorMap Cmap,int keyval){
    auto imap=Cmap.find(keyval);
    double energy=0.0;
    if (imap !=Cmap.end()) {
        vector <int> inVect = (*imap).second;
        for (int ij=0; ij<inVect.size(); ij++){
            int mapCol=inVect[ij];
            energy=energy+cal_ener_branch->at(mapCol);
        }
    }
    return energy;
}
