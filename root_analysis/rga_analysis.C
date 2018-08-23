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
    gStyle->SetLineWidth(2);
    gStyle->SetHistLineWidth(2);
    gStyle->SetTitleXOffset(0.98);
    gStyle->SetTitleYOffset(1.3);
    gStyle->SetTitleSize(0.04,"X");
    gStyle->SetTitleSize(0.04,"Y");
    gStyle->SetLabelSize(0.04,"X");
    gStyle->SetLabelSize(0.04,"Y");
    gPad->SetLeftMargin(0.16);
    gPad->SetRightMargin(0.16);
    gROOT->ForceStyle(true);
    
    //TH1D *h_rec_charge = new TH1D("h_rec_charge","rec_charge",100,0,0);
    //h_rec_charge->SetXTitle("Rec.Charge");
    //h_rec_charge->SetLineWidth(2);
    //h_rec_charge->SetCanExtend(TH1::kXaxis);
    
    //*********** Creating Histograms ************************
    
    TH1D *h_rec_momen_elec = new TH1D("h_rec_momen_elec","Rec.Elec.Momentum",300,0.0,12.0);
    TH1D *h_rec_momen_posit = new TH1D("h_rec_momen_posit","Rec.Posit.Momentum",250,0.0,10.0);
    TH1D *h_rec_momen_proton = new TH1D("h_rec_momen_proton","Rec.Proton.Momentum",200,0.0,8.0);
    
    
    TH1D *h_rec_theta_elec = new TH1D("h_rec_theta_elec","Rec.Elec.theta",90,0.0,45.0);
    TH1D *h_rec_theta_posit = new TH1D("h_rec_theta_posit","Rec.Posit.theta",90,0.0,45.0);
    TH1D *h_rec_theta_proton = new TH1D("h_rec_theta_proton","Rec.Proton.theta",320,0.0,160.0);
    
    TH1D *h_rec_phi_elec = new TH1D("h_rec_phi_elec","Rec.Elec.phi",360,-180.0,180.0);
    TH1D *h_rec_phi_posit = new TH1D("h_rec_phi_posit","Rec.Posit.phi",360,-180.0,180.0);
    TH1D *h_rec_phi_proton = new TH1D("h_rec_phi_proton","Rec.Proton.phi",360,-180.0,180.0);
    
    TH1D *h_elec_sampl_frac= new TH1D("h_elec_sampl_frac","Elec.Sampling.fraction",100,0.1,0.35);
    TH1D *h_posit_sampl_frac= new TH1D("h_posit_sampl_frac","Posit.Sampling.fraction",100,0.1,0.35);
    TH2D *h_elec_samplfrac_vs_P= new TH2D("h_elec_samplfrac_vs_P","Elec.Sampling.fraction.vs.P",300,0.0,12.0,100,0.1,0.35);
    TH2D *h_posit_samplfrac_vs_P= new TH2D("h_posit_samplfrac_vs_P","Posit.Sampling.fraction.vs.P",300,0.0,12.0,100,0.1,0.35);
    
    TH2D *h_theta_phi_elec= new TH2D("h_theta_phi_elec","Rec.Theta.vs.Phi.elec",360,-180.0,180.0,100,0.0,50.0);
    TH2D *h_theta_phi_posit= new TH2D("h_theta_phi_posit","Rec.Theta.vs.Phi.posit",360,-180.0,180.0,100,0.0,50.0);
    TH2D *h_theta_phi_proton= new TH2D("h_theta_phi_proton","Rec.Theta.vs.Phi.proton",360,-180.0,180.0,320,0.0,160.0);
    TH2D *h_theta_mom_elec= new TH2D("h_theta_mom_elec","Rec.Theta.vs.Mom.elec",300,0.0,12.0,120,0.0,60.0);
    TH2D *h_theta_mom_posit= new TH2D("h_theta_mom_posit","Rec.Theta.vs.Mom.posit",275,0.0,11.0,120,0.0,60.0);
    TH2D *h_theta_mom_proton= new TH2D("h_theta_mom_proton","Rec.Theta.vs.Mom.proton",225,0.0,9.0,360,0.0,180.0);
    
    TH2D *h_ECin_vs_ECout_neg_particle= new TH2D("h_ECin_vs_ECout_neg_particle","REC.ECin.vs.ECout.neg.particle",200,0.0,1.2,200,0,1.0);
    TH2D *h_ECin_vs_ECout_pos_particle= new TH2D("h_ECin_vs_ECout_pos_particle","REC.ECin.vs.ECout.pos.particle",200,0.0,1.2,200,0,1.0);
    TH2D *h_beta_vs_mom_neg_particle= new TH2D("h_beta_vs_mom_neg_particle","REC.beta.vs.Mom.neg.particle",275,0.0,0.0,200,0.0,0.0);
    TH2D *h_beta_vs_mom_pos_particle= new TH2D("h_beta_vs_mom_pos_particle","REC.beta.vs.Mom.pos.particle",275,0.0,0.0,200,0.0,0.0);
    //***********************************************
    
    /*
    // add to get default plot
    h_rec_momen_elec->SetCanExtend(TH1::kXaxis);
    h_rec_momen_elec->SetBinsLength(-1);
    h_rec_momen_posit->SetCanExtend(TH1::kXaxis);
    h_rec_momen_posit->SetBinsLength(-1);
    h_rec_momen_proton->SetCanExtend(TH1::kXaxis);
    h_rec_momen_proton->SetBinsLength(-1);
    */
    
    //----------- Functions -----------------------
    vectorMap loadMapIndex(vector<int>* branch_pindex);
    void Calorimeter_energy(vector<double>* cal_ener_branch,vector<int>* cal_layer_branch, vectorMap Cmap,int keyval, double &Cal_energy, double &Pcal_energy, double &ECout_energy, double &ECin_energy);
    //---------------------------------------------
    
    if (fChain == 0) return;
    
    Long64_t nentries = fChain->GetEntriesFast();
    
    Long64_t nbytes = 0, nb = 0;
    
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //start loop over all the events/jentry
    for (Long64_t jentry=0; jentry<nentries;jentry++) {
        
        //if(jentry>500) break; //for debugging
        
        double elec_cal_energy=0.0, posit_cal_energy=0.0;
        
        Long64_t ientry = LoadTree(jentry);
        if (ientry < 0) break;
        nb = fChain->GetEntry(jentry);   nbytes += nb;
        
        //-----------------------------------------------------
        //------load map for each entry/event for pindex-------
        vectorMap calomap = loadMapIndex(REC_Calorimeter_pindex);
        //-----------------------------------------------------
        
        //*********************************************************
        // loop over vectors elements of "particle" branch/bank
        for(int vecEntry=0; vecEntry< REC_Particle_pid->size();++vecEntry){
            
            //h_rec_charge->Fill(REC_Particle_charge->at(vecEntry));
            
            //------------------------------------------------------
            //-------- Looking the calorimeter map -----------------
            double Cal_energy=0.0, Pcal_energy=0.0, ECout_energy=0.0, ECin_energy=0.0;
            Calorimeter_energy(REC_Calorimeter_energy,REC_Calorimeter_layer,calomap,vecEntry,Cal_energy, Pcal_energy,ECout_energy,ECin_energy);
            //------------------------------------------------------
            
            
            //#########################################
            // Choosing different particle ID and getting 4momentum
            
            auto PID = REC_Particle_pid->at(vecEntry);
            auto px = REC_Particle_px->at(vecEntry);
            auto py = REC_Particle_py->at(vecEntry);
            auto pz = REC_Particle_pz->at(vecEntry);
            TVector3 Mom3Vector(px,py,pz);
            auto rec_momentum = Mom3Vector.Mag();
            
            auto beta= REC_Particle_beta->at(vecEntry);
            if(REC_Particle_charge->at(vecEntry)==-1 && rec_momentum<11.0 && rec_momentum>=1.0 && beta>=-1 && beta<=1){
                h_beta_vs_mom_neg_particle->Fill(rec_momentum,beta);
            }else if(REC_Particle_charge->at(vecEntry)==1 && rec_momentum<11.0 && rec_momentum>=1.0 && beta>=-1 && beta<=1){
                h_beta_vs_mom_pos_particle->Fill(rec_momentum,beta);
            }
            
            // Calorimeter plots
            if(calomap.find(vecEntry)!=calomap.end()){
                if(ECout_energy!=0 && ECin_energy!=0 && REC_Particle_charge->at(vecEntry)==-1){
                    h_ECin_vs_ECout_neg_particle->Fill(ECin_energy,ECout_energy);
                }else if(ECout_energy!=0 && ECin_energy!=0 && REC_Particle_charge->at(vecEntry)==1){
                    h_ECin_vs_ECout_pos_particle->Fill(ECin_energy,ECout_energy);
                }
            }
            
            //.........................................................
            // selecting particular PID
            if(PID == 11){
                TVector3 elecMom3Vector(Mom3Vector);
                h_rec_momen_elec->Fill(elecMom3Vector.Mag());
                h_rec_theta_elec->Fill(elecMom3Vector.Theta()*TMath::RadToDeg());
                h_rec_phi_elec->Fill(elecMom3Vector.Phi()*TMath::RadToDeg());
                //TLorentzVector elecMom4Vector(elecMom3Vector,sqrt(elecMom3Vector.Mag2()+0.0005*0.0005));
                h_theta_phi_elec->Fill(elecMom3Vector.Phi()*TMath::RadToDeg(),elecMom3Vector.Theta()*TMath::RadToDeg());
                h_theta_mom_elec->Fill(elecMom3Vector.Mag(),elecMom3Vector.Theta()*TMath::RadToDeg());
                
                if(calomap.find(vecEntry)!=calomap.end()){
                    h_elec_sampl_frac->Fill(Cal_energy/elecMom3Vector.Mag());
                    h_elec_samplfrac_vs_P->Fill(elecMom3Vector.Mag(),Cal_energy/elecMom3Vector.Mag());
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
                    h_posit_sampl_frac->Fill(Cal_energy/posMom3Vector.Mag());
                    h_posit_samplfrac_vs_P->Fill(posMom3Vector.Mag(),Cal_energy/posMom3Vector.Mag());
                }
            }else if(PID== 2212){
                TVector3 protMom3Vector(Mom3Vector);
                h_rec_momen_proton->Fill(protMom3Vector.Mag());
                h_rec_theta_proton->Fill(protMom3Vector.Theta()*TMath::RadToDeg());
                h_rec_phi_proton->Fill(protMom3Vector.Phi()*TMath::RadToDeg());
                //TLorentzVector protMom4Vector(protMom3Vector,sqrt(protMom3Vector.Mag2()+0.938*0.938));
                h_theta_phi_proton->Fill(protMom3Vector.Phi()*TMath::RadToDeg(),protMom3Vector.Theta()*TMath::RadToDeg());
                h_theta_mom_proton->Fill(protMom3Vector.Mag(),protMom3Vector.Theta()*TMath::RadToDeg());
            }else if(PID== 22){
                TVector3 photMom3Vector(Mom3Vector);
                
                //TLorentzVector photMom4Vector(photMom3Vector,sqrt(photMom3Vector.Mag2()));
            }else{
            }
            //...................................................
            
            //####################################################
            
        }//end of "for" loop over vecEntry of particle branch/bank
        //********************************************************
        
    } // end of the "for" loop over all events/entry "jentry<=nentries" of chain
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++
    
    //c1->cd(); h_rec_charge->Draw();
    //c1->Print("figs/rec_charge.png"); c1->Clear();
    ///*
    c1->Divide(3,1);
    c1->cd(1); h_rec_momen_elec->Draw();
    h_rec_momen_elec->GetXaxis()->SetTitle("e- mom [GeV]");
    h_rec_momen_elec->GetYaxis()->SetTitle("counts");
    c1->cd(2); h_rec_momen_posit->Draw();
    h_rec_momen_posit->GetXaxis()->SetTitle("e+ mom [GeV]");
    h_rec_momen_posit->GetYaxis()->SetTitle("counts");
    c1->cd(3); h_rec_momen_proton->Draw();
    h_rec_momen_proton->GetXaxis()->SetTitle("e- mom [GeV]");
    h_rec_momen_proton->GetYaxis()->SetTitle("counts");
    c1->Print("figs/rec_mom.png"); //c1->Clear();
    
    c2->Divide(3,2);
    c2->cd(1); h_rec_theta_elec->Draw();
    h_rec_theta_elec->GetXaxis()->SetTitle("theta [degree]");
    h_rec_theta_elec->GetYaxis()->SetTitle("counts");
    c2->cd(2); h_rec_theta_posit->Draw();
    h_rec_theta_posit->GetXaxis()->SetTitle("theta [degree]");
    h_rec_theta_posit->GetYaxis()->SetTitle("counts");
    c2->cd(3); h_rec_theta_proton->Draw();
    h_rec_theta_proton->GetXaxis()->SetTitle("theta [degree]");
    h_rec_theta_proton->GetYaxis()->SetTitle("counts");
    c2->cd(4); h_rec_phi_elec->Draw();
    h_rec_phi_elec->GetXaxis()->SetTitle("phi [degree]");
    h_rec_phi_elec->GetYaxis()->SetTitle("counts");
    c2->cd(5); h_rec_phi_posit->Draw();
    h_rec_phi_posit->GetXaxis()->SetTitle("phi [degree]");
    h_rec_phi_posit->GetYaxis()->SetTitle("counts");
    c2->cd(6); h_rec_phi_proton->Draw();
    h_rec_phi_proton->GetXaxis()->SetTitle("phi [degree]");
    h_rec_phi_proton->GetYaxis()->SetTitle("counts");
    c2->Print("figs/rec_theta_phi.png");
    
    c3->Divide(2,2);
    c3->cd(1); h_elec_sampl_frac->Draw();
    h_elec_sampl_frac->GetXaxis()->SetTitle("e- E_Cal/P");
    h_elec_sampl_frac->GetYaxis()->SetTitle("counts");
    c3->cd(2); h_posit_sampl_frac->Draw();
    h_posit_sampl_frac->GetXaxis()->SetTitle("e+ E_Cal/P");
    h_posit_sampl_frac->GetYaxis()->SetTitle("counts");
    c3->cd(3); h_elec_samplfrac_vs_P->Draw("COLZ");
    h_elec_samplfrac_vs_P->SetXTitle("e- P [GeV]");
    h_elec_samplfrac_vs_P->SetYTitle("sampl.frac");
    c3->cd(4); h_posit_samplfrac_vs_P->Draw("COLZ");
    h_posit_samplfrac_vs_P->SetXTitle("e+ P [GeV]");
    h_posit_samplfrac_vs_P->SetYTitle("sampl.frac");
    c3->Print("figs/sampling_frac_&mom.png");
    
    c4->Divide(3,2);
    c4->cd(1); h_theta_phi_elec->Draw("COLZ");
    h_theta_phi_elec->GetXaxis()->SetTitle("phi [degree]");
    h_theta_phi_elec->GetYaxis()->SetTitle("theta [degree]");
    c4->cd(2); h_theta_phi_posit->Draw("COLZ");
    h_theta_phi_posit->GetXaxis()->SetTitle("phi [degree]");
    h_theta_phi_posit->GetYaxis()->SetTitle("theta [degree]");
    c4->cd(3); h_theta_phi_proton->Draw("COLZ");
    h_theta_phi_proton->GetXaxis()->SetTitle("phi [degree]");
    h_theta_phi_proton->GetYaxis()->SetTitle("theta [degree]");
    c4->cd(4); h_theta_mom_elec->Draw("COLZ");
    h_theta_mom_elec->GetXaxis()->SetTitle("e- mom [GeV]");
    h_theta_mom_elec->GetYaxis()->SetTitle("theta [degree]");
    c4->cd(5); h_theta_mom_posit->Draw("COLZ");
    h_theta_mom_posit->GetXaxis()->SetTitle("e+ mom [GeV]");
    h_theta_mom_posit->GetYaxis()->SetTitle("theta [degree]");
    c4->cd(6); h_theta_mom_proton->Draw("COLZ");
    h_theta_mom_proton->GetXaxis()->SetTitle("e+ mom [GeV]");
    h_theta_mom_proton->GetYaxis()->SetTitle("theta [degree]");
    c4->Print("figs/theta_phi_mom.png");
    //*/
    c5->Divide(2,2);
    c5->cd(1);
    gPad->SetLogz();
    h_ECin_vs_ECout_neg_particle->Draw("COLZ");
    h_ECin_vs_ECout_neg_particle->SetXTitle("ECin [GeV]");
    h_ECin_vs_ECout_neg_particle->SetYTitle("ECout [GeV]");
    c5->cd(2);
    gPad->SetLogz();
    h_ECin_vs_ECout_pos_particle->Draw("COLZ");
    h_ECin_vs_ECout_pos_particle->SetXTitle("ECin [GeV]");
    h_ECin_vs_ECout_pos_particle->SetYTitle("ECout [GeV]");
    c5->cd(3);
    gPad->SetLogz();
    h_beta_vs_mom_neg_particle->Draw("COLZ");
    h_beta_vs_mom_neg_particle->SetXTitle("P [GeV]");
    h_beta_vs_mom_neg_particle->SetYTitle("beta");
    c5->cd(4);
    gPad->SetLogz();
    h_beta_vs_mom_pos_particle->Draw("COLZ");
    h_beta_vs_mom_pos_particle->SetXTitle("P [GeV]");
    h_beta_vs_mom_pos_particle->SetYTitle("beta");
    c5->Print("figs/ECin_ECout_&betaP.png");
    
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
void Calorimeter_energy(vector<double>* cal_ener_branch,vector<int>* cal_layer_branch, vectorMap Cmap,int keyval, double &Cal_energy, double &Pcal_energy, double &ECout_energy, double &ECin_energy){
    auto imap=Cmap.find(keyval);
    if (imap !=Cmap.end()) {
        vector <int> inVect = (*imap).second;
        for (int ix=0; ix<inVect.size(); ix++){
            int mapCol=inVect[ix];
            if(cal_layer_branch->at(mapCol)==1){
                Pcal_energy=cal_ener_branch->at(mapCol);
            }else if(cal_layer_branch->at(mapCol)==4){
                ECin_energy=cal_ener_branch->at(mapCol);
            }else if(cal_layer_branch->at(mapCol)==7){
                ECout_energy=cal_ener_branch->at(mapCol);
            }
        }
        Cal_energy=Pcal_energy+ECin_energy+ECout_energy;
    }
    return;
}
