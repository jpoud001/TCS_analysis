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
    
    //*********** Creating Histograms ************************
    
    TH1D *h_rec_momen_elec = new TH1D("h_rec_momen_elec","Rec.Elec.Momentum",300,0.0,12.0);
    TH1D *h_rec_momen_posit = new TH1D("h_rec_momen_posit","Rec.Posit.Momentum",250,0.0,10.0);
    TH1D *h_rec_momen_proton = new TH1D("h_rec_momen_proton","Rec.Proton.Momentum",200,0.0,8.0);
    
    TH1I *h_rec_charge = new TH1I("h_rec_charge","Rec.Particle.Charge",6,-3,3);
    TH1I *h_rec_pid= new TH1I("h_rec_pid", "Rec.Particle.PID",1000,0,0);
    TH1D *h_rec_beta= new TH1D("h_rec_beta","Rec.Particle.beta",1000,-2.0,2.0);
    TH1D *h_rec_chi2pid_elec = new TH1D("h_rec_chi2pid_elec","Rec.Elec.chi2pid",30,0.0,5.0);
    TH1D *h_rec_chi2pid_posit = new TH1D("h_rec_chi2pid_posit","Rec.Posit.chi2pid",60,0.0,10.0);
    TH1D *h_rec_chi2pid_proton = new TH1D("h_rec_chi2pid_proton","Rec.Proton.chi2pid",300,0.0,50.0);
    
    
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
    TH2D *h_ECal_vs_Pcal_neg_particle= new TH2D("h_ECal_vs_Pcal_neg_particle","REC.ECal.vs.Pcal.neg.particle",200,0.0,1.2,200,0,1.0);
    TH2D *h_ECal_vs_Pcal_pos_particle= new TH2D("h_ECal_vs_Pcal_pos_particle","REC.ECal.vs.Pcal.pos.particle",200,0.0,1.2,200,0,1.0);
    
    TH1D *h_rec_lu_neg=new TH1D("h_rec_lu_neg","REC.lu.neg.particle",1000,0.0,0.0);
    TH1D *h_rec_lv_neg=new TH1D("h_rec_lv_neg","REC.lv.neg.particle",1000,0.0,0.0);
    TH1D *h_rec_lw_neg=new TH1D("h_rec_lw_neg","REC.lw.neg.particle",1000,0.0,0.0);
    TH1D *h_rec_lu_pos=new TH1D("h_rec_lu_pos","REC.lu.pos.particle",1000,0.0,0.0);
    TH1D *h_rec_lv_pos=new TH1D("h_rec_lv_pos","REC.lv.pos.particle",1000,0.0,0.0);
    TH1D *h_rec_lw_pos=new TH1D("h_rec_lw_pos","REC.lw.pos.particle",1000,0.0,0.0);
    
    TH2D *h_beta_vs_mom_neg_particle= new TH2D("h_beta_vs_mom_neg_particle","REC.beta.vs.Mom.neg.particle",275,0.0,0.0,200,0.0,0.0);
    TH2D *h_beta_vs_mom_pos_particle= new TH2D("h_beta_vs_mom_pos_particle","REC.beta.vs.Mom.pos.particle",275,0.0,0.0,200,0.0,0.0);
    TH2D *h_rec_pcalXY_neg_particle=new TH2D("h_rec_pcalXY_neg_particle","REC.PCalXY.neg.particle",800,0,0,800,0,0);
    TH2D *h_rec_pcalXY_pos_particle=new TH2D("h_rec_pcalXY_pos_particle","REC.PCalXY.pos.particle",800,0,0,800,0,0);
    
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
    void Calorimeter_bank(vector<double>* cal_ener_branch,vector<int>* cal_layer_branch,vector<double>* cal_lu_branch,vector<double>* cal_lv_branch,vector<double>* cal_lw_branch,vector<double>* cal_X_branch,vector<double>* cal_Y_branch,vector<double>* cal_Z_branch, vectorMap Cmap,int keyval,vector<double> &Cal_energy, vector<double> &lu, vector<double> &lv, vector<double> &lw, vector<double> &Xcal,vector<double> &Ycal,vector<double> &Zcal);
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
            
            //#########################################
            // Choosing different particle ID and getting 4momentum
            
            auto PID = REC_Particle_pid->at(vecEntry);
            auto beta= REC_Particle_beta->at(vecEntry);
            auto particle_charge= REC_Particle_charge->at(vecEntry);
            auto px = REC_Particle_px->at(vecEntry);
            auto py = REC_Particle_py->at(vecEntry);
            auto pz = REC_Particle_pz->at(vecEntry);
            auto chi2 = REC_Particle_chi2pid->at(vecEntry);
            TVector3 Mom3Vector(px,py,pz);
            auto rec_momentum = Mom3Vector.Mag();
            
            h_rec_charge->Fill(particle_charge);
            h_rec_pid->Fill(PID);
            h_rec_beta->Fill(beta);
            
            if(particle_charge==-1 && rec_momentum<11.0 && rec_momentum>=0.25){
                h_beta_vs_mom_neg_particle->Fill(rec_momentum,beta);
            }else if(particle_charge==1 && rec_momentum<11.0 && rec_momentum>=0.25){
                h_beta_vs_mom_pos_particle->Fill(rec_momentum,beta);
            }
            
            
            //------------------------------------------------------
            //-------- Looking the calorimeter map -----------------
            vector<double> Cal_energy; vector<double> lu; vector<double> lv; vector<double> lw; vector<double> Xcal; vector<double> Ycal; vector<double> Zcal;
            Cal_energy.clear(); lu.clear(); lv.clear(); lw.clear(); Xcal.clear(); Ycal.clear(); Zcal.clear();
           //Implementation of Calorimeter bank function
            Calorimeter_bank(REC_Calorimeter_energy,REC_Calorimeter_layer,REC_Calorimeter_lu,REC_Calorimeter_lv,REC_Calorimeter_lw,REC_Calorimeter_x,REC_Calorimeter_y,REC_Calorimeter_z,calomap,vecEntry, Cal_energy,lu,lv,lw,Xcal,Ycal,Zcal);
            //------------------------------------------------------
            
            
            // Calorimeter plots
            auto Cal_total_energy=0.0;
            if(calomap.find(vecEntry)!=calomap.end()){
                for(int jj=0; jj<Cal_energy.size();++jj){
                    auto lu_value= lu.at(jj);
                    auto lv_value= lv.at(jj);
                    auto lw_value= lw.at(jj);
                    auto cal_x_position=Xcal.at(jj);
                    auto cal_y_position=Ycal.at(jj);
                    Cal_total_energy=Cal_total_energy+Cal_energy.at(jj);
                    
                    if(particle_charge==-1){
                        h_rec_lu_neg->Fill(lu_value);
                        h_rec_lv_neg->Fill(lv_value);
                        h_rec_lw_neg->Fill(lw_value);
                        h_rec_pcalXY_neg_particle->Fill(cal_x_position,cal_y_position);
                        if(jj==2){
                            //auto PCal_energy=Cal_energy.at(0);
                            //auto ECin_energy=Cal_energy.at(1);
                            //auto ECout_energy=Cal_energy.at(2);
                            h_ECin_vs_ECout_neg_particle->Fill(Cal_energy.at(1),Cal_energy.at(2));
                            h_ECal_vs_Pcal_neg_particle->Fill(Cal_energy.at(1)+Cal_energy.at(2),Cal_energy.at(0));
                        }
                    }else if(particle_charge==1){
                        h_rec_lu_pos->Fill(lu_value);
                        h_rec_lv_pos->Fill(lv_value);
                        h_rec_lw_pos->Fill(lw_value);
                        h_rec_pcalXY_pos_particle->Fill(cal_x_position,cal_y_position);
                        if(jj==2){
                            h_ECin_vs_ECout_pos_particle->Fill(Cal_energy.at(1),Cal_energy.at(2));
                            h_ECal_vs_Pcal_pos_particle->Fill(Cal_energy.at(1)+Cal_energy.at(2),Cal_energy.at(0));
                        }
                    }
                }
            }
            
            
            //.........................................................
            // selecting particular PID
            if(PID == 11){
                TVector3 elecMom3Vector(Mom3Vector);
                h_rec_momen_elec->Fill(elecMom3Vector.Mag());
                h_rec_chi2pid_elec->Fill(chi2);
                h_rec_theta_elec->Fill(elecMom3Vector.Theta()*TMath::RadToDeg());
                h_rec_phi_elec->Fill(elecMom3Vector.Phi()*TMath::RadToDeg());
                //TLorentzVector elecMom4Vector(elecMom3Vector,sqrt(elecMom3Vector.Mag2()+0.0005*0.0005));
                h_theta_phi_elec->Fill(elecMom3Vector.Phi()*TMath::RadToDeg(),elecMom3Vector.Theta()*TMath::RadToDeg());
                h_theta_mom_elec->Fill(elecMom3Vector.Mag(),elecMom3Vector.Theta()*TMath::RadToDeg());
                
                if(calomap.find(vecEntry)!=calomap.end()){
                    h_elec_sampl_frac->Fill(Cal_total_energy/elecMom3Vector.Mag());
                    h_elec_samplfrac_vs_P->Fill(elecMom3Vector.Mag(),Cal_total_energy/elecMom3Vector.Mag());
                }
            } else if(PID== -11){
                TVector3 postrnMom3Vector(Mom3Vector);
                h_rec_momen_posit->Fill(postrnMom3Vector.Mag());
                h_rec_chi2pid_posit->Fill(chi2);
                h_rec_theta_posit->Fill(postrnMom3Vector.Theta()*TMath::RadToDeg());
                h_rec_phi_posit->Fill(postrnMom3Vector.Phi()*TMath::RadToDeg());
                //TLorentzVector posMom4Vector(postrnMom3Vector,sqrt(postrnMom3Vector.Mag2()+0.0005*0.0005));
                h_theta_phi_posit->Fill(postrnMom3Vector.Phi()*TMath::RadToDeg(),postrnMom3Vector.Theta()*TMath::RadToDeg());
                h_theta_mom_posit->Fill(postrnMom3Vector.Mag(),postrnMom3Vector.Theta()*TMath::RadToDeg());
                
                if(calomap.find(vecEntry)!=calomap.end()){
                    h_posit_sampl_frac->Fill(Cal_total_energy/postrnMom3Vector.Mag());
                    h_posit_samplfrac_vs_P->Fill(postrnMom3Vector.Mag(),Cal_total_energy/postrnMom3Vector.Mag());
                }
            }else if(PID== 2212){
                TVector3 protMom3Vector(Mom3Vector);
                h_rec_momen_proton->Fill(protMom3Vector.Mag());
                h_rec_chi2pid_proton->Fill(chi2);
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
    
    c2->cd(1); h_theta_phi_elec->Draw("COLZ");
    h_theta_phi_elec->GetXaxis()->SetTitle("phi [degree]");
    h_theta_phi_elec->GetYaxis()->SetTitle("theta [degree]");
    c2->cd(2); h_theta_phi_posit->Draw("COLZ");
    h_theta_phi_posit->GetXaxis()->SetTitle("phi [degree]");
    h_theta_phi_posit->GetYaxis()->SetTitle("theta [degree]");
    c2->cd(3); h_theta_phi_proton->Draw("COLZ");
    h_theta_phi_proton->GetXaxis()->SetTitle("phi [degree]");
    h_theta_phi_proton->GetYaxis()->SetTitle("theta [degree]");
    c2->cd(4); h_theta_mom_elec->Draw("COLZ");
    h_theta_mom_elec->GetXaxis()->SetTitle("e- momemtum [GeV]");
    h_theta_mom_elec->GetYaxis()->SetTitle("theta [degree]");
    c2->cd(5); h_theta_mom_posit->Draw("COLZ");
    h_theta_mom_posit->GetXaxis()->SetTitle("e+ momemtum [GeV]");
    h_theta_mom_posit->GetYaxis()->SetTitle("theta [degree]");
    c2->cd(6); h_theta_mom_proton->Draw("COLZ");
    h_theta_mom_proton->GetXaxis()->SetTitle("p momentum [GeV]");
    h_theta_mom_proton->GetYaxis()->SetTitle("theta [degree]");
    c2->Print("figs/theta_phi_mom.png");
    
    c2->cd(1);
    h_rec_lu_neg->Draw();
    h_rec_lu_neg->SetXTitle("lu");
    h_rec_lu_neg->SetYTitle("counts");
    c2->cd(2);
    h_rec_lv_neg->Draw();
    h_rec_lv_neg->SetXTitle("lv");
    h_rec_lv_neg->SetYTitle("counts");
    c2->cd(3);
    h_rec_lw_neg->Draw();
    h_rec_lw_neg->SetXTitle("lw");
    h_rec_lw_neg->SetYTitle("counts");
    c2->cd(4);
    h_rec_lu_pos->Draw();
    h_rec_lu_pos->SetXTitle("lu");
    h_rec_lu_pos->SetYTitle("counts");
    c2->cd(5);
    h_rec_lv_pos->Draw();
    h_rec_lv_pos->SetXTitle("lv");
    h_rec_lv_pos->SetYTitle("counts");
    c2->cd(6);
    h_rec_lw_pos->Draw();
    h_rec_lw_pos->SetXTitle("lw");
    h_rec_lw_pos->SetYTitle("counts");
    c2->Print("figs/lu_lv_lw.png");
    
    
    c2->cd(1); h_rec_charge->Draw();
    h_rec_charge->SetXTitle("charge");
    c2->cd(2); h_rec_pid->Draw();
    h_rec_pid->SetXTitle("pid");
    c2->cd(3); h_rec_beta->Draw();
    h_rec_beta->SetXTitle("beta");
    c2->cd(4); h_rec_chi2pid_elec->Draw();
    h_rec_chi2pid_elec->SetXTitle("chi2pid");
    c2->cd(5); h_rec_chi2pid_posit->Draw();
    h_rec_chi2pid_posit->SetXTitle("chi2pid");
    c2->cd(6); h_rec_chi2pid_proton->Draw();
    h_rec_chi2pid_proton->SetXTitle("chi2pid");
    c2->Print("figs/Rec_particle.png");
    
    
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
    
    c3->cd(1);
    gPad->SetLogz();
    h_ECin_vs_ECout_neg_particle->Draw("COLZ");
    h_ECin_vs_ECout_neg_particle->SetXTitle("ECin [GeV]");
    h_ECin_vs_ECout_neg_particle->SetYTitle("ECout [GeV]");
    c3->cd(2);
    gPad->SetLogz();
    h_ECin_vs_ECout_pos_particle->Draw("COLZ");
    h_ECin_vs_ECout_pos_particle->SetXTitle("ECin [GeV]");
    h_ECin_vs_ECout_pos_particle->SetYTitle("ECout [GeV]");
    c3->cd(3);
    gPad->SetLogz();
    h_ECal_vs_Pcal_neg_particle->Draw("COLZ");
    h_ECal_vs_Pcal_neg_particle->SetXTitle("ECal [GeV]");
    h_ECal_vs_Pcal_neg_particle->SetYTitle("PCal [GeV]");
    c3->cd(4);
    gPad->SetLogz();
    h_ECal_vs_Pcal_pos_particle->Draw("COLZ");
    h_ECal_vs_Pcal_pos_particle->SetXTitle("ECal [GeV]");
    h_ECal_vs_Pcal_pos_particle->SetYTitle("PCal [GeV]");
    c3->Print("figs/ECin_ECout_Pcal.png");
    
    c3->cd(1);
    h_rec_pcalXY_neg_particle->Draw("COLZ");
    h_rec_pcalXY_neg_particle->SetXTitle("Cal hit-X position");
    h_rec_pcalXY_neg_particle->SetYTitle("Cal hit-Y position");
    c3->cd(2);
    h_rec_pcalXY_pos_particle->Draw("COLZ");
    h_rec_pcalXY_pos_particle->SetXTitle("Cal hit-X position");
    h_rec_pcalXY_pos_particle->SetYTitle("Cal hit-Y position");
    c3->cd(3);
    gPad->SetLogz();
    h_beta_vs_mom_neg_particle->Draw("COLZ");
    h_beta_vs_mom_neg_particle->SetXTitle("P [GeV]");
    h_beta_vs_mom_neg_particle->SetYTitle("beta");
    c3->cd(4);
    gPad->SetLogz();
    h_beta_vs_mom_pos_particle->Draw("COLZ");
    h_beta_vs_mom_pos_particle->SetXTitle("P [GeV]");
    h_beta_vs_mom_pos_particle->SetYTitle("beta");
    c3->Print("figs/PCal_XY_betaP.png");
    
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
void Calorimeter_bank(vector<double>* cal_ener_branch,vector<int>* cal_layer_branch,vector<double>* cal_lu_branch,vector<double>* cal_lv_branch,vector<double>* cal_lw_branch,vector<double>* cal_X_branch,vector<double>* cal_Y_branch,vector<double>* cal_Z_branch, vectorMap Cmap,int keyval, vector<double> &Cal_energy, vector<double> &lu, vector<double> &lv, vector<double> &lw, vector<double> &Xcal,vector<double> &Ycal,vector<double> &Zcal){
    auto imap=Cmap.find(keyval);
    if (imap !=Cmap.end()) {
        vector <int> inVect = (*imap).second;
        for (int ix=0; ix<inVect.size(); ix++){
            int mapCol=inVect[ix];
            if(cal_layer_branch->at(mapCol)==1){
                Cal_energy.push_back(cal_ener_branch->at(mapCol));
                lu.push_back(cal_lu_branch->at(mapCol));
                lv.push_back(cal_lv_branch->at(mapCol));
                lw.push_back(cal_lw_branch->at(mapCol));
                Xcal.push_back(cal_X_branch->at(mapCol));
                Ycal.push_back(cal_Y_branch->at(mapCol));
                Zcal.push_back(cal_Z_branch->at(mapCol));
            }else if(cal_layer_branch->at(mapCol)==4){
                Cal_energy.push_back(cal_ener_branch->at(mapCol));
                lu.push_back(cal_lu_branch->at(mapCol));
                lv.push_back(cal_lv_branch->at(mapCol));
                lw.push_back(cal_lw_branch->at(mapCol));
                Xcal.push_back(cal_X_branch->at(mapCol));
                Ycal.push_back(cal_Y_branch->at(mapCol));
                Zcal.push_back(cal_Z_branch->at(mapCol));
            }else if(cal_layer_branch->at(mapCol)==7){
                Cal_energy.push_back(cal_ener_branch->at(mapCol));
                lu.push_back(cal_lu_branch->at(mapCol));
                lv.push_back(cal_lv_branch->at(mapCol));
                lw.push_back(cal_lw_branch->at(mapCol));
                Xcal.push_back(cal_X_branch->at(mapCol));
                Ycal.push_back(cal_Y_branch->at(mapCol));
                Zcal.push_back(cal_Z_branch->at(mapCol));
            }
        }
    }
    return;
}
