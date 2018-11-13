void asymmetry(){//plot the beam spin asymmetry from helicity_phi.txt file
    gROOT->Reset();
    
    ifstream fin("helicity_phi.txt");
    
    TCanvas *c1 = new TCanvas("asymmetry","BSA",2500,2500);
    TH1D *h1=new TH1D("h1","positive",9,0,360);
    TH1D *h2=new TH1D("h2","negative",9,0,360);
    TH1D *h3=new TH1D("h3","sum",9,0,360);
    TH1D *h4=new TH1D("h4","sub",9,0,360);
    TH1D *h=new TH1D("h","asym",9,0,360);
    
    TF1 *f1= new TF1("f1","[0]*sin(x*TMath::DegToRad())",0,360);
    
    double helic,phi_angle;
    
    while(!fin.eof()){
        fin>>helic>>phi_angle;
        if(helic==1){
            h1->Fill(phi_angle);
        }else if(helic==0){
            h2->Fill(phi_angle);
        }
    }
    c1->Divide(3,1);
    
    c1->cd(1);
    h1->SetMarkerStyle(8);
    h1->SetMarkerSize(1.15);
    h1->SetMarkerColor(4);
    h1->Draw("P");
    
    c1->cd(2);
    h2->SetMarkerStyle(20);
    h2->SetMarkerSize(1.15);
    h2->SetMarkerColor(4);
    h2->Draw("P");
    
    h3->Add(h1,h2,1,1);
    h4->Add(h1,h2,1,-1);
    
    c1->cd(3);
    h->Divide(h4,h3,1,1);
    h->SetMarkerStyle(22);
    h->SetMarkerSize(1.15);
    h->SetMarkerColor(4);
    h->SetStats(0000);
    h->Draw("P");
    h->Fit(f1);
    
    
    
}
