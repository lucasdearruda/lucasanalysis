void gspline() {
   TCanvas *c1 = new TCanvas("c1","gerrors2",200,10,700,500);
   const Int_t n = 10;
   Double_t x[]  = {-0.22, 0.05, 0.25, 0.35, 0.5, 0.61,0.7,0.85,0.89,0.95};
   Double_t y[]  = {1,2.9,5.6,7.4,9,9.6,8.7,6.3,4.5,1};

   TGraph *gr = new TGraph(n,x,y);
   gr->SetMarkerStyle(21);
   TSpline3 *s3 = new TSpline3("s3",gr->GetX(),gr->GetY(),gr->GetN());
   s3->SetLineColor(kRed);
   gr->Draw("alp");
   s3->Draw("l same");
}
