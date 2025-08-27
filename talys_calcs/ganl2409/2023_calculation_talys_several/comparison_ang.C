{

Int_t mycolor[5]={kBlack,kRed,kGreen+2,kMagenta,kPink};

TGraph *grFe=new TGraph("Fe56/pp0020.000ang.L00","%lg %lg");
TGraph *grC=new TGraph("C12/pp0020.000ang.L00","%lg %lg");
//TGraph *grH=new TGraph("./H1/");


grFe->SetLineColor(mycolor[1]);
grC->SetLineColor(mycolor[2]);
//grH->SetLineColor(mycolor[3]);

grFe->Draw("al");
grC->Draw("l");
//grH->Draw("l");

TLegend *legend=new TLegend(0.6,0.88,0.7,0.88);
legend->SetFillColor(0);
legend->SetFillStyle(0);

//legend->AddEntry(grH,"H1","l");
legend->AddEntry(grC,"C12","l");
legend->AddEntry(grFe,"Fe56","l");
legend->Draw();

}
