{
//========= Macro generated from object: deuterons/Graph
//========= by ROOT version6.24/06
   
   cutg = new TCutG("deuteronsT3",33);
   cutg->SetVarX("Medley_1_dE2+Medley_1_Eres");
   cutg->SetVarY("Medley_1_dE1");
   cutg->SetTitle("Graph");
   cutg->SetFillStyle(1000);

   ci = TColor::GetColor("#00ff00");
   cutg->SetLineColor(ci);
   cutg->SetLineWidth(2);
   cutg->SetPoint(0,0.0722176,3.1734);
   cutg->SetPoint(1,0.911163,2.66279);
   cutg->SetPoint(2,1.77325,2.23213);
   cutg->SetPoint(3,2.49334,1.95903);
   cutg->SetPoint(4,3.92773,1.6215);
   cutg->SetPoint(5,5.54191,1.3351);
   cutg->SetPoint(6,7.68267,1.11518);
   cutg->SetPoint(7,9.67878,0.953285);
   cutg->SetPoint(8,12.9767,0.791395);
   cutg->SetPoint(9,16.94,0.642783);
   cutg->SetPoint(10,20.7008,0.547241);
   cutg->SetPoint(11,25.3584,0.471604);
   cutg->SetPoint(12,30.7392,0.398621);
   cutg->SetPoint(13,39.823,0.322984);
   cutg->SetPoint(14,39.9097,0.229362);
   cutg->SetPoint(15,34.6157,0.259907);
   cutg->SetPoint(16,27.2966,0.30878);
   cutg->SetPoint(17,22.8415,0.363761);
   cutg->SetPoint(18,18.4732,0.421797);
   cutg->SetPoint(19,15.8696,0.482888);
   cutg->SetPoint(20,12.7453,0.568415);
   cutg->SetPoint(21,9.07127,0.721141);
   cutg->SetPoint(22,6.69908,0.916631);
   cutg->SetPoint(23,4.90547,1.15454);
   cutg->SetPoint(24,4.05201,1.27939);
   cutg->SetPoint(25,3.18205,1.47641);
   cutg->SetPoint(26,2.65904,1.61387);
   cutg->SetPoint(27,2.35352,1.7223);
   cutg->SetPoint(28,1.88318,1.83241);
   cutg->SetPoint(29,1.15417,2.15476);
   cutg->SetPoint(30,0.778089,2.36365);
   cutg->SetPoint(31,0.106933,2.72468);
   cutg->SetPoint(32,0.0722176,3.1734);
   //cutg->Draw("");
}
