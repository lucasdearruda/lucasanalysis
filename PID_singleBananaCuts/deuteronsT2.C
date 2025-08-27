{
//========= Macro generated from object: deuterons/Graph
//========= by ROOT version6.24/06
   
   cutg = new TCutG("deuteronsT2",33);
   cutg->SetVarX("Medley_1_dE2+Medley_1_Eres");
   cutg->SetVarY("Medley_1_dE1");
   cutg->SetTitle("Graph");
   cutg->SetFillStyle(1000);

   ci = TColor::GetColor("#00ff00");
   cutg->SetLineColor(ci);
   cutg->SetLineWidth(2);
   cutg->SetPoint(0,0.0579221,2.7663);
   cutg->SetPoint(1,0.855388,2.30402);
   cutg->SetPoint(2,1.77435,1.88844);
   cutg->SetPoint(3,2.73578,1.59894);
   cutg->SetPoint(4,4.1698,1.30477);
   cutg->SetPoint(5,6.23621,1.03394);
   cutg->SetPoint(6,8.46375,0.847164);
   cutg->SetPoint(7,10.8649,0.725759);
   cutg->SetPoint(8,14.2207,0.595015);
   cutg->SetPoint(9,17.9525,0.496957);
   cutg->SetPoint(10,21.5687,0.431585);
   cutg->SetPoint(11,27.0652,0.352205);
   cutg->SetPoint(12,31.9253,0.31018);
   cutg->SetPoint(13,39.8519,0.249478);
   cutg->SetPoint(14,39.8808,0.165428);
   cutg->SetPoint(15,32.7642,0.198114);
   cutg->SetPoint(16,26.1973,0.240139);
   cutg->SetPoint(17,21.829,0.282164);
   cutg->SetPoint(18,17.9525,0.350509);
   cutg->SetPoint(19,15.725,0.38839);
   cutg->SetPoint(20,13.1503,0.448554);
   cutg->SetPoint(21,9.76557,0.571668);
   cutg->SetPoint(22,7.68267,0.693073);
   cutg->SetPoint(23,5.88906,0.833155);
   cutg->SetPoint(24,4.76082,0.95923);
   cutg->SetPoint(25,3.86402,1.09931);
   cutg->SetPoint(26,2.99614,1.27675);
   cutg->SetPoint(27,2.44649,1.41683);
   cutg->SetPoint(28,1.93467,1.55225);
   cutg->SetPoint(29,1.17077,1.83241);
   cutg->SetPoint(30,0.784102,2.01919);
   cutg->SetPoint(31,0.0862148,2.38807);
   cutg->SetPoint(32,0.0579221,2.7663);
   //cutg->Draw("");
}
