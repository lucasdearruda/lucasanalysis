{
//========= Macro generated from object: protons/Graph
//========= by ROOT version6.24/06
   
   TCutG *cutg = new TCutG("protonsT2",33);
   cutg->SetVarX("Medley_1_dE2+Medley_1_Eres");
   cutg->SetVarY("Medley_1_dE1");
   cutg->SetTitle("Graph");
   cutg->SetFillStyle(1000);

   Int_t ci;      // for color index setting
   TColor *color; // for color definition with alpha
   ci = TColor::GetColor("#ff0000");
   cutg->SetLineColor(ci);
   cutg->SetLineWidth(2);
   cutg->SetPoint(0,0.0703755,2.16861);
   cutg->SetPoint(1,0.737597,1.81373);
   cutg->SetPoint(2,1.46221,1.48221);
   cutg->SetPoint(3,2.54555,1.16935);
   cutg->SetPoint(4,3.62889,0.968569);
   cutg->SetPoint(5,4.87007,0.833155);
   cutg->SetPoint(6,7.07262,0.651048);
   cutg->SetPoint(7,9.76303,0.515635);
   cutg->SetPoint(8,11.9931,0.450263);
   cutg->SetPoint(9,14.8282,0.382819);
   cutg->SetPoint(10,17.9814,0.328226);
   cutg->SetPoint(11,20.9611,0.281432);
   cutg->SetPoint(12,25.4452,0.240208);
   cutg->SetPoint(13,30.5946,0.204555);
   cutg->SetPoint(14,37.9715,0.16556);
   cutg->SetPoint(15,39.823,0.158875);
   cutg->SetPoint(16,39.8741,0.0813785);
   cutg->SetPoint(17,34.0341,0.0907174);
   cutg->SetPoint(18,27.7278,0.100056);
   cutg->SetPoint(19,21.3641,0.123403);
   cutg->SetPoint(20,18.1786,0.14675);
   cutg->SetPoint(21,12.31,0.216792);
   cutg->SetPoint(22,9.47605,0.291503);
   cutg->SetPoint(23,6.55606,0.370883);
   cutg->SetPoint(24,4.56157,0.520304);
   cutg->SetPoint(25,3.07646,0.71642);
   cutg->SetPoint(26,2.26575,0.907866);
   cutg->SetPoint(27,1.83529,1.01993);
   cutg->SetPoint(28,1.21828,1.20204);
   cutg->SetPoint(29,0.981528,1.30944);
   cutg->SetPoint(30,0.421922,1.58493);
   cutg->SetPoint(31,0.0703755,1.81373);
   cutg->SetPoint(32,0.0703755,2.16861);
   //cutg->Draw("");
}
