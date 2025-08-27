{
//========= Macro generated from object: protons/Graph
//========= by ROOT version6.30/04
   
   TCutG *cutg = new TCutG("protons",20);
   cutg->SetVarX("si2_1+csi_1");
   cutg->SetVarY("si1_1");
   cutg->SetTitle("Graph");
   cutg->SetFillStyle(1000);
   cutg->SetPoint(0,0.132438,2.20183);
   cutg->SetPoint(1,0.984564,1.80272);
   cutg->SetPoint(2,1.94178,1.32379);
   cutg->SetPoint(3,3.49155,0.977902);
   cutg->SetPoint(4,5.26924,0.751741);
   cutg->SetPoint(5,8.50554,0.585446);
   cutg->SetPoint(6,13.0637,0.452411);
   cutg->SetPoint(7,17.4395,0.405848);
   cutg->SetPoint(8,29.9745,0.29942);
   cutg->SetPoint(9,32.208,0.259509);
   cutg->SetPoint(10,32.2992,0.0932143);
   cutg->SetPoint(11,25.2796,0.126473);
   cutg->SetPoint(12,17.1477,0.106518);
   cutg->SetPoint(13,9.73624,0.232902);
   cutg->SetPoint(14,6.91018,0.372589);
   cutg->SetPoint(15,4.49435,0.532232);
   cutg->SetPoint(16,1.71387,0.911384);
   cutg->SetPoint(17,0.528747,1.29054);
   cutg->SetPoint(18,0.12349,1.77612);
   cutg->SetPoint(19,0.132438,2.20183);
   cutg->Draw("same");
}
