{
//Run the test
Double_t par[10];

par[0]=0.72;         //Gamma
par[1]=0.17;         //DGrate
par[2]=0.14;         //Rt
par[3]=0.64;         //Rp
par[4]=0.000;        //tphase1
par[5]=TMath::Pi();  //tphase2
par[6]=-0.04;        //w phase
par[7]=20.00;        //dms

RandomDists *exp = new RandomDists(100);
exp->RunExperiment(par,8);


}
