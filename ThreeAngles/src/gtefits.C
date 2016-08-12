{
  
GetResults *tst = new GetResults();

tst->addDataList("results_list_xxDG_fitRt2.dat");
tst->addDataSpec("../exe/inputParameters.txt");
tst->getFitResults();

//tst->addHistogram("phi_s",50,-0.2,0.1);
//tst->addHistogram("phi_s_error",50,-0.02,0.01);
//tst->addHistogram("phi_s_pull",100,-10.0,10.0);

//tst->addHistogram("Gbar",50,0.00,0.80);
//tst->addHistogram("Gbar_error",50,0.001,0.010);
//tst->addHistogram("Gbar_pull",100,-10.0,10.0);

//tst->addHistogram("DeltaG",50,0.00,0.30);
//tst->addHistogram("DeltaG_error",50,0.001,0.010);
//tst->addHistogram("DeltaG_pull",100,-10.0,10.0);

tst->addHistogram("Rt",50,0.10,0.40);
tst->addHistogram("Rt_error",50,0.005,0.008);
tst->addHistogram("Rt_pull",100,0.0,55.0);

//tst->addHistogram("R0",50,0.50,0.55);
//tst->addHistogram("R0_error",50,0.005,0.008);
//tst->addHistogram("R0_pull",100,0.0,14.0);

tst->plotFittedValue(2,0);
tst->plotFittedError(2,1);
tst->plotPull(2,2);

}
