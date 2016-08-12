{
  
GetResults *tst = new GetResults();

tst->addDataList("list_xxxS1.dat");
tst->addDataSpec("../exe/inputParameters.txt");
tst->getFitResults();


tst->addSelection( 4 , "GBar"   , 0.50 , 0.70);
tst->addSelection( 3 , "DeltaG" , -0.28 , 0.28);
tst->addSelection( 9 , "Rt"     , 0.16 , 0.24);
tst->addSelection( 2 , "phi_s"  , -1.0 , 1.0);

tst->drawHistograms();


}
