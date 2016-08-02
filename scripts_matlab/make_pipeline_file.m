function make_pipeline_file( filename )

fin = fopen(filename,'wt');
fprintf(fin,'MOTCOR=[0,1]\nCENSOR=[0,1]\nRETROICOR=[0,1]\nTIMECOR=[0,1]\nSMOOTH=[6]\nDETREND=[0,1,2,3,4,5]\nMOTREG=[0,1]\nTASK=[0,1]\nGSPC1=[0,1]\nPHYPLUS=[0,1]\nCUSTOMREG=[0]\nLOWPASS=[0]\n');
fclose(fin);
