recover(const char * dir)
{

  using namespace std;

   ofstream f("out-r.txt");

    TH1F * h_pairs[5][15][5][15];


  for(int w1=0;w1<5;w1++)
   for(int w2=0;w2<5;w2++)
    for(int s1=1;s1<15;s1++)
     for(int s2=1;s2<15;s2++)
       {
          std::stringstream s;
  
        s<< dir << "/W" << w1-2 <<"S" << s1 << "_VS_" << "W" << w2-2 <<"S" << s2  ;
          h_pairs[w1][s1][w2][s2]=(TH1F *)_file0->Get(s.str().c_str());
          if(h_pairs[w1][s1][w2][s2]->GetEntries() > 0)
            {
             f << w1 << " " << s1 << " " << w2 << " " << s2 << " " << h_pairs[w1][s1][w2][s2]->GetEntries() << " " <<  h_pairs[w1][s1][w2][s2]->GetMean() << " " <<  h_pairs[w1][s1][w2][s2]->GetRMS() << endl;
            }
      }




}
