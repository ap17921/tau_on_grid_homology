#include "hfk_class_graph2.cpp"

int Progress(const char *message, int percent){
  switch(percent) {
  case 0:
    cout << message << "       ";
    break;
  case -1:
    cout << "\n" << message << "\n";
    break;
  default:
    cout << "\b\b\b\b\b\b\b" << setw(4) << percent << "%  " << flush;
  }
  return 0;
}

int main(int argc, char *argv[]){
  int m, a;

  const int gridsize = 8;
  int black[8] = {0, 1, 2, 3, 4, 5, 6, 7};
  int white[8] = {3, 6, 0, 1, 7, 2, 4, 5}; //T(2,3)#T(2,3)
  //const int gridsize = 10;
  //int black[10] = {0, 1, 2, 8, 9, 5, 6, 7, 3, 4};
  //int white[10] = {6, 7, 9, 0, 3, 1, 2, 4, 5, 8}; //12n403
  //const int gridsize = 12;
  //int black[12] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
  //int white[12] = {8, 6, 5, 10, 7, 9, 2, 0, 11, 4, 1, 3}; //k14
  //const int gridsize = 10;
  //int black[10] = {0,1,2,9,8,5,6,7,4,3};
  //int white[10] = {6,8,0,7,4,1,3,5,2,9}; //10_155
    

  if(gridsize % 2 != 0) {
    cout << "Edges2 requires even gridsize.\n";
    return 0;
  }

  if(!ValidGrid(gridsize, black, white)) {
    cout << "Invalid grid!!\n";
    return 0;
  }

  Link link(gridsize, black, white, &Progress);
  if ( link.aborted ) {
    cout << "\n*****   Computation aborted!  *****\n";
    return -1;
  }

  cout << "\nHFK^ ranks:\n";
  for(a = link.HFK_maxA; a >= -link.HFK_maxA; a--) {
    for(m = link.HFK_minM; m <= link.HFK_maxM; m++)
      cout << setw(3) << link.HFK_Rank(m,a);
    cout << "\n";
  }
  cout << "\n";

  link.PrintReducedEdges2(cout);
  link.PrintEdges2HomologyRanks(cout);
  return 0;
}
