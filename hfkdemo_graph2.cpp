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

  // 12n403
  //const int gridsize = 10;
  //int black[10] = {0, 1, 2, 8, 9, 5, 6, 7, 3, 4};
  //int white[10] = {6, 7, 9, 0, 3, 1, 2, 4, 5, 8};

  const int gridsize = 8;
  int black[8] = {0, 1, 2, 3, 4, 5, 6, 7};
  int white[8] = {3, 6, 0, 1, 7, 2, 4, 5};
    
    
  if(gridsize % 2 != 0) {
    cout << "Edges2 requires even gridsize.\n";
    return 0;
  }

  if(!ValidGrid(gridsize, black, white)) {
    cout << "Invalid grid!!\n";
    return 0;
  }

  cout << "Running 12n403 only.\n";
  cout << "black = {0, 1, 2, 8, 9, 5, 6, 7, 3, 4}\n";
  cout << "white = {6, 7, 9, 0, 3, 1, 2, 4, 5, 8}\n\n";

  Link link(gridsize, black, white, &Progress);
  if ( link.aborted ) {
    cout << "\n*****   Computation aborted!  *****\n";
    return -1;
  }

  cout << "\nHFK^ ranks:\n";
  for(a = link.HFK_maxA; a >= -link.HFK_maxA; a--) {
    cout << "A=" << setw(3) << a << ":";
    for(m = link.HFK_minM; m <= link.HFK_maxM; m++)
      cout << setw(3) << link.HFK_Rank(m,a);
    cout << "\n";
  }
  cout << "\n";

  link.PrintReducedEdges2(cout);
  link.PrintEdges2HomologyRanks(cout);
  return 0;
}
