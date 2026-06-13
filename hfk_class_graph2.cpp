#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <vector>
#include <list>
using namespace std;

#if defined(__GNUC__)
typedef long long int Int64;
#elif defined(_MSC_VER)
typedef __int64 Int64;
#endif

// Class Declarations

class Vertex {
public:
  Int64 perm;
  list<int> out;
  list<int> in;
  bool alive;
  Vertex(Int64 p);
  ~Vertex();
};

class Generator {
public:
  int m;
  int a;
  Generator(int m, int a);
  ~Generator();
};

class Homology {
  friend class Link;
public:
  Homology();
  ~Homology();
private:
  list < Generator > gens;
  void add(Generator g);
  int maxM;
  int maxA;
};

// Parallel vertex data for the second operator Edges2.  Graph stores the
// ordinary differential; Graph2 stores only Edges2, using the same indices
// as Graph in each grading level.
class Vertex2 {
public:
  list<int> out;
  list<int> in;
  bool alive;
  Vertex2();
  ~Vertex2();
};

// Permanent storage for the reduced Edges2 graph, after the temporary
// Graph/Graph2 level is cleared.
class ReducedVertex2 {
public:
  Int64 perm;
  int m;
  int a;
  vector<int> out;
  ReducedVertex2(Int64 p, int m, int a);
  ~ReducedVertex2();
};

class BiArray {
public:
  int M;
  int A;
  int& operator()(int m, int a);
  int operator()(int m, int a) const;
  BiArray(int M, int A);
  ~BiArray();
private:
  int *data;
};

class Link {
public:
  int gridsize;
  int *black;
  int *white;
  int Ashift;
  int HFK_maxA;
  int HFK_minM;
  int HFK_maxM;
  int aborted;
  int WindingNumber(int x, int y);
  int MaslovGrading(int g[]);
  int MOS_Rank(int m, int a);
  int HFK_Rank(int m, int a);
  void PrintReducedEdges2(ostream &os=cout);
  void PrintEdges2HomologyRanks(ostream &os=cout);
  Link(int gridsize,
       int *black, 
       int *white, 
       int (*Progress)(const char *, int),
       bool quiet=1);
  ~Link();
private:
  int graphsize;
  vector < vector <Int64> > Generators;
  vector < vector <Vertex> > Graph;
  vector < vector <Vertex2> > Graph2;
  vector < vector <ReducedVertex2> > ReducedGraph2;
  Homology MOS;
  unsigned char Rectangles[16*16*16*16];
  BiArray *MOS_Array;
  BiArray *HFK_Array;
  int (*Progress)(const char *message, int percent);
  bool quiet;
  short *counter;
  int *g;
  int *gij;
  int *g2;
  int  AlexanderShift();
  int  RectDotFree(int xll, int yll, int xur, int yur, int which);
  int  BuildVertices();
  int  BuildEdges(int MM, const char *msg, int *count);
  int  BuildEdges2(int MM, const char *msg, int *count);
  int  Reduce(int MM, const char *msg, int *count);
  int  FindHomology();
  void BuildRectangles();
  void ShiftHalfTurn(int src[], int dest[]);
  void ToggleEdge2(int level, int from, int to);
  void AddColumn2(int level, int dest, int src);
  void AddRow2(int level, int dest, int src);
  void RemoveAllEdges2(int level, int v);
  void StoreReducedEdges2(int level);
  int RankMod2(vector < vector <unsigned long long> > &rows, int ncols);
  void ComputeMOSRanks();
  void ComputeHFKRanks();
};

// Function Prototypes
void      Int2Perm(Int64 k, int h [], int size);
Int64     Perm2Int(int y [], int size);
void      NextPerm(short *counter, int *h, int size);
bool      ValidGrid(int gridsize, int black[], int white[]);
int       Find(vector <Vertex> &V, Int64 x);

// Global Data

#if defined(__GNUC__)
Int64 Factorial[16] = {
  1LL,1LL,2LL,6LL,24LL,120LL,720LL,5040LL,40320LL,362880LL,3628800LL,
  39916800LL,479001600LL,6227020800LL,87178291200LL,1307674368000LL};
#elif defined(_MSC_VER)
Int64 Factorial[16] = {
  1,1,2,6,24,120,720,5040,40320,362880,3628800,
  39916800,479001600,6227020800,87178291200,1307674368000};
#endif

#define Binomial(n,k) (Factorial[n]/(Factorial[k]*Factorial[n-k]))

// Class Methods
Vertex::Vertex(Int64 p): perm(p){
  alive = 1;
}

Vertex::~Vertex(){
}

Generator::Generator(int m, int a): m(m), a(a){
}

Generator::~Generator(){
}

Vertex2::Vertex2(){
  alive = 1;
}

Vertex2::~Vertex2(){
}

ReducedVertex2::ReducedVertex2(Int64 p, int m, int a): perm(p), m(m), a(a) {
}

ReducedVertex2::~ReducedVertex2(){
}

Homology::Homology(){
  maxM = 0;
  maxA = 0;
}

Homology::~Homology(){};

void Homology::add(Generator g){
    if (abs(g.m) > maxM) maxM = abs(g.m);
    if (abs(g.a) > maxA) maxA = abs(g.a);
    gens.push_back(g);
}

BiArray::BiArray(int M, int A): M(M), A(A) {
    int i;
    data = new int[(2*M+1)*(2*A+1)];
    for (i=0; i<(2*M+1)*(2*A+1); i++) data[i]=0;
}

BiArray::~BiArray(){
  delete[] data;
}

inline
int& BiArray::operator()(int m, int a) {
  if ( (m < - M) || (m > M) || (a < -A) || (a > A) ) {
    cout <<"Out of bounds ("<< m <<","<< a <<") in ("<< M <<","<< A <<")\n";
  }
  return data[(2*A+1)*(M+m) + (A+a)];
}

inline
int BiArray::operator()(int m, int a) const {
  if ( (m < - M) || (m > M) || (a < -A) || (a > A) ) {
    cout <<"Out of bounds ("<< m <<","<< a <<") in ("<< M <<","<< A <<")\n";
  }
  return data[(2*A+1)*(M+m) + (A+a)];
}

Link::Link(int gridsize, 
	   int *black, 
	   int *white, 
	   int (*Progress)(const char *, int),
	   bool quiet) :
  gridsize(gridsize),
  black(black),
  white(white),
  Progress(Progress),
  quiet(quiet){
    MOS_Array = NULL;
    HFK_Array = NULL;
    Graph.resize(4*gridsize);
    Graph2.resize(4*gridsize);
    ReducedGraph2.resize(4*gridsize);
    Generators.resize(4*gridsize);
    Ashift = AlexanderShift();
    BuildRectangles();
    counter = new short [gridsize-1];
    g = new int [gridsize];
    gij = new int [gridsize];
    g2 = new int [gridsize];
    aborted = ( gridsize > 16 || BuildVertices() || FindHomology() );
    if (aborted )return;
    MOS_Array = new BiArray(MOS.maxM, MOS.maxA);
    ComputeMOSRanks();
    HFK_Array = new BiArray(MOS.maxM+gridsize, MOS.maxA+gridsize);
    ComputeHFKRanks();
  };

Link::~Link(){
  if (MOS_Array) delete MOS_Array;
  if (HFK_Array) delete HFK_Array;
  delete [] counter;
  delete [] g;
  delete [] gij;
  delete [] g2;
};

// Compute winding number around (x,y)
int Link::WindingNumber(int x, int y)
{
  int i, ret=0;
  for(i=0; i<x; i++) {
    if ((black[i] >= y) && (white[i] < y)) ret++;
    if ((white[i] >= y) && (black[i] < y)) ret--;
  }
  return ret;
}

// Compute the Alexander grading shift.
int Link::AlexanderShift(){
  int i, shift, WN=0;
  for(i=0; i<gridsize; i++) {
    WN += WindingNumber(i, black[i]);
    WN += WindingNumber(i, (black[i]+1)%gridsize);
    WN += WindingNumber((i+1)%gridsize, black[i]);
    WN += WindingNumber((i+1)%gridsize, (black[i]+1)%gridsize);
  } 
  for(i=0; i<gridsize; i++) {
    WN += WindingNumber(i , white[i]);
    WN += WindingNumber(i ,(white[i]+1)%gridsize);
    WN += WindingNumber((i+1)%gridsize, white[i]);
    WN += WindingNumber((i+1)%gridsize, (white[i]+1)%gridsize);
  }
  /*
  if ( WN > 0 ) {
    int *temp;
    temp = black;
    black = white;
    white = temp;
    WN = -WN;
  }
  */
  shift = (WN - 4 * gridsize + 4)/8;
  if (!quiet)
    cout << "Alexander Grading Shift: " << shift << "\n";
  return shift;
}

// Returns 1 if the given rectangle has no dot; 0 if it has a dot;
// or -1 on error.
int Link::RectDotFree(int xll, int yll, int xur, int yur, int which){
  int dotfree = 1;
  switch (which) {
  case 0: 
    for(int x=xll; x<xur && dotfree; x++) {
      if (white[x] >= yll && white[x] < yur) dotfree = 0;
      if (black[x] >= yll && black[x] < yur) dotfree = 0;
    }
    return dotfree;
  case 1:
    for(int x=0; x<xll && dotfree; x++) {
      if (white[x] < yll || white[x] >= yur) dotfree = 0;
      if (black[x] < yll || black[x] >= yur) dotfree = 0;
    }
    for(int x=xur; x<gridsize && dotfree; x++) {
      if (white[x] < yll || white[x] >= yur) dotfree = 0;
      if (black[x] < yll || black[x] >= yur) dotfree = 0;
    }
    return dotfree;
  case 2:
    for(int x=xll; x<xur && dotfree; x++) {
      if (white[x] < yll || white[x] >= yur) dotfree = 0;
      if (black[x] < yll || black[x] >= yur) dotfree = 0;
    }
    return dotfree;
  case 3:
    for(int x=0; x<xll && dotfree; x++) {
      if (white[x] >= yll && white[x] < yur) dotfree = 0;
      if (black[x] >= yll && black[x] < yur) dotfree = 0;
    }
    for(int x=xur; x<gridsize && dotfree; x++) {
      if (white[x] >= yll && white[x] < yur) dotfree = 0;
      if (black[x] >= yll && black[x] < yur) dotfree = 0;
    }
    return dotfree;
  }
  return -1; //Error!
}

 // Compute Maslov grading via the formula:
 // 4M(y)
 // = 4M(white)+4P_y(R_{y, white})+4P_{white}(R_{y, white})-8W(R_{y, white})
 // = 4-4*gridsize+4P_y(R_{y, white})+4P_{x_0}(R_{y, white})-8W[j](R_{y, white})

int Link::MaslovGrading(int g[]) {
  int i, j, k, wi, Wi, Wj, gi, Gi, Gj;
  int P = 4 - 4*gridsize; // Four times the Maslov grading of x_0
  for(i=0; i<gridsize; i++) {
    Wi = white[i];
    wi = Wi - 1;
    if (wi < 0) wi += gridsize;
    Gi = g[i];
    gi = Gi - 1;
    if (gi < 0) gi += gridsize;

    // Calculate incidence number R_{y x_0}.S for each of the four
    // squares S having (i,white[i]) as a corner and each of the 
    // four squares having (i,g[i]) as a corner and shift P appropriately


    for(j=0; j<=i; j++) {
      Wj = white[j];
      Gj = g[j];
      // Squares whose BL corners are (i,white[i]) and (i,g[i])
      // multiply by 7 because of the -8W(R_{y, white}) contribution
      P -= 7*((Wj > Wi) & (Gj <= Wi));
      P += 7*((Gj > Wi) & (Wj <= Wi)); 
      P +=   ((Wj > Gi) & (Gj <= Gi));
      P -=   ((Gj > Gi) & (Wj <= Gi)); 
      // Squares whose TL corners are (i,white[i]) and (i,g[i]) (mod gridsize)
      P += ((Wj > wi) & (Gj <= wi));
      P -= ((Gj > wi) & (Wj <= wi)); 
      P += ((Wj > gi) & (Gj <= gi));
      P -= ((Gj > gi) & (Wj <= gi)); 
    }
    k = i-1;
    if (k < 0) k += gridsize;
    // Squares whose BR corners are ...
    for(j=0; j<= k; j++) {
      Wj = white[j];
      Gj = g[j];
      P += ((Wj > Wi) & (Gj <= Wi));
      P -= ((Gj > Wi) & (Wj <= Wi)); 
      P += ((Wj > Gi) & (Gj <= Gi));
      P -= ((Gj > Gi) & (Wj <= Gi)); 
    // Squares whose TR corners are ...
      P += ((Wj > wi) & (Gj <= wi));
      P -= ((Gj > wi) & (Wj <= wi)); 
      P += ((Wj > gi) & (Gj <= gi));
      P -= ((Gj > gi) & (Wj <= gi)); 
    }
  }
  return (P/4);
}

// Add a vertex for each permutation in every Alexander grading.
int Link::BuildVertices() {
  int i, AGrading, percent;
  Int64 count = 0, end = 0, step;
  stringstream msg;
  int winding_numbers[16][16];

  msg << "Constructing generators ... ";

  // Build a table of winding numbers.
  for(int x=0; x < gridsize; x++) {
    for(int y=0; y < gridsize; y++) {
      winding_numbers[x][y] = WindingNumber(x, y);
    }
  }
  if (!quiet) {
    cout << "Matrix of winding numbers:\n";
    for(int y=gridsize - 1; y>=0; y--) {
      for(int x=0; x<gridsize; x++) {
	cout << setw(2) << winding_numbers[x][y];
      }
      cout << "\n";
    }
  }

  for(i=0; i < gridsize-1; i++)
    counter[i] = 0;
  
  step = 100000 < Factorial[gridsize] ? 100000 : Factorial[gridsize];
  Progress(msg.str().c_str(), 0);
  while (end < Factorial[gridsize]) {
    end = count + step;
    if (end > Factorial[gridsize])
      end = Factorial[gridsize];
    percent = (int)(100.0*(count/(float)Factorial[gridsize]));
    if ( Progress(msg.str().c_str(), max(percent,1)) )
      return -1;
     //This loop accounts for most of the computation.
    for(; count < end; count++) {
      NextPerm(counter, g, gridsize);
      AGrading = Ashift;
      for(i=0; i<gridsize; i++)
	AGrading -= winding_numbers[i][g[i]];
      // Keep generators in all Alexander gradings.  The old version
      // kept only A >= 0 and later recovered A < 0 by symmetry; for the
      // equivariant/Edges2 computation we need the actual negative-A
      // generators present in the complex.
      Generators[gridsize + gridsize + MaslovGrading(g)].push_back(count);
    }
  }
  Progress(msg.str().c_str(), 100);
  msg.str("");
  graphsize = 0;
  for(i=0; i < 4*gridsize; i++)
    graphsize += Generators[i].size();
  msg << "Number of generators: " << graphsize << "\n";
  Progress(msg.str().c_str(), -1);
  return 0;
}

void Link::BuildRectangles(){
  // Build a table showing whether a rectangle on the torus contains a dot.
  for(int xll=0; xll < gridsize; xll++) {
    for(int xur=xll+1; xur < gridsize; xur++) {
      for(int yll=0; yll < gridsize; yll++) {
	for(int yur=yll+1; yur < gridsize; yur++) {
	  Rectangles[xll<<12 | yll<<8 | xur<<4 | yur] = 
	    RectDotFree(xll,yll,xur,yur,0)    |
	    RectDotFree(xll,yll,xur,yur,1)<<1 |
	    RectDotFree(xll,yll,xur,yur,2)<<2 |
	    RectDotFree(xll,yll,xur,yur,3)<<3;
	}
      }
    }
  }
}

// Add edges describing the boundary map.
int Link::BuildEdges(int MM, const char *msg, int *count) {
  int i, j, k, Gi, Gj, indexgij, rectinfo; 
  int index = 0, edges = 0, end = 0, step, percent;
  bool firstrect = 0;
  bool secondrect = 0;

  step = min(20000, (int)Graph[MM].size());

  // Add the edges.
  while (end < (int)Graph[MM].size()) {
    end = min(index + step, (int)Graph[MM].size());
    percent = (int)((*count + index)/(1 + float(graphsize)/50));
    if ( Progress(msg, max(1, percent)) )
      return -1;
    for(; index < end; index++) {
      Int2Perm(Graph[MM][index].perm, g, gridsize);
      for(i=0; i<gridsize; i++) {
	Gi = g[i];
	for(j=i+1; j<gridsize; j++) {
	  Gj = g[j];
	  if(Gi < Gj) {
	    rectinfo = Rectangles[i<<12 | Gi<<8 | j<<4 | Gj];
	    firstrect = (bool)(rectinfo & 1);
	    for(k=i+1; k<j && firstrect; k++) {
	      firstrect = !(bool)((Gi < g[k]) & (g[k] < Gj));
	    }
	    secondrect = (bool)(rectinfo & 2);
	    for(k=0; k<i && secondrect; k++) {
	      secondrect = !(bool)((g[k] < Gi) | (g[k] > Gj));
	    }
	    for(k=j+1; k<gridsize && secondrect; k++) {
	      secondrect = !(bool)((g[k] < Gi) | (g[k] > Gj));
	    }
	  }
	  if(Gj < Gi) {
	    rectinfo = Rectangles[i<<12 | Gj<<8 | j<<4 | Gi];
	    firstrect = rectinfo & 4;
	    for(k=i+1; k<j && firstrect; k++) {
	      firstrect = !((g[k] < Gj) | (g[k] > Gi));
	    }
	    secondrect = rectinfo & 8;
	    for(k=0; k<i && secondrect; k++) {
	      secondrect = !((g[k] > Gj) & (g[k] < Gi));
	    }
	    for(k=j+1; k<gridsize && secondrect; k++) {
	      secondrect = !((g[k] > Gj) & (g[k] < Gi));
	    }
	  }
	  if(firstrect != secondrect) { // Exactly one rectangle is a boundary
	    for(k=0; k<gridsize; k++) {
	      gij[k] = g[k];
	    }
	    gij[i] = Gj;
	    gij[j] = Gi;
	    indexgij = Find(Graph[MM-1], Perm2Int(gij, gridsize));
	    if(indexgij==-1) {
	      Progress("Error with Alexander grading!!", -1);
	      return -1;
	    }
	    Graph[MM][index].out.push_back( indexgij );
	    Graph[MM-1][indexgij].in.push_back( index );     
	    edges++;
	  }
	}
      }
    }
  }
  *count += end;
  return 0;
}

// Translate a grid state by (gridsize/2, gridsize/2) modulo gridsize.
// If src[x] is the selected y-coordinate in column x, then dest is the
// permutation representing the translated state.
void Link::ShiftHalfTurn(int src[], int dest[]) {
  int h = gridsize/2;
  for(int x=0; x<gridsize; x++) {
    dest[(x+h)%gridsize] = (src[x]+h)%gridsize;
  }
}

// Toggle one Edges2 arrow in Graph2, keeping out and in synchronized.
void Link::ToggleEdge2(int level, int from, int to) {
  list<int>::iterator search;

  for(search = Graph2[level][from].out.begin();
      search != Graph2[level][from].out.end(); search++) {
    if (*search == to) {
      Graph2[level][from].out.erase(search);
      Graph2[level][to].in.remove(from);
      return;
    }
  }

  Graph2[level][from].out.push_back(to);
  Graph2[level][to].in.push_back(from);
}

// Add column src to column dest in the Edges2 matrix for a fixed level.
// In terms of the operator T represented by Graph2, this implements
// T(dest) <- T(dest) + T(src) over F_2.
void Link::AddColumn2(int level, int dest, int src) {
  vector<int> terms(Graph2[level][src].out.begin(),
                    Graph2[level][src].out.end());

  for(int i=0; i<(int)terms.size(); i++) {
    ToggleEdge2(level, dest, terms[i]);
  }
}

// Add row src to row dest in the Edges2 matrix for a fixed level.
// Equivalently, every vertex that points to src is toggled so that it
// also points to dest.
void Link::AddRow2(int level, int dest, int src) {
  vector<int> terms(Graph2[level][src].in.begin(),
                    Graph2[level][src].in.end());

  for(int i=0; i<(int)terms.size(); i++) {
    ToggleEdge2(level, terms[i], dest);
  }
}

// Remove every Edges2 arrow into or out of v.  This is used when the
// ordinary reduction kills a vertex, so no surviving vertex keeps a stale
// Edges2 reference to it.
void Link::RemoveAllEdges2(int level, int v) {
  vector<int> outs(Graph2[level][v].out.begin(),
                   Graph2[level][v].out.end());
  vector<int> ins(Graph2[level][v].in.begin(),
                  Graph2[level][v].in.end());

  for(int i=0; i<(int)outs.size(); i++) {
    Graph2[level][outs[i]].in.remove(v);
  }

  for(int i=0; i<(int)ins.size(); i++) {
    Graph2[level][ins[i]].out.remove(v);
  }

  Graph2[level][v].out.clear();
  Graph2[level][v].in.clear();
}

// Store the surviving vertices in this grading level, together with the
// induced Edges2 arrows between the survivors.  The stored arrow indices
// are reindexed to the reduced graph for this level, not the temporary
// Graph/Graph2 indices.
void Link::StoreReducedEdges2(int level) {
  vector<int> oldToNew(Graph[level].size(), -1);

  ReducedGraph2[level].clear();

  for(int i=0; i<(int)Graph[level].size(); i++) {
    if(Graph[level][i].alive) {
      Int2Perm(Graph[level][i].perm, g, gridsize);
      int AGrading = Ashift;
      for(int j=0; j<gridsize; j++)
        AGrading -= WindingNumber(j, g[j]);
      oldToNew[i] = ReducedGraph2[level].size();
      ReducedGraph2[level].push_back(ReducedVertex2(Graph[level][i].perm,
                                                   MaslovGrading(g),
                                                   AGrading));
    }
  }

  for(int i=0; i<(int)Graph2[level].size(); i++) {
    if(!Graph2[level][i].alive) continue;
    if(i >= (int)oldToNew.size() || oldToNew[i] == -1) continue;
    int from = oldToNew[i];
    for(list<int>::iterator it = Graph2[level][i].out.begin();
        it != Graph2[level][i].out.end(); it++) {
      if((*it >= 0) && (*it < (int)oldToNew.size()) && oldToNew[*it] != -1) {
        ReducedGraph2[level][from].out.push_back(oldToNew[*it]);
      }
    }
  }
}

// Add Edges2 arrows describing the half-grid translation operator on one
// fixed Maslov grading level.  For the grids considered here, this operator
// is assumed to preserve both Maslov and Alexander gradings.
int Link::BuildEdges2(int MM, const char *msg, int *count) {
  int index = 0, end = 0, step, percent;
  int target;

  if (gridsize % 2 != 0) {
    Progress("Error building Edges2: gridsize must be even!!", -1);
    return -1;
  }

  if (Graph[MM].size() != Graph2[MM].size()) {
    Progress("Error building Edges2: Graph and Graph2 sizes differ!!", -1);
    return -1;
  }

  step = min(20000, (int)Graph[MM].size());

  while (end < (int)Graph[MM].size()) {
    end = min(index + step, (int)Graph[MM].size());
    percent = (int)((*count + index)/(1 + float(graphsize)/50));
    if ( Progress(msg, max(1, percent)) )
      return -1;

    for(; index < end; index++) {
      Int2Perm(Graph[MM][index].perm, g, gridsize);
      ShiftHalfTurn(g, g2);
      target = Find(Graph[MM], Perm2Int(g2, gridsize));

      if(target == -1) {
        Progress("Error building Edges2: shifted generator not found!!", -1);
        return -1;
      }

      Graph2[MM][index].out.push_back(target);
      Graph2[MM][target].in.push_back(index);
    }
  }

  return 0;
}
  
// Eliminate edges to find a homology basis.
int Link::Reduce(int MM, const char *msg, int *count){
  int source = 0, target;
  int end = 0, step, percent;

  step = min(20000, (int)Graph[MM].size());

  while (end < (int)Graph[MM].size()) {
    end = min(source + step, (int) Graph[MM].size());
    percent = (int)((*count + source)/(1 + float(graphsize)/50));
    if ( Progress(msg, max(1,percent)) )
      return -1;
    for(; source<end; source++) {
      if( (!Graph[MM][source].alive) || Graph[MM][source].out.size()==0)
	continue;

      target = Graph[MM][source].out.front();
      Graph[MM][source].alive = 0;
      Graph2[MM][source].alive = 0;

      // Store the rows/columns involved in this cancellation before the
      // ordinary differential edges are modified.
      vector<int> Js;
      vector<int> Ks;
      for(list<int>::iterator jj = Graph[MM-1][target].in.begin();
          jj != Graph[MM-1][target].in.end(); jj++) {
        if((*jj != source) && Graph[MM][*jj].alive) Js.push_back(*jj);
      }
      for(list<int>::iterator kk = Graph[MM][source].out.begin();
          kk != Graph[MM][source].out.end(); kk++) {
        Ks.push_back(*kk);
      }

      // Edges2 update in degree MM: each surviving j with j->target is
      // changed to j+source, so T(j) gets toggled by T(source).
      for(int jj=0; jj<(int)Js.size(); jj++) {
        AddColumn2(MM, Js[jj], source);
      }

      // Edges2 update in degree MM-1: since target will be killed and
      // source points to target plus the other k's, every Edges2 arrow
      // into target is replaced by arrows into those other k's.
      for(int kk=0; kk<(int)Ks.size(); kk++) {
        if(Ks[kk] != target) AddRow2(MM-1, Ks[kk], target);
      }
      
      // For every j->target, j != source
      for(list<int>::iterator j = Graph[MM-1][target].in.begin();
	  j != Graph[MM-1][target].in.end(); j++) {
	    if( !Graph[MM][*j].alive ) continue;
	    // For every source->k
	    for(list<int>::iterator k = Graph[MM][source].out.begin();
		k != Graph[MM][source].out.end(); k++) {
	      // Search for j->k. If found, remove it; if not, add it.
	      list<int>::iterator search;
	      for(search = Graph[MM][*j].out.begin();
		  search != Graph[MM][*j].out.end(); search++ )
		if (*search == *k) {
		  Graph[MM][*j].out.erase( search );
		  if( *k != target) Graph[MM-1][*k].in.remove(*j);
		  break;
		}
	      if (search == Graph[MM][*j].out.end()){
	        Graph[MM][*j].out.push_back(*k);
	        Graph[MM-1][*k].in.push_back(*j);
	      } 
	    }
	  }
      
      // For each source->j, remove source from j.in
      for(list<int>::iterator j = Graph[MM][source].out.begin();
	  j != Graph[MM][source].out.end(); j++)
	Graph[MM-1][*j].in.remove(source);
      
      // Remove all Edges2 arrows incident to the killed vertices.
      RemoveAllEdges2(MM, source);
      RemoveAllEdges2(MM-1, target);

      //Kill source and target
      Graph[MM-1][target].alive = 0;
      Graph2[MM-1][target].alive = 0;
      Graph[MM-1][target].out.clear();
      Graph[MM-1][target].in.clear();
      Graph[MM][source].out.clear();
      Graph[MM][source].in.clear();
    }
  }
  *count += end;
  return 0;
}

int Link::FindHomology() {
  int i, j, count=0, MM;
  const char *msg="Computing homology ...";

  Progress(msg, 0);

  MM = 4*gridsize;
  for(i=0; i < (int)Generators[MM-1].size(); i++) {
    Graph[MM-1].push_back(Vertex(Generators[MM-1][i]));
    Graph2[MM-1].push_back(Vertex2());
  }
  if (BuildEdges2(MM-1, msg, &count)) return -1;

  for (MM = 4*gridsize -1; MM > 0; MM--) {
    for(i=0; i < (int)Generators[MM-1].size(); i++) {
      Graph[MM-1].push_back(Vertex( Generators[MM-1][i]));
      Graph2[MM-1].push_back(Vertex2());
    }
    Generators[MM-1].clear();

    if (BuildEdges2(MM-1, msg, &count)) return -1;

    aborted = BuildEdges(MM, msg, &count) || Reduce(MM, msg, &count);
    if ( aborted ) return aborted;

    StoreReducedEdges2(MM);
    for(i=0; i < (int)ReducedGraph2[MM].size(); i++) {
      MOS.add(Generator(ReducedGraph2[MM][i].m, ReducedGraph2[MM][i].a));
    }

    Graph[MM].clear();
    Graph2[MM].clear();
  }

  Progress(msg, 100);
  Progress("", -1);
  return 0;
}

// Print the induced Edges2 operator on the reduced graph.
// Vertices are grouped by Maslov grading.  The displayed vertex labels are
// reduced-basis labels; the original grid-state tuples are intentionally
// suppressed.
void Link::PrintReducedEdges2(ostream &os) {
  os << "\nEdges2:\n";

  for(int level=0; level<(int)ReducedGraph2.size(); level++) {
    if(ReducedGraph2[level].size()==0) continue;

    os << "\nMaslov " << level - 2*gridsize << "\n";

    for(int i=0; i<(int)ReducedGraph2[level].size(); i++) {
      os << "v" << i
         << " M=" << ReducedGraph2[level][i].m
         << " A=" << ReducedGraph2[level][i].a
         << " -> ";

      if(ReducedGraph2[level][i].out.size()==0) {
        os << "{}";
      } else {
        os << "{";
        for(int j=0; j<(int)ReducedGraph2[level][i].out.size(); j++) {
          if(j) os << ",";
          os << "v" << ReducedGraph2[level][i].out[j];
        }
        os << "}";
      }
      os << "\n";
    }
  }
}


// Compute the rank over F_2 of a matrix represented by its rows.
// The bit in column c of row r is stored in rows[r][c/64].
int Link::RankMod2(vector < vector <unsigned long long> > &rows, int ncols) {
  int nrows = rows.size();
  int rank = 0;

  for(int col=0; col<ncols && rank<nrows; col++) {
    int word = col/64;
    unsigned long long mask = 1ULL << (col%64);
    int pivot = -1;

    for(int r=rank; r<nrows; r++) {
      if(rows[r][word] & mask) {
        pivot = r;
        break;
      }
    }

    if(pivot == -1) continue;

    if(pivot != rank) {
      vector<unsigned long long> temp = rows[pivot];
      rows[pivot] = rows[rank];
      rows[rank] = temp;
    }

    for(int r=0; r<nrows; r++) {
      if(r == rank) continue;
      if(rows[r][word] & mask) {
        for(int w=0; w<(int)rows[r].size(); w++)
          rows[r][w] ^= rows[rank][w];
      }
    }

    rank++;
  }

  return rank;
}

// Print the ranks of the homology of the reduced MOS complex with respect
// to the differential D = 1 + Edges2. Everything is over F_2.
void Link::PrintEdges2HomologyRanks(ostream &os) {
  os << "\nE2 page in all Alexander gradings:\n\n";

  for(int level=0; level<(int)ReducedGraph2.size(); level++) {
    int n = ReducedGraph2[level].size();
    if(n == 0) continue;

    vector<int> used(n, 0);

    for(int start=0; start<n; start++) {
      if(used[start]) continue;

      int m = ReducedGraph2[level][start].m;
      int a = ReducedGraph2[level][start].a;

      vector<int> oldToLocal(n, -1);
      vector<int> localToOld;

      for(int i=0; i<n; i++) {
        if(ReducedGraph2[level][i].m == m && ReducedGraph2[level][i].a == a) {
          oldToLocal[i] = localToOld.size();
          localToOld.push_back(i);
          used[i] = 1;
        }
      }

      int dim = localToOld.size();
      int words = (dim + 63)/64;
      vector < vector <unsigned long long> > rows(dim, vector<unsigned long long>(words, 0ULL));

      // Build the matrix of D = I + T.
      // Columns are source generators, rows are target generators.
      // The printed out-list gives T(source).
      for(int localSource=0; localSource<dim; localSource++) {
        int source = localToOld[localSource];

        // Identity term in I+T.
        rows[localSource][localSource/64] ^= (1ULL << (localSource%64));

        // T term.
        for(int e=0; e<(int)ReducedGraph2[level][source].out.size(); e++) {
          int target = ReducedGraph2[level][source].out[e];

          if(target < 0 || target >= n) {
            os << "Warning: Edges2 target out of range at Maslov "
               << level - 2*gridsize << "\n";
            continue;
          }

          int localTarget = oldToLocal[target];

          if(localTarget == -1) {
            os << "Warning: Edges2 does not preserve the (M,A) block at Maslov "
               << level - 2*gridsize
               << ", source v" << source << " -> v" << target << "\n";
            continue;
          }

          rows[localTarget][localSource/64] ^= (1ULL << (localSource%64));
        }
      }

      int rankD = RankMod2(rows, dim);
      int homRank = dim - 2*rankD;

      os << " M=" << m
         << " A=" << a
         << " rank=" << homRank
         << "\n";
    }
  }
}

// Compute MOS ranks.
void Link::ComputeMOSRanks(){
  for (list<Generator>::iterator 
	 G=MOS.gens.begin(); G != MOS.gens.end(); G++)
    (*MOS_Array)((*G).m, (*G).a)++;
  if (!quiet) {
    cout << "MOS ranks for all Alexander gradings:\n";
    for(int a = MOS.maxA; a >= -MOS.maxA; a--) {
      for(int m = -MOS.maxM; m <= MOS.maxM; m++) {
	cout << setw(5) << (*MOS_Array)(m,a);
      }
      cout << "\n";
    }
  }
}

// Compute HFK^ ranks from MOS ranks.
void Link::ComputeHFKRanks(){
  int a, m, i;

  // MOS has now been computed in every Alexander grading, so initialize
  // HFK_Array from all MOS ranks rather than only from A >= 0.
  for(a = -MOS.maxA; a <= MOS.maxA; a++)
    for(m = -MOS.maxM; m <= MOS.maxM; m++)
      (*HFK_Array)(m,a) = MOS_Rank(m,a);

  // MOS =  \hat HFK \otimes K^{gridsize-1}.
  for(a = MOS.maxA + gridsize ; a >= -MOS.maxA - gridsize; a--) {
    for(m = MOS.maxM + gridsize; m >= -MOS.maxM - gridsize; m--) {
      if( (*HFK_Array)(m,a) > 0) {
	for(i=1; i <= gridsize-1; i++)
	  (*HFK_Array)(m-i,a-i) -= (*HFK_Array)(m,a)*Binomial(gridsize-1,i);
      }
    }
  }

  HFK_maxA = MOS.maxA;
  HFK_minM = 0;
  HFK_maxM = 0;

  for (m = -MOS.maxM - gridsize; m <= MOS.maxM + gridsize; m++) {
    for (a = -MOS.maxA; a <= MOS.maxA; a++) {
      if ( (*HFK_Array)(m,a) ) {
        HFK_minM = m;
        goto found_min_m;
      }
    }
  }
found_min_m:

  for (m = MOS.maxM + gridsize; m >= -MOS.maxM - gridsize; m--) {
    for (a = -MOS.maxA; a <= MOS.maxA; a++) {
      if ( (*HFK_Array)(m,a) ) {
        HFK_maxM = m;
        goto found_max_m;
      }
    }
  }
found_max_m:
  ;
}

int Link::MOS_Rank(int m, int a) {
  return (*MOS_Array)(m,a);
}

int Link::HFK_Rank(int m, int a) {
  return (*HFK_Array)(m,a);
}

// Returns i if V[i]=x or -1 if x isn't a member of V.
int Find(vector <Vertex> &V, Int64 x) {
  int above=V.size()-1;
  int below=0;
  while(above - below > 1) {
    if (x >= V[below+(above-below)/2].perm) below += (above-below)/2;
    else above = below+(above-below)/2;
  }
  if (V[below].perm == x) return below;
  if (V[above].perm == x) return above;
  return -1;
}

// Check that we have exactly 2 dots of each color in each row and column.
bool ValidGrid(int gridsize, int black[], int white[]) {
 int numwhite=0;
 int numblack=0;
 for(int i=0; i<gridsize; i++) {
  for(int j=0; j<gridsize; j++) {
   if (white[j]==i) numwhite++;
   if (black[j]==i) numblack++;
  }
  if (numwhite != 1 || numblack != 1)
   return 0;
  numwhite=0;
  numblack=0;
 }
 return 1;
}

// Maps a permutation of size n to an integer < n!
// See: Knuth, Volume 2, Section 3.3.2, Algorithm P
Int64 Perm2Int(int Perm[], int size) {
 Int64 Int=0;
 int temp, m, r = size;
 while (r > 0){
   for (m=0; m < r; m++){
     if (r - Perm[m] == 1)
       break;
   }
   Int = Int*r + m;
   r -= 1;
   temp = Perm[r];
   Perm[r] = Perm[m];
   Perm[m] = temp;
 }
 return Int;
}

// Inverse mapping, from integers < n! to permutations of size n
// Writes the permutation corresponding to Int into the array Perm.
void Int2Perm(Int64 Int, int Perm[], int size) {
  int r, m, temp;
  for(int i=0; i<size; i++) Perm[i]=i;
  r = 1;
  while (r < size)  {
    m = Int%(r+1);
    Int = Int/(r+1);
    temp = Perm[r];
    Perm[r] = Perm[m];
    Perm[m] = temp;
    r += 1;
  }
  return;
}

// Generator for permutations.  Inputs a factorial based counter and
// an array.  Writes the permutation indexed by the counter into the
// array, and then increments the counter.
void NextPerm(short counter[], int P[], int size) {
  int i, m, temp;
  for(i=0; i<size; i++) P[i]=i;
  for (i=1; i < size; i++) {
    m = counter[i-1];
    temp = P[i];
    P[i] = P[m];
    P[m] = temp;
  }
  for (i=0; i<size-1; i++) {
    counter[i] += 1;
    if (counter[i] == i+2)
      counter[i] = 0;
    else
      break;
  }
  return;
}
