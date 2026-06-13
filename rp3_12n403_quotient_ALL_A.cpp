#include <algorithm>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

using namespace std;

struct Gen {
  int permRank;
  int pcoordCode;
  vector<int> perm;     // residues 0,...,n-1, one for each alpha-circle
  vector<int> pcoord;   // in {0,...,p-1}
  vector<int> abscoord; // length p*n, periodic extension of the first n entries
  int spin;
  int M4;               // four times Maslov grading
  int A4;               // four times Alexander grading
};

struct MatrixKey {
  int spin;
  int A4;
  int M4; // source Maslov grading, differential goes M4 -> M4-4
  bool operator<(const MatrixKey &o) const {
    if (spin != o.spin) return spin < o.spin;
    if (A4 != o.A4) return A4 < o.A4;
    return M4 < o.M4;
  }
};

struct GroupKey {
  int spin;
  int A4;
  int M4;
  bool operator<(const GroupKey &o) const {
    if (spin != o.spin) return spin < o.spin;
    if (A4 != o.A4) return A4 < o.A4;
    return M4 < o.M4;
  }
};

static long long factorial_int(int n) {
  long long r = 1;
  for (int i = 2; i <= n; ++i) r *= i;
  return r;
}

static int mod_int(int a, int m) {
  int r = a % m;
  if (r < 0) r += m;
  return r;
}

static int gcd_int(int a, int b) {
  if (a < 0) a = -a;
  if (b < 0) b = -b;
  while (b) { int r = a % b; a = b; b = r; }
  return a ? a : 1;
}

static string grade4_to_string(int x4) {
  if (x4 % 4 == 0) return to_string(x4 / 4);
  string sign = "";
  int a = x4;
  if (a < 0) { sign = "-"; a = -a; }
  int g = gcd_int(a, 4);
  return sign + to_string(a / g) + "/" + to_string(4 / g);
}

static long long perm_code(const vector<int> &perm, int n) {
  long long code = 0;
  long long pow = 1;
  for (int i = 0; i < n; ++i) {
    code += pow * perm[i];
    pow *= n;
  }
  return code;
}

// Celoria's IAB function, specialized to integer arrays.  The half-integer
// comparisons are rewritten using integer inequalities.
static int IAB(const vector<int> &A, const vector<int> &B,
               int n, int p, int q, int paramA, int paramB) {
  (void)n; (void)q;
  int N = (int)A.size(); // should be p*n
  int count = 0;

  if (((paramA + paramB) % 2) == 0) {
    for (int i = 0; i < N; ++i)
      for (int j = 0; j < N; ++j)
        if (A[i] < B[j] && i < j) count++;
  } else if (paramA == 0 && paramB == 1) {
    for (int i = 0; i < N; ++i)
      for (int j = 0; j < N; ++j)
        if (A[i] <= B[j] && i <= j) count++;
  } else if (paramA == 1 && paramB == 0) {
    for (int i = 0; i < N; ++i)
      for (int j = 0; j < N; ++j)
        if (A[i] < B[j] && i < j) count++;
  }

  return count;
}

static int component_count_lens(const vector<int> &X, const vector<int> &O,
                                int n, int p, int q) {
  (void)p; (void)q;
  vector<int> A(n), B(n), Binv(n);
  int equal_count = 0;
  for (int i = 0; i < n; ++i) {
    A[i] = mod_int(X[i], n);
    B[i] = mod_int(O[i], n);
    Binv[B[i]] = i;
    if (A[i] == B[i]) equal_count++;
  }

  // VAR = PO * PX^{-1}; in 0-index notation this is A o B^{-1}.
  vector<int> f(n);
  for (int r = 0; r < n; ++r) f[r] = A[Binv[r]];

  vector<int> seen(n, 0);
  int cycles = 0;
  for (int i = 0; i < n; ++i) {
    if (seen[i]) continue;
    cycles++;
    int x = i;
    while (!seen[x]) {
      seen[x] = 1;
      x = f[x];
    }
  }
  return cycles + equal_count;
}

static vector<int> extend_periodic(const vector<int> &first, int n, int p, int q) {
  vector<int> out(p * n);
  for (int i = 0; i < n; ++i) out[i] = first[i];
  for (int i = n; i < p * n; ++i) {
    out[i] = mod_int(out[i - n] + q * n, p * n);
  }
  return out;
}

static vector<int> abs_from_perm_pcoord(const vector<int> &perm,
                                        const vector<int> &pcoord,
                                        int n, int p, int q) {
  vector<int> abscoord(p * n);
  for (int i = 0; i < n; ++i) abscoord[i] = n * pcoord[i] + perm[i];
  for (int i = n; i < p * n; ++i) {
    abscoord[i] = mod_int(abscoord[i - n] + q * n, p * n);
  }
  return abscoord;
}

// Celoria's controllarett.  Returns 1 if the corresponding parallelogram is
// empty of X/O markings and generator points; otherwise returns 0.
static int empty_parallelogram(const vector<int> &A, const vector<int> &B,
                               int n, int p, int q,
                               const vector<int> &O, const vector<int> &X,
                               int parametro) {
  (void)q;
  int N = n * p;
  int base[2] = {0, 0};
  int altezza[2] = {0, 0};
  int trasposizione[2] = {0, 0};
  int flag = 0;

  for (int i = 0; i < n; ++i) {
    if (A[i] != B[i]) {
      if (flag >= 2) return 0;
      if (A[i] < B[i]) base[flag] = B[i] - A[i];
      if (A[i] > B[i]) base[flag] = N - A[i] + B[i];
      trasposizione[flag] = i;
      flag++;
    }
  }
  if (flag != 2) return 0;

  for (int j = 0; j < N; ++j) {
    if (B[j] == A[trasposizione[0]]) {
      if (j > trasposizione[0]) altezza[0] = j - trasposizione[0];
      if (j < trasposizione[0]) altezza[0] = N + j - trasposizione[0];
    }
    if (B[j] == A[trasposizione[1]]) {
      if (j > trasposizione[1]) altezza[1] = j - trasposizione[1];
      if (j < trasposizione[1]) altezza[1] = N + j - trasposizione[1];
    }
  }

  int r = parametro;
  int flag_inside = 0;
  int flag_dot = 0;
  for (int ii = 0; ii < base[r]; ++ii) {
    for (int jj = 0; jj < altezza[r]; ++jj) {
      int col = mod_int(trasposizione[r] + jj, N);
      int value = mod_int(A[trasposizione[r]] + ii, N);

      if (mod_int(B[col], N) == value) flag_inside = 1;
      if (mod_int(A[col], N) == value && ii != 0) flag_inside = 1;
      if (X[col] == value) flag_dot = 1;
      if (O[col] == value) flag_dot = 1;
    }
  }

  if (flag_dot) return 0;
  if (flag_inside) return 0;
  return 1;
}

static vector<int> xor_sorted_vectors(const vector<int> &a, const vector<int> &b) {
  vector<int> out;
  out.reserve(a.size() + b.size());
  size_t i = 0, j = 0;
  while (i < a.size() || j < b.size()) {
    if (j == b.size() || (i < a.size() && a[i] < b[j])) {
      out.push_back(a[i++]);
    } else if (i == a.size() || b[j] < a[i]) {
      out.push_back(b[j++]);
    } else {
      i++; j++; // same entry cancels over F_2
    }
  }
  return out;
}

static int sparse_rank_mod2(vector<vector<int> > columns) {
  unordered_map<int, vector<int> > pivot;
  int rank = 0;

  for (size_t c = 0; c < columns.size(); ++c) {
    vector<int> v = columns[c];
    sort(v.begin(), v.end());
    v.erase(unique(v.begin(), v.end()), v.end());

    while (!v.empty()) {
      int lead = v.back();
      unordered_map<int, vector<int> >::iterator it = pivot.find(lead);
      if (it == pivot.end()) {
        pivot[lead] = v;
        rank++;
        break;
      }
      v = xor_sorted_vectors(v, it->second);
    }
  }

  return rank;
}

class RP3GridHomologyF2 {
public:
  RP3GridHomologyF2(int n_, int p_, int q_, const vector<int> &X_input,
                    const vector<int> &O_input)
      : n(n_), p(p_), q(q_) {
    if (p != 2 || q != 1) {
      throw runtime_error("This stripped-down calculator is intended for L(2,1)=RP3 only.");
    }
    if ((int)X_input.size() != n || (int)O_input.size() != n) {
      throw runtime_error("X and O must have length n.");
    }
    X = extend_periodic(X_input, n, p, q);
    O = extend_periodic(O_input, n, p, q);
    num_components = component_count_lens(X, O, n, p, q);
    valO = IAB(O, O, n, p, q, 1, 1);
    valX = IAB(X, X, n, p, q, 1, 1);
  }

  void run() {
    generate_permutations();
    generate_generators();
    compute_local_indices();
    build_differential_and_ranks();
    print_output();
  }

private:
  int n, p, q;
  int num_components;
  int valO, valX;
  vector<int> X, O;
  vector<vector<int> > perms;
  unordered_map<long long, int> permRankByCode;
  vector<Gen> gens;
  vector<int> localIndex;
  unordered_map<long long, int> keptIndexByKey;

  map<GroupKey, int> dimByGroup;
  map<MatrixKey, vector<vector<int> > > columnsByMatrix;
  map<MatrixKey, int> rankByMatrix;
  map<GroupKey, int> homRankByGroup;
  map<pair<int,int>, int> rankBySpinAlex;
  map<int, int> rankByAlex;
  int totalRank = 0;

  long long key_from_perm_rank_pcoord(int permRank, int pcoordCode) const {
    return (long long)pcoordCode * (long long)perms.size() + (long long)permRank;
  }

  void generate_permutations() {
    vector<int> perm(n);
    for (int i = 0; i < n; ++i) perm[i] = i;
    do {
      int rank = (int)perms.size();
      permRankByCode[perm_code(perm, n)] = rank;
      perms.push_back(perm);
    } while (next_permutation(perm.begin(), perm.end()));
  }

  void generate_generators() {
    int generatorespin = 0;
    for (int l = 0; l < n; ++l) generatorespin += O[l] / n;

    int ppower = 1;
    for (int i = 0; i < n; ++i) ppower *= p;

    for (int pc = 0; pc < ppower; ++pc) {
      vector<int> pcoord(n, 0);
      int temp = pc;
      int spinprov = 0;
      for (int s = 0; s < n; ++s) {
        pcoord[s] = temp % p;
        temp /= p;
        spinprov += pcoord[s];
      }

      for (int pr = 0; pr < (int)perms.size(); ++pr) {
        Gen g;
        g.permRank = pr;
        g.pcoordCode = pc;
        g.perm = perms[pr];
        g.pcoord = pcoord;
        g.abscoord = abs_from_perm_pcoord(g.perm, g.pcoord, n, p, q);
        g.spin = mod_int(q - 1 + spinprov - generatorespin, p);

        int maslov_num = IAB(g.abscoord, g.abscoord, n, p, q, 0, 0)
                       - IAB(g.abscoord, O, n, p, q, 0, 1)
                       - IAB(O, g.abscoord, n, p, q, 1, 0)
                       + valO + 1;
        // For L(2,1), correction(2,1,0)+(p-1)/p = 1/4.
        g.M4 = 2 * maslov_num + 1;

        int alex_num = valO
                     + IAB(X, g.abscoord, n, p, q, 1, 0)
                     + IAB(g.abscoord, X, n, p, q, 0, 1)
                     - valX
                     - IAB(g.abscoord, O, n, p, q, 0, 1)
                     - IAB(O, g.abscoord, n, p, q, 1, 0);
        // For p=2, A = alex_num/4 + (components-n)/2.
        g.A4 = alex_num + 2 * (num_components - n);

        {
          long long key = key_from_perm_rank_pcoord(pr, pc);
          int idx = (int)gens.size();
          gens.push_back(g);
          keptIndexByKey[key] = idx;
          GroupKey group = {g.spin, g.A4, g.M4};
          dimByGroup[group]++;
        }
      }
    }
  }

  void compute_local_indices() {
    localIndex.assign(gens.size(), 0);
    map<GroupKey, int> count;
    for (int i = 0; i < (int)gens.size(); ++i) {
      GroupKey gk = {gens[i].spin, gens[i].A4, gens[i].M4};
      localIndex[i] = count[gk]++;
    }
  }

  void add_column_for_source(const Gen &src, int srcIndex) {
    (void)srcIndex;
    MatrixKey mk = {src.spin, src.A4, src.M4};
    vector<int> col;

    for (int a = 0; a < n; ++a) {
      for (int b = a + 1; b < n; ++b) {
        vector<int> bperm = src.perm;
        swap(bperm[a], bperm[b]);
        for (int pa = 0; pa < p; ++pa) {
          for (int pb = 0; pb < p; ++pb) {
            vector<int> bpcoord = src.pcoord;
            bpcoord[a] = pa;
            bpcoord[b] = pb;
            int bpc = 0;
            int mult = 1;
            for (int i = 0; i < n; ++i) {
              bpc += bpcoord[i] * mult;
              mult *= p;
            }

            long long bcode = perm_code(bperm, n);
            unordered_map<long long, int>::iterator prit = permRankByCode.find(bcode);
            if (prit == permRankByCode.end()) continue;
            long long key = key_from_perm_rank_pcoord(prit->second, bpc);
            unordered_map<long long, int>::iterator kit = keptIndexByKey.find(key);
            if (kit == keptIndexByKey.end()) continue;

            const Gen &tar = gens[kit->second];
            if (tar.spin != src.spin) continue;
            if (tar.A4 != src.A4) continue;
            if (tar.M4 != src.M4 - 4) continue;

            int r0 = empty_parallelogram(src.abscoord, tar.abscoord, n, p, q, O, X, 0);
            int r1 = empty_parallelogram(src.abscoord, tar.abscoord, n, p, q, O, X, 1);
            int coeff = (r0 + r1) & 1;
            if (!coeff) continue;

            GroupKey targetGroup = {tar.spin, tar.A4, tar.M4};
            (void)targetGroup;
            int targetLocal = localIndex[kit->second];
            col.push_back(targetLocal);
          }
        }
      }
    }

    if (!col.empty()) columnsByMatrix[mk].push_back(col);
    else columnsByMatrix[mk].push_back(vector<int>());
  }

  void build_differential_and_ranks() {
    cerr << "Kept generators in all Alexander gradings: " << gens.size() << "\n";

    // Make one column, possibly empty, for each source generator.  Empty columns
    // are harmless and keep the matrix rank calculation honest.
    for (int i = 0; i < (int)gens.size(); ++i) {
      if (i % 50000 == 0) cerr << "Building differential columns: " << i << "/" << gens.size() << "\n";
      add_column_for_source(gens[i], i);
    }

    cerr << "Computing matrix ranks...\n";
    for (map<MatrixKey, vector<vector<int> > >::iterator it = columnsByMatrix.begin();
         it != columnsByMatrix.end(); ++it) {
      rankByMatrix[it->first] = sparse_rank_mod2(it->second);
    }

    for (map<GroupKey, int>::iterator it = dimByGroup.begin(); it != dimByGroup.end(); ++it) {
      GroupKey gk = it->first;
      int dim = it->second;
      MatrixKey d_down = {gk.spin, gk.A4, gk.M4};
      MatrixKey d_up = {gk.spin, gk.A4, gk.M4 + 4};
      int r_down = rankByMatrix.count(d_down) ? rankByMatrix[d_down] : 0;
      int r_up = rankByMatrix.count(d_up) ? rankByMatrix[d_up] : 0;
      int h = dim - r_down - r_up;
      if (h < 0) {
        cerr << "Warning: negative homology rank at spin=" << gk.spin
             << " A4=" << gk.A4 << " M4=" << gk.M4
             << "; dim=" << dim << " r_down=" << r_down << " r_up=" << r_up << "\n";
      }
      if (h != 0) {
        homRankByGroup[gk] = h;
        rankBySpinAlex[make_pair(gk.spin, gk.A4)] += h;
        rankByAlex[gk.A4] += h;
        totalRank += h;
      }
    }
  }

  void print_output() {
    cout << "Raw RP3 grid homology ranks over F2 in all Alexander gradings\n";
    cout << "No V-factor cancellation has been performed.\n";
    cout << "n=" << n << " p=" << p << " q=" << q << " components=" << num_components << "\n\n";

    for (int spin = 0; spin < p; ++spin) {
      cout << "Spin^c = " << spin << "\n";
      int spinTotal = 0;
      for (map<GroupKey, int>::iterator it = homRankByGroup.begin(); it != homRankByGroup.end(); ++it) {
        const GroupKey &gk = it->first;
        if (gk.spin != spin) continue;
        cout << "  F2^" << it->second
             << " [ M=" << grade4_to_string(gk.M4)
             << ", A=" << grade4_to_string(gk.A4) << " ]\n";
        spinTotal += it->second;
      }
      cout << "  rank = " << spinTotal << "\n\n";
    }

    cout << "Ranks by Alexander grading, summed over Spin^c and Maslov gradings:\n";
    for (map<int,int>::iterator it = rankByAlex.begin(); it != rankByAlex.end(); ++it) {
      cout << "  A=" << grade4_to_string(it->first) << " rank=" << it->second << "\n";
    }
    cout << "Total rank = " << totalRank << "\n";
  }
};

static vector<int> parse_list(const string &s) {
  vector<int> out;
  string token;
  stringstream ss(s);
  while (getline(ss, token, ',')) {
    if (!token.empty()) out.push_back(atoi(token.c_str()));
  }
  return out;
}

int main() {
  const int n = 5;
  const int p = 2;
  const int q = 1;

  vector<int> X = {0, 1, 2, 8, 9};
  vector<int> O = {3, 5, 6, 4, 7};

  cout << "\n============================================================\n";
  cout << "12n403_quotient\n";
  cout << "============================================================\n";
  cout << "RP3 quotient input: n=5 p=2 q=1\n";

  cout << "X=";
  for (size_t i = 0; i < X.size(); ++i) {
    if (i) cout << ",";
    cout << X[i];
  }
  cout << "\n";

  cout << "O=";
  for (size_t i = 0; i < O.size(); ++i) {
    if (i) cout << ",";
    cout << O[i];
  }
  cout << "\n\n";

  try {
    RP3GridHomologyF2 calc(n, p, q, X, O);
    calc.run();
  } catch (const exception &e) {
    cerr << "Error: " << e.what() << "\n";
    cout << "ERROR: " << e.what() << "\n";
  }

  cout << "\nFinished 12n403_quotient\n";
  cout << flush;
  cerr << flush;

  return 0;
}
