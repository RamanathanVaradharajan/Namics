
#include <iostream>
#include <vector>
#include "../alias.h"
#include "../input.h"
#include "../lattice.h"
#include "../molecule.h"
#include "../namics.h"
#include "../mesodyn.h"
#include "../cleng/cleng.h"
#include "../teng.h"
//#include "../newton.h"
#include "../output.h"
#include "../segment.h"
#include "../state.h"
#include "../reaction.h"
#include "../system.h"
#include "../tools.h"
#include "../variate.h"
#include "../sfnewton.h"
#include "../solve_scf.h"


using namespace std;

class NamicsConfig {
public:
    vector<Input*> In;     // Inputs read from file
    vector<Output*> Out;   // Outputs written to file
    vector<Lattice*> Lat;  // Properties of the lattice
    vector<Molecule*> Mol; // Properties of entire molecule
    vector<Segment*> Seg;  // Properties of molecule segments
    vector<State*> Sta;
    vector<Reaction*> Rea;
    vector<Solve_scf*> New;   // Solvers and iteration schemes
    vector<System*> Sys;
    vector<Variate*> Var;
    vector<Mesodyn*> Mes;
    vector<Cleng*> Cle; // engine for clampled molecules
    vector<Teng*> Ten;  // engine for clampled molecules


    NamicsConfig();

    bool createInClass(string filename);
    bool createLatClass(int &start);
    bool createSysClass(int &start,  bool &cuda);
    bool createSegClass();
    bool createStaClass(int &start);
    bool createReaClass(int &start);
    bool createMolClass(int &start);
    bool createVarClass(
            int &start,
            int &search_nr, int &scan_nr, int &target_nr, int &ets_nr, int &bm_nr, int &etm_nr,
            int &n_search, int &n_bm, int &n_scan, int &n_etm, int &n_ets, int &n_target);
    bool createNewClass(int &start);
    bool createOutClass(int &start, int &ii, int &n_out);
    bool guess_geometry(
            int &start,
            Real* &X, string &METHOD, vector<string> MONLIST, vector<string> STATELIST, bool &CHARGED,
            int &MX, int &MY, int &MZ, int &fjc_old
            );

    bool initCleng(
            int &start,
            Real* &X, string &METHOD, vector<string> MONLIST, vector<string> STATELIST, bool &CHARGED,
            int &MX, int &MY, int &MZ, int &fjc_old,
            bool save_vector = false
            );

    bool initMesodyn(
            int &start,
            Real* &X, string &METHOD, vector<string> MONLIST, vector<string> STATELIST, bool &CHARGED,
            int &MX, int &MY, int &MZ, int &fjc_old
            );

    bool initTeng(
            int &start,
            Real* &X, string &METHOD, vector<string> MONLIST, vector<string> STATELIST, bool &CHARGED,
            int &MX, int &MY, int &MZ, int &fjc_old
    );

    bool initSCF(
            int &start,
            int &subloop, int &substart, int &scan_nr, int &search_nr, int &ets_nr, int &etm_nr,
            int &target_nr, int &bm_nr,
            Real *&X, string &METHOD, const vector<string>& MONLIST, const vector<string>& STATELIST,
            bool &CHARGED, int &MX, int &MY, int &MZ, int &fjc_old
            );


    bool testCaseCleng(string &filename, bool save_vector = false);

    bool clear_(int &n_out);

    ~NamicsConfig();

    string s1;


};


