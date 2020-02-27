#define MAINxH
#include "namics_config/namics_config.h"
// Tests
#define CATCH_CONFIG_RUNNER
#include "test_unit/catch.hpp"
// Cleng test
#include "test_unit/test_cleng.h"
// You can put your header with tests.
// #include "namics_config/test_***.h"

string version = "2.1.1.1.1.1.1";
// meaning:
// newton version number =2
// system version number =1
// lattice version number =1
// molecule version number =1
// segment version number =1
// alias version number =1
// output version number =1
Real check = 0.4534345;
Real e = 1.60217e-19;
Real T = 298.15;
Real k_B = 1.38065e-23;
Real k_BT = k_B * T;
Real eps0 = 8.85418e-12;
Real PIE = 3.14159265;
int DEBUG_BREAK = 1;
//Used for command line switches
bool debug = false;

//Output when the user malforms input. Update when adding new command line switches.
void help() {
    cerr << "Usage: namics [-options/--options] [filename]." << endl;

    cerr << "Options available:" << endl;
    cerr << endl;
    cerr << "-t [--tests] [show list of tests by tags. " << endl;
    cerr << "              Choose all by <-#> or by specific tag name [name_of_test]" << endl;
    cerr << "              Use additional -l flag for showing all available test cases. [-?] for details.]" << endl;
    cerr << "-d [Enables debugging mode]" << endl;
    cerr << "-GPU [index] Sets the GPU to be used in multi-GPU systems." << endl;
}

/**
 * Main
 * @param argc
 * @param argv
 * @return
 *
 * Main is slightly elaborated.
 * New:
 *      * introduced NamicsConfig class
 *          ** NamicsConfig create the rest classes such as InClass, LatClass, SysClass and so on.
 *          ** The NamicsConfig is needed because of test_unit. Look TEST_CASE procedure.
 *      * slightly changed reading flags from console
 */

int main(int argc, char *argv[]) {

    int errors_ = 0;  // counting errors for executing namics
    bool GPU = false;
    bool has_filename = false;
    string filename_;

    if (argc == 1) { help(); return 1; } // if namics without any arguments

    for( int i = 1; i < argc; ++i ) {

        if (strcmp( argv[i], "-d") == 0)   debug = true;
        if (strcmp( argv[i], "-GPU") == 0) GPU   = true;
        if (*argv[i] != '-') { filename_ = argv[i]; has_filename = true; }

        if (strcmp( argv[i], "--tests") == 0 or (strcmp( argv[i], "-t") == 0 )) {
            std::vector<char*> argv_; // valid parameters for testing
            for( int jj = 0; jj < argc; ++jj ) {
                if (strcmp(argv[jj], "-GPU") != 0 and
                    strcmp(argv[jj], "-d") != 0 and
                    strcmp(argv[jj], "-t") != 0 and
                    strcmp(argv[jj], "--tests") != 0)
                { argv_.push_back(argv[jj]); }}

            int argc_ = (int) argv_.size();
            if (argc_ == 1) { argv_.push_back((char*)"-t"); argc_ ++; } // because -t flag is used Catch2
            Catch::Session session;
            int returnCode = session.applyCommandLine(argc_, argv_.data());
            // Indicates a command line error
            if( returnCode != 0 ) return returnCode;
            int result = session.run();
            return result;
        }
    }

    // check last time
    if (errors_ != 0 or !has_filename) { help(); return 1; }

// If the specified filename has no extension: add the extension specified below.
    string extension = "in";
    ostringstream filename;
    filename << filename_;
    bool hasNoExtension = (filename.str().substr(filename.str().find_last_of('.') + 1) != extension);
    if (hasNoExtension) filename << "." << extension;

    int cudaDeviceIndex = 0;
    bool cuda = false;
    // If the switch -GPU is given, select GPU.
    if (GPU) {
        vector<string> args(argv, argv+argc);
        try { cudaDeviceIndex = load_argument_value(args, "-GPU", cudaDeviceIndex); }
        catch (int) { help(); exit(0); }
    }

    int start = 0;
    int n_starts = 0;

    string initial_guess;
    string final_guess;
    string METHOD;
    Real *X = nullptr;
    int MX = 0, MY = 0, MZ = 0;
    int fjc_old = 0;
    bool CHARGED = false;
    vector<string> MONLIST;
    vector<string> STATELIST;

#ifdef CUDA
	GPU_present(cudaDeviceIndex);
	cuda = true;
	SUM_RESULT = (Real*)AllOnDev(1);
#else
	cuda = false;
	SUM_RESULT = new Real;
#endif

// Now NamicsConfig knows all about class instances
    bool success;
    NamicsConfig config;

// Create input class instance and handle errors
    success = config.createInClass(filename.str());
    if (!success) return 0;

    n_starts = config.In[0]->GetNumStarts();
    if (n_starts == 0) n_starts++; // Default to 1 start..

/******** This while loop basically contains the rest of main, initializes all classes and performs  ********/
/******** calculations for a given number (start) of cycles ********/

    while (start < n_starts)
	{

        start++;
        config.In[0]->MakeLists(start);
        cout << "Problem nr " << start << " out of " << n_starts << endl;

/******** Class creation starts here ********/

// Create lattice class instance and check inputs
        success = config.createLatClass(start); if (!success) return 0;
// Create segment class instance and check inputs
        success = config.createSegClass();      if (!success) return 0;
//Create state class instance and check inputs
        success = config.createStaClass(start); if (!success) return 0;
//Create reaction class instance and check inputs
        success = config.createReaClass(start); if (!success) return 0;
// Create molecule class instance and check inputs
        success = config.createMolClass(start); if (!success) return 0;
// Create system class instance and check inputs
        success = config.createSysClass(start, cuda); if (!success) return 0;

//// all variables below go to createVarClass as links
        int search_nr = -1, scan_nr = -1, target_nr = -1, ets_nr = -1, bm_nr = -1, etm_nr = -1;
        int n_search = 0,   n_scan = 0,   n_ets = 0,      n_bm = 0,    n_etm = 0,  n_target = 0;
//// until here
// Create variate class instance and check inputs (reference above)
        success = config.createVarClass(
                start,
                search_nr, scan_nr, target_nr, ets_nr, bm_nr, etm_nr,
                n_search, n_bm, n_scan, n_etm, n_ets, n_target
        ); if (!success) return 0;

// Error code for faulty variate class creation
        if (n_etm > 1) { cout << "too many equate_to_mu's in var statements. The limit is 1 " << endl;      return 0;}
        if (n_ets > 1) { cout << "too many equate_to_solvent's in var statements. The limit is 1 " << endl; return 0;}
        if (n_bm > 1) { cout << "too many balance membrane's in var statements. The limit is 1 " << endl; return 0;}
        if (n_search > 1) { cout << "too many 'searches' in var statements. The limit is 1 " << endl; return 0;}
        if (n_scan > 1) { cout << "too many 'scan's in var statements. The limit is 1 " << endl; return 0; }
        if (n_target > 1) { cout << "too many 'target's in var statements. The limit is 1 " << endl; return 0;}
        if (n_search > 0 && n_target == 0) { cout << "lonely search. Please specify in 'var' a target function, " <<
                 " e.g. 'var : sys-NN : grand_potential : 0'" << endl; return 0;}
        if (n_target > 0 && n_search == 0) { cout << "lonely target. Please specify in 'var' a search function, " <<
                 " e.g. 'var : mol-lipid : search : theta'" << endl; return 0; }

// Create newton class instance and check inputs (reference above)
        success = config.createNewClass(start); if (!success) return 0;

//Guesses geometry
        success = config.guess_geometry(start, X, METHOD, MONLIST, STATELIST, CHARGED, MX, MY, MZ, fjc_old);
        if (!success) return 0;

        int substart = 0, subloop = 0;
        if (scan_nr < 0) substart = 0;
        else substart = config.Var[scan_nr]->num_of_cals;
        if (substart < 1) substart = 1; // Default to 1 substart

        EngineType TheEngine;
        TheEngine = SCF;
        if (!config.In[0]->MesodynList.empty()) { TheEngine = MESODYN;
        }
		if (!config.In[0]->ClengList.empty()) { TheEngine = CLENG;
        }
		if (!config.In[0]->TengList.empty()) { TheEngine = TENG; }

        int n_out = 0;
        switch (TheEngine) {
            case SCF:
                success = config.initSCF(
                        start,
                        subloop, substart, scan_nr, search_nr, ets_nr, etm_nr, target_nr, bm_nr,
                        X, METHOD, MONLIST, STATELIST, CHARGED, MX, MY, MZ, fjc_old);
                if (!success) return success;
                break;

            case MESODYN:
                success = config.initMesodyn(start, X, METHOD, MONLIST, STATELIST, CHARGED, MX, MY, MZ, fjc_old);
                if (!success) return success;
                break;

            case CLENG:
                success = config.initCleng(start, X, METHOD, MONLIST, STATELIST, CHARGED, MX, MY, MZ, fjc_old);
                if (!success) return success;
                break;

            case TENG:
                success = config.initTeng(start, X, METHOD, MONLIST, STATELIST, CHARGED, MX, MY, MZ, fjc_old);
                if (!success) return success;
                break;

                // this unreachable part of code! Do we need it?
            default:
                cout << "TheEngine is unknown. Programming error " << endl;
                return 0;
                break;
                // end of unreachable part
        }

		for (auto all_segments : config.Seg) all_segments->prepared = false;
        if (scan_nr > -1) config.Var[scan_nr]->ResetScanValue();

        // TODO: What does it do?
        if (config.Sys[0]->initial_guess == "previous_result") {
            METHOD = config.New[0]->SCF_method; //check this..
            MX = config.Lat[0]->MX;
            MY = config.Lat[0]->MY;
            MZ = config.Lat[0]->MZ;
            CHARGED = config.Sys[0]->charged;
            int IV_new = config.New[0]->iv; //check this
            if (start > 1 || (start == 1 && config.Sys[0]->initial_guess == "file"))
#ifdef CUDA
                cudaFree(X);
#else
                free(X);
#endif
#ifdef CUDA
            X = (Real *)AllOnDev(IV_new);
#else
            X = (Real *)malloc(IV_new * sizeof(Real));
#endif
            //for (int i = 0; i < IV_new; i++) X[i] = config.New[0]->xx[i];
            Cp(X, config.New[0]->xx, IV_new);
            fjc_old = config.Lat[0]->fjc;
            int mon_length = config.Sys[0]->ItMonList.size();
            int state_length = config.Sys[0]->ItStateList.size();
            MONLIST.clear();
            STATELIST.clear();
            for (int i = 0; i < mon_length; i++) MONLIST.push_back(config.Seg[config.Sys[0]->ItMonList[i]]->name);
            for (int i = 0; i < state_length; i++) STATELIST.push_back(config.Sta[config.Sys[0]->ItStateList[i]]->name);
        }

        // TODO: What does it do (2)?
        if (config.Sys[0]->final_guess == "file") {
            MONLIST.clear();
            STATELIST.clear();
            int mon_length = config.Sys[0]->ItMonList.size();
            int state_length = config.Sys[0]->ItStateList.size();
            for (int i = 0; i < mon_length; i++) MONLIST.push_back(config.Seg[config.Sys[0]->ItMonList[i]]->name);
            for (int i = 0; i < state_length; i++) STATELIST.push_back(config.Sta[config.Sys[0]->ItStateList[i]]->name);
            config.Lat[0]->StoreGuess(
                    config.Sys[0]->guess_outputfile, config.New[0]->xx, config.New[0]->SCF_method, MONLIST, STATELIST,
                    config.Sys[0]->charged, start);
        }

/******** Clear all class instances ********/
        config.clear_(n_out);
    } //loop over starts.
    #ifdef CUDA
	cudaFree(X);
#else
	free(X);
#endif
    delete config.In[0];
    config.In.clear();
    return 0;
}
