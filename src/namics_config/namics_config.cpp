#include <utility>
#include "namics_config.h"


NamicsConfig::NamicsConfig() {
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

}

NamicsConfig::~NamicsConfig() = default;

bool NamicsConfig::createInClass(string filename) {
    bool success;
    In.push_back(new Input(std::move(filename)) );
    if (In[0]->Input_error) {
        if (debug) cout << "Input_error raises error!" << endl;
        success = false;
    } else success = true;
    return success;
}

bool NamicsConfig::createLatClass(int &start) {
    bool success;
    Lat.push_back(new Lattice(In, In[0]->LatList[0]));
    if (!Lat[0]->CheckInput(start)) {
        if (debug) cout << "CheckInput in Lat raises error!" << endl;
        success = false;
    } else success = true;
    return success;
}

bool NamicsConfig::createSegClass() {
    bool success = true;
    int n_seg = In[0]->MonList.size();  // unsigned long to int !?!
    for (int i = 0; i < n_seg; i++) {
        Seg.push_back(new Segment(In, Lat, In[0]->MonList[i], i, n_seg));
    }
    return success;
}

bool NamicsConfig::createStaClass(int &start) {
    bool success = true;
    int n_stat = In[0]->StateList.size();
    int n_seg = In[0]->MonList.size();  // unsigned long to int !?!

    for (int i = 0; i < n_stat; i++) Sta.push_back(new State(In, Seg, In[0]->StateList[i]));

    for (int i = 0; i < n_seg; i++) {
        for (int k = 0; k < n_seg; k++)  Seg[i]->PutChiKEY(Seg[k]->name);
        for (int k = 0; k < n_stat; k++) Seg[i]->PutChiKEY(Sta[k]->name);
        if (!Seg[i]->CheckInput(start)){
            if (debug) cout << "CheckInput in Seg raises error!" << endl;
            success = false;
            return success;
        }
    }
    for (int i = 0; i < n_stat; i++) {
        for (int k = 0; k < n_seg; k++)  Sta[i]->PutChiKEY(Seg[k]->name);
        for (int k = 0; k < n_stat; k++) Sta[i]->PutChiKEY(Sta[k]->name);
        if (!Sta[i]->CheckInput(start)){
            if (debug) cout << "CheckInput in Sta raise error!" << endl;
            success = false;
            return success;
        }
    }
    return success;
}

bool NamicsConfig::createReaClass(int &start) {
    bool success = true;
    int n_rea = In[0]->ReactionList.size();
    for (int i = 0; i < n_rea; i++) {
        Rea.push_back(new Reaction(In, Seg, Sta, In[0]->ReactionList[i]));
        if (!Rea[i]->CheckInput(start)) {
            if (debug) cout << "CheckInput in Rea raises error!" << endl;
            success = false;
            return success;
        }
    }
    return success;
}

bool NamicsConfig::createMolClass(int &start) {
    bool success = true;
    int n_mol = In[0]->MolList.size();
    for (int i = 0; i < n_mol; i++) {
        Mol.push_back(new Molecule(In, Lat, Seg, In[0]->MolList[i]));
        if (!Mol[i]->CheckInput(start)) {
            if (debug) cout << "CheckInput in Mol raises error!" << endl;
            success = false;
            return success;
        }
    }
    return success;
}

bool NamicsConfig::createVarClass(
        int &start,
        int &search_nr, int &scan_nr, int &target_nr, int &ets_nr, int &bm_nr, int &etm_nr,
        int &n_search, int &n_bm, int &n_scan, int &n_etm, int &n_ets, int &n_target\
        ) {
    bool success = true;
    int n_var = In[0]->VarList.size();

    for (int k = 0; k < n_var; k++) {
        Var.push_back(new Variate(In, Lat, Seg, Sta, Rea, Mol, Sys, In[0]->VarList[k]));
        if (!Var[k]->CheckInput(start)) {
            if (debug) cout << "CheckInput in Mol raises error!" << endl;
            success = false;
            return success;
        }

        if (Var[k]->scanning > -1) {  scan_nr   = k; n_scan++;   }
        if (Var[k]->searching > -1) { search_nr = k; n_search++; }
        if (Var[k]->targeting > -1) { target_nr = k; n_target++; }
        if (Var[k]->eq_to_solvating > -1) { ets_nr = k; n_ets++; }
        if (Var[k]->balance_membraning > -1) { bm_nr = k; n_bm++;}
        if (Var[k]->eq_to_mu > -1) { etm_nr = k; n_etm++; }
    }

    return success;
}

bool NamicsConfig::createNewClass(int &start) {
    bool success = true;
    New.push_back(new Solve_scf(In, Lat, Seg, Sta, Rea, Mol, Sys, Var, In[0]->NewtonList[0]));
    if (!New[0]->CheckInput(start))  {
        if (debug) cout << "CheckInput in New raises error!" << endl;
        success = false;
        return success;
    }

    return success;
}

bool NamicsConfig::createOutClass(int &start, int &ii, int &n_out) {
    bool success = true;
    Out.push_back(new Output(In, Lat, Seg, Sta, Rea, Mol, Sys, New, In[0]->OutputList[ii], ii, n_out));
    if (!Out[ii]->CheckInput(start)) {
        cout << "input_error in output " << endl;
        success = false;
        return success;
    }
    return success;
}

bool NamicsConfig::guess_geometry(
        int &start,
        Real* &X, string &METHOD, vector<string> MONLIST, vector<string> STATELIST, bool &CHARGED,
        int &MX, int &MY, int &MZ, int &fjc_old
        ) {
    bool success = true;
    if (Sys[0]->initial_guess == "file") {
        MONLIST.clear();
        STATELIST.clear();
        if (!Lat[0]->ReadGuess(Sys[0]->guess_inputfile, X, METHOD, MONLIST, STATELIST, CHARGED, MX, MY, MZ, fjc_old,
                               0))
        {
            // last argument 0 is to first checkout sizes of system.
            success = false;
            return success;
        }
        int nummon = MONLIST.size();
        int numstate = STATELIST.size();
        int m;
        if (MY == 0) m = MX + 2;
        else {
            if (MZ == 0) m = (MX + 2) * (MY + 2);
            else m = (MX + 2) * (MY + 2) * (MZ + 2);
        }
        int IV = (nummon + numstate) * m;

        if (CHARGED) IV += m;
        if (start > 1) {
#ifdef CUDA
            cudaFree(X);
            X = (Real *)AllOnDev(IV);
#else
            free(X);
            X = (Real *) malloc(IV * sizeof(Real));
#endif
        }
        MONLIST.clear();
        STATELIST.clear();
        Lat[0]->ReadGuess(Sys[0]->guess_inputfile, X, METHOD, MONLIST, STATELIST, CHARGED, MX, MY, MZ, fjc_old, 1);
// last argument 1 is to read guess in X.
    }

    return success;
}

bool NamicsConfig::createSysClass(int &start, bool &cuda) {
    bool success = true;

    int n_seg = In[0]->MonList.size();  // unsigned long to int !?!

    Sys.push_back(new System(In, Lat, Seg, Sta, Rea, Mol, In[0]->SysList[0]));
    Sys[0]->cuda = cuda;
    if (!Sys[0]->CheckInput(start)) {
        cout << "input_error in output " << endl;
        success = false;
        return success;
    }
    if (!Sys[0]->CheckChi_values(n_seg)) {
        cout << "input_error in output " << endl;
        success = false;
        return success;
    }

    return success;
}

bool NamicsConfig::initCleng(
        int &start,
        Real* &X, string &METHOD, vector<string> MONLIST, vector<string> STATELIST, bool &CHARGED,
        int &MX, int &MY, int &MZ, int &fjc_old,
        bool save_vector
        ) {
    bool success = true;

    New[0]->AllocateMemory();
    New[0]->Guess(X, METHOD, std::move(MONLIST), std::move(STATELIST), CHARGED, MX, MY, MZ, fjc_old);
    if (!debug) cout << "Creating Cleng module" << endl;
    Cle.push_back(new Cleng(In, Lat, Seg, Sta, Rea, Mol, Sys, New, In[0]->ClengList[0]));
    if (!Cle[0]->CheckInput(start, save_vector)) {
        if (debug) cout << "CheckInput in Cle raises error! " << endl;
        success = false;
        return success;
    }
//    delete Cle[0];
//    Cle.clear();
    return success;
}

bool NamicsConfig::initMesodyn(
        int &start,
        Real *&X, string &METHOD, vector<string> MONLIST, vector<string> STATELIST,
        bool &CHARGED, int &MX, int &MY, int &MZ, int &fjc_old
        ) {
    bool success = true;
    New[0]->mesodyn = true;
    New[0]->AllocateMemory();
    New[0]->Guess(X, METHOD, std::move(MONLIST), std::move(STATELIST), CHARGED, MX, MY, MZ, fjc_old);
    if (debug) cout << "Creating mesodyn" << endl;
    Mes.push_back(new Mesodyn(start, In, Lat, Seg, Sta, Rea, Mol, Sys, New, In[0]->MesodynList[0]));
    if (!Mes[0]->CheckInput()) {
        if (debug) cout << "CheckInput in Mes raises error! " << endl;
        success = false;
        return success;
    }
    try { Mes[0]->mesodyn();}
    catch (error RC) {
        cout << "Exiting Mesodyn with error: " << RC << "."
             << "Please check the documentation for details." << endl;
        exit(RC);
    }
    delete Mes[0];
    Mes.clear();
    return success;
}

bool NamicsConfig::initTeng(
        int &start,
        Real *&X, string &METHOD, vector<string> MONLIST, vector<string> STATELIST,
        bool &CHARGED, int &MX, int &MY, int &MZ, int &fjc_old
        ) {
    bool success = true;

    New[0]->AllocateMemory();
    New[0]->Guess(X, METHOD, std::move(MONLIST), std::move(STATELIST), CHARGED, MX, MY, MZ, fjc_old);
    cout << "Solving Teng problem" << endl;
    Ten.push_back(new Teng(In, Lat, Seg, Sta, Rea, Mol, Sys, New, In[0]->TengList[0]));
    if (!Ten[0]->CheckInput(start)) {
        if (debug) cout << "CheckInput in Ten raises error! " << endl;
        success = false;
        return success;
    }

    return  success;
}

bool NamicsConfig::initSCF(
        int &start,
        int &subloop, int &substart, int &scan_nr, int &search_nr, int &ets_nr, int &etm_nr,
        int &target_nr, int &bm_nr,
        Real *&X, string &METHOD, const vector<string>& MONLIST, const vector<string>& STATELIST,
        bool &CHARGED, int &MX, int &MY, int &MZ, int &fjc_old
        ) {
    bool success = true;

    // Prepare, catch errors for output class creation
    int n_out = In[0]->OutputList.size();
    if (n_out == 0) cout << "Warning: no output defined!" << endl;

    // Create output class instance and check inputs (reference above)
    for (int ii = 0; ii < n_out; ii++) {
        success = createOutClass(start, ii, n_out);
        if (!success) return success;
    }

    while (subloop < substart) {
        if (scan_nr > -1) Var[scan_nr]->PutVarScan(subloop);
        New[0]->AllocateMemory();
        New[0]->Guess(X, METHOD, MONLIST, STATELIST, CHARGED, MX, MY, MZ, fjc_old);
        if (search_nr < 0 && ets_nr < 0 && etm_nr < 0) { New[0]->Solve(true);
        } else {
            if (!debug) cout << "Solve towards superiteration " << endl;
            New[0]->SuperIterate(search_nr, target_nr, ets_nr, etm_nr, bm_nr);
        }
        New[0]->PushOutput();

        for (int ii = 0; ii < n_out; ii++) Out[ii]->WriteOutput(subloop);
        subloop++;
    }

    return success;
}



bool NamicsConfig::clear_(int &n_out) {
    bool success = true;

    int n_var  = In[0]->VarList.size();
    int n_seg  = In[0]->MonList.size();
    int n_stat = In[0]->StateList.size();
    int n_rea  = In[0]->ReactionList.size();
    int n_mol  = In[0]->MolList.size();

    for (int i = 0; i < n_out; i++)  {delete Out[i];} Out.clear();
    for (int i = 0; i < n_var; i++)  {delete Var[i];} Var.clear();
    for (int i = 0; i < n_mol; i++)  {delete Mol[i];} Mol.clear();
    for (int i = 0; i < n_seg; i++)  {delete Seg[i];} Seg.clear();
    for (int i = 0; i < n_stat; i++) {delete Sta[i];} Sta.clear();
    for (int i = 0; i < n_rea; i++)  {delete Rea[i];} Rea.clear();
    delete New[0]; New.clear();
    delete Sys[0]; Sys.clear();
    delete Lat[0]; Lat.clear();

    return success;
}


/**
 * testCaseCleng
 *
 * How it currently works:
 *      1) Initially you should create a script for testing.
 *  (this case it is linear chain terminated by 2 clamped node)
 *  the testing script should lie in reachiable directory. (currently it is "inputs/test_scripts")
 *      2) For testing you should provide "filename" of the testing script.
 *      3) For testing it is necessary to create all Namics classes such as In, Lat, Seg, Rea and so on.
 *
 * @param filename is a name of script for testing case
 * @param save_vector is a flag for saving your data in some variable for asking later
 * (e.g. could be used for comparison expected "free_energy_vector" and current "free_energy_vector" (from Cleng)
 * @return success - bool variable of classes success
 */

bool NamicsConfig::testCaseCleng(string &filename, bool save_vector) {

// some compulsory part of testing;
    bool success;
    bool cuda = false;
    string METHOD;
    Real* X = nullptr;
    int MX = 0, MY = 0, MZ = 0;
    int fjc_old = 0;
    bool CHARGED = false;
    vector<string> MONLIST;
    vector<string> STATELIST;

    // Create input class instance and handle errors (more details in namics.cpp)
    success = createInClass(filename);
    if (!success) {cout << "Problem is in In" << endl; return success;}

    int n_starts = 0;
    n_starts = In[0]->GetNumStarts();
    if (n_starts == 0) n_starts++; // Default to 1 start..
    int start = 0;
    while (start < n_starts) {
        start++;
        In[0]->MakeLists(start);

        success = createLatClass(start);
        if (!success) {if (debug) cout << "Problem is in Lat" << endl; return success;}
        success = createSegClass();
        if (!success) {if (debug) cout << "Problem is in Seg" << endl; return success;}
        success = createStaClass(start);
        if (!success) {if (debug) cout << "Problem is in Sta" << endl; return success;}
        success = createReaClass(start);
        if (!success) {if (debug) cout << "Problem is in Rea" << endl; return success;}
        success = createMolClass(start);
        if (!success) {if (debug) cout << "Problem is in Mol" << endl; return success;}
        success = createSysClass(start, cuda);
        if (!success) {if (debug) cout << "Problem is in Sys" << endl; return success;}

        int search_nr = -1, scan_nr = -1, target_nr = -1, ets_nr = -1, bm_nr = -1, etm_nr = -1;
        int n_search = 0; int n_scan = 0; int n_ets = 0;
        int n_bm = 0;     int n_etm = 0;  int n_target = 0;
        success = createVarClass(
                start,
                search_nr, scan_nr, target_nr, ets_nr,bm_nr, etm_nr,
                n_search, n_bm, n_scan, n_etm, n_ets, n_target
        );
        if (!success) {if (debug) cout << "Problem is in Var" << endl; return success;}

        success = createNewClass(start);
        if (!success) {if (debug) cout << "Problem is in New" << endl; return success;}

// end of compulsory part of testing;

        // Can be exchanged by any others engine. /*
        success = initCleng(start, X, METHOD, MONLIST, STATELIST, CHARGED, MX, MY, MZ, fjc_old, save_vector);
        if (!success) {if (debug) cout << "Problem is in initCleng" << endl; return success;}
        // */
    }

    return success;
}
