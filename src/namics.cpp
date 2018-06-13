#define MAINxH
#include "alias.h"
#include "input.h"
#include "lattice.h"
#include "molecule.h"
#include "namics.h"
#include "mesodyn.h"
#include "cleng.h"
#include "teng.h"
//#include "newton.h"
#include "output.h"
#include "segment.h"
#include "system.h"
#include "tools.h"
#include "variate.h"
#include "sfnewton.h"
#include "solve_scf.h"

string version = "2.1.1.1.1.1.1";
// meaning:
// newton version number =1
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

//Used for command line switches
bool debug = false;
bool suppress = false;

//Output when the user malforms input. Update when adding new command line switches.
void improperInput() {
  cerr << "Improper usage: namics [filename] [-options]." << endl << "Options available:" << endl;
  cerr << "-d Enables debugging mode." << endl;
  cerr << "-s Suppresses newton's extra output." << endl;
}

int main(int argc, char* argv[]) {
	vector<string> args(argv, argv+argc);
//Output error if no filename has been specified.
 	if (argc == 1) {
    		improperInput();
    		return 1;
  	}
//Output error if user starts with a commandline switch. (Also catches combination with forgotten filename)
  	if ( (args[1])[0] == '-' ) {
    		improperInput();
    		return 1;
  	}

// If the specified filename has no extension: add the extension specified below.
	string extension = "in";
	ostringstream filename;
	filename << argv[1];
	bool hasNoExtension = (filename.str().substr(filename.str().find_last_of(".")+1) != extension);
 	if (hasNoExtension) filename << "." << extension;

//If the switch -d is given, enable debug. Add new switches by copying and replacing -d and debug = true.
  	if ( find(args.begin(), args.end(), "-d") != args.end() ) {
    		debug = true;
 	 }

  	if ( find(args.begin(), args.end(), "-s") != args.end() ) {
    		suppress = true;
  	}

  	bool cuda;
  	int start = 0;
  	int n_starts = 0;

  	string initial_guess;
  	string final_guess;
  	string METHOD = "";
  	Real* X = NULL;
  	int MX = 0, MY = 0, MZ = 0;
  	int fjc_old;
  	bool CHARGED = false;
  	vector<string> MONLIST;

#ifdef CUDA
  GPU_present();
  cuda = true;
#else
  cuda = false;
#endif

// All class instances are stored in the following vectors
  	vector<Input*> In;     // Inputs read from file
  	vector<Output*> Out;   // Outputs written to file
  	vector<Lattice*> Lat;  // Properties of the lattice
  	vector<Molecule*> Mol; // Properties of entire molecule
  	vector<Segment*> Seg;  // Properties of molecule segments
  	vector<Solve_scf*> New;   // Solvers and iteration schemes
  	vector<System*> Sys;
  	vector<Variate*> Var;
  	vector<Mesodyn*> Mes;
  	vector<Cleng*> Cle; //enginge for clampled molecules
  	vector<Teng*> Ten; //enginge for clampled molecules

// Create input class instance and handle errors(reference above)
  	In.push_back(new Input(filename.str()) );
  	if (In[0]->Input_error) {
    		return 0;
  	}
  	n_starts = In[0]->GetNumStarts();
  	if (n_starts == 0)
    		n_starts++; // Default to 1 start..

/******** This while loop basically contains the rest of main, initializes all classes and performs  ********/
/******** calculations for a given number (start) of cycles ********/

	while (start < n_starts) {

    		start++;
   		In[0]->MakeLists(start);
    		cout << "Problem nr " << start << " out of " << n_starts << endl;

/******** Class creation starts here ********/

// Create lattice class instance and check inputs (reference above)
    		Lat.push_back(new Lattice(In, In[0]->LatList[0]));
    		if (!Lat[0]->CheckInput(start)) {
      			return 0;
    		}

// Create segment class instance and check inputs (reference above)
    		int n_seg = In[0]->MonList.size();
    		for (int i = 0; i < n_seg; i++)
      		Seg.push_back(new Segment(In, Lat, In[0]->MonList[i], i, n_seg));
    		for (int i = 0; i < n_seg; i++) {
      			for (int k = 0; k < n_seg; k++) {
        			Seg[i]->PutChiKEY(Seg[k]->name);
     			}
      			if (!Seg[i]->CheckInput(start))
        			return 0;
    		}

// Create segment class instance and check inputs (reference above)
    		int n_mol = In[0]->MolList.size();
    		for (int i = 0; i < n_mol; i++) {
      			Mol.push_back(new Molecule(In, Lat, Seg, In[0]->MolList[i]));
      			if (!Mol[i]->CheckInput(start))
        			return 0;
   		}

// Create system class instance and check inputs (reference above)
    		Sys.push_back(new System(In, Lat, Seg, Mol, In[0]->SysList[0]));
    		Sys[0]->cuda = cuda;
    		if (!Sys[0]->CheckInput(start)) {
      			return 0;
    		}
    		if (!Sys[0]->CheckChi_values(n_seg))
      			return 0;

// Prepare variables used in variate class creation
// TODO: What do these variables do?
    		int n_var = In[0]->VarList.size();
    		int n_search = 0;
    		int n_scan = 0;
    		int n_ets = 0;
    		int n_etm = 0;
    		int n_target = 0;
    		int search_nr = -1, scan_nr = -1, target_nr = -1, ets_nr = -1, etm_nr = -1;

// Create variate class instance and check inputs (reference above)
    		for (int k = 0; k < n_var; k++) {
      			Var.push_back(new Variate(In, Lat, Seg, Mol, Sys, In[0]->VarList[k]));
      			if (!Var[k]->CheckInput(start)) {
        		return 0;
      			}

      			if (Var[k]->scanning > -1) {
        			scan_nr = k;
        			n_scan++;
     			}
      			if (Var[k]->searching > -1) {
        			search_nr = k;
        			n_search++;
      			}
      			if (Var[k]->targeting > -1) {
        			target_nr = k;
        			n_target++;
      			}
      			if (Var[k]->eq_to_solvating > -1) {
        			ets_nr = k;
        			n_ets++;
      			}
      			if (Var[k]->eq_to_mu > -1) {
        			etm_nr = k;
        			n_etm++;
      			}
		}

// Error code for faulty variate class creation
		if (n_etm > 1) {
      			cout << "too many equate_to_mu's in var statements. The limit is 1 " << endl;
      			return 0;
    		}
    		if (n_ets > 1) {
      			cout << "too many equate_to_solvent's in var statements. The limit is 1 " << endl;
      			return 0;
    		}
    		if (n_search > 1) {
      			cout << "too many 'search'es in var statements. The limit is 1 " << endl;
      			return 0;
    		}
    		if (n_scan > 1) {
      			cout << "too many 'scan's in var statements. The limit is 1 " << endl;
      			return 0;
    		}
    		if (n_target > 1) {
      			cout << "too many 'target's in var statements. The limit is 1 " << endl;
      			return 0;
    		}
    		if (n_search > 0 && n_target == 0) {
      			cout << "lonely search. Please specify in 'var' a target function, e.g. 'var : sys-NN : grand_potential : 0'" << endl;
      			return 0;
    		}
    		if (n_target > 0 && n_search == 0) {
      			cout << "lonely target. Please specify in 'var' a search function, e.g. 'var : mol-lipid : search : theta'" << endl;
      			return 0;
    		}

// Create newton class instance and check inputs (reference above)
    		New.push_back(new Solve_scf(In, Lat, Seg, Mol, Sys, Var, In[0]->NewtonList[0]));
    		if (!New[0]->CheckInput(start)) {
      			return 0;
    		}
    		if (suppress == true) New[0]->e_info = false;


//Guesses geometry
    		if (Sys[0]->initial_guess == "file") {
      			MONLIST.clear();
      			if (!Lat[0]->ReadGuess(Sys[0]->guess_inputfile, X, METHOD, MONLIST, CHARGED, MX, MY, MZ, fjc_old, 0)) {
// last argument 0 is to first checkout sizes of system.
        			return 0;
      			}

			int nummon = MONLIST.size();
      			int m;
      			if (MY == 0) {
      				m = MX + 2;
      			} else {
        			if (MZ == 0) {
         				m = (MX + 2) * (MY + 2);
        			} else {
          				m = (MX + 2) * (MY + 2) * (MZ + 2);
        			}
      			}
      			int IV = nummon * m;

      			if (CHARGED)
        			IV += m;
      			if (start > 1)
        			free(X);
X = (Real*)malloc(IV * sizeof(Real));
      			MONLIST.clear();
      			Lat[0]->ReadGuess(Sys[0]->guess_inputfile, X, METHOD, MONLIST, CHARGED, MX, MY, MZ, fjc_old, 1);
// last argument 1 is to read guess in X.
		}

    		int substart = 0;
    		int subloop = 0;
    		if (scan_nr < 0)
      			substart = 0;
   		else
      			substart = Var[scan_nr]->num_of_cals;
    		if (substart < 1)
      			substart = 1; // Default to 1 substart

		EngineType TheEngine;
		TheEngine=SCF;
		if (In[0]->MesodynList.size()>0) {TheEngine=MESODYN;}
		if (In[0]->ClengList.size()>0) {TheEngine=CLENG;}
		if (In[0]->TengList.size()>0) {TheEngine=TENG;}

		int ii,kk,length,length_al;
		int n_out=0;
		switch(TheEngine) {
			case SCF:
				// Prepare, catch errors for output class creation
    				n_out = In[0]->OutputList.size();
    				if (n_out == 0)
      					cout << "Warning: no output defined!" << endl;

				// Create output class instance and check inputs (reference above)
    				for (int ii = 0; ii < n_out; ii++) {
      					Out.push_back(new Output(In, Lat, Seg, Mol, Sys, New, In[0]->OutputList[ii], ii, n_out));
      					if (!Out[ii]->CheckInput(start)) {
        					cout << "input_error in output " << endl;
        					return 0;
      					}
    				}

				while (subloop < substart) {
      					if (scan_nr > -1)
        					Var[scan_nr]->PutVarScan(subloop);
      					New[0]->AllocateMemory();
      					New[0]->Guess(X, METHOD, MONLIST, CHARGED, MX, MY, MZ,fjc_old);
					if (search_nr < 0 && ets_nr < 0 && etm_nr < 0) {
       	  				New[0]->Solve(true);
          				} else {
            					if (!debug) cout << "Solve towards superiteration " << endl;
            					New[0]->SuperIterate(search_nr, target_nr, ets_nr, etm_nr);
         				}
					Lat[0]->PushOutput();
      					New[0]->PushOutput();
      					length = In[0]->MonList.size();
      					for (ii = 0; ii < length; ii++)
        					Seg[ii]->PushOutput();
      					length = In[0]->MolList.size();
      					for (ii = 0; ii < length; ii++) {
        				  length_al = Mol[ii]->MolAlList.size();
        				for (kk = 0; kk < length_al; kk++) {
          				Mol[ii]->Al[kk]->PushOutput();
        					}
       			 		Mol[ii]->PushOutput();
      					}
      					// length = In[0]->AliasList.size();
      					// for (ii=0; ii<length; ii++) Al[i]->PushOutput();
     					Sys[0]->PushOutput(); // needs to be after pushing output for seg.

      					for (ii = 0; ii < n_out; ii++) {
        					Out[ii]->WriteOutput(subloop);
      					}
      					subloop++;
    				}

				break;
			case MESODYN:
				New[0]->AllocateMemory();
      				New[0]->Guess(X, METHOD, MONLIST, CHARGED, MX, MY, MZ,fjc_old);
				if (debug) cout << "Creating mesodyn" << endl;
        			Mes.push_back(new Mesodyn(In, Lat, Seg, Mol, Sys, New, In[0]->MesodynList[0]));
        			if (!Mes[0]->CheckInput(start)) {
          				return 0;
        			}
        			Mes[start-1]->mesodyn();
				break;
			case CLENG:
				New[0]->AllocateMemory();
      				New[0]->Guess(X, METHOD, MONLIST, CHARGED, MX, MY, MZ,fjc_old);
				if (!debug) cout << "Creating Cleng module" << endl;
				Cle.push_back(new Cleng(In, Lat, Seg, Mol, Sys, New, In[0]->ClengList[0]));
				if (!Cle[0]->CheckInput(start)) {return 0;}
				break;
			case TENG:
				New[0]->AllocateMemory();
      				New[0]->Guess(X, METHOD, MONLIST, CHARGED, MX, MY, MZ,fjc_old);
				cout << "Solving Teng problem" << endl;
				Ten.push_back(new Teng(In, Lat, Seg, Mol, Sys, New, In[0]->TengList[0]));
				if (!Ten[0]->CheckInput(start)) {return 0;}
				break;
			default:
				cout <<"TheEngine is unknown. Programming error " << endl; return 0;
				break;
<<<<<<< HEAD
		}



=======
    }
>>>>>>> 07db4c70110a77e032cc3c011780530ee6de9888
    		if (scan_nr > -1)
      			Var[scan_nr]->ResetScanValue();
    		if (Sys[0]->initial_guess == "previous_result") {
     			int IV_new=New[0]->iv; //check this
      			METHOD = New[0]->SCF_method; //check this..
      			MX = Lat[0]->MX;
      			MY = Lat[0]->MY;
      			MZ = Lat[0]->MZ;
      			CHARGED = Sys[0]->charged;
      			if (start > 1 || (start == 1 && Sys[0]->initial_guess == "file")) free(X);
X = (Real *)malloc(IV_new * sizeof(Real));
      			for (int i = 0; i < IV_new; i++) X[i] = New[0]->xx[i];
      			fjc_old=Lat[0]->fjc;
      			int length = Sys[0]->SysMonList.size();
      			MONLIST.clear();
      			for (int i = 0; i < length; i++) {
       			MONLIST.push_back(Seg[Sys[0]->SysMonList[i]]->name);
      			}
    		}
    		if (Sys[0]->final_guess == "file") {
      			MONLIST.clear();
      			int length = Sys[0]->SysMonList.size();
      			for (int i = 0; i < length; i++) {
        			MONLIST.push_back(Seg[Sys[0]->SysMonList[i]]->name);
      			}
      			Lat[0]->StoreGuess(Sys[0]->guess_outputfile, New[0]->xx, New[0]->SCF_method, MONLIST, Sys[0]->charged, start);
   		}

/******** Clear all class instances ********/
    		//for (int i = 0; i < n_out; i++) delete Out[i];
		Out.clear();
   	 	//for (int i = 0; i < n_var; i++)
      			//delete Var[i];
    		Var.clear();
   	 	//delete New[0];
    		New.clear();
   	 	//delete Sys[0];
    		Sys.clear();
   	 	//for (int i = 0; i < n_mol; i++) {
      			//delete Mol[i];
    		//}
    		Mol.clear();
    		//for (int i = 0; i < n_seg; i++) {
			//delete Seg[i];
		//}
   	 	Seg.clear();
    		//delete Lat[0];
    		Lat.clear();
	} //loop over starts.
return 0;
}
