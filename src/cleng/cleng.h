#ifndef CLENGxH
#define CLENGxH

#include "../input.h"
#include "../namics.h"
#include "../solve_scf.h"
#include "../system.h"
#include "../output.h"
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <memory>
#include <ostream>
#include <map>
#include <cassert>
#include <fstream>
#include <utility>
#include <chrono>
// features
#include "nodes/simple_node.h"
#include "nodes/monolit.h"
#include "nodes/point.h"
//#include "matrix.h"
#include "checkpoint/checkpoint.h"
#include "random/random.h"

#ifdef CLENG_EXPERIMENTAL
#include "cwriter/cwriter.h"
#endif

#include <csignal>

class Cleng {
private:
    const string name;
    const vector<Input *> In;
    const vector<Lattice *> Lat;
    const vector<Segment *> Seg;
    const vector<State *> Sta;
    const vector<Reaction *> Rea;
    const vector<Molecule *> Mol;
    const vector<System *> Sys;
    const vector<Solve_scf *> New;
    const string brand;

    int cmajor_version = 1;
    int cminor_version = 4;
    int cpatch = 1;
    int cversion = 129;
    int pseed{};
    string internal_name = "[Cleng] ";
    string metropolis_name = "[Metropolis] ";
    string analysis_name = "[Analysis] ";

    void fillXYZ();

public:
    Cleng(
            vector<Input *>,
            vector<Lattice *>,
            vector<Segment *>,
            vector<State *>,
            vector<Reaction *>,
            vector<Molecule *>,
            vector<System *>,
            vector<Solve_scf *>,
            string
    );

    ~Cleng();

    vector<string> KEYS;
    vector<string> PARAMETERS;
    vector<string> VALUES;

    string filename;
    Point box{Lat[0]->MX, Lat[0]->MY, Lat[0]->MZ};
    Point J{Lat[0]->JX, Lat[0]->JY, 1};
    vector<int> dims_vtk{box.x * box.y * box.z, 1};
    vector<int> dims_3 = {1, 3};

    int clamp_seg{};
    int clp_mol{};
    int n_boxes{};
    int n_out{};
    Point sub_box_size;
    int MCS{};
    int delta_step{};
    int pivot_move{};
    int pivot_axis{};
    vector<int> pivot_node_ids;
    int pivot_arms{};
    bool pivot_one_node{};
    bool pivot_one_bond{};

    map<int, Point> nodeIDs_clampedMove;
    map<int, vector<int>> pivot_arm_nodes;

// TODO: expand to template
    Matrix<int> rotation_matrix;
    Matrix<Real> scaling_matrix;

    int t{};
    int delta_save{};
    bool checkpoint_save{};
    bool checkpoint_load{};
    bool simultaneously{};
    bool metropolis{};
    int axis{};
    Real prefactor_kT{};
    Real n_times_mu{};
    string sign_move;
    bool cleng_pos{};
    bool cleng_dis{};
    bool two_ends_extension{};
    vector<int> ids_node4move;
    vector<int> ids_node4fix;
    std::chrono::steady_clock::time_point begin_simulation;
    std::chrono::steady_clock::time_point end_simulation;
    int MCS_checkpoint = 0;
    bool loaded = false;
    string save_filename;
    Point BC;
    Random rand;

    vector<Output *> Out;
    vector<int> P;
    vector<shared_ptr<SimpleNode>> simpleNodeList;
    map<int, vector<std::shared_ptr<Node>>> nodes_map;


    // temporary arrays for keep nodes_map coordinates for output
    int *xs = nullptr;
    int *ys = nullptr;
    int *zs = nullptr;

    ofstream out;
    Real free_energy_current{};
    Real free_energy_trial{};

    Real accepted{};
    Real rejected{};
    Real cleng_rejected{};
    int MC_attempt{};

    vector<Real> test_vector;

    bool MonteCarlo(bool save_vector);

    bool CP(transfer);

    bool MakeMove(bool back);

    bool MakeChecks(int id_node_for_move);

    void _moveClampedNode(bool back, int id_node_for_move, const Point &clamped_move);

    bool _pivotMoveClampedNode(const bool &back);

    bool _pivotMoveOneBond(const bool &back);

    bool _oneNodeMoveClampedNode(const bool &back);

    Point preparePivotClampedMove(int id_node_for_move);

    bool Checks(int id_node_for_move);

    bool InSubBoxRange(int id_node_for_move);

    bool IsCommensuratable();

    bool NotCollapsing(int id_node_for_move);

    bool InRange(int id_node_for_move);

    void PushOutput();

    void WriteOutput();

    void WriteClampedNodeDistance();

    void WriteClampedNodePosition();

    void make_BC();

    int getLastMCS();

    void update_ids_node4move();

    Point prepareMove(const string& type_move);

    int prepareIdNode();

    void prepareIdsNode();

    template<class T>
    Matrix<T> prepareRotationMatrix();

    template<class T>
    Matrix<T> prepareScalingMatrix();

    template<class T>
    Matrix<T> _create_rotational_matrix(int axis_rotation, int grad);

    template<class T>
    Matrix<T> _create_scaling_matrix(Real scaling_coef);

    vector<Real> prepare_vtk();

    bool CheckInput(int start, bool save_vector);

    string GetValue(const string&);

    Real GetN_times_mu();
};

#endif
