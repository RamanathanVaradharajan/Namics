#include "cleng.h"
#include "iterator/EnumerateIterator.h"

volatile sig_atomic_t cleng_flag_termination = 0;

using Fp_info = numeric_limits<double>;

inline auto is_ieee754_nan(double const x) -> bool {
    static constexpr bool is_claimed_ieee754 = Fp_info::is_iec559;
    static constexpr int n_bits_per_byte = CHAR_BIT;
    using Byte = unsigned char;

    static_assert(is_claimed_ieee754, "!");
    static_assert(n_bits_per_byte == 8, "!");
    static_assert(sizeof(x) == sizeof(uint64_t), "!");

#ifdef _MSC_VER
    uint64_t const bits = reinterpret_cast<uint64_t const&>( x );
#else
    Byte bytes[sizeof(x)];
    memcpy(bytes, &x, sizeof(x));
    uint64_t int_value;
    memcpy(&int_value, bytes, sizeof(x));
    uint64_t const &bits = int_value;
#endif

    static constexpr uint64_t sign_mask = 0x8000000000000000;
    static constexpr uint64_t exp_mask = 0x7FF0000000000000;
    static constexpr uint64_t mantissa_mask = 0x000FFFFFFFFFFFFF;

    (void) sign_mask;
    return (bits & exp_mask) == exp_mask and (bits & mantissa_mask) != 0;
}

shared_ptr<SimpleNode> fromSystemToNode(int x, int y, int z, int id, const Point &box) {
    return make_shared<SimpleNode>(Point(x, y, z), id, box);
}

map<int, vector<shared_ptr<Node>>> createNodes(const vector<shared_ptr<SimpleNode>> &simple_nodes, map<int, vector<int>> &pivot_arm_nodes, int &arms) {
    map<int, vector<shared_ptr<Node>>> result_map;
    int index = 0;
    arms = 0;
    map<SimpleNode, vector<shared_ptr<SimpleNode>>> m;
    for (auto &&n : simple_nodes) m[*n].push_back(n);
    for (auto &&entry : m) {
        if (entry.second.size() == 1) result_map[index] = {entry.second[0]};
        else result_map[index] = {make_shared<Monolit>(entry.second)};
        index++;
    }
    // pivot_arm_nodes
    Point P_first = simple_nodes[0]->point();
    for (auto &&n : simple_nodes) {
        int ID_test = n->get_ID();
//        cout << "ID_test: " << ID_test << endl;
        if (n->point() == P_first) arms++;
        for (auto &&node : result_map) {
            if (node.second.data()->get()->isIdInside(ID_test)) {pivot_arm_nodes[arms].push_back(node.first);}
        }
    }
    // clean up the repeating elements
    for (auto &&pair : pivot_arm_nodes){
        int arm_ = pair.first;
        auto last = std::unique(pivot_arm_nodes[arm_].begin(), pivot_arm_nodes[arm_].end());
        pivot_arm_nodes[arm_].erase(last, pivot_arm_nodes[arm_].end());
    }

    return result_map;
}

vector<int> makeExcludedArray(int step) {
    vector<int> result;
    for (int i = -step + 1; i < step; i++) result.push_back(i);
    return result;
}

bool Cleng::InSubBoxRange(int id_node_for_move) {
    bool success = true;
    vector<int> id_nodes;
    Point sub_box = {sub_box_size.x, sub_box_size.y, sub_box_size.z};  // change it in future

    if (!nodes_map[id_node_for_move].data()->get()->inSubBoxRange(sub_box)) {
        id_nodes.push_back(id_node_for_move);
        success = false;
    }
    return success;
}

bool Cleng::IsCommensuratable() {
    bool success = true;
    int length = (int) In[0]->MolList.size();

    // TODO: procedure goes throw all Molecules in the System...
    // TODO: Need to take only Polymer molecules...
    for (int i = 0; i < length; i++) {
        int chain_length = Mol[i]->chainlength;
        if (chain_length < 10) continue;

        for (auto &&SN : Enumerate(simpleNodeList)) {
            auto p3 = SN.second->get_cnode()->get_system_point() - SN.second->get_system_point();
            int path_length = abs(p3.x) + abs(p3.y) + abs(p3.z);

            int path_length_even = path_length % 2;
            int chain_length_even = chain_length % 2;

//            cout << "path_length:  "  << path_length << endl;
//            cout << "chain_length: " << chain_length << endl;
//            if (path_length > 20) cout << "Ooowww!! path_len > 20! "<< endl;
//            if (chain_length > 20) cout << "Ooowww!! chain_len > 20! "<< endl;

            if (path_length_even == chain_length_even) success = false;
            if (path_length >= chain_length) success = false;
        }
    }
    if (!success) {
        cout << endl;
        cout << "[WARNING] paths between clamps is not commensurate with the chain length!" << endl;
    }
    return success;
}

bool Cleng::NotCollapsing(int id_node_for_move) {
    bool not_collapsing = true;
    Point shifted_point(nodes_map[id_node_for_move].data()->get()->point());
    double min_dist = 0; // minimal distance between nodes_map. 0 for a while...

    for (auto &&n : nodes_map) {
        Point test_point = n.second.data()->get()->point();
        if (shifted_point != test_point) {
            Real distance = shifted_point.distance(n.second.data()->get()->point());
            if (distance <= min_dist) {
                cout << endl;
                cout << "Nodes are too close to each other." << endl;
                cout << "Shifted point from nodes_map: "     << shifted_point.to_string() << endl;
                cout << "Test point from nodes_map:    "     << test_point.to_string()    << endl;
                not_collapsing = false;
            }
        }
    }
    return not_collapsing;
}

bool Cleng::InRange(int id_node_for_move) {
    bool in_range = false;
    Point down_boundary = {1, 1, 1};
    Point up_boundary = box - down_boundary;
    Point shifted_point(nodes_map[id_node_for_move].data()->get()->point());
    if ((down_boundary.all_elements_less_than(shifted_point)) and
        (shifted_point.all_elements_less_than(up_boundary)))
        in_range = true;
    return in_range;
}

void signalHandler(int signum) {
    cout << "---> Termination..." << endl;
    cleng_flag_termination++;
    if (cleng_flag_termination > 1) exit(signum);
}

void Cleng::make_BC() {
//    make boundary condition point => BC = {1,1,1} => minor; BC = {0,0,0} => periodic
    BC = {
            Lat[0]->BC[0] == "mirror",
            Lat[0]->BC[2] == "mirror",
            Lat[0]->BC[4] == "mirror",
    };
}

vector<Real> Cleng::prepare_vtk() {
    int Size = 0;
    vector<Real> vtk;
    Real *X = Out[0]->GetPointer(Out[0]->OUT_key[0], Out[0]->OUT_name[0], Out[0]->OUT_prop[0], Size);
    for (int i = 1; i < box.x + 1; i++) {
        for (int j = 1; j < box.y + 1; j++) {
            for (int k = 1; k < box.z + 1; k++) {
                vtk.push_back(X[i * J.x + j * J.y + k]);
            }
        }
    }
    return vtk;
}

int Cleng::prepareIdNode() {
    if (debug) cout << "prepareIdNode in Cleng" << endl;
    int id_node_for_move = 0;
    if (ids_node4move[0] != -1) { id_node_for_move = ids_node4move[rand.getInt(0, (int) ids_node4move.size() - 1)]; }
    else { id_node_for_move = rand.getInt(0, (int) nodes_map.size() - 1); }

    if (id_node_for_move > (int) nodes_map.size() - 1) {
        cout << "###" << endl;
        cout << "Node id is too high. Trying to move node id: " << id_node_for_move << "." << endl;
        cout << "Available nodes_map: " << endl;
        for (auto &&n: nodes_map)
            cout << "id: " << n.first << " " << n.second.data()->get()->to_string() << endl;
        cout << "Termination..." << endl;
        exit(0);
    }

    return id_node_for_move;
}

void Cleng::prepareIdsNode() {
    if (debug) cout << "prepareIdsNode in Cleng" << endl;
    pivot_node_ids.clear();
    int _arm = rand.getInt(1, pivot_arms);
    int _start_index = rand.getInt(0, pivot_arm_nodes[_arm].size()-2);
    for (size_t i = _start_index; i<pivot_arm_nodes[_arm].size(); i++)
        pivot_node_ids.push_back(pivot_arm_nodes[_arm][i]);
//    for (auto &&iD: pivot_node_ids) cout << "id: " << iD << endl;
}

Point Cleng::prepareMove(const string& type_move) {
    if (debug) cout << "prepareMove in Cleng" << endl;
    Point clamped_move;

    if (type_move == "pivot_move") {
        prepareIdsNode();
        prepareRotationMatrix<int>();
//        prepareScalingMatrix<Real>();

    } else {
        if (type_move == "pivot_one_bond_move") {prepareIdsNode();}
        if (axis) {
            int c1 = 1;  // sing of movement +1 or -1
            int c2 = 2;  // MC step by default 2; could be 1 in two_ends_extension mode

            if (two_ends_extension) c2 = 1;
            if (sign_move == "-") c1 = -1;
            // currently it is implemented for step
            if (axis == 1) clamped_move = {c1 * c2, 0, 0};
            if (axis == 2) clamped_move = {0, c1 * c2, 0};
            if (axis == 3) clamped_move = {0, 0, c1 * c2};
        } else {
            clamped_move = {
                    rand.getIntExcludeArray(-delta_step, delta_step, makeExcludedArray(delta_step)),
                    rand.getIntExcludeArray(-delta_step, delta_step, makeExcludedArray(delta_step)),
                    rand.getIntExcludeArray(-delta_step, delta_step, makeExcludedArray(delta_step)),
            };

            if ((delta_step % 2) == 1) {  // 1 3 5
                int id4rm = rand.getInt(0, 2);
                if (id4rm == 0) clamped_move.x = 0;
                else { if (id4rm == 1) clamped_move.y = 0; else clamped_move.z = 0; }
            } else {  // 2 4 6
                int id4rm1 = rand.getInt(0, 2);
                int id4rm2 = rand.getInt(0, 2);

                if (id4rm1 == 0) clamped_move.x = 0;
                else { if (id4rm1 == 1) clamped_move.y = 0; else clamped_move.z = 0; }

                if (id4rm2 == 0) clamped_move.x = 0;
                else { if (id4rm2 == 1) clamped_move.y = 0; else clamped_move.z = 0; }
            }
        }
    }
    return clamped_move;
}

template<class T>
Matrix<T> Cleng::prepareRotationMatrix() {
    if (debug) cout << "prepareRotationMatrix in Cleng" << endl;
    int pivot_coef = rand.getIntExcludeArray(-1, 1, {0});
    if (pivot_axis == -1) {
        int pivot_axis_current = rand.getInt(1, 3);
        rotation_matrix = _create_rotational_matrix<T>(pivot_axis_current, pivot_move*pivot_coef);
    } else rotation_matrix = _create_rotational_matrix<T>(pivot_axis, pivot_move*pivot_coef);
    return rotation_matrix;
}

template<class T>
Matrix<T> Cleng::prepareScalingMatrix() {
    if (debug) cout << "prepareScalingMatrix in Cleng" << endl;
    vector<Real> scaling_coeff = {0.5, 1.0, 2.0};
    int scaling_index = rand.getInt(0, 2);
    scaling_matrix = _create_scaling_matrix<Real>(scaling_coeff[scaling_index]);

    cout << "Scaling matrix" << endl;
    cout << scaling_matrix << endl;
    return scaling_matrix;
}

template<class T>
Matrix<T> Cleng::_create_rotational_matrix(int axis_rotation, int grad) {
    Matrix<T> rotation_matrix_(3, 3);
    switch (axis_rotation) {
        case 1:
            rotation_matrix_.put(0, 0, 1);
            rotation_matrix_.put(0, 1, 0);
            rotation_matrix_.put(0, 2, 0);

            rotation_matrix_.put(1, 0, 0);
            rotation_matrix_.put(1, 1, cos(grad * PIE / 180.0));
            rotation_matrix_.put(1, 2, -sin(grad * PIE / 180.0));

            rotation_matrix_.put(2, 0, 0);
            rotation_matrix_.put(2, 1, sin(grad * PIE / 180.0));
            rotation_matrix_.put(2, 2, cos(grad * PIE / 180.0));

            break;
        case 2:
            rotation_matrix_.put(0, 0, cos(grad * PIE / 180.0));
            rotation_matrix_.put(0, 1, 0);
            rotation_matrix_.put(0, 2, sin(grad * PIE / 180.0));

            rotation_matrix_.put(1, 0, 0);
            rotation_matrix_.put(1, 1, 1);
            rotation_matrix_.put(1, 2, 0);

            rotation_matrix_.put(2, 0, -sin(grad * PIE / 180.0));
            rotation_matrix_.put(2, 1, 0);
            rotation_matrix_.put(2, 2, cos(grad * PIE / 180.0));

            break;
        case 3:
            rotation_matrix_.put(0, 0, cos(grad * PIE / 180.0));
            rotation_matrix_.put(0, 1, -sin(grad * PIE / 180.0));
            rotation_matrix_.put(0, 2, 0);

            rotation_matrix_.put(1, 0, sin(grad * PIE / 180.0));
            rotation_matrix_.put(1, 1, cos(grad * PIE / 180.0));
            rotation_matrix_.put(1, 2, 0);

            rotation_matrix_.put(2, 0, 0);
            rotation_matrix_.put(2, 1, 0);
            rotation_matrix_.put(2, 2, 1);
            break;

        default:
            cout << "[Warning] Strange axis: " << axis_rotation << endl;
    }
    return rotation_matrix_;
}

template<class T>
Matrix<T> Cleng::_create_scaling_matrix(Real scaling_coef) {
    Matrix<T> scaling_matrix_(3, 3);
    scaling_matrix_.put(0, 0, 1*scaling_coef);
    scaling_matrix_.put(0, 1, 0);
    scaling_matrix_.put(0, 2, 0);

    scaling_matrix_.put(1, 0, 0);
    scaling_matrix_.put(1, 1, 1*scaling_coef);
    scaling_matrix_.put(1, 2, 0);

    scaling_matrix_.put(2, 0, 0);
    scaling_matrix_.put(2, 1, 0);
    scaling_matrix_.put(2, 2, 1*scaling_coef);
    return scaling_matrix_;
}

int Cleng::getLastMCS() {

    int MC_attemp = 0;

    //read kal file
    ifstream infile(filename + ".kal");

    if (infile) {
        string line;
        while (infile >> ws && getline(infile, line)); // skip empty lines

        istringstream iss(line);
        vector<string> results(istream_iterator<string>{iss}, istream_iterator<string>());
        MC_attemp = stoi(results[0]);
    } else cout << "Unable to open kal file.\n";

    return MC_attemp;
}

void Cleng::WriteOutput() {
    if (debug) cout << "WriteOutput in Cleng" << endl;
    cout << "[Cleng] Saving... " << MC_attempt + MCS_checkpoint << endl;
    PushOutput();
    New[0]->PushOutput();
    for (int i = 0; i < n_out; i++) Out[i]->WriteOutput(MC_attempt + MCS_checkpoint);
}

void Cleng::PushOutput() {

    for (int i = 0; i < n_out; i++) {
        Out[i]->PointerVectorInt.clear();
        Out[i]->PointerVectorReal.clear();
        Out[i]->SizeVectorInt.clear();
        Out[i]->SizeVectorReal.clear();
        Out[i]->strings.clear();
        Out[i]->strings_value.clear();
        Out[i]->bools.clear();
        Out[i]->bools_value.clear();
        Out[i]->Reals.clear();
        Out[i]->Reals_value.clear();
        Out[i]->ints.clear();
        Out[i]->ints_value.clear();

        if (Out[i]->name == "ana" || Out[i]->name == "kal") {
            Out[i]->push("MC_attempt", MC_attempt + MCS_checkpoint);
            Out[i]->push("MCS", MCS);
            Out[i]->push("free_energy", free_energy_current);
        }

        if (Out[i]->name == "ana" ||
            Out[i]->name == "vec") { //example for putting an array of Reals of arbitrary length to output
            string s = "vector;0"; //the keyword 'vector' is used for Reals; the value 0 is the first vector, use 1 for the next etc,
            Out[i]->push("gn", s); //"gn" is the name that will appear in output file
            Out[i]->PointerVectorReal.push_back(
                    Mol[clp_mol]->gn); //this is the pointer to the start of the 'vector' that is reported to output.
            Out[i]->SizeVectorReal.push_back(
                    sizeof(Mol[clp_mol]->gn)); //this is the size of the 'vector' that is reported to output
        }

        if (Out[i]->name == "ana" ||
            Out[i]->name == "pos") { //example for putting 3 array's of Integers of arbitrary length to output
            string s = "array;0";
            Out[i]->push("X", s);

//			 TODO using normal point in output and remove xs, ys, zs
            fillXYZ();
//			point=X.data();
            Out[i]->PointerVectorInt.push_back(xs);
            Out[i]->SizeVectorInt.push_back((int) nodes_map.size());
            s = "array;1";
            Out[i]->push("Y", s);
//			point = Y.data();
            Out[i]->PointerVectorInt.push_back(ys);
            Out[i]->SizeVectorInt.push_back((int) nodes_map.size());
            s = "array;2";
            Out[i]->push("Z", s);
//			point=Z.data();
            Out[i]->PointerVectorInt.push_back(zs);
            Out[i]->SizeVectorInt.push_back((int) nodes_map.size());
        }
    }
}

void Cleng::fillXYZ() {
    delete[] xs;
    delete[] ys;
    delete[] zs;

    int n = nodes_map.size();
    xs = new int[n];
    ys = new int[n];
    zs = new int[n];
    for (int i = 0; i < n; i++) {
        xs[i] = nodes_map[i].data()->get()->point().x;
        ys[i] = nodes_map[i].data()->get()->point().y;
        zs[i] = nodes_map[i].data()->get()->point().z;
    }
}

void Cleng::WriteClampedNodeDistance() {
    vector<Real> distPerMC;
    int i = 0;
    for (auto &&SN :  simpleNodeList) {
        if (!(i % 2)) distPerMC.push_back(SN->distance(SN->get_cnode()->get_system_point()));
        i++;
    }

    ofstream outfile;
    outfile.open(filename + ".cdis", std::ios_base::app);
    outfile << MC_attempt + MCS_checkpoint << " ";
    for (auto n : distPerMC)outfile << n << " ";
    outfile << endl;
}

void Cleng::WriteClampedNodePosition() {

    ofstream outfile;
    outfile.open(filename + ".cpos", std::ios_base::app);

//    cout
    outfile
            << "#step " << MC_attempt + MCS_checkpoint << " {X, Y, Z} #" << endl;

    int i = 0;
    for (auto &&SN :  simpleNodeList) {
        if (!(i % 2)) {
//            cout
            outfile
                    << SN->get_system_point().to_string() << SN->get_cnode()->get_system_point().to_string() << endl;
//            << SN->get_system_point().to_string() << endl;
        }
        i++;
    }
}


Real Cleng::GetN_times_mu() {
    int n_mol = (int) In[0]->MolList.size();
    Real n_times_mu = 0;
    for (int i = 0; i < n_mol; i++) {
        Real Mu = Mol[i]->Mu;
        Real n = Mol[i]->n;
        if (Mol[i]->IsClamped()) n = Mol[i]->n_box;
        n_times_mu += n * Mu;
    }
    return n_times_mu;
}