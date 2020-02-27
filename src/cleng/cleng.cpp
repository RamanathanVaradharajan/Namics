#include "cleng.h"
#include "cleng_tools.h"

using namespace std;

Cleng::Cleng(
        vector<Input *> In_,
        vector<Lattice *> Lat_,
        vector<Segment *> Seg_,
        vector<State *> Sta_,
        vector<Reaction *> Rea_,
        vector<Molecule *> Mol_,
        vector<System *> Sys_,
        vector<Solve_scf *> New_,
        string name_
) : name(std::move(name_)),
    In(std::move(In_)),
    Lat(std::move(Lat_)),
    Seg(std::move(Seg_)),
    Sta(std::move(Sta_)),
    Rea(std::move(Rea_)),
    Mol(std::move(Mol_)),
    Sys(std::move(Sys_)),
    New(std::move(New_)) {

    if (debug) cout << "Cleng initialized" << endl;
    KEYS.emplace_back("MCS");
    KEYS.emplace_back("delta_step");
    KEYS.emplace_back("delta_save");
    KEYS.emplace_back("save_filename");
    KEYS.emplace_back("seed");
    KEYS.emplace_back("checkpoint_save");
    KEYS.emplace_back("checkpoint_load");
    KEYS.emplace_back("cleng_pos");
    KEYS.emplace_back("cleng_dis");
    KEYS.emplace_back("simultaneously");
    KEYS.emplace_back("movement_along");
    KEYS.emplace_back("sign_move");
    KEYS.emplace_back("user_node_id_move");
    KEYS.emplace_back("two_ends_extension");
    KEYS.emplace_back("metropolis");
    KEYS.emplace_back("prefactor_kT");
    KEYS.emplace_back("pivot_move");
    KEYS.emplace_back("pivot_axis");
    KEYS.emplace_back("pivot+one_node");
    KEYS.emplace_back("pivot+one_bond");

    // Debug.log
    //out.open("debug.out", ios_base::out);
}

Cleng::~Cleng() {
    delete[] xs;
    delete[] ys;
    delete[] zs;
    // Debug.log closing
    //out.close();
}

string Cleng::GetValue(const string& parameter) {
    int i = 0;
    int length = (int) PARAMETERS.size();
    while (i < length) {
        if (parameter == PARAMETERS[i]) {
            return VALUES[i];
        }
        i++;
    }
    return "";
}

bool Cleng::CheckInput(int start, bool save_vector) {
    if (debug) cout << "CheckInput in Cleng" << endl;
    bool success;

    success = In[0]->CheckParameters("cleng", name, start, KEYS, PARAMETERS, VALUES);
    if (success) {

        // MCS
        if (!GetValue("MCS").empty()) {
            success = In[0]->Get_int(GetValue("MCS"), MCS, 1, 10000000, "The number of Monte Carlo steps should be between 1 and 10000000");
            if (!success) { cout << "MCS will be equal to 5" << endl; MCS = 1; }
        } else MCS = 5;
        if (debug) cout << "MCS is " << MCS << endl;

        // seed
        if (!GetValue("seed").empty()) {
            success = In[0]->Get_int(GetValue("seed"), pseed, 1, 1000, "The seed should be between 1 and 1000");
            if (!success) { cout << "The seed will be equal 1" << endl; pseed = 1; }
        } else pseed = 0;
        rand = pseed==0 ? Random() : Random(pseed);
        if (debug) cout << "seed is " << pseed << endl;

        // delta_step
        if (!GetValue("delta_step").empty()) {
            success = In[0]->Get_int(GetValue("delta_step"), delta_step, 1, 5, "The number of delta_step should be between 1 and 5");
            if (!success) { cout << "The delta_step will be equal 1" << endl; delta_step = 1;}
        } else delta_step = 0;
        if (debug) cout << "delta_step is " << delta_step << endl;

        // pivot_move
        if (!GetValue("pivot_move").empty()) {
            success = In[0]->Get_int(GetValue("pivot_move"), pivot_move, 1, 360, "The angle of pivot_move should be between 1 and 360");
            if (!success) { cout << "The pivot_move will be disable." << endl; pivot_move = 0;}
        } else pivot_move = 0;
        if (debug) cout << "pivot_move is " << pivot_move << endl;

        // pivot_axis
        if (pivot_move) {
            if (!GetValue("pivot_axis").empty()) {
                success = In[0]->Get_int(GetValue("pivot_axis"), pivot_axis, 1, 3,
                                         "The axis of pivot_move should be between 1 and 3");
                if (!success) {
                    cout << "The pivot_axis will be all axis." << endl;
                    pivot_axis = -1;
                }
            } else pivot_axis = -1;
            if (debug) cout << "pivot_axis is " << pivot_axis << endl;
        } else {pivot_axis = 0;} // pivot_axis = 0 means disable pivot.

        // pivot+one_node
        if (!GetValue("pivot+one_node").empty()) {
            pivot_one_node = In[0]->Get_bool(GetValue("pivot+one_node"), false);
        } else pivot_one_node = false;
        if (debug) cout << "pivot+one_node " << pivot_one_node << endl;
        if ((pivot_one_node) and (!pivot_move)) {
            cout << "Sorry, if you like to use combination pivot+one_node move, please, specify pivot movement "
                    "by parameter 'pivot_move'"
                    "\n Termination..." << endl;
            exit(0);
        }

        // pivot+one_bond
        if (!GetValue("pivot+one_bond").empty()) {
            pivot_one_bond = In[0]->Get_bool(GetValue("pivot+one_bond"), false);
        } else pivot_one_bond = false;
        if (debug) cout << "pivot+one_node " << pivot_one_bond << endl;
        if ((pivot_one_bond) and (!pivot_move)) {
            cout << "Sorry, if you like to use combination pivot+one_bond move, please, specify pivot movement "
                    "by parameter 'pivot_move'"
                    "\n Termination..." << endl;
            exit(0);
        }

        if (pivot_one_bond) pivot_one_node = false;
        // TODO: think about user
//        if ((delta_step) and (pivot_move) and ()) {
//            cout << "Sorry, but you have to choose either one_node_move by `delta_step` parameter or `pivot_move`"
//                    "\n Termination..." << endl;
//            exit(0);
//        }

        // delta_save
        if (!GetValue("delta_save").empty()) {
            success = In[0]->Get_int(GetValue("delta_save"), delta_save, 1, MCS, "The delta_save interval should be between 1 and " + to_string(MCS));
        } else delta_save = 1;
        if (debug) cout << "delta_save_interval " << delta_save << endl;

        // Cleng molecules
        if (Sys[0]->SysClampList.empty()) {
            cout << "Cleng needs to have clamped molecules in the system" << endl;
            success = false;
        } else {
            clamp_seg = Sys[0]->SysClampList[0];
            if (Sys[0]->SysClampList.size() > 1) {
                success = false;
                cout << "Currently the clamping is limited to one molecule per system. " << endl;
            }
        }

        // checkpoint save
        if (!GetValue("checkpoint_save").empty()) {checkpoint_save = In[0]->Get_bool(GetValue("checkpoint_save"), false);}
        else checkpoint_save = false;
        if (debug) cout << "checkpoint_save " << checkpoint_save << endl;

        // checkpoint load
        if (!GetValue("checkpoint_load").empty()) {checkpoint_load = In[0]->Get_bool(GetValue("checkpoint_load"), false);}
        else checkpoint_load = false;
        if (debug) cout << "checkpoint_load " << checkpoint_load << endl;

        // saving cleng position of nodes_map
        if (!GetValue("cleng_pos").empty()) {cleng_pos = In[0]->Get_bool(GetValue("cleng_pos"), false);}
        else cleng_pos = false;
        if (debug) cout << "cleng_pos " << cleng_pos << endl;

        // saving distance between of nodes_map
        if (!GetValue("cleng_dis").empty()) {cleng_dis = In[0]->Get_bool(GetValue("cleng_dis"), false);}
        else cleng_dis = false;
        if (debug) cout << "cleng_dis " << cleng_dis << endl;

        // simultaneously
        if (!GetValue("simultaneously").empty()) simultaneously = In[0]->Get_bool(GetValue("simultaneously"), false);
        else simultaneously = false;
        if (debug) cout << "simultaneously move " << simultaneously << endl;

        // movement_along
        if (!GetValue("movement_along").empty()) {
            cout << "Warning!!! In movement_along mode delta step will be ignored! Monte Carlo step will be {sign_move*2} depending on axis " << endl;
            success = In[0]->Get_int(GetValue("movement_along"), axis, 1, 3, "The number of delta_step should be between 1 and 3");
            if (!success) {
                cout << "Sorry, you provide incorrect axis number. The movement_along will be disable" << endl;
                axis = 0;
                success = true;
            }
        } else axis = 0;
        if (debug) cout << "movement_along axis " << axis << endl;

        // sign_move
        vector<string> options {"+", "-"};
        if (axis) {
            if (!GetValue("sign_move").empty()) {
                success = In[0]->Get_string(GetValue("sign_move"), sign_move, options, "The sigh could be ether + or -");
            } else { sign_move = "+";}
        } else {
            if (!GetValue("sign_move").empty()) {
                cout << "Sorry, but you cannot use sign_move without movement_along axis flag" << endl;
                success = false;
                return success;
            }
        }
        if (debug) cout << "sign_move " << sign_move << endl;

        // user_node_move_id
        string struser_node_move_id;
        if (!GetValue("user_node_id_move").empty()) {
            success = In[0]->Get_string(GetValue("user_node_id_move"), struser_node_move_id, "");
            string node_id;
            string _ ;
            stringstream stream(struser_node_move_id);
            while( getline(stream, node_id, ',') ) {
                if (node_id[0] == '~') {
                    _ = node_id.erase(0,1);
                    try { ids_node4fix.push_back(stoi(_)); }
                    catch (invalid_argument &e) {
                        success = false;
                        cout << "Check yours node id structure" << endl;
                        return success;
                    }
                } else {
                    try { ids_node4move.push_back(stoi(node_id)); }
                    catch (invalid_argument &e) {
                        success = false;
                        cout << "Check yours node id structure" << endl;
                        return success;
                    }
                }
            }
        }
        else ids_node4move = {-1};
        if (debug) for (auto &&id:ids_node4move) cout << "user_node_id_move: " << id << endl;
        if (debug) for (auto &&id:ids_node4fix)  cout << "user_node_id_fix:  "  << id << endl;

        // 2 end extension mode
        if (!GetValue("two_ends_extension").empty()) {two_ends_extension = In[0]->Get_bool(GetValue("two_ends_extension"), false);}
        else two_ends_extension = false;
        if (debug) cout << "two_ends_extension " << two_ends_extension << endl;

        // checkpoint save
        if (!GetValue("metropolis").empty()) {metropolis = In[0]->Get_bool(GetValue("metropolis"), false);}
        else metropolis = true;
        if (debug) cout << "metropolis " << metropolis << endl;

        // prefactor kT constant
        if (!GetValue("prefactor_kT").empty()) {
            success = In[0]->Get_Real(GetValue("prefactor_kT"), prefactor_kT, 0, 10,
                    "Prefactor_kT is a constant in Metropolis algorithm (-1/C1) * (delta Fs/kT), where C1 - prefactor.\n"
                    "Currently available range is from 0 to 10.");
            if (!success) {return success;}
        } else prefactor_kT = 1;
        if (debug) cout << "prefactor_kT is " << prefactor_kT << endl;

        // TODO: EXTEND CLENG
        if (success) {
            n_boxes = Seg[clamp_seg]->n_box;
            sub_box_size = {Seg[clamp_seg]->mx, Seg[clamp_seg]->my, Seg[clamp_seg]->mz};
        }
        clp_mol = -1;
        int length = (int) In[0]->MolList.size();
        for (int i = 0; i < length; i++) if (Mol[i]->freedom == "clamped") clp_mol = i;
    }

    if (success) {
        n_out = (int) In[0]->OutputList.size();
        if (n_out == 0) cout << "Warning: no output defined!" << endl;

        for (int i = 0; i < n_out; i++) {
            Out.push_back(new Output(In, Lat, Seg, Sta, Rea, Mol, Sys, New, In[0]->OutputList[i], i, n_out));
            if (!Out[i]->CheckInput(start)) {
                cout << "input_error in output " << endl;
                success = false;
                return success;
            }
        }

        vector<string> sub;
        In[0]->split(In[0]->name, '.', sub);
        filename = In[0]->output_info.getOutputPath() + sub[0];

        cout << cmajor_version << "." << cminor_version << "." << cpatch << " v." << cversion << endl;
        begin_simulation = std::chrono::steady_clock::now();
        success = MonteCarlo(save_vector);
    }
    return success;
}

void Cleng::update_ids_node4move() {
    vector<int> _ = ids_node4move;
    if (!ids_node4fix.empty()) {
        ids_node4move.clear();
        for (auto const &n_: nodes_map) {
            bool exists = any_of(begin(ids_node4fix), end(ids_node4fix), [&](int i) { return i == n_.first; });
            if (!exists) ids_node4move.push_back(n_.first);
        }
    }
}

bool Cleng::CP(transfer tofrom) {
    if (debug) cout << "CP in Cleng" << endl;

    bool success = true;

    Segment *clamped = Seg[clamp_seg];
    map<int, Point> system_points;
    switch (tofrom) {
        case to_cleng:
            simpleNodeList.clear();
            for (int i = 0; i < n_boxes; i++) {
                auto first_node  = fromSystemToNode(clamped->px1[i], clamped->py1[i], clamped->pz1[i], 2 * i, box);
                auto second_node = fromSystemToNode(clamped->px2[i], clamped->py2[i], clamped->pz2[i], 2 * i + 1, box);
                first_node ->set_cnode(second_node);
                second_node->set_cnode(first_node);
                //
                simpleNodeList.push_back(first_node);
                simpleNodeList.push_back(second_node);}
            nodes_map = createNodes(simpleNodeList, pivot_arm_nodes, pivot_arms);
            break;

        case to_segment:
            //Zero(clamped->H_MASK, Lat[0]->M);  [CLENG]
            // merged:
            std::fill(clamped->H_MASK, clamped->H_MASK+Lat[0]->M, 0);
            //std::fill(Seg[clamp_seg]->H_MASK, Seg[clamp_seg]->H_MASK+M, 0); [Master]

            for (auto &&SN : Enumerate(simpleNodeList)) {
                size_t index = SN.first;    //
                if (index % 2 == 1) continue;
                size_t n_box = index / 2;   // n_boxes
                Point p1 = SN.second->point();                // all time in box [0:box_l]
                Point p2 = SN.second->get_cnode()->point();   // all time in box [0:box_l]

                Point p1sys  = SN.second->get_system_point();
                Point p2sys  = SN.second->get_cnode()->get_system_point();

//                cout << " n_box: " << n_box << endl;
//                cout << " index: " << index << " v: " << p1.x << " " << "("<<p1sys.x<<")" << " " << p1.y << " "  << "("<<p1sys.y<<")" << " " << p1.z << " "  << "("<<p1sys.z<<")" << endl;
//                cout << " index: " << index << " v: " << p2.x << " " << "("<<p2sys.x<<")" << " " << p2.y << " "  << "("<<p2sys.y<<")" << " " << p2.z << " "  << "("<<p2sys.z<<")" << endl;

                if ((p1.x - p2.x) != (p1sys.x - p2sys.x)) {if (p1.x < p2.x) p1.x += box.x; else p2.x += box.x;}
                if ((p1.y - p2.y) != (p1sys.y - p2sys.y)) {if (p1.y < p2.y) p1.y += box.y; else p2.y += box.y;}
                if ((p1.z - p2.z) != (p1sys.z - p2sys.z)) {if (p1.z < p2.z) p1.z += box.z; else p2.z += box.z;}

                clamped->px1[n_box] = p1.x;
                clamped->py1[n_box] = p1.y;
                clamped->pz1[n_box] = p1.z;

                clamped->px2[n_box] = p2.x;
                clamped->py2[n_box] = p2.y;
                clamped->pz2[n_box] = p2.z;
// box
                clamped->bx[n_box] = (p2.x + p1.x - sub_box_size.x) / 2;
                clamped->by[n_box] = (p2.y + p1.y - sub_box_size.y) / 2;
                clamped->bz[n_box] = (p2.z + p1.z - sub_box_size.z) / 2;

                if (clamped->bx[n_box] < 1) { clamped->bx[n_box] += box.x; clamped->px1[n_box] += box.x; clamped->px2[n_box] += box.x;}
                if (clamped->by[n_box] < 1) { clamped->by[n_box] += box.y; clamped->py1[n_box] += box.y; clamped->py2[n_box] += box.y;}
                if (clamped->bz[n_box] < 1) { clamped->bz[n_box] += box.z; clamped->pz1[n_box] += box.z; clamped->pz2[n_box] += box.z;}
// clearing
                auto hp1x = ((clamped->px1[n_box] - 1) % box.x + 1) * J.x;
                auto hp1y = ((clamped->py1[n_box] - 1) % box.y + 1) * J.y;
                auto hp1z =  (clamped->pz1[n_box] - 1) % box.z + 1;

                auto hp2x = ((clamped->px2[n_box] - 1) % box.x + 1) * J.x;
                auto hp2y = ((clamped->py2[n_box] - 1) % box.y + 1) * J.y;
                auto hp2z = (clamped->pz2[n_box] - 1) % box.z + 1;

                clamped->H_MASK[hp1x + hp1y + hp1z] = 1;
                clamped->H_MASK[hp2x + hp2y + hp2z] = 1;
            }
            break;

        default:
            success = false;
            cout << "error in transfer" << endl;
            break;
    }
    return success;
}

bool Cleng::Checks(int id_node_for_move) {
    bool in_range;

    bool in_subbox_range = InSubBoxRange(id_node_for_move);

    bool not_collapsing = NotCollapsing(id_node_for_move);

    // check distance between all nodes_map and constrains (walls)
    // BC.x, BC.y, BC.z = mirror
    if (BC.x and BC.y and BC.z) in_range = InRange(id_node_for_move);
    // BC.x and/or BC.y and/or BC.z != mirror
    else in_range = true;

    // check distances between all nodes_map => checking possibility to achieve clamped nodes
    bool commensurate = IsCommensuratable();

    bool result = not_collapsing and in_range and in_subbox_range and commensurate;
//    cout << "[Checks] result: " << result << endl;
    return result;
}

bool Cleng::MakeChecks(int id_node_for_move) {
    bool success = true;
    while (!Checks(id_node_for_move)) {
        cout << internal_name << "Prepared MC step for the node ["<< id_node_for_move << "] does not pass checks. Rejected." << endl;
        MakeMove(true);
        rejected++;
        cleng_rejected++;
        success = false;
        return success; // TODO: CHECK THIS RETURN
    }
    return success;
}

void Cleng::_moveClampedNode(bool back, int id_node_for_move, const Point& clamped_move) {
    if (back) {
        Point _clamped_move = clamped_move.negate();
        if (!simultaneously) nodes_map[id_node_for_move].data()->get()->shift(_clamped_move);
        else for (auto &node : nodes_map) node.second.data()->get()->shift(_clamped_move);
    } else {
        if (!simultaneously) {
            cout << "[Prepared node] id: " << id_node_for_move << " clamped_move: " << clamped_move.to_string() << "... ";
            nodes_map[id_node_for_move].data()->get()->shift(clamped_move);
            cout << "[Moved]" << endl;
        } else {
            // simultaneously procedure (+ two_ends_extension)
            if (two_ends_extension) {
                // TODO: deprecate in future...
                for (auto &&node : nodes_map) {
                    int index = node.first;
                    if (index % 2 == 0) node.second.data()->get()->shift(clamped_move.negate());
                    else node.second.data()->get()->shift(clamped_move);
                }
                cout << "[Moved] " << "*All* "
                << "MC step: " << clamped_move.to_string() << " and " << clamped_move.negate().to_string() << endl;
            } else {
                for (auto &&node : nodes_map) node.second.data()->get()->shift(clamped_move);
                cout << "[Moved] " << "*All* " << "MC step: " << clamped_move.to_string() << endl;
            }
        }
    }
}

Point Cleng::preparePivotClampedMove(int id_node_for_move) {
    Point center_of_rotation = nodes_map[pivot_node_ids[0]].data()->get()->point();
    // moving center of rotation to the center
    for (auto &&node : nodes_map) node.second.data()->get()->shift(center_of_rotation.negate());
    Point current_Point = nodes_map[id_node_for_move].data()->get()->_returnSystemPoint();
    Point point_shifted_by_matrix = rotation_matrix.dot(current_Point);;
    Point clamped_move = point_shifted_by_matrix - current_Point;
    // center +1
    if (clamped_move.x < -(box.x / 2) +1) clamped_move.x += box.x;
    if (clamped_move.y < -(box.y / 2) +1) clamped_move.y += box.y;
    if (clamped_move.z < -(box.z / 2) +1) clamped_move.z += box.z;

    if (clamped_move.x > (box.x / 2)+1) clamped_move.x -= box.x;
    if (clamped_move.y > (box.y / 2)+1) clamped_move.y -= box.y;
    if (clamped_move.z > (box.z / 2)+1) clamped_move.z -= box.z;
    // moving center of rotation to initial place
    for (auto &&node : nodes_map) node.second.data()->get()->shift(center_of_rotation);
    return clamped_move;
}

bool Cleng::MakeMove(bool back) {
    if (debug) cout << "MakeMove in Cleng" << endl;
    bool success = true;

    if (back) {
        cout << internal_name << "Moving back...";
        for (auto &&nodeID_clampedMove : nodeIDs_clampedMove) {
            _moveClampedNode(back, nodeID_clampedMove.first, nodeID_clampedMove.second);
        }
        cout << "[Moved back]" << endl;
    } else {
        nodeIDs_clampedMove.clear();
        if (pivot_move) {
            if (pivot_one_bond) {
                int type_move = rand.getInt(0, 1);
                if (type_move) success = _pivotMoveClampedNode(false);
                else success = _pivotMoveOneBond(false);
            }
            if (pivot_one_node) {
                int type_move = rand.getInt(0, 1);
                if (type_move) success = _pivotMoveClampedNode(false);
                else success = _oneNodeMoveClampedNode(false);
            }
            if (!pivot_one_node and !pivot_one_bond) {
                success = _pivotMoveClampedNode(false);
            }
        } else {
            success = _oneNodeMoveClampedNode(false);
        }
    }
    return success;
}

bool Cleng::_pivotMoveClampedNode(const bool &back) {
    Point clamped_move;
    bool success = true;
    cout << internal_name << "[pivot_move]" << endl;
    prepareMove("pivot_move");
    for (size_t node_position_in_vector=1; node_position_in_vector != pivot_node_ids.size(); node_position_in_vector++) {
        clamped_move = preparePivotClampedMove(pivot_node_ids[node_position_in_vector]);  // need to create for each node_id!
        _moveClampedNode(back, pivot_node_ids[node_position_in_vector], clamped_move);
        nodeIDs_clampedMove[pivot_node_ids[node_position_in_vector]] = clamped_move;
    }
    cout << internal_name << "checking positions...";
    for (auto &&node_id : pivot_node_ids) {
        success = MakeChecks(node_id);
        if (!success) return success;
    }
    cout << " OK" << endl;
    return  success;
}

bool Cleng::_pivotMoveOneBond(const bool &back) {
    Point clamped_move;
    bool success = true;
    cout << internal_name << "[pivot_one_bond]" << endl;
    clamped_move = prepareMove("pivot_one_bond_move");  // need to create at once and move all
    for (size_t node_position_in_vector=1; node_position_in_vector != pivot_node_ids.size(); node_position_in_vector++) {
        _moveClampedNode(back, pivot_node_ids[node_position_in_vector], clamped_move);
        nodeIDs_clampedMove[pivot_node_ids[node_position_in_vector]] = clamped_move;
    }
    cout << internal_name << "checking positions...";
    for (auto &&node_id : pivot_node_ids) {
        success = MakeChecks(node_id);
        if (!success) return success;
    }
    cout << " OK" << endl;
    return  success;
}

bool Cleng::_oneNodeMoveClampedNode(const bool &back) {
    int id_node_for_move;
    Point clamped_move;
    bool success;

    cout << internal_name << "[one_node_move]" << endl;
    id_node_for_move = prepareIdNode();
    clamped_move = prepareMove("one_node_move");
    _moveClampedNode(back, id_node_for_move, clamped_move);

    cout << internal_name << "checking positions...";
    success = MakeChecks(id_node_for_move);
    if (!success) return success;
    nodeIDs_clampedMove[id_node_for_move] = clamped_move;
    cout << " OK" << endl;
    return success;
}


// Main procedure
bool Cleng::MonteCarlo(bool save_vector) {
    if (debug) cout << "Monte Carlo in Cleng" << endl;
    bool success = true;

    signal(SIGINT, signalHandler);
#ifdef CLENG_EXPERIMENTAL
    CWriter cleng_writer(filename);
#endif
    Checkpoint checkpoint;
    if (checkpoint_load) {
        if (checkpoint.isLoadable()) {
            simpleNodeList = checkpoint.loadCheckpoint(simpleNodeList, box);
            nodes_map = createNodes(simpleNodeList, pivot_arm_nodes, pivot_arms);
            cout << internal_name << "From checkpoint next nodes are available: " << endl;
            for (auto &&n : nodes_map) cout << "id: " << n.first << " " << n.second.data()->get()->to_string() << endl;
            CP(to_segment);
            if (getLastMCS() != 0) MCS_checkpoint = getLastMCS() + 1;
            loaded = true;
            for (auto &&pair_pivot : pivot_arm_nodes) {
                cout << "---> arm: " << pair_pivot.first << " ids: ";
                for (auto &&ids: pair_pivot.second) cout << ids << " ";
                cout << endl;
            }
            if (pivot_move) {ids_node4fix.clear();ids_node4fix.push_back(pivot_arm_nodes[1][0]);}
        }
    }

// Analysis MC
    accepted = 0.0;
    rejected = 0.0;
    cleng_rejected = 0.0;
    MC_attempt = 0;
    make_BC();


// init system outlook
    if (!loaded) {
        CP(to_cleng);
        if (pivot_move) {ids_node4fix.clear();ids_node4fix.push_back(pivot_arm_nodes[1][0]);}
    }
    if (!Checks(0)) {cout << internal_name << "Some checks are not passed. Termination..." << endl; exit(1); }
    if (checkpoint_save) {checkpoint.saveCheckpoint(simpleNodeList);}

    bool success_iteration = New[0]->Solve(true);
    free_energy_current = Sys[0]->GetFreeEnergy();
    if (is_ieee754_nan(free_energy_current)) {
        cout << "%?% Sorry, Free Energy is NaN. Termination..." << endl;
        cout << "%?% Solver output " << success_iteration       << endl;
        exit(1);
    }
    if (save_vector) test_vector.push_back(Sys[0]->GetFreeEnergy());

// init save
    WriteOutput();
    if (cleng_dis) WriteClampedNodeDistance();
    if (cleng_pos) WriteClampedNodePosition();

#ifdef CLENG_EXPERIMENTAL
    vector<Real>vtk = prepare_vtk();
    cleng_writer.write("/VTK_data", "vtk"+to_string(MC_attempt+MCS_checkpoint), dims_vtk, vtk);

    n_times_mu = GetN_times_mu();
    vector<Real> MC_free_energy = {static_cast<Real>(MC_attempt+MCS_checkpoint), free_energy_current, free_energy_current-n_times_mu};
    cleng_writer.append("/Free_energy", "free_energy", dims_3, MC_free_energy);
#endif

    cout << "Initialization done.\n" << endl;
    update_ids_node4move();
    cout << internal_name << "System: " << endl;
    for (auto &&n : nodes_map) cout << "Node id: " << n.first << " " << n.second.data()->get()->to_string() << endl;
    cout << endl;
    if (pivot_move) {
        cout << internal_name << "System [pivot -> stars only]: " << endl;
        for (auto &&pair_pivot : pivot_arm_nodes) {
            cout << "---> arm: " << pair_pivot.first << " ids: ";
            for (auto &&ids: pair_pivot.second) cout << ids << " ";
            cout << endl;
        }
        cout << endl;
    }
    cout << internal_name <<  "Here we go..." << endl;
    bool success_;
    for (MC_attempt = 1; MC_attempt <= MCS; MC_attempt++) { // main loop for trials
        success_ = MakeMove(false);
        if (success_) {
            CP(to_segment);

            cout << endl;
            cout << internal_name << "System for calculation: " << endl;
            for (auto &&n : nodes_map) cout << "Node id: " << n.first << " " << n.second.data()->get()->to_string() << endl;
            cout << endl;
            if (pivot_move) {
                cout << internal_name << "System [pivot -> stars only]: " << endl;
                for (auto &&pair_pivot : pivot_arm_nodes) {
                    cout << "---> arm: " << pair_pivot.first << " ids: ";
                    for (auto &&ids: pair_pivot.second) cout << ids << " ";
                    cout << endl;
                }
                cout << endl;
            }

            success_iteration = New[0]->Solve(true);
            //// Simulation without rescue procedure ---> TODO: [1nd solution]
            if (is_ieee754_nan(Sys[0]->GetFreeEnergy())) {
                cout << "%?% Sorry, Free Energy is NaN. " << endl;
                cout << "%?% Here is result from solver: " << success_iteration << endl;

                cout << "%?% The step will be rejected! "
                        "Probably your system is too dense! "
                        "Simulation will continue... " << endl;
                cout << internal_name << "[CRASH STATE] " << "returning back the system configuration... " << endl;
                MakeMove(true);
                rejected++;
                cleng_rejected++;
                continue;
            } else {free_energy_trial = Sys[0]->GetFreeEnergy();}
            //// Simulation without rescue procedure <---

            //// Simulation with rescue procedure --->
//            if (is_ieee754_nan(Sys[0]->GetFreeEnergy())) {
//                cout << "%?% Sorry, Free Energy is NaN.  " << endl;
//                cout << "%?% Here is result from solver: " << success_iteration << endl;
////                New[0]->rescue_status = NONE;  # TODO need to rescue?
//                New[0]->attempt_DIIS_rescue();
//                cout << "%?% Restarting iteration." << endl;
//                success_iteration = New[0]->Solve(true);
//
//                if (is_ieee754_nan(Sys[0]->GetFreeEnergy())) {
//                    cout << "%?% Sorry, Free Energy is still NaN. " << endl;
//                    cout << "%?% Here is result from solver: " << success_iteration << endl;
//
//                    cout << "%?% The step will be rejected! "
//                            "Probably your system is too dense! "
//                            "Simulation will continue... " << endl;
//                    cout << internal_name << "[CRASH STATE] " << "returning back the system configuration... " << endl;
//                    MakeMove(true);
//                    rejected++;
//                    cleng_rejected++;
//                    continue;
//                }
//            } else {free_energy_trial = Sys[0]->GetFreeEnergy();}

            //// Simulation with rescue procedure <---
            if (save_vector) test_vector.push_back(Sys[0]->GetFreeEnergy());

            cout << "Free Energy (c): " << free_energy_current << endl;
            cout << "            (t): " << free_energy_trial   << endl;
            cout << "   prefactor kT: " << prefactor_kT        << endl;
            cout << endl;

            if (!metropolis) {
                cout << "Metropolis is disabled. " << endl;
                free_energy_current = free_energy_trial;
                accepted++;
            } else {
                if (free_energy_trial - free_energy_current <= 0.0) {
                    cout << internal_name << metropolis_name << "Accepted" << endl;
                    n_times_mu = GetN_times_mu();
                    free_energy_current = free_energy_trial;
                    accepted++;
                } else {
                    Real acceptance = rand.getReal(0, 1);

                    if (acceptance < exp((-1.0/prefactor_kT) * (free_energy_trial - free_energy_current))) {
                        cout << internal_name << metropolis_name << "Accepted with probability" << endl;
                        free_energy_current = free_energy_trial;
                        n_times_mu = GetN_times_mu();
                        accepted++;
                    } else {
                        cout << internal_name << metropolis_name << "Rejected" << endl;
                        MakeMove(true);
                        CP(to_segment);
                        cout << internal_name << metropolis_name << "returning back the system configuration... " << endl;
                        New[0]->Solve(true);
                        // TODO check cleng returned to normal numbers... [2nd solution]
                        rejected++;
                        cout << "... [Done]" << endl;
                    }
                }
            }
        }
        if (metropolis) {
            cout << internal_name << analysis_name << "Monte Carlo attempt: " << MC_attempt << endl;
            cout << internal_name << analysis_name << "Accepted: # " << accepted << " | " << 100 * (accepted / MC_attempt) << "%" << endl;
            cout << internal_name << analysis_name << "Rejected: # " << rejected << " | " << 100 * (rejected / MC_attempt) << "%" << endl;
            cout << internal_name << analysis_name << "Cleng_rejected: # " << cleng_rejected << endl;
        }
#ifdef CLENG_EXPERIMENTAL
        vtk = prepare_vtk();
        cleng_writer.write("/VTK_data", "vtk" + to_string(MC_attempt+MCS_checkpoint), dims_vtk, vtk);

        MC_free_energy = {static_cast<Real>(MC_attempt+MCS_checkpoint), free_energy_current, free_energy_current-n_times_mu};
        cleng_writer.append("/Free_energy", "free_energy", dims_3, MC_free_energy);
#endif
        if (((MC_attempt + MCS_checkpoint) % delta_save) == 0) WriteOutput();
        if (checkpoint_save) checkpoint.updateCheckpoint(simpleNodeList);
        if (cleng_dis) WriteClampedNodeDistance();
        if (cleng_pos) WriteClampedNodePosition();
        if (cleng_flag_termination) break;

        cout << endl;
    }

    cout << endl;
    cout << "Finally:" << endl;
    cout << internal_name << analysis_name << "Monte Carlo attempts: " << MC_attempt-1 << endl;
    cout << internal_name << analysis_name << "Accepted: # " << accepted << " | " << 100 * (accepted / (MC_attempt-1)) << "%" << endl;
    cout << internal_name << analysis_name << "Rejected: # " << rejected << " | " << 100 * (rejected / (MC_attempt-1)) << "%" << endl;
    cout << internal_name << analysis_name << "Cleng_rejected: # " << cleng_rejected << endl;

    end_simulation= std::chrono::steady_clock::now();
    std::cout << "# It took (s) = " << std::chrono::duration_cast<std::chrono::seconds> (end_simulation - begin_simulation).count() <<
    " or (m) = " << std::chrono::duration_cast<std::chrono::minutes> (end_simulation - begin_simulation).count() <<
    " or (h) = " << std::chrono::duration_cast<std::chrono::hours> (end_simulation - begin_simulation).count() << std::endl;

    cout << "Have a fun. " << endl;
    return success;
}
