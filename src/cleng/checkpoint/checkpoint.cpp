#include "checkpoint.h"
#include "../iterator/EnumerateIterator.h"

using namespace std;

const string Checkpoint::IN_CLASS_NAME = ".checkpoint";
const string Checkpoint::DEFAULT_CHECKPOINT_PATH = "./checkpoints/";
const string Checkpoint::DEFAULT_CHECKPOINT_NAME = "0";

Checkpoint::Checkpoint() {
    checkpoint_path = DEFAULT_CHECKPOINT_PATH;
    checkpoint_name = DEFAULT_CHECKPOINT_NAME;
}

bool Checkpoint::isLoadable() {
    bool success;
    updateCheckpointName();
    success = isCheckpointExists(checkpoint_path + checkpoint_name + IN_CLASS_NAME);
    return  success;
}

bool Checkpoint::isCheckpointExists(const std::string &name) const {
    struct stat info{};
    return (stat(name.c_str(), &info) == 0);
}

void Checkpoint::updateCheckpointName(const bool plusOne) {
    int index = 0;
    while (isCheckpointExists(checkpoint_path + to_string(index) + IN_CLASS_NAME)) {index++;}
    if (!plusOne and index > 0) index--;
    checkpoint_name = to_string(index);
}

void Checkpoint::saveCheckpoint(const vector<std::shared_ptr<SimpleNode>>& simpleNodeList) {
    ifstream checkpointDir(checkpoint_path);
    if (checkpointDir.fail()) { mkdir(checkpoint_path.c_str(), 0777); }
    updateCheckpointName(true);
    updateCheckpoint(simpleNodeList);
}

void Checkpoint::updateCheckpoint(const vector<std::shared_ptr<SimpleNode>>& simpleNodeList) {
    string filename;
    ofstream outfile;
    outfile.open(checkpoint_path + checkpoint_name + IN_CLASS_NAME);
    for (auto &&n : Enumerate(simpleNodeList)) {
        if (n.first % 2 == 1) continue;
        auto p  = n.second->to_string(); auto pc = n.second->get_cnode()->to_string();
        outfile << p << pc << endl;}
}

shared_ptr<SimpleNode> fromFileToNode(int x, int y, int z, int id, const Point &box) {
    return make_shared<SimpleNode>(Point(x, y, z), id, box);
}

vector<shared_ptr<SimpleNode>> Checkpoint::loadCheckpoint(vector<shared_ptr<SimpleNode>> simpleNodeList, const Point& box) {
    if (isCheckpointExists(checkpoint_path + checkpoint_name + IN_CLASS_NAME)) {
        cout << "Loading checkpoint ..." << endl;
        updateCheckpointName();
        cout << "checkpoint: " << checkpoint_name + IN_CLASS_NAME << endl;
        ifstream infile(checkpoint_path + checkpoint_name + IN_CLASS_NAME);

        simpleNodeList.clear();
        string line;
        cout << "Reading file..." << endl;
        while (getline(infile, line)) {
            istringstream iss(line);
            int id1, x1, y1, z1, id2, x2, y2, z2;
            string a;
            if (!(iss >> a >> id1 >> a >> x1 >> a >> y1 >> a >> z1 >> a >> id2 >> a >> x2 >> a >> y2 >> a >> z2)) {
                cout << "Error! The checkpoint has error inside! " << endl; break;
            }

            auto first_node  = fromFileToNode(x1, y1, z1, id1, box);
            auto second_node = fromFileToNode(x2, y2, z2, id2, box);
            first_node ->set_cnode(second_node);
            second_node->set_cnode(first_node);
            //
            simpleNodeList.push_back(first_node);
            simpleNodeList.push_back(second_node);
        }
        //
    } else cout << "Sorry! Unable to open file!" << endl;
    return simpleNodeList;
}
