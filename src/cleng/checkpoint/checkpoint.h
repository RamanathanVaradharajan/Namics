#pragma once

#include <string>
#include <vector>
#include <memory>
#include "../nodes/simple_node.h"
#include "../nodes/monolit.h"

#include <sys/types.h>
#include <sys/stat.h>
#include <algorithm>

#include <fstream>

using namespace std;

class Checkpoint {
public:
    Checkpoint();

    void updateCheckpoint(const vector<std::shared_ptr<SimpleNode>> &simpleNodeList);

    void saveCheckpoint(const vector<std::shared_ptr<SimpleNode>> &simpleNodeList);

    vector<shared_ptr<SimpleNode>> loadCheckpoint(vector<shared_ptr<SimpleNode>> simpleNodeList, const Point &box);

    bool isCheckpointExists(const std::string &name) const;

    void updateCheckpointName(bool plusOne = false);

    bool isLoadable();

    static const string IN_CLASS_NAME;
    static const string DEFAULT_CHECKPOINT_PATH;
    static const string DEFAULT_CHECKPOINT_NAME;


private:

    string checkpoint_name;
    string checkpoint_path;
};

