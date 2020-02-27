#include "../../namics.h"
#include <memory>
#include <H5Cpp.h>

using namespace H5;
using namespace std;

class CWriter {
public:
    int RANK = 2;
    shared_ptr<H5File> file;

    shared_ptr<DataSpace> dataspace_vtk;
    shared_ptr<DataSpace> dataspace_kal;
    shared_ptr<DataSet> dataset_vtk;
    shared_ptr<DataSet> dataset_kal;
    shared_ptr<DSetCreatPropList> plist_vtk;
    shared_ptr<DSetCreatPropList> plist_kal;

    hsize_t maxdims[2] = {H5S_UNLIMITED, H5S_UNLIMITED};

    // Variables used in extending and writing to the extended portion of dataset
    hsize_t size[2] = {1, 3};
    hsize_t offset[2] = {1, 0};
    hsize_t dimsext[2] = {1, 3}; // extend dimensions
    shared_ptr<DataSpace> filespace;
    shared_ptr<DataSpace> memspace;

    explicit CWriter(const string &filename);

    bool write(const string &group_name, const string &dataset_name, vector<int> &dims, vector<Real> data);

    bool append(const string &group_name, const string &dataset_name, vector<int> &dims, vector<Real> data);

private:
    string getFreeFileName(const string &filename);

};
