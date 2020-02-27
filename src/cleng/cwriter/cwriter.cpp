#ifdef CLENG_EXPERIMENTAL  // experimental

#include "cwriter.h"

string CWriter::getFreeFileName(const string &filename) {
    string FreeFileName = filename;
    ifstream myFile(FreeFileName + ".h5");
    int index = 0;

    if (myFile.good()) {
        // exist
        while (true) {
            FreeFileName = filename + "_" + to_string(index);
            ifstream myFile1(FreeFileName + ".h5");
            if (myFile1.good()) index++;
            else break;
        }
    }
    return FreeFileName + ".h5";
}

CWriter::CWriter(const string &filename) {

    string FREE_FILE_NAME = getFreeFileName(filename);
    // Create a new file using the default property lists.

    file = make_shared<H5File>(H5File(FREE_FILE_NAME, H5F_ACC_TRUNC));
    Exception::dontPrint();
}

bool CWriter::append(const string &group_name, const string &dataset_name, vector<int> &dims, vector<Real> data) {
    bool success = true;

    hsize_t dims_[2] = {static_cast<hsize_t>(dims[0]), static_cast<hsize_t>(dims[1])}; // dataset dimensions

    // Modify dataset creation property to enable chunking
    plist_kal = make_shared<DSetCreatPropList>(DSetCreatPropList());
    hsize_t chunk_dims[2] = {1, 2};
    plist_kal->setChunk(2, chunk_dims);
    plist_kal->setDeflate(6);

    dataspace_kal = make_shared<DataSpace>(DataSpace(RANK, dims_, maxdims));
    try {
        Group group = file->openGroup(group_name);
    }
    catch (...) {
        Group group = file->createGroup(group_name);
        cerr << "Created: " << group_name << " group" << endl;
    }
    Group group = file->openGroup(group_name);
    try {
        // Create the dataset.
        dataset_kal = make_shared<DataSet>(
                DataSet(group.createDataSet(dataset_name, PredType::IEEE_F64LE, *dataspace_kal, *plist_kal)));
        // Write the data to the dataset using default memory space, file space, and transfer properties.
        dataset_kal->write(&data[0], PredType::NATIVE_DOUBLE);
    }
    catch (...) {

        size[0] = size[0] + dimsext[0];
        size[1] = size[1];
        dataset_kal->extend(size);

        // Select a hyperslab in extended portion of the dataset.
        filespace = make_shared<DataSpace>(DataSpace(dataset_kal->getSpace()));
        filespace->selectHyperslab(H5S_SELECT_SET, dimsext, offset);
        // Define memory space.
        memspace = make_shared<DataSpace>(DataSpace(RANK, dimsext, nullptr));

        // Write data to the extended portion of the dataset.
        dataset_kal->write(&data[0], PredType::NATIVE_DOUBLE, *memspace, *filespace);

        offset[0] = size[0];
    }

    return success;
}

bool CWriter::write(const string &group_name, const string &dataset_name, vector<int> &dims, vector<Real> data) {
    bool success = true;
    hsize_t dims_[2] = {static_cast<hsize_t>(dims[0]), static_cast<hsize_t>(dims[1])}; // dataset dimensions
    hsize_t chunk_dims[2] = {static_cast<hsize_t>((int) data.size() / 2), 1};

    // Modify dataset creation property to enable chunking
    plist_vtk = make_shared<DSetCreatPropList>(DSetCreatPropList());
    plist_vtk->setChunk(2, chunk_dims);
    plist_vtk->setDeflate(6);

    dataspace_vtk = make_shared<DataSpace>(DataSpace(RANK, dims_, maxdims));
    try { Group group = file->openGroup(group_name); }
    catch (...) {
        Group group = file->createGroup(group_name);
        cerr << "Created: " << group_name << " group\n ";
    }
    Group group = file->openGroup(group_name); // for debugging

    // Create the dataset.
    dataset_vtk = make_shared<DataSet>(
            DataSet(group.createDataSet(dataset_name, PredType::IEEE_F64LE, *dataspace_vtk, *plist_vtk)));
    // Write the data to the dataset using default memory space, file space, and transfer properties.
    dataset_vtk->write(&data[0], PredType::NATIVE_DOUBLE);

    return success;
}

#endif
