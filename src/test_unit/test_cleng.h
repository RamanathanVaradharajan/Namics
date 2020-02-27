#include "catch.hpp"
#include "../namics.h"

/**
 * Test 1
 */

TEST_CASE("Free energy after 10 steps", "[chain10steps], [cleng], [working], [short]") {

// create a configuration of Namics
    NamicsConfig config;
    string filename = "chain_after10step.in";

// providing filename to testCaseCleng
    bool success = config.testCaseCleng(filename);
    if (!success) exit(0);

// some useful logic
    Real threshold = 1e-5; // set some precision
    Real f = -38.3316149995;   // expected number
    Real Free_energy = config.Cle[0]->free_energy_trial;  // asking free_energy value from the engine
    REQUIRE( abs(Free_energy - f)  < threshold );
}

/**
 * Test 2
 */

TEST_CASE("Free energy vector after 20 steps", "[chain20steps_vector], [cleng], [working], [short]") {

// create a configuration of Namics
    NamicsConfig config;
    string filename = "chain_after20step.in";

    const bool save_vector = true;  // if I would like to save something in my engine (optional)
    bool success = config.testCaseCleng(filename, save_vector);
    if (!success) exit(0);

// useful logic
    Real threshold = 1e-5; // set some precision
    vector<Real> e;
    int lenght = 21;
    for (int i = 0; i < lenght; i++) { e.push_back(-38.3316149995); } // expected vector
    vector<Real> v = config.Cle[0]->test_vector; // if save_vector was provided I can ask variable (optional)
    for (std::size_t i = 0; i < v.size(); ++i) REQUIRE( abs(e[i] - v[i]) < threshold );
}

/**
 * Test 3
 */

TEST_CASE("Free energy of chain along 1 axis", "[chain_along1_axis], [cleng], [working], [long], [chain_along_axises]") {

// create a configuration of Namics
    NamicsConfig config;
    string filename = "chain_go_along_axis1.in";

    const bool save_vector = true;  // if I would like to save something in my engine (optional)
    bool success = config.testCaseCleng(filename, save_vector);
    if (!success) exit(0);

// useful logic
    Real threshold = 1e-5;          // set some precision
    vector<Real> e;
    int lenght = 201;
    for (int i = 0; i < lenght; i++) { e.push_back(-38.3316149995); } // expected vector
    vector<Real> v = config.Cle[0]->test_vector; // if save_vector was provided I can ask variable (optional)
    for (size_t i = 0; i < v.size(); ++i) { REQUIRE( abs(e[i] - v[i]) < threshold);}
}

/**
 * Test 4
 */

TEST_CASE("Free energy of chain along 2 axis", "[chain_along2_axis], [cleng], [working], [long], [chain_along_axises]") {

// create a configuration of Namics
NamicsConfig config;
string filename = "chain_go_along_axis2.in";

const bool save_vector = true;  // if I would like to save something in my engine (optional)
bool success = config.testCaseCleng(filename, save_vector);
if (!success) exit(0);

// useful logic
Real threshold = 1e-5;          // set some precision
vector<Real> e;
int lenght = 201;
for (int i = 0; i < lenght; i++) { e.push_back(-38.3316149995); } // expected vector
vector<Real> v = config.Cle[0]->test_vector; // if save_vector was provided I can ask variable (optional)
for (size_t i = 0; i < v.size(); ++i) { REQUIRE( abs(e[i] - v[i]) < threshold);}
}

/**
 * Test 5
 */

TEST_CASE("Free energy of chain along 3 axis", "[chain_along3_axis], [cleng], [working], [long], [chain_along_axises]") {

// create a configuration of Namics
NamicsConfig config;
string filename = "chain_go_along_axis3.in";

const bool save_vector = true;  // if I would like to save something in my engine (optional)
bool success = config.testCaseCleng(filename, save_vector);
if (!success) exit(0);

// useful logic
Real threshold = 1e-5;          // set some precision
vector<Real> e;
int lenght = 201;
for (int i = 0; i < lenght; i++) { e.push_back(-38.3316149995); } // expected vector
vector<Real> v = config.Cle[0]->test_vector; // if save_vector was provided I can ask variable (optional)
for (size_t i = 0; i < v.size(); ++i) { REQUIRE( abs(e[i] - v[i]) < threshold);}
}

/**
 * Test 6
 */

TEST_CASE("Free energy of 3 segments along 1 axis", "[3segments_along1_axis], [cleng], [working], [long], [3segments_along_axises]") {

// create a configuration of Namics
NamicsConfig config;
string filename = "3segments_go_along_axis1.in";

const bool save_vector = true;  // if I would like to save something in my engine (optional)
bool success = config.testCaseCleng(filename, save_vector);
if (!success) exit(0);

// useful logic
Real threshold = 1e-5;          // set some precision
vector<Real> e;
int lenght = 201;
for (int i = 0; i < lenght; i++) { e.push_back(-40.8498); } // expected vector
vector<Real> v = config.Cle[0]->test_vector; // if save_vector was provided I can ask variable (optional)
for (size_t i = 0; i < v.size(); ++i) { REQUIRE( abs(e[i] - v[i]) < threshold);}
}

/**
 * Test 6
 */

TEST_CASE("Free energy of 3 segments along 2 axis", "[3segments_along2_axis], [cleng], [working], [long], [3segments_along_axises]") {

// create a configuration of Namics
NamicsConfig config;
string filename = "3segments_go_along_axis2.in";

const bool save_vector = true;  // if I would like to save something in my engine (optional)
bool success = config.testCaseCleng(filename, save_vector);
if (!success) exit(0);

// useful logic
Real threshold = 1e-5;          // set some precision
vector<Real> e;
int lenght = 201;
for (int i = 0; i < lenght; i++) { e.push_back(-40.8498); } // expected vector
vector<Real> v = config.Cle[0]->test_vector; // if save_vector was provided I can ask variable (optional)
for (size_t i = 0; i < v.size(); ++i) { REQUIRE( abs(e[i] - v[i]) < threshold);}
}

/**
 * Test 7
 */

TEST_CASE("Free energy of 3 segments along 3 axis", "[3segments_along3_axis], [cleng], [working], [long], [3segments_along_axises]") {

// create a configuration of Namics
NamicsConfig config;
string filename = "3segments_go_along_axis3.in";

const bool save_vector = true;  // if I would like to save something in my engine (optional)
bool success = config.testCaseCleng(filename, save_vector);
if (!success) exit(0);

// useful logic
Real threshold = 1e-5;          // set some precision
vector<Real> e;
int lenght = 201;
for (int i = 0; i < lenght; i++) { e.push_back(-40.8498); } // expected vector
vector<Real> v = config.Cle[0]->test_vector; // if save_vector was provided I can ask variable (optional)
for (size_t i = 0; i < v.size(); ++i) { REQUIRE( abs(e[i] - v[i]) < threshold);}
}


/**
 * Test 8
 */

TEST_CASE("Free energy of chain extension ", "[chain_extension], [cleng], [working], [short]") {

// create a configuration of Namics
NamicsConfig config;
string filename = "chain_extension.in";

const bool save_vector = true;  // if I would like to save something in my engine (optional)
bool success = config.testCaseCleng(filename, save_vector);
if (!success) exit(0);
// useful logic
Real threshold = 1e-4; // set some precision
vector<Real> e = {-38.6145, -38.7202, -38.6307, -38.2885, -37.6829, -36.8088, -35.6629, -34.2294};
vector<Real> v = config.Cle[0]->test_vector; // if save_vector was provided I can ask variable (optional)
for (size_t i = 0; i < v.size(); ++i) { REQUIRE( abs(e[i] - v[i]) < threshold );}
}

/**
 * Test 9
 */

TEST_CASE("Free energy of chain compression ", "[chain_compression], [cleng], [working], [short]") {

// create a configuration of Namics
NamicsConfig config;
string filename = "chain_compression.in";

const bool save_vector = true;  // if I would like to save something in my engine (optional)
bool success = config.testCaseCleng(filename, save_vector);
if (!success) exit(0);
// useful logic
Real threshold = 1e-4; // set some precision
vector<Real> e = {-34.2294, -35.6629, -36.8088, -37.6829, -38.2886, -38.6307, -38.7202, -38.6145};
vector<Real> v = config.Cle[0]->test_vector;
for (size_t i = v.size()-1; i > 0; --i) { REQUIRE( abs(e[i] - v[i]) < threshold);}
}