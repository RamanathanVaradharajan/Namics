#include "output_info.h"
#include <string>

using std::string;

const std::string OutputInfo::IN_CLASS_NAME = "out_info";
const std::string OutputInfo::DEFAULT_OUTPUT_PATH = "./output/";

OutputInfo::OutputInfo() {
	output_path = DEFAULT_OUTPUT_PATH;
}

void OutputInfo::addProperty(const std::string &elem, const std::string &param, const std::string &value) {
	if (elem == "folder") {
		addFolderInfo(param, value);
	}
}

void OutputInfo::addFolderInfo(const std::string &param, const std::string &value) {
	if (param == "path") {
		output_path = processOutputPath(value);
	}
}

const string &OutputInfo::getOutputPath() const {
	return output_path;
}

string OutputInfo::processOutputPath(string path) {
	if (path.back() != '/') {
		path.push_back('/');
	}
	return path;
}
