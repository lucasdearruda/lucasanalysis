#include <vector>
#include <string>
#include <sstream>
#include <iostream>
#include <TNamed.h>

std::string makeRunsString(const std::vector<int>& runs) {
    std::stringstream ss;
    for (size_t i = 0; i < runs.size(); i += 2) {
        if (i + 1 < runs.size() && runs[i] != runs[i+1]) {
            // caso seja um intervalo
            ss << runs[i] << "-" << runs[i+1];
        } else {
            // caso seja um único run
            ss << runs[i];
        }
        if (i + 2 < runs.size()) ss << ", ";
    }
    return ss.str();
}

int test() {
    std::vector<int> runsForward = {397,405,408,408,410,410,425,438};
    std::vector<int> runsBackward = {410,412};

    std::string forwardStr = makeRunsString(runsForward);
    std::string backwardStr = makeRunsString(runsBackward);

    TNamed *runs_infoF = new TNamed("Forward runs:", forwardStr.c_str());
    TNamed *runs_infoB = new TNamed("Backward runs:", backwardStr.c_str());

    std::cout << "Forward: " << forwardStr << std::endl;
    std::cout << "Backward: " << backwardStr << std::endl;

    return 0;
}