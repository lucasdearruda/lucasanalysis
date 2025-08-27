#include <iostream>
#include <fstream>
#include <sstream>
#include <map>

int main() {
    // Define a map for storing energy values
    std::map<char, std::map<int, float>> energy_map;

    // Open the text file
    std::ifstream file("your_data.txt");
    if (!file.is_open()) {
        std::cerr << "Failed to open the file." << std::endl;
        return 1;
    }

    // Read line by line
    std::string line;
    while (std::getline(file, line)) {
        if (line.empty() || line[0] == '#') {
            continue; // Skip empty lines and comments
        }

        // Parse the line
        char particle;
        int telescope;
        float energy;
        std::istringstream iss(line);
        iss >> particle >> telescope >> energy;

        // Store the energy value in the map
        energy_map[particle][telescope] = energy;
    }

    // Close the file
    file.close();

    // Example: Accessing energy values
    char particle = 'p';  // Example particle 'p'
    int telescope = 1;    // Example telescope 1

    // Retrieve the energy value
    float energy_value = energy_map[particle][telescope];

    // Print the retrieved energy value
    std::cout << "Energy value for particle " << particle << " at telescope " << telescope << " is " << energy_value << std::endl;

    return 0;
}

