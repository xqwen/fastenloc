#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <zlib.h>

std::vector<std::string> readLines(const std::string& filename) {
    std::vector<std::string> lines;

    // Check if file is gzip compressed
    std::ifstream file(filename, std::ios::binary);
    if (!file) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return lines;
    }

    unsigned char magic[2];
    file.read(reinterpret_cast<char*>(magic), 2);
    file.close();

    if (magic[0] == 0x1F && magic[1] == 0x8B) { // gzip magic numbers
        gzFile gz = gzopen(filename.c_str(), "rb");
        if (!gz) {
            std::cerr << "Error opening gzip file: " << filename << std::endl;
            return lines;
        }

        char buffer[4096];
        while (true) {
            if (gzgets(gz, buffer, sizeof(buffer)) == NULL) break;
            lines.push_back(std::string(buffer));
        }
        gzclose(gz);
    } else {
        std::ifstream txtFile(filename);
        if (!txtFile) {
            std::cerr << "Error opening text file: " << filename << std::endl;
            return lines;
        }

        std::string line;
        while (std::getline(txtFile, line)) {
            lines.push_back(line);
        }
        txtFile.close();
    }

    return lines;
}

int main(int argc, char **argv ) {
    std::string filename = argv[1]; // Replace with your filename
    std::vector<std::string> lines = readLines(filename);
    for (const auto& line : lines) {
        std::cout << line << std::endl;
    }
    return 0;
}

