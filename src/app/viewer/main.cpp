#include <iostream>
#include <string>
#include "math/math.hpp"
#include "util/strings.hpp"

int main(int argc, char** argv) {
    std::string mode = (argc > 1) ? argv[1] : "demo";
    std::cout << "[viewer] mode=" << mode << "\n";
    std::cout << "2+2=" << mathx::add(2,2) << ", 3*5=" << mathx::mul(3,5) << "\n";
    std::cout << "title=" << util::to_lower(std::string("APHELION 3D VIEWER")) << "\n";
    return 0;
}
