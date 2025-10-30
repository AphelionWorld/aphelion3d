#include <iostream>
#include "math/math.hpp"
#include "util/strings.hpp"

int main() {
    using std::cout;
    using util::to_upper;


    int a = 7, b = 6;
    cout << "add: " << mathx::add(a,b) << "\n";
    cout << "mul: " << mathx::mul(a,b) << "\n";
    cout << "upper: " << to_upper(std::string("Hello MinGW + CMake + VS Code")) << "\n";
    return 0;
}
