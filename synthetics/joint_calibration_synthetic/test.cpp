#include <ql/quantlib.hpp>
#include <iostream>

int main() {
#ifdef QL_ENABLE_SESSIONS
    std::cout << "Sessions enabled\n";
#else
    std::cout << "Sessions NOT enabled\n";
#endif
}