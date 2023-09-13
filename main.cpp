#include "dithernumbers.hpp"

int main() {
    const double r3 = calculateConstantForDesiredDNumberSingleDigitPeriod('3', "10", 500, 0, true);
    DNumber dr3 = ditherNumberFromRight("123423423423",r3);
    printDNumber(dr3);
    return 0;
}