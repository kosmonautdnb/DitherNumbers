// ......
#include <string>
// ......
// some string based math "lib", nothing to be proud of..
// needed because my compiler was massively broken even the simplest math functions failed
// ......
typedef std::string Number; // the actual numbers just as a string
typedef std::pair<std::string,std::string> DNumber; // dither numbers consisting of a base and fractional
// ......
int dotPos(const Number& n); // returns the position of the decimal point
Number cleanUp(const Number& n); // cleans front zeroes like "000010" gets "10" and "00000.20" gets "0.20"
Number cleanUpAll(const Number& n); // cleans front + back like "0000010.01003200" gets "10.010032"
Number fixAll(const Number& n); // cleans front + back + some rather odd cases like "0000010." gets "10.0" and "0000112" gets "112.0"
Number concatNumberTo(const Number& n, int vorkomma, int nachkomma); // expands digits before point and expands digits behind point for (4,5) and "1.2" this would give "0001.20000"
Number toNumber(double n); // this would make a string number out of the double value (rather imprecise so to say)
double fromNumber(Number n); // this would make a double number out of a string number (rather imprecise so to say)
double cleanDouble(double n); // this converts the double to string and back to double to get a somehow cleaned double/float representation (not sure if that makes much sense)
Number makeWhole(const Number& n); // actually this does "124234" to "000000124234" not sure why have to check again (same numbers of zeroes prepended like the actual digits of the source number)
int makeSameBase(Number& a, Number& b); // this makes both numbers to the actual same string length so to say "1.2" and "100.20" gets "001.20"
Number add(const Number& n1, const Number& n2); // addition of two numbers
Number multiplyWithTen(const Number& number); // multiplication with 10 (actually just "change of digit length")
Number divideByTen(const Number& number); // division by 10 (actually just "change of digit length")
Number multiplyWithTwo(const Number& number); // multiplication with 2 (additions)
Number multiplyWithThree(const Number& number); // multiplication with 3 (additions)
Number multiplyWithFour(const Number& number); // multiplication with 4 (additions)
Number multiplyWithFive(const Number& number); // multiplication with 5 (additions)
Number multiplyWithSix(const Number& number); // multiplication with 6 (additions)
Number multiplyWithSeven(const Number& number); // multiplication with 7 (additions)
Number multiplyWithEight(const Number& number); // multiplication with 8 (additions)
Number multiplyWithNine(const Number& number); // multiplication with 9 (additions)
Number divideByTwo(const Number& number); // divide by two somehow derived from 6502 divisions
Number divideByLessThan2Test(const Number& number, double divider = 2.0); // maybe not fully working
Number divideByFour(const Number& number); // divide by four through two divide by twos
Number multiplyBy2Point5(const Number& number); // multiply with 2.5 (addition and addition of halve)
Number divideBy2Point5(const Number& number); // * 4 / 10
Number divideBy1Point5(const Number& number); // * 3 / 10
int countSameDigits(Number& a, int digit = -1); // counts all digits which are like "digit" e.g "33334.311" would give 5 if looking for "3" and 1 if looking for "4" and 2 if looking for "1"
// ......
void printNumber2(const Number& n);
void printNumber3(const Number& n);
void printDNumber(const DNumber& n);
// ......
Number multiplyWith0to9(const Number &n, int number);
Number multiply(Number a, Number b);
// ......