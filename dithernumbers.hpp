#include "string_math1.hpp"

// dither numbers actually are some sort of a step below addition
// in 6502 assembly additions consist of bit operations, one step is the carry bit which just transfers the last operations overflow onto the next bit(s)
// dithernumbers just represent the overflow but for any sort of scalings in case of dither numbers its the range "0".."9" so to say it's the overflow of a range of 10
// not sure about the implications of that, however it's a rather big field of maths and logics involved in that or may be derived from it
// ...
// the idea just was to reconstruct a stream of digits out of any possible source string. So to say a fail safe error correction that can correct anything into the right shape again.
// if you have "1213112.3432" and want to get a "5555555555" as a resulting dither pattern out of that then just apply a dither factor to the stream of "1213112.3432" and it should result in a dither number of "5555555555" as overflow pattern of the dithersequence of "1213112.3432"
// have to check now whats all the code does..
// actually I wanted to have a set of functions to construct any error correction picture out of any given source number
// comment: there are no irational numbers in dither space :)
// ...
// this is not at all final just as a starting point for others to think about that
// actually our sun system somehow also can be seen as some sort of error correction on a rather broad sense (you know I am with angles and so on, so sorry for that type of view)
// i think the universe is f..... . and someone tries to correct that. (sorry... That is nothing from angels just my view on that. I hope my soul was not responsible for that f...up, at least I just have seen that this time (in the sun system) and not being producing it.)
// ...

// ...
// some test cases to receive a appropiate sequence of bits for like 500 (or more) digits of "1" like "11111111"("1" 500 times) by dithering almost any number
// ...
// no final values just nearings to good numbers for these cases
// ...
// 500 digits
// 1 = 0.95672855568135196
// 2 = 1.0913546777421734
// 3 = 1.0528880595330141
// 4 = 1.0576995003539853
// 5 = 1.9278915557977128
// 6 = 4.0528933678267345
// 7 = 4.1442357267503036
// 8 = 4.2355936087987844
// 9 = 4.3267487105901559
// ...
// 1500 digits
// 1 = 0.95673040869297221
// 2 = 1.0913486688439284
// 3 = 1.0528866809025061
// 4 = 1.0576933106362807
// 5 = 1.9278823867782122
// 6 = 4.0528863849627479
// 7 = 4.1442332077704469
// 8 = 4.2355790494594112
// 9 = 4.3267487107561902
// ...
// 5000 digits
// 1 = 0.95673162883552898
// 2 = 1.0913467011377385
// 3 = 1.0528849576131016
// 4 = 1.0576929667625439
// 5 = 1.9278841697173710
// 6 = 4.0528854297591872
// 7 = 4.1442313185350512
// 8 = 4.2355775995653229
// 9 = 4.3267480306792905
// ...

#define DITHERNUMBERS1TO9 \
const double r1 = calculateConstantForDesiredDNumberSingleDigitPeriod('1', "10", precision, 0, meanSquares); \
const double r2 = calculateConstantForDesiredDNumberSingleDigitPeriod('2', "10", precision, 0, meanSquares); \
const double r3 = calculateConstantForDesiredDNumberSingleDigitPeriod('3', "10", precision, 0, meanSquares); \
const double r4 = calculateConstantForDesiredDNumberSingleDigitPeriod('4', "10", precision, 0, meanSquares); \
const double r5 = calculateConstantForDesiredDNumberSingleDigitPeriod('5', "10", precision, 0, meanSquares); \
const double r6 = calculateConstantForDesiredDNumberSingleDigitPeriod('6', "10", precision, 0, meanSquares); \
const double r7 = calculateConstantForDesiredDNumberSingleDigitPeriod('7', "10", precision, 0, meanSquares); \
const double r8 = calculateConstantForDesiredDNumberSingleDigitPeriod('8', "10", precision, 0, meanSquares); \
const double r9 = calculateConstantForDesiredDNumberSingleDigitPeriod('9', "10", precision, 0, meanSquares);

// this generates a dithered version of the  number for the given number "v" and the ditherfactor "ditherFactor"
DNumber ditherNumberFromRight(Number v, double ditherFactor = 1.0, int clipRestPeriodAt = 1000, bool flippedDitherFactor = false);
double calculateConstantForDesiredDNumberSingleDigitPeriod(char desiredDigit, const std::string &base = "10", int maxNumberOfDigitsToCountIn = 10000, int EPSILON = 1, bool useMeanSquares = true, int antiConvergenceBreakThreshold = 10);
double calculateConstantForDesiredDNumberInitialStringPeriod(const std::string &desiredInitialPeriodString, const std::string& base = "10", int maxNumberOfDigitsToCountIn = -1, int EPSILON = 1);
double calculateConstantForDesiredDNumberByPeriodSizes(const std::string& desiredInitialPeriodString, const std::string& base = "10", int maxNumberOfDigitsToCountIn = 100000, int EPSILON = 1);
