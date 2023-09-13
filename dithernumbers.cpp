#include "dithernumbers.hpp"
#include <set>
#include <map>
#include <functional>
#include <math.h>

std::map<Number, int> antiConvergence;

// Newton like error correction function
double errorReduce(double rangeLow, double rangeHi, double threshold, std::function<double(double v, double threshold)> checky) {
    double rv = 0;
    int steps = 10000000;
    double lastError = 10000000000000000.0;
    double lowerBound = rangeLow;
    double upperBound = rangeHi;
    double calcConstant = 0.0;
    while (lastError > 0.0001) { // kilometers
        double k = (lowerBound + upperBound) / 2.0;
        double error1 = checky((k + upperBound) / 2.0, threshold);
        double error2 = checky((k + lowerBound) / 2.0, threshold);
        error1 *= error1; // least squares
        error2 *= error2; // least squares
        if (error1 > error2) {
            upperBound = k;
            lastError = error1;
        }
        else {
            lowerBound = k;
            lastError = error2;
        }
        double lastErrorB = lastError;
        rv = k;
    }
    return rv;
}

DNumber ditherNumberFromRight(Number v, double ditherFactor, int clipRestPeriodAt, bool flippedDitherFactor) {
    Number errorValue = "0.0";
    Number errorValueHigh = "0.0";
    double EPSILON = 0.0001;
    Number result;
    Number resultRestPeriod;
    std::set<Number> alreadyEncounteredRests;
    for (int i = v.size() - 1; i == i; --i) {
        unsigned char c = i >= 0 ? v[i] : 0;
        if (c != '.') {
            c -= '0';
            Number ditherVal = add(add(toNumber(c * ditherFactor), errorValue), errorValueHigh);
            if (flippedDitherFactor) { // dunno if that makes sense
                ditherVal = add(add(toNumber(c / ditherFactor), errorValue), errorValueHigh);
            }
            int snapped = (int)trunc(fromNumber(ditherVal));
            errorValue = toNumber(fromNumber(ditherVal) - snapped);
            errorValueHigh = toNumber(snapped/10);
            snapped = snapped % 10;
            if (i >= 0) {
                result = std::to_string(snapped) + result;
            }
            else {
                if (alreadyEncounteredRests.find(errorValue) != alreadyEncounteredRests.end())
                    break;
                alreadyEncounteredRests.insert(errorValue);
                resultRestPeriod = std::to_string(snapped) + resultRestPeriod;
                if (alreadyEncounteredRests.size() > clipRestPeriodAt) {
                    resultRestPeriod = "*" + resultRestPeriod; // period too long
                    break;
                }
            }
        }
    }
    Number before = fixAll(v);
    Number after = fixAll(result);
    Number resultRest = fixAll(resultRestPeriod);
    return std::make_pair(result, resultRest);
}


double calculateConstantForDesiredDNumberSingleDigitPeriod(char desiredDigit, const std::string &base, int maxNumberOfDigitsToCountIn, int EPSILON, bool useMeanSquares, int antiConvergenceBreakThreshold) {
    antiConvergence.clear();
    int iterations = 0;
    double checkkky = errorReduce(0.1, 9.0, 0.0, [&](double v, double thresh) -> double {
        iterations++;
        DNumber d = ditherNumberFromRight(base, pow(v, 1.0 / 1.5), maxNumberOfDigitsToCountIn);
        if (useMeanSquares) {
            int i = 0;
            double error = 0;
            if (i < d.second.size())
                if (d.second[i] == '*')
                    i++;
            for (; i < d.second.size() && i < maxNumberOfDigitsToCountIn; ++i) {
                int delta = 0;
                delta = (int)d.second[i] - (int)desiredDigit;
                error += delta * delta;
            }
            antiConvergence[d.second]++;
            if (antiConvergence[d.second] > antiConvergenceBreakThreshold)
                return 0;
            error = sqrt(error);
            return error;
        }
        else {
            int count = 0;
            int i = 0;
            if (i < d.second.size())
                if (d.second[i] == '*')
                    i++;
            int lastNotMet = 0;
            //int k = 0;
            for (; i < d.second.size(); ++i) {
                if (d.second[i] == desiredDigit) {
                    count++;
                }
                else {
                    lastNotMet = i;
                }
                //k++;
                //if (k >= maxNumberOfDigitsToCountIn)
                //    break;
            }
            int notAlright = maxNumberOfDigitsToCountIn - count;
            // the numbers don't follow the same pattern as rational numbers
            // it's not possible to move the error more to the back, because this would destroy the long period
            //if (notAlright == 1) {
            //    return (double)(d.second.size()-1- lastNotMet) / (d.second.size()+1);
            //}
            antiConvergence[d.second]++;
            if (antiConvergence[d.second] > antiConvergenceBreakThreshold)
                return 0;
            return notAlright;
        }
    });
    double n = pow(checkkky, 1.0 / 1.5);
    return n;
}

double calculateConstantForDesiredDNumberInitialStringPeriod(const std::string &desiredInitialPeriodString, const std::string& base, int maxNumberOfDigitsToCountIn, int EPSILON) {
    if (maxNumberOfDigitsToCountIn < 0)
        maxNumberOfDigitsToCountIn = desiredInitialPeriodString.size();
    int iterations = 0;
    double checkkky = errorReduce(0.1, 9.0, 0.0, [&](double v, double thresh) -> double {
        iterations++;
        DNumber d = ditherNumberFromRight(base, pow(v, 1.0 / 1.5), maxNumberOfDigitsToCountIn);
        double error = 0;
        for (int i = 0; i < desiredInitialPeriodString.size(); ++i) {
            int delta = 0;
            int i2 = i;
            if (i2 < d.second.size())
                if (d.second[i2] == '*')
                    i2++;
            if (i2 < d.second.size())
                delta = (int)d.second[i2] - (int)desiredInitialPeriodString[i];
            error += delta * delta;
        }
        if ((iterations % 20) == 0) {
            printf("lolo:%d\n", iterations);
            printNumber2(d.second);
        }
        error = sqrt(error);
        if (error < EPSILON)
            error = 0;
        return error;
    });
    double n = pow(checkkky, 1.0 / 1.5);
    return n;
}

// dunno whats that for..
double calculateConstantForDesiredDNumberByPeriodSizes(const std::string& desiredInitialPeriodString, const std::string& base, int maxNumberOfDigitsToCountIn, int EPSILON) {
    int iterations = 0;
    double checkkky = errorReduce(0.1, 9.0, 0.0, [&](double v, double thresh) -> double {
        iterations++;
        DNumber d = ditherNumberFromRight(base, pow(v, 1.0 / 1.5), maxNumberOfDigitsToCountIn);
        std::string sectionString = toNumber(v/1.5);
        sectionString.erase(sectionString.begin() + dotPos(sectionString));
        double value = 0;
        int i = 0;
        if (i < d.second.size())
            if (d.second[i] == '*')
                i++;
        int thisSectionCount = 0;
        int thisSectionChar = sectionString.size() > i ? sectionString[i] : -1;
        Number b;
        for (; i < sectionString.size(); ++i) {
            if (sectionString[i] == thisSectionChar) {
                thisSectionCount++;
            }
            else {
                b = multiplyWithTen(b);
                b = add(b, toNumber(thisSectionCount+1));
                thisSectionCount = 0;
                thisSectionChar = sectionString[i];
            }
        }
        b = multiplyWithTen(b);
        b = add(b, toNumber(thisSectionCount + 1));
        std::string b2 = b.substr(0, dotPos(b));
        double error = 0;
        for (int i = 0; i < desiredInitialPeriodString.size() && i < b2.size(); ++i) {
            int delta = 0;
            delta = (int)b2[b2.size()-1-i] - (int)desiredInitialPeriodString[desiredInitialPeriodString.size()-1-i];
            error += delta * delta;
        }
        error = sqrt(error);
        if ((iterations % 50) == 0) {
            printf("lolo:%d, %f, %f\n", iterations, error, v);
            printNumber2(b2);
        }
        return error;
    });
    double n = pow(checkkky, 1.0 / 1.5);
    return n;
}


void basicConstants() {
    std::string kr = "31415926535897932";
    //std::string kr = "3333333333333333333333333333333333333333333333333333333333333333333333333333333333333333";
    //std::string kr = flipString("31415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679821480865132823066470938446095505822317253594081284811174502841027019385211055596446229489549303819644288109756659334461284756482337867831652712019091456485669234603486104543266482133936072602491412737245870066063155881748815209209628292540917153643678925903600113305305488204665213841469519415116094330572703657595919530921861173819326117931051185480744623799627495673518857527248912279381830119491298336733624406566430860213949463952247371907021798609437027705392171762931767523846748184676694051320005681271452635608277857713427577896091736371787214684409012249534301465495853710507922796892589235420199561121290219608640344181598136297747713099605187072113499999983729780499510597317328160963185950244594553469083026425223082533446850352619311881710100031378387528865875332083814206171776691473035982534904287554687311595628638823537875937519577818577805321712268066130019278766111959092164201989380952572010654858632788659361533818279682303019520353018529689957736225994138912497217752834791315155748572424541506959508295331168617278558890750983817546374649393192550604009277016711390098488240128583616035637076601047101819429555961989467678374494482553797747268471040475346462080466842590694912933136770289891521047521620569660240580381501935112533824300355876402474964732639141992726042699227967823547816360093417216412199245863150302861829745557067498385054945885869269956909272107975093029553211653449872027559602364806654991198818347977535663698074265425278625518184175746728909777727938000816470600161452491921732172147723501414419735685481613611573525521334757418494684385233239073941433345477624168625189835694855620992192221842725502542568876717904946016534668049886272327917860857843838279679766814541009538837863609506800642251252051173929848960841284886269456042419652850222106611863067442786220391949450471237137869609563643719172874677646");
    //std::string kr = "135";
    const bool meanSquares = true;
    const int precision = 5000; // meanSquares ? kr.size() : 600;
    DITHERNUMBERS1TO9
    const int precisiona = 5000;
    const double r1const5000mean = 0.95673162883552897994832164840772747993469238281250;
    const double r5const5000mean = 1.92788416971737097682648709451314061880111694335938;
    const double r5r1const5000mean = r5const5000mean - r1const5000mean;

    // without mean squares there has at least to be an epsilon of 1 for most cases (i think 3 got through with an epsilon of 0 digits being other than "3" in the final dither number for any length of "3" digits)
    //const double r1 = calculateConstantForDesiredDNumberSingleDigitPeriod('1', precision, 1);
    //const double r2 = calculateConstantForDesiredDNumberSingleDigitPeriod('2', precision, 1);
    //const double r4 = calculateConstantForDesiredDNumberSingleDigitPeriod('4', precision, 1);
    //const double r5 = calculateConstantForDesiredDNumberSingleDigitPeriod('5', precision, 1);
    //const double r6 = calculateConstantForDesiredDNumberSingleDigitPeriod('6', precision, 1);
    //const double r8 = calculateConstantForDesiredDNumberSingleDigitPeriod('8', precision, 1);
    //const double r9 = calculateConstantForDesiredDNumberSingleDigitPeriod('9', precision, 1);
    const double r1_prec = calculateConstantForDesiredDNumberSingleDigitPeriod('1', "10", precisiona, 0, true);
    const double r5_prec = calculateConstantForDesiredDNumberSingleDigitPeriod('5', "10", precisiona, 0, true);
    const double point5 = fabs((r5_prec-trunc(r5_prec)) - (r1_prec-trunc(r1_prec))); // not sure maybe division by 0.5 or something (pattern would be interessting here)
    const double r57 = calculateConstantForDesiredDNumberByPeriodSizes(kr, "10", 40, 1);

    // just some test cases for some arbitrary choosen numbers to check the result values
    DNumber k_1 = ditherNumberFromRight("12345", r1);
    DNumber k_2 = ditherNumberFromRight("12345", r2);
    DNumber k_3 = ditherNumberFromRight("12345", r3);
    DNumber k_4 = ditherNumberFromRight("12345", r4);
    DNumber k_5 = ditherNumberFromRight("12345", r5);
    DNumber k_6 = ditherNumberFromRight("12345", r6);
    DNumber k_7 = ditherNumberFromRight("12345", r7);
    DNumber k_8 = ditherNumberFromRight("12345", r8);
    DNumber k_9 = ditherNumberFromRight("12345", r9);
    DNumber kb10_1 = ditherNumberFromRight("12345", r1);
    DNumber kb10_2 = ditherNumberFromRight("12345", r2);
    DNumber kb10_3 = ditherNumberFromRight("12345", r3);
    DNumber kb10_4 = ditherNumberFromRight("12345", r4);
    DNumber kb10_5 = ditherNumberFromRight("12345", r5);
    DNumber kb10_6 = ditherNumberFromRight("12345", r6);
    DNumber kb10_7 = ditherNumberFromRight("12345", r7);
    DNumber kb10_8 = ditherNumberFromRight("12345", r8);
    DNumber kb10_9 = ditherNumberFromRight("12345", r9);
    DNumber k3 = ditherNumberFromRight("1234", r3);
    DNumber k7_1 = ditherNumberFromRight("1", r7);
    DNumber k7_2 = ditherNumberFromRight("2", r7);
    DNumber k7_3 = ditherNumberFromRight("3", r7);
    DNumber k7_4 = ditherNumberFromRight("4", r7);
    DNumber k7_5 = ditherNumberFromRight("5", r7);
    DNumber k7_6 = ditherNumberFromRight("6", r7);
    DNumber k7_7 = ditherNumberFromRight("7", r7);
    DNumber k7_8 = ditherNumberFromRight("8", r7);
    DNumber k7_9 = ditherNumberFromRight("9", r7);
    DNumber k7_10 = ditherNumberFromRight("10", r7);
    DNumber k7_20 = ditherNumberFromRight("20", r7);
    DNumber k7_50 = ditherNumberFromRight("50", r7);
    DNumber k7_100 = ditherNumberFromRight("100", r7);
    DNumber k7_666 = ditherNumberFromRight("666", r7);
    DNumber k7_1234 = ditherNumberFromRight("1234", r7);
    DNumber k7_21234 = ditherNumberFromRight("21234", r7);
    DNumber k3_1 = ditherNumberFromRight("1", r3);
    DNumber k3_2 = ditherNumberFromRight("2", r3);
    DNumber k3_3 = ditherNumberFromRight("3", r3);
    DNumber k3_4 = ditherNumberFromRight("4", r3);
    DNumber k3_5 = ditherNumberFromRight("5", r3);
    DNumber k3_6 = ditherNumberFromRight("6", r3);
    DNumber k3_7 = ditherNumberFromRight("7", r3);
    DNumber k3_8 = ditherNumberFromRight("8", r3);
    DNumber k3_9 = ditherNumberFromRight("9", r3);
    DNumber k3_10 = ditherNumberFromRight("10", r3);
    DNumber k3_20 = ditherNumberFromRight("20", r3);
    DNumber k3_50 = ditherNumberFromRight("50", r3);
    DNumber k3_100 = ditherNumberFromRight("100", r3);
    DNumber k3_1234 = ditherNumberFromRight("1234", r3);
    DNumber k3_21234 = ditherNumberFromRight("21234", r3);
}

static void checkycheck() {
    // just to look for the count of iterations in this case (not 100% sure why that)
    // ...
    int iterations = 0;
    const int EPSILON = 1; // dunno why this is needed // 42 iterations here
    double checkkky = errorReduce(0.1, 9.0, 0.0, [&](double v, double thresh) -> double {
        iterations++;
        int fracPart = 1000;
        int count = 0;
        DNumber d = ditherNumberFromRight("10.0", pow(v, 1.0 / 1.5), fracPart);
        for (int i = 0; i < d.second.size(); ++i) {
            if (d.second[i] == '3')
                count++;
        }
        return fabs(fracPart - count - EPSILON);
    });
    double n = pow(checkkky, 1.0 / 1.5);
    DNumber b = ditherNumberFromRight("100.0", n);
    // ...
}
