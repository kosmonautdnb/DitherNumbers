//#define _CRT_SECURE_NO_WARNINGS
//#include <stdio.h>
//#include <string>
//#include <vector>
//#include "png/png.h"
//#include "windows.h"
//#include "string_math1.hpp"
//
//void interesstingValues() {
//    // pow(3, 1.0 / 1.5) is an interessting "dithering value"
//    // 1.0/1.5 always leads to same patterns
//    DNumber b = ditherNumberFromRight("99.0", pow(3, 1.0 / 1.5));
//}
//
//const double pi = 3.1415927;
//double standardPi2 = 3.141592653589793238;
//double keralaPi2 = 104348.0 / 33215.0;
//#include <functional>
//
//double errorReduce(double rangeLow, double rangeHi, double threshold, std::function<double(double v, double threshold)> checky);
//
//#include <map>
//std::map<Number, int> antiConvergence;
//
//double calculateConstantForDesiredDNumberSingleDigitPeriod(char desiredDigit, const std::string &base = "10", int maxNumberOfDigitsToCountIn = 10000, int EPSILON = 1, bool useMeanSquares = true, int antiConvergenceBreakThreshold = 10) {
//    antiConvergence.clear();
//    int iterations = 0;
//    double checkkky = errorReduce(0.1, 9.0, 0.0, [&](double v, double thresh) -> double {
//        iterations++;
//        DNumber d = ditherNumberFromRight(base, pow(v, 1.0 / 1.5), maxNumberOfDigitsToCountIn);
//        if (useMeanSquares) {
//            int i = 0;
//            double error = 0;
//            if (i < d.second.size())
//                if (d.second[i] == '*')
//                    i++;
//            for (; i < d.second.size() && i < maxNumberOfDigitsToCountIn; ++i) {
//                int delta = 0;
//                delta = (int)d.second[i] - (int)desiredDigit;
//                error += delta * delta;
//            }
//            antiConvergence[d.second]++;
//            if (antiConvergence[d.second] > antiConvergenceBreakThreshold)
//                return 0;
//            error = sqrt(error);
//            return error;
//        }
//        else {
//            int count = 0;
//            int i = 0;
//            if (i < d.second.size())
//                if (d.second[i] == '*')
//                    i++;
//            int lastNotMet = 0;
//            //int k = 0;
//            for (; i < d.second.size(); ++i) {
//                if (d.second[i] == desiredDigit) {
//                    count++;
//                }
//                else {
//                    lastNotMet = i;
//                }
//                //k++;
//                //if (k >= maxNumberOfDigitsToCountIn)
//                //    break;
//            }
//            int notAlright = maxNumberOfDigitsToCountIn - count;
//            // the numbers don't follow the same pattern as rational numbers
//            // it's not possible to move the error more to the back, because this would destroy the long period
//            //if (notAlright == 1) {
//            //    return (double)(d.second.size()-1- lastNotMet) / (d.second.size()+1);
//            //}
//            antiConvergence[d.second]++;
//            if (antiConvergence[d.second] > antiConvergenceBreakThreshold)
//                return 0;
//            return notAlright;
//        }
//    });
//    double n = pow(checkkky, 1.0 / 1.5);
//    return n;
//}
//
//double calculateConstantForDesiredDNumberInitialStringPeriod(const std::string &desiredInitialPeriodString, const std::string& base = "10", int maxNumberOfDigitsToCountIn = -1, int EPSILON = 1) {
//    if (maxNumberOfDigitsToCountIn < 0)
//        maxNumberOfDigitsToCountIn = desiredInitialPeriodString.size();
//    int iterations = 0;
//    double checkkky = errorReduce(0.1, 9.0, 0.0, [&](double v, double thresh) -> double {
//        iterations++;
//        DNumber d = ditherNumberFromRight(base, pow(v, 1.0 / 1.5), maxNumberOfDigitsToCountIn);
//        double error = 0;
//        for (int i = 0; i < desiredInitialPeriodString.size(); ++i) {
//            int delta = 0;
//            int i2 = i;
//            if (i2 < d.second.size())
//                if (d.second[i2] == '*')
//                    i2++;
//            if (i2 < d.second.size())
//                delta = (int)d.second[i2] - (int)desiredInitialPeriodString[i];
//            error += delta * delta;
//        }
//        if ((iterations % 20) == 0) {
//            printf("lolo:%d\n", iterations);
//            printNumber2(d.second);
//        }
//        error = sqrt(error);
//        if (error < EPSILON)
//            error = 0;
//        return error;
//        });
//    double n = pow(checkkky, 1.0 / 1.5);
//    return n;
//}
//
//// dunno whats that for..
//double calculateConstantForDesiredDNumberByPeriodSizes(const std::string& desiredInitialPeriodString, const std::string& base = "10", int maxNumberOfDigitsToCountIn = 100000, int EPSILON = 1) {
//    int iterations = 0;
//    double checkkky = errorReduce(0.1, 9.0, 0.0, [&](double v, double thresh) -> double {
//        iterations++;
//        DNumber d = ditherNumberFromRight(base, pow(v, 1.0 / 1.5), maxNumberOfDigitsToCountIn);
//        std::string sectionString = toNumber(v/1.5);
//        sectionString.erase(sectionString.begin() + dotPos(sectionString));
//        double value = 0;
//        int i = 0;
//        if (i < d.second.size())
//            if (d.second[i] == '*')
//                i++;
//        int thisSectionCount = 0;
//        int thisSectionChar = sectionString.size() > i ? sectionString[i] : -1;
//        Number b;
//        for (; i < sectionString.size(); ++i) {
//            if (sectionString[i] == thisSectionChar) {
//                thisSectionCount++;
//            }
//            else {
//                b = multiplyWithTen(b);
//                b = add(b, toNumber(thisSectionCount+1));
//                thisSectionCount = 0;
//                thisSectionChar = sectionString[i];
//            }
//        }
//        b = multiplyWithTen(b);
//        b = add(b, toNumber(thisSectionCount + 1));
//        std::string b2 = b.substr(0, dotPos(b));
//        double error = 0;
//        for (int i = 0; i < desiredInitialPeriodString.size() && i < b2.size(); ++i) {
//            int delta = 0;
//            delta = (int)b2[b2.size()-1-i] - (int)desiredInitialPeriodString[desiredInitialPeriodString.size()-1-i];
//            error += delta * delta;
//        }
//        error = sqrt(error);
//        if ((iterations % 50) == 0) {
//            printf("lolo:%d, %f, %f\n", iterations, error, v);
//            printNumber2(b2);
//        }
//        return error;
//        });
//    double n = pow(checkkky, 1.0 / 1.5);
//    return n;
//}
//
////r3 = 3.4759626170161191 the old ones
////r7 = 3.0625111171076562 the old ones
////     3.4759617871283073
////     3.0620907600389233
//
////r3 = calculateConstantForDesiredDNumberSingleDigitPeriod('3', "10", 60, 1);
////r7 = calculateConstantForDesiredDNumberSingleDigitPeriod('7', "10", 60, 1);
//
//std::string flipString(const std::string& s) {
//    std::string r;
//    for (int i = 0; i < s.size(); ++i) {
//        r += s[s.size() - 1 - i];
//    }
//    return r;
//}
//
//double approximatePi() {
//    //int iterations = 100000000;
//    //int p = 0;
//    //srand(1);
//    //for (int i = 0; i < iterations; ++i) {
//    //    double dx = (double)rand() / RAND_MAX;
//    //    double dy = (double)rand() / RAND_MAX;
//    //    double dist = sqrt(dx * dx + dy * dy);
//    //    if (dist <= 1.0)
//    //        p++;
//    //}
//    //double pi = (double)p / iterations * 4.0;
//
//    double p = 0;
//    //int SEGMENTS = 1000;
//    //for (int x = 0; x < SEGMENTS; ++x) {
//    //    for (int y = 0; y < SEGMENTS; ++y) {
//    //        double dx = (double)x / SEGMENTS;
//    //        double dy = (double)y / SEGMENTS;
//    //        double dist = sqrt(dx * dx + dy * dy);
//    //        if (dist <= 1.0)
//    //            p++;
//    //    }
//    //}
//
//    int count = 10000000;
//    for (int x = 0; x < count; ++x) {
//        double dx = (double)x / (count-1);
//        double dy = (double)x / (count-1);
//        dy *= dx;
//        double dist = sqrt(dx * dx + dy * dy);
//        if (dist <= 1.0)
//            p++;
//    }
//
//    double pi = asin(1.0) * 2.0;
//
//    int k = 0;
//    return pi;
//}
//
//Number multiplyWith0to9(const Number &n, int number) {
//    switch (number) {
//    case 0:return toNumber(0);
//    case 1:return n;
//    case 2:return multiplyWithTwo(n);
//    case 3:return multiplyWithThree(n);
//    case 4:return multiplyWithFour(n);
//    case 5:return multiplyWithFive(n);
//    case 6:return multiplyWithSix(n);
//    case 7:return multiplyWithSeven(n);
//    case 8:return multiplyWithEight(n);
//    case 9:return multiplyWithNine(n);
//    }
//    return toNumber(-1);
//}
//
//// the binary way but as tenary :)
//Number multiply(Number a, Number b) {
//    Number result;
//    Number adder = fixAll(b);
//    Number number = fixAll(a);
//    makeSameBase(adder, number);
//    int dp = number.length() - dotPos(number);
//    number.erase(number.length() - dp,1);
//    adder.erase(adder.length() - dp,1);
//    dp--;
//    adder = fixAll(adder);
//    for (int i = 0; i < dp*2; ++i) {
//        adder = divideByTen(adder);
//    }
//    for (int i = 0; i < number.size(); ++i) {
//        int letterHere = number[number.size() - 1 - i];
//        if (letterHere != '.') {
//            int multiplyBy = letterHere - '0';
//            result = add(result, multiplyWith0to9(adder, multiplyBy));
//        }
//        adder = multiplyWithTen(adder);
//    }
//    double res = fromNumber(a) * fromNumber(b);
//    fixAll(result);
//    return result;
//}
//
//Number pin = "3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679821480865132823066470938446095505822317253594081284811174502841027019385211055596446229489549303819644288109756659334461284756482337867831652712019091456485669234603486104543266482133936072602491412737245870066063155881748815209209628292540917153643678925903600113305305488204665213841469519415116094330572703657595919530921861173819326117931051185480744623799627495673518857527248912279381830119491298336733624406566430860213949463952247371907021798609437027705392171762931767523846748184676694051320005681271452635608277857713427577896091736371787214684409012249534301465495853710507922796892589235420199561121290219608640344181598136297747713099605187072113499999983729780499510597317328160963185950244594553469083026425223082533446850352619311881710100031378387528865875332083814206171776691473035982534904287554687311595628638823537875937519577818577805321712268066130019278766111959092164201989380952572010654858632788659361533818279682303019520353018529689957736225994138912497217752834791315155748572424541506959508295331168617278558890750983817546374649393192550604009277016711390098488240128583616035637076601047101819429555961989467678374494482553797747268471040475346462080466842590694912933136770289891521047521620569660240580381501935112533824300355876402474964732639141992726042699227967823547816360093417216412199245863150302861829745557067498385054945885869269956909272107975093029553211653449872027559602364806654991198818347977535663698074265425278625518184175746728909777727938000816470600161452491921732172147723501414419735685481613611573525521334757418494684385233239073941433345477624168625189835694855620992192221842725502542568876717904946016534668049886272327917860857843838279679766814541009538837863609506800642251252051173929848960841284886269456042419652850222106611863067442786220391949450471237137869609563643719172874677646";
//
//Number expm1(int exponent, int taylorLength) {
//    Number result;
//    Number quadrippler = toNumber(exponent);
//    for (int i = 0; i < taylorLength; ++i) {
//    }
//    return result;
//}
//
//Number ehochx(const Number &k, int precission = 10) {
//    Number result;
//    Number up = toNumber(1);
//    Number down = toNumber(1);
//    for (int i = 1; i < precission; ++i) {
//        up = multiply(up, k);
//        down = multiply(down, toNumber(i));
//    }
//    // geht noch nicht keine division (reciproke multiplication would be cool)
//    return result;
//}
//
//// todo: implement that for subtraction :)
//Number tennerComplement(const Number& n) { // not completed yet
//    Number k;
//    Number addi;
//    Number einsi = n;
//    int lastEinsi = -1;
//    for (int i = 0; i < einsi.size(); ++i) {
//        if (einsi[i] >= '0' && einsi[i] <= '9') {
//            einsi[i] = '0';
//            lastEinsi = i;
//        }
//    }
//    if (lastEinsi != -1) {
//        einsi[lastEinsi] = '1';
//    }
//    for (int i = n.size() - 1; i >= 0; --i) {
//        if (n[i] != '.') {
//            char c = (((9 - (n[i] - '0')) % 10) + '0');
//            k = c + k;
//        }
//        else {
//            k = '.' + k;
//        }
//    }
//    k = add(k,einsi);
//    return k;
//}
//
//// very funny way to divide not really a good way at all // not completed yet
//Number divideByTrying(Number dividend, const Number& divisor, int maxIterations = 1000) {
//    int iterations = 0;
//    // errorReduce would be needed as Number version (so it's not exact enough)
//    // and we somehow need a subtraction or something to compare the values
//    double v = errorReduce(0, 100, 0, [&](double v, double thresh) -> double {
//        iterations++;
//        if (iterations > maxIterations)
//            return 0;
//        Number k = multiply(divisor, toNumber(v));
//        makeSameBase(k, dividend);
//        if (k.empty())
//            return 0;
//        int dp = dotPos(k);
//        dp--;
//        Number error = toNumber(0);
//        Number err = toNumber(1.0);
//        for (int i = 0; i < dp; ++i)
//            err = multiplyWithTen(err);
//        for (int i = 0; i < k.size(); ++i) {
//            if (k[i] != '.' && dividend[i] != '.') {
//                Number delta = toNumber((k[i]-'0') ^ (dividend[i]-'0'));
//                delta = multiply(delta,err);
//                error = add(error,multiply(delta,delta));
//            }
//            err = divideByTen(err);
//        }
//        //for (int i = 0; i < k.size(); ++i) {
//        //    if (k[i] != dividend[i]) {
//        //        int d = (dividend[i] - '0') - (k[i] - '0');
//        //        return fabs(d);
//        //    }
//        //}
//        return fromNumber(error); // this would need direct evaluation not a floating point number here
//    });
//    Number k = toNumber(v);
//    return k;
//}
//
////const double r1_prec = calculateConstantForDesiredDNumberSingleDigitPeriod('1', "10", precisiona, 0, true);
////const double r5_prec = calculateConstantForDesiredDNumberSingleDigitPeriod('5', "10", precisiona, 0, true);
//const double r1const5000mean = 0.95673162883552897994832164840772747993469238281250;
//const double r5const5000mean = 1.92788416971737097682648709451314061880111694335938;
//const double r5r1const5000mean = r5const5000mean - r1const5000mean;
//
//std::string intensity(const std::string& n) {
//    const char digitsToIntensities[] = {
//        ' ',
//        '.',
//        ',',
//        ';',
//        'i',
//        'I',
//        'O',
//        'Q',
//        '0',
//        'B',
//        '&',
//        '$',
//    };
//    std::string r;
//    for (int i = 0; i < n.size(); ++i) {
//        char c = n[i];
//        if (c >= '0' && c <= '9')
//            c = digitsToIntensities[c - '0'];
//        else if (c == '.') c = '|';
//        r += c;
//    }
//    return r;
//}
//
////#include "nativefiledialog/nfd.h"
//
//#define MAXWIDTH 1000
//#define MAXHEIGHT 4096
//
//unsigned int pictureS[MAXWIDTH * MAXHEIGHT];
//unsigned int pictureWriteOut[MAXWIDTH * MAXHEIGHT];
//int pictureWidth;
//int pictureHeight;
//
//void abort_(const char* s, ...)
//{
//    va_list args;
//    va_start(args, s);
//    vfprintf(stderr, s, args);
//    fprintf(stderr, "\n");
//    va_end(args);
//    abort();
//}
//
//void SavePNG(const char* file_name, bool flipped)
//{
//    png_structp png_ptr;
//    png_infop info_ptr;
//
//    /* create file */
//    FILE* fp = fopen(file_name, "wb");
//    if (!fp)
//        abort_("[write_png_file] File %s could not be opened for writing", file_name);
//
//
//    /* initialize stuff */
//    png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
//
//    if (!png_ptr) {
//        fclose(fp);
//        abort_("[write_png_file] png_create_write_struct failed");
//    }
//
//    info_ptr = png_create_info_struct(png_ptr);
//    if (!info_ptr) {
//        fclose(fp);
//        abort_("[write_png_file] png_create_info_struct failed");
//    }
//
//    if (setjmp(png_jmpbuf(png_ptr))) {
//        fclose(fp);
//        abort_("[write_png_file] Error during init_io");
//    }
//
//    png_init_io(png_ptr, fp);
//
//    /* write header */
//    if (setjmp(png_jmpbuf(png_ptr))) {
//        fclose(fp);
//        abort_("[write_png_file] Error during writing header");
//    }
//
//    png_set_IHDR(png_ptr, info_ptr, pictureWidth, pictureHeight,
//        8, PNG_COLOR_TYPE_RGBA, PNG_INTERLACE_NONE,
//        PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);
//
//    png_write_info(png_ptr, info_ptr);
//
//
//    /* write bytes */
//    if (setjmp(png_jmpbuf(png_ptr))) {
//        fclose(fp);
//        abort_("[write_png_file] Error during writing bytes");
//    }
//
//    png_bytep* row_pointers;
//
//    unsigned char* array = new unsigned char[pictureWidth * pictureHeight * 4 + 4];
//    for (int i = pictureWidth * pictureHeight; i >= 0; --i)
//    {
//        array[i * 4 + 0] = ((unsigned char*)pictureWriteOut)[i * 4 + (flipped ? 2 : 0)];
//        array[i * 4 + 1] = ((unsigned char*)pictureWriteOut)[i * 4 + 1];
//        array[i * 4 + 2] = ((unsigned char*)pictureWriteOut)[i * 4 + (flipped ? 0 : 2)];
//        array[i * 4 + 3] = ((unsigned char*)pictureWriteOut)[i * 4 + 3];
//    }
//
//
//    // Allocate pointers...
//    row_pointers = new png_bytep[pictureHeight];
//
//    for (int i = 0; i < pictureHeight; i++)
//        row_pointers[i] = (png_bytep)(array + (i)*pictureWidth * 4); // we flip it
//
//    png_write_image(png_ptr, row_pointers);
//
//    /* end write */
//    if (setjmp(png_jmpbuf(png_ptr)))
//        abort_("[write_png_file] Error during end of write");
//
//    png_write_end(png_ptr, NULL);
//
//    /* cleanup heap allocation */
//    delete[] row_pointers;
//
//    fclose(fp);
//}
//
//int yline = 0; // actually x line due to "fast programming mode"
//
//void saveImage(const std::string& name) {
//    memcpy(pictureWriteOut, pictureS, 4 * pictureWidth * pictureHeight);
//    SavePNG(name.c_str(), false);
//}
//
//void newImage() {
//    yline = 0;
//    pictureWidth = MAXWIDTH;
//    pictureHeight = 1;
//    memset(pictureS, 0, MAXWIDTH * MAXHEIGHT * 4);
//}
//
//void putPixel(int x, int y, unsigned int color) {
//    pictureS[x + y * pictureWidth] = color;
//}
//
//std::string pictureName = "hello.png";
//void addIntensityToImage(const std::string &nm) {
//    int oneTabSize = 12;
//    double xscale = oneTabSize * 0.9;
//    for (int i = 0; i < nm.size(); ++i) {
//        double x = (double)i * oneTabSize;
//        char c = nm[i];
//        if (c >= '0' && c <= '9') {
//            c -= '0';
//            x += xscale * c / 9.0;
//            putPixel(yline, x, 0xffffffff);
//        }
//        else if (c == '.') {
//            putPixel(yline, x, 0xff0000ff);
//        }
//        else if (c == ' ') {
//            //putPixel(yline, x, 0xff00ff00);
//        }
//        if (x >= pictureHeight) {
//            pictureHeight = x;
//        }
//    }
//    yline++;
//    if (pictureWidth > 10 && pictureHeight > 10) {
//        saveImage(pictureName);
//    }
//}
//
//std::string intensityString(const double n) {
//    Number extend = toNumber(1000.00000000001);
//    Number nm = toNumber(n);
//    makeSameBase(extend, nm);
//    return intensity(nm);
//}
//
//void outputPeriodRepeatToImage(const Number &vorkommaA, const Number& vorkommaB, const double n, int mostDigitsOf = 0) {
//    Number extend = toNumber(100.00000000001);
//    Number nff = vorkommaA; // resultlast.first r5r1constant*verhaeltnis
//    Number nfff = vorkommaB; // resulthere.first startverhaeltnis*verhaeltnis
//    Number nf = std::to_string(mostDigitsOf);
//    Number nn = toNumber(n);
//    makeSameBase(extend, nn);
//    makeSameBase(extend, nff);
//    makeSameBase(extend, nfff);
//    std::string combined = nff + "     " + nf + "  " + nn + "   " + nfff;
//    addIntensityToImage(combined);
//}
//
//double checkVerhaeltnis(const Number &baseNumber, double ditherValHere, double ditherValLast, int digits) {
//    DNumber resultHere = ditherNumberFromRight(baseNumber, ditherValHere, digits);
//    DNumber resultLast = ditherNumberFromRight(baseNumber, ditherValLast, digits);
//    if (!resultHere.second.empty() && resultHere.second[0] == '*')
//        resultHere.second = resultHere.second.substr(1, resultHere.second.size() - 1);
//    if (!resultLast.second.empty() && resultLast.second[0] == '*')
//        resultLast.second = resultLast.second.substr(1, resultLast.second.size() - 1);
//    int minSize = resultHere.second.size() < resultLast.second.size() ? resultHere.second.size() : resultLast.second.size();
//    int validDigits = minSize;
//    const int convergingSize = 10;
//    validDigits -= convergingSize; // this is a fix, so that the initial converging part doesn't influences the result too much
//    if (validDigits > 2) {
//        resultLast.second = resultLast.second.substr(0, validDigits); // this takes the later part of the period, maybe it's better to check that anyway
//        resultHere.second = resultHere.second.substr(0, validDigits);
//
//        double approxPeriodLengthHere = validDigits;
//        if (validDigits != countSameDigits(resultHere.second))
//            approxPeriodLengthHere = (double)validDigits / ((double)validDigits - (double)countSameDigits(resultHere.second));
//
//        double approxPeriodLengthLast = validDigits;
//        if (validDigits != countSameDigits(resultLast.second))
//            approxPeriodLengthLast = (double)validDigits / ((double)validDigits - (double)countSameDigits(resultLast.second));
//
//        double approxPeriodDistance = (double)approxPeriodLengthHere / approxPeriodLengthLast;
//        outputPeriodRepeatToImage(resultLast.first, resultHere.first, approxPeriodDistance, maxI);
//        printf("periodverhaeltnis:[%s]%f period1:%f period2:%f ditherNumber1:%.10f ditherNumber2:%.10f\n", intensityString(approxPeriodDistance).c_str(), approxPeriodDistance, approxPeriodLengthHere, approxPeriodLengthLast, ditherValHere, ditherValLast);
//        return approxPeriodDistance;
//    }
//    return -1;
//}
//
//void printNumber3(const Number& n) {
//    printf("%s\n", intensity(n).c_str());
//}
//
//void testi() {
//    newImage();
//    int a = 2442;
//    int bx = 3443;
//    Number tc = tennerComplement(toNumber(bx));
//    Number r = add(toNumber(a),tc); // tobe implemented, hopefully :)
//    int ab = a + bx;
//    //int debug = 1;
//
//    //Number resa = divideByTrying(toNumber(33), toNumber(11), 2500);
//    divideByLessThan2Test("20.0", 1.5);
//    approximatePi();
//
//    //for (int j = 1; j < 10; ++j) {
//    //    double factor = 1.0 / j;
//    //    printf("-------------------------------------\n");
//    //    for (int i = 0; i < 10; ++i) {
//    //        double factorHere = factor * i;
//    //        printf("i,j,%f->\n", i,j,factorHere);
//    //        Number r = toNumber(rand());
//    //        DNumber k_1x = ditherNumberFromRight(r.substr(0,dotPos(r)), 0.95673162883552898 + factor * i);
//    //        int debugk_1x = countSameDigits(k_1x.second);
//    //        printf("%d\n", debugk_1x);
//    //        printDNumber(k_1x);
//    //    }
//    //}
//    printf("-------------------------------------\n");
//    printf("-------------------------------------\n");
//    printf("-------------------------------------\n");
//    printf("-------------------------------------\n");
//    printf("-------------------------------------\n");
//
//    int precisiona = 5000;
//    //const double r1_prec = calculateConstantForDesiredDNumberSingleDigitPeriod('1', "10", precisiona, 0, true);
//    //const double r5_prec = calculateConstantForDesiredDNumberSingleDigitPeriod('5', "10", precisiona, 0, true);
//    //const double point5 = fabs((r5_prec-trunc(r5_prec)) - (r1_prec-trunc(r1_prec)));
//
//    //ditherVal += r5r1const5000mean*200; // period of 18
//
//    srand(1);
//    for (int j = 0; j < 4; ++j) {
//        const int digits = 2000;
//        Number baseNumber = toNumber(23747935);
//        double ditherVal = r1const5000mean;
//        printf("baseNumber: %s kommaStellen:%d\n", baseNumber.c_str(), digits);
//        double verhaeltnis = 8;
//        int verhaeltnisCount = 125;
//
//
//        //double checkycheckor = 0.001;
//        //double multiply = 1.05;
//        //double add = 0.0;
//        //double power = 1.0;
//
//        int maxihj = 50;
//        for (int i = 0; i < maxihj; ++i) {
//            double checkycheckor = 0.01;
//            double baseFrequency = 1;
//            double power = 1.001 * baseFrequency;// 1.0 / 1.001;
//            double multiply = 1.002 * baseFrequency;
//            double add = 0.0;// 1.0 / 1.001;
//            double startVerhaeltnis = (double)i / maxihj * 2.33;
//
//            pictureName = "anim_" + std::to_string(i) + ".png";
//            newImage();
//            while (true) {
//                verhaeltnis = checkycheckor;
//                checkycheckor = pow(checkycheckor, power);
//                checkycheckor *= multiply;
//                checkycheckor += add;
//                if (yline > 500)
//                    break;
//                // r5r1const5000mean ist ne constante die die periode in der länge modifiziert und selbst wenn man sie erhöht ein schwingungspattern in der Periodenlänge erzeugt
//                checkVerhaeltnis(baseNumber, ditherVal + verhaeltnis * startVerhaeltnis, ditherVal + verhaeltnis * r5r1const5000mean, digits);
//                //verhaeltnis *= (faktorialCounter + 1);
//                verhaeltnis *= 10;
//            }
//        }
//        //saveImage("hello.png");
//        checkVerhaeltnis(baseNumber, ditherVal, ditherVal + 8 * 9 * 10 * r5r1const5000mean, digits);
//        printf("------------------------\n");
//        for (int i = 0; i <= 100; ++i) {
//            double ditherValHere = ditherVal + r5r1const5000mean * 10 * i;
//            double ditherValLast = ditherVal + r5r1const5000mean * 0;
//            printf("%d: ", i);
//            double verhaeltniss = checkVerhaeltnis(baseNumber, ditherValHere, ditherValLast, digits);
//        }
//    }
//
//    //printDNumber(result);
//
//    printf("r1: %.50f\n", r1const5000mean);
//    printf("r5: %.50f\n", r5const5000mean);
//    printf("r5-r1: %.50f\n", r5r1const5000mean);
//
//    //pow(3.0, 1.0 / 1.5); // looks interessting
//    int lastOptimumWas = 0;
//    for (int i = 0; i < 400; ++i) {
//        //double here = i*0.025 + 0.95673162883552898;
//        double here = r1const5000mean + r5r1const5000mean*i;
//        DNumber k_1x = ditherNumberFromRight("10", here, 10000);
//        int debugk_1x = countSameDigits(k_1x.second);
//        printf("%d\n", i);
//        //if (debugk_1x > lastOptimumWas) {
//            lastOptimumWas = debugk_1x;
//            printDNumber(k_1x);
//            printf("up there is: %i,%f->\n", i, here);
//            printf("%d\n", debugk_1x);
//        //}
//    }
//
//    int debug = 1;
//
//    // there are no irational numbers in dither space :)
//
//    //std::string kr = flipString("31415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679821480865132823066470938446095505822317253594081284811174502841027019385211055596446229489549303819644288109756659334461284756482337867831652712019091456485669234603486104543266482133936072602491412737245870066063155881748815209209628292540917153643678925903600113305305488204665213841469519415116094330572703657595919530921861173819326117931051185480744623799627495673518857527248912279381830119491298336733624406566430860213949463952247371907021798609437027705392171762931767523846748184676694051320005681271452635608277857713427577896091736371787214684409012249534301465495853710507922796892589235420199561121290219608640344181598136297747713099605187072113499999983729780499510597317328160963185950244594553469083026425223082533446850352619311881710100031378387528865875332083814206171776691473035982534904287554687311595628638823537875937519577818577805321712268066130019278766111959092164201989380952572010654858632788659361533818279682303019520353018529689957736225994138912497217752834791315155748572424541506959508295331168617278558890750983817546374649393192550604009277016711390098488240128583616035637076601047101819429555961989467678374494482553797747268471040475346462080466842590694912933136770289891521047521620569660240580381501935112533824300355876402474964732639141992726042699227967823547816360093417216412199245863150302861829745557067498385054945885869269956909272107975093029553211653449872027559602364806654991198818347977535663698074265425278625518184175746728909777727938000816470600161452491921732172147723501414419735685481613611573525521334757418494684385233239073941433345477624168625189835694855620992192221842725502542568876717904946016534668049886272327917860857843838279679766814541009538837863609506800642251252051173929848960841284886269456042419652850222106611863067442786220391949450471237137869609563643719172874677646");
//    //std::string kr = "135";
//    //const double r57 = calculateConstantForDesiredDNumberByPeriodSizes(kr, "10", 40, 1);
//
//    //std::string kr = "31415926535897932";
//    //std::string kr = "3333333333333333333333333333333333333333333333333333333333333333333333333333333333333333";
//    //const double r5 = calculateConstantForDesiredDNumberInitialStringPeriod(kr, "10", -1, 1);
//    bool meanSquares = true;
//        // 500 digits
//        // 1 = 0.95672855568135196
//        // 2 = 1.0913546777421734
//        // 3 = 1.0528880595330141
//        // 4 = 1.0576995003539853
//        // 5 = 1.9278915557977128
//        // 6 = 4.0528933678267345
//        // 7 = 4.1442357267503036
//        // 8 = 4.2355936087987844
//        // 9 = 4.3267487105901559
//
//        // 1500 digits
//        //0.95673040869297221
//        //1.0913486688439284
//        //1.0528866809025061
//        //1.0576933106362807
//        //1.9278823867782122
//        //4.0528863849627479
//        //4.1442332077704469
//        //4.2355790494594112
//        //4.3267487107561902
//
//        // 5000 digits
//        //0.95673162883552898
//        //1.0913467011377385
//        //1.0528849576131016
//        //1.0576929667625439
//        //1.9278841697173710
//        //4.0528854297591872
//        //4.1442313185350512
//        //4.2355775995653229
//        //4.3267480306792905
//
//    {
//        int precision = 5000;
//        const double rb3_1 = calculateConstantForDesiredDNumberSingleDigitPeriod('3', "10", precision, 0, true);
//        const double rb3_1b = calculateConstantForDesiredDNumberSingleDigitPeriod('3', "10", precision, 1, false);
//        precision /= 10;
//        const double rb3_2 = calculateConstantForDesiredDNumberSingleDigitPeriod('3', "10", precision, 0, true);
//        const double rb3_2b = calculateConstantForDesiredDNumberSingleDigitPeriod('3', "10", precision, 1, false);
//        precision /= 10;
//        const double rb3_3 = calculateConstantForDesiredDNumberSingleDigitPeriod('3', "10", precision, 0, true);
//        const double rb3_3b = calculateConstantForDesiredDNumberSingleDigitPeriod('3', "10", precision, 1, false);
//    }
//
//    int precision = 5000; // meanSquares ? kr.size() : 600;
//    const double r1 = calculateConstantForDesiredDNumberSingleDigitPeriod('1', "10", precision, 0, meanSquares);
//    const double r2 = calculateConstantForDesiredDNumberSingleDigitPeriod('2', "10", precision, 0, meanSquares);
//    const double r3 = calculateConstantForDesiredDNumberSingleDigitPeriod('3', "10", precision, 0, meanSquares);
//    const double r4 = calculateConstantForDesiredDNumberSingleDigitPeriod('4', "10", precision, 0, meanSquares);
//    const double r5 = calculateConstantForDesiredDNumberSingleDigitPeriod('5', "10", precision, 0, meanSquares);
//    const double r6 = calculateConstantForDesiredDNumberSingleDigitPeriod('6', "10", precision, 0, meanSquares);
//    const double r7 = calculateConstantForDesiredDNumberSingleDigitPeriod('7', "10", precision, 0, meanSquares);
//    const double r8 = calculateConstantForDesiredDNumberSingleDigitPeriod('8', "10", precision, 0, meanSquares);
//    const double r9 = calculateConstantForDesiredDNumberSingleDigitPeriod('9', "10", precision, 0, meanSquares);
//
//    DNumber k_1 = ditherNumberFromRight("12345", r1);
//    DNumber k_2 = ditherNumberFromRight("12345", r2);
//    DNumber k_3 = ditherNumberFromRight("12345", r3);
//    DNumber k_4 = ditherNumberFromRight("12345", r4);
//    DNumber k_5 = ditherNumberFromRight("12345", r5);
//    DNumber k_6 = ditherNumberFromRight("12345", r6);
//    DNumber k_7 = ditherNumberFromRight("12345", r7);
//    DNumber k_8 = ditherNumberFromRight("12345", r8);
//    DNumber k_9 = ditherNumberFromRight("12345", r9);
//
//    DNumber kb10_1 = ditherNumberFromRight("12345", r1);
//    DNumber kb10_2 = ditherNumberFromRight("12345", r2);
//    DNumber kb10_3 = ditherNumberFromRight("12345", r3);
//    DNumber kb10_4 = ditherNumberFromRight("12345", r4);
//    DNumber kb10_5 = ditherNumberFromRight("12345", r5);
//    DNumber kb10_6 = ditherNumberFromRight("12345", r6);
//    DNumber kb10_7 = ditherNumberFromRight("12345", r7);
//    DNumber kb10_8 = ditherNumberFromRight("12345", r8);
//    DNumber kb10_9 = ditherNumberFromRight("12345", r9);
//
//
//    //DNumber k3 = ditherNumberFromRight("1234", r3);
//    DNumber k7_1 = ditherNumberFromRight("1", r7);
//    DNumber k7_2 = ditherNumberFromRight("2", r7);
//    DNumber k7_3 = ditherNumberFromRight("3", r7);
//    DNumber k7_4 = ditherNumberFromRight("4", r7);
//    DNumber k7_5 = ditherNumberFromRight("5", r7);
//    DNumber k7_6 = ditherNumberFromRight("6", r7);
//    DNumber k7_7 = ditherNumberFromRight("7", r7);
//    DNumber k7_8 = ditherNumberFromRight("8", r7);
//    DNumber k7_9 = ditherNumberFromRight("9", r7);
//    DNumber k7_10 = ditherNumberFromRight("10", r7);
//    DNumber k7_20 = ditherNumberFromRight("20", r7);
//    DNumber k7_50 = ditherNumberFromRight("50", r7);
//    DNumber k7_100 = ditherNumberFromRight("100", r7);
//    DNumber k7_666 = ditherNumberFromRight("666", r7);
//    DNumber k7_1234 = ditherNumberFromRight("1234", r7);
//    DNumber k7_21234 = ditherNumberFromRight("21234", r7);
//    DNumber k3_1 = ditherNumberFromRight("1", r3);
//    DNumber k3_2 = ditherNumberFromRight("2", r3);
//    DNumber k3_3 = ditherNumberFromRight("3", r3);
//    DNumber k3_4 = ditherNumberFromRight("4", r3);
//    DNumber k3_5 = ditherNumberFromRight("5", r3);
//    DNumber k3_6 = ditherNumberFromRight("6", r3);
//    DNumber k3_7 = ditherNumberFromRight("7", r3);
//    DNumber k3_8 = ditherNumberFromRight("8", r3);
//    DNumber k3_9 = ditherNumberFromRight("9", r3);
//    DNumber k3_10 = ditherNumberFromRight("10", r3);
//    DNumber k3_20 = ditherNumberFromRight("20", r3);
//    DNumber k3_50 = ditherNumberFromRight("50", r3);
//    DNumber k3_100 = ditherNumberFromRight("100", r3);
//    DNumber k3_1234 = ditherNumberFromRight("1234", r3);
//    DNumber k3_21234 = ditherNumberFromRight("21234", r3);
//    int debug2 = 1;
//
//    //double r1 = calculateConstantForDesiredDNumberSingleDigitPeriod('1', precision, 1);
//    //double r2 = calculateConstantForDesiredDNumberSingleDigitPeriod('2', precision, 1);
//    //double r4 = calculateConstantForDesiredDNumberSingleDigitPeriod('4', precision, 1);
//    //double r5 = calculateConstantForDesiredDNumberSingleDigitPeriod('5', precision, 1);
//    //double r6 = calculateConstantForDesiredDNumberSingleDigitPeriod('6', precision, 1);
//    //double r8 = calculateConstantForDesiredDNumberSingleDigitPeriod('8', precision, 1);
//    //double r9 = calculateConstantForDesiredDNumberSingleDigitPeriod('9', precision, 1);
//
//    int iterations = 0;
//    const int EPSILON = 1; // dunno why this is needed // 42 iterations here
//    double checkkky = errorReduce(0.1, 9.0, 0.0, [&](double v, double thresh) -> double {
//        iterations++;
//        int fracPart = 1000;
//        int count = 0;
//        DNumber d = ditherNumberFromRight("10.0", pow(v, 1.0 / 1.5), fracPart);
//        for (int i = 0; i < d.second.size(); ++i) {
//            if (d.second[i] == '3')
//                count++;
//        }
//        return fabs(fracPart - count - EPSILON);
//    });
//    double n = pow(checkkky, 1.0 / 1.5);
//    DNumber b = ditherNumberFromRight("100.0", n);
//    int checky = 1;
//}
//
//
//
//const double rotationOfMoonAroundEarthInDays = 27.32;
//const double moonAngleToEarthAfterOneSolarYear = 360*366/(27.32) ; // 142.840409956 degree phase shift per year
//const double moonQuantizationToReach360Degrees = 360/142.840409956; // 2.52029520295
//const double oneMoonCycle = (360*2.52029520295)*142.840409956; //129600 years??
//const double enochConstantDaysFor8Years = 2912; // would give 364 days per year
//const double sunConstantDividedBy364 = 74381.7069363;
////const double sunConstantDividedBy365 = 74177.9214378;
////const double sunConstantDividedBy366 = 73975.2495213;
//
//// maybe it's ok or needed, if the moon is 0.5 shifted in a day after a perfection yuga, or computers are not sufficient to compute this
//// maybe this is needed to step through the numbers one letter by another (like the execution of a programm, consisting of 4 steps till the next programm step)
//// this means the programm would be defined by these full "yuga" steps but underlying this would be the actual "constant" the program this iterates over
//// maybe the yuga constant actually defines the program itself
//// this would be the program definition: hindu yuga year of perfection is 1728000 (would define 7 major program steps if no zero steps follow)
//// this program would be not so perfect: 1296000 (dunno if this is the right number for it)
//// this lesser perfect: 864000 (dunno if this is the right number for it)
//// this would be actually ?destructive?: 432000
//const double afterOnePerfectionYugaTheMoonWouldBeAtThisPhaseOnTheSameDayOfASolarYear = 1593955660.5166094; // these values may also be rounded. // a perfection yuga is 1728000 ?hindu? years whilst it could be that a year is not a solar year of 366 days there
//const double afterSecondOfTheseTheMoonWouldBeAtThisPhaseOnTheSameDayOfASolarYear  = 1195466745.3874571; // if this is actually 3/4 of the perfection yuga
//const double afterThisOfTheseTheMoonWouldBeAtThisPhaseOnTheSameDayOfASolarYear = 796977830.25830472; // if this is actually 2/4 of the perfection yuga
//const double afterOneKaliYugaTheMoonWouldBeAtThisPhaseOnTheSameDayOfASolarYear = 398488915.12915236; // if this is actually 1/4 of the perfection yuga
//const double afterCompleteYugaCycleTheMoonWouldBeAtThisPhaseOnTheSameDayOfASolarYear = 3984889151.2915239; // dunno what this is actually it's a multiplication by 10 of the perfection value
//// 1cycles + 0.75cycles + 0.5cycle + 0.25cycles = (0.25cycles * 10) or 2.5cycles :D sorry.. (moon would be at the same position of day if the number would consist of only the same ziffern)
//// after a complete ?yuga? cycle the fractional part is one cypher more of the current ?kali? "program"one is "integered" / "defractionalized" / "defloated" :) or a zero added to this kali
//// actually it's lastPerfection/4*10 = newPerfection = ; newPerfection = lastPerfection * 2.5 is actually the full program execution step and shifting the decimal dot by that
//// so 2.5 is the best constant there I assume by the data I got.. And somehow part of the programs execution unit / "program counter advancement".
//
//// first(1) hindu yuga year of perfection is 1728000 // earth rotations: 632448000
//// second?? seems to be 1296000 // earth rotations: 474336000
//// third?? seems to be  864000 // earth rotations: 316224000
//// last(4) hindu yug is 432000 // earth rotations: 158112000
//
//const double assumedSunConstant = 27074941.324810438; // whatever that is
//// truncated
//const double strange_estimated_NumberOfSunYears_NeededForAlignmentOf_SunMoonAndEarth_OnAYearStart = 8374.9552154541052; // 2827457443.0799184; // if physics is 100% perfect and without friction etc.. this would be some sort of "age" in years (solar years with 366 days each)
//const double strange_estimated22_NumberOfSunYears_NeededForAlignmentOf_SunMoonAndEarth_OnADayStart = 620.21098632812505; // 572100.00715989841; // period of 42 days and 1563 solar years // 0 stelle found at 22 iterations
//const double strange_estimated72_NumberOfSunYears_NeededForAlignmentOf_SunMoonAndEarth_OnADayStart = 35156.896203803350; // 32429707.008244935; // period of 277 days and 88605 solar year // 0 stelle found at 72 iterations
///// rounded (would add up as periods every cycle 0.5 more?) maybe yugas correct themselves by the 4 different parts?
//const double strange_estimatedy_NumberOfSunYears_NeededForAlignmentOf_SunMoonAndEarth_OnAYearStart = 8013.5472580671321; // 7391920.7250077892 // if physics is 100% perfect and without friction etc.. this would be some sort of "age" in years (solar years with 366 days each)
//const double strange_estimatedy22_NumberOfSunYears_NeededForAlignmentOf_SunMoonAndEarth_OnADayStart = 7207.0591796874996; // 6647993.5041325707 // period of 42 days and 1563 solar years // 0 stelle found at 22 iterations
//const double strange_estimatedy72_NumberOfSunYears_NeededForAlignmentOf_SunMoonAndEarth_OnADayStart = 18238.122314453125; // 16823355.497868545 // period of 277 days and 88605 solar year // 0 stelle found at 72 iterations
//
////Formulas(consequtive)
////0. ? sin() ? ; seems to incorporate two dimensions(positive / negative)
////1. deSign()->vminus1to1 * 2 - 2 (2 or -2); or +2 respective, dunno maybe fast way of binary "communicating"
////2. doubling / amplifying / contrast(2 + 2) or (-2 + -2) respective; *2
////3. ? stabilizing through ? amplifying the same via loop(2 * 2 * 2 * 2); potenzierung; ^ 1.07798593225 (stable circular orbit)
////4. going back in time(for instance from 2 * 2 * 2 * 2 to 2 * 2 * 2), from before looped of the alltime same; wurzel(1 / 1.09337413665)
////5. update to before impossible(error correction), (2 + 2 + 2) + 2; multiplyand add
////6. ? stabilizing through ? amplifying the same via loop);
///////6. going back into the before unknown corrected time steps (for instance from (2*2*2*2+2) to wurzel(2*2*2*2+2)) (somehow with an mercury aspect of instable y axis of planet?)
////7. doubling / amplifying / contrast; *2
////8. ? stabilizing through ? amplifying the same via loop); (stable circular orbit) potenzierung with different time timeframe(seems to be a smaller "loop") though the "y axis deviation" differs now largely from earth
////9. Pretty sure this is some sort of subtraction..Maybe a sign() function... No clue, yet..Lets Check
////Colors(consequtive) (*means not 100 % clear to me)
////0. * full spektrum plasma(which peaks in green according to science)
////1. grau
////2. * initialer farbton
////3. * partielle inversion(earth seems "to be" for reducing ? specific ? yellow parts ? ? )
////4. * recolor to ? red ? ; somehow combining all channels into in red whilst diminishing greenand blue(removing / diminishing of a channel(s))
////5. Lab->a in the center->seems to add yellow again to form orange and brown and white..whyever
////6. Lightness: HSB->S value + HSB->L value -> ? brightness and contrast increase ?
////7. Lab->b value
////8. HSB->S value increase(from uranus to neptune)
////9. Lab->a value Pluto somehow seems lila(some sort of brown / lila with some blue and white) to me.
////Mars to Neptune is Lab->b(somehow uranus seems to be in charge of that)
////Pluto may be an inverse sun colorwise ? ? (it's brown, Nasa picture also covers some blue)
//
//
//// For long, scientists suspected the life cycle of Earth's magnetic field to be 200 million years. // 332 and 416 million years ago
//// The magnetic poles flip approximately every 200,000 to 300,000 years
//// They found that between 332 and 416 million years ago, the strength of the geomagnetic field preserved in these rocks was less than quarter of what it is today, and similar to a previously identified period of low // // magnetic field strength that started around 120 million years ago. The researchers have coined this period “the Mid-Palaeozoic Dipole low (MPDL).”
//
//// from the current 366 earth rotations per year
//// first(1) hindu yuga year of perfection is 1,728,000 // earth rotations: 632448000
//// second?? seems to be 1,296,000 // earth rotations: 474336000
//// third?? seems to be 864,000 // earth rotations: 316224000
//// last(4) hindu yug is 432,000 // earth rotations: 158112000
//
//// 6000 days in the bible :D (this strange sun constant divided by 432000 also gives something with a 6er part 62.6734752889)
//// 27074941.324810438/6000.66666666666 = 4511.9888888888888
//// 45119.888
//// 0.666666 as daytime would be = 15.59999 (not interessting?:/)
//// if one year is a solar year than the earth would have rotated 241.56 times at these 0.66 "days"
//
////2: maxEsoteric:1
////15 : maxEsoteric : 2
////36 : maxEsoteric : 3
////42 : maxEsoteric : 4
////84 : maxEsoteric : 6
////126 : maxEsoteric : 8
////252 : maxEsoteric : 12
////576 : maxEsoteric : 13
////756 : maxEsoteric : 18
////1584 : maxEsoteric : 22
////1848 : maxEsoteric : 23
////2268 : maxEsoteric : 24
////2772 : maxEsoteric : 26
////3168 : maxEsoteric : 28
////4284 : maxEsoteric : 30
////4752 : maxEsoteric : 32
////5544 : maxEsoteric : 38
////9504 : maxEsoteric : 40
////11088 : maxEsoteric : 45
////16632 : maxEsoteric : 54
////31416 : maxEsoteric : 55
////33264 : maxEsoteric : 65
////49896 : maxEsoteric : 70
////66528 : maxEsoteric : 81
////94248 : maxEsoteric : 86
////133056 : maxEsoteric : 97
////199584 : maxEsoteric : 105
////266112 : maxEsoteric : 113
////282744 : maxEsoteric : 118
////399168 : maxEsoteric : 125
////465696 : maxEsoteric : 129
////576576 : maxEsoteric : 138
//
//const int rotationsOfTheEarthInOneYear = 366;
//const int secondsForOneEarthRotation = 86400; // 432 is halve (hindu yuga decrease time)
//// int earthTimeHours = 23;
//// int earthTimeMinutes = 56;
//// int earthTimeSeconds = 0;
//
//struct checkySetup {
//
//    // seems they formed only in a range of 100000000 years, so no clear distinction can be made here, yet
//    // there is some differency in notations from the net
//    // clearly stated is that jupiter was the first and planets after jupiter where formed next
//    // theory is that there was an astroid belt which split the solar system in two parts
//    // so these numbers here may be nonsense except that jupiter was first
//    const double ageFromSun_mercury = 4.503;
//    const double ageFromSun_venus = 4.503;
//    const double ageFromSun_earth = 4.543; // (4.53 to 4.58) from wikipedia // 2. potencing
//    const double ageFromSun_mars = 4.603; // 1. potencing
//    const double ageFromSun_jupiter = 4.603; // 1. correcting
//    const double ageFromSun_saturn = 4.503;
//    const double ageFromSun_uranus = 4.503;
//    double assumedSunConstant             = 27074941.324810438; // whatever that is
//    const double distanceToSun_mercury    = 35000000; // correction function
//    const double distanceToSun_venus      = 67000000; // scaling
//    const double distanceToSun_earth      = 93000000; // potencing
//    const double distanceToSun_mars       = 142000000; // potencing
//    const double distanceToSun_jupiter    = 484000000; // correction function
//    const double distanceToSun_saturn     = 889000000; // potencing
//    const double distanceToSun_uranus     = 1790000000; // scaling
//    const double distanceToSun_neptune    = 2880000000; // potencing
//    const double distanceToSun_pluto      = 3670000000; // correction function
//    const double distanceMercuySun = 350;
//    const double distanceVenusMercuy = 320; // decrease of distance
//    const double distanceEarthVenus = 260; // decrease of distance
//    const double distanceMarsEarth = 490;
//    const double distanceJupiterMars = 3420;
//    const double distanceSaturnJupiter = 4050;
//    const double distanceUranusSaturn = 9010;
//    const double distanceNeptuneUranus = 10900;
//    const double distancePlutoNeptune = 7900; // decrease of distance
//
//                                                        // these values should be adequate if planet distances above are right
//    double potenzKonstante1 = log(distanceToSun_earth) / log(distanceToSun_venus);      // earth 1.0181966246737080
//    double potenzKonstante2 = log(distanceToSun_mars) / log(distanceToSun_earth);       // mars 1.0230665481768368
//    double potenzKonstante3 = log(distanceToSun_saturn) / log(distanceToSun_jupiter);   // saturn 1.0304042718471391
//    double potenzKonstante4 = log(distanceToSun_neptune) / log(distanceToSun_uranus);   // neptune 1.0223217051095284
//    double multiplikationsKonstante1 = distanceToSun_venus / distanceToSun_mercury;     // venus 1.9142857142857144
//    double multiplikationsKonstante2 = distanceToSun_uranus / distanceToSun_saturn;     // uranus 2.0134983127109112
//
//    // this is some sort of error correction
//    // I have no clue on how to calculate the right values here
//    // mathematically correct and perfect seem to be these ones (that would be clean)
//    // but actually maybe would correct nothing by these
//    // these error corrections would be optimal for an "sunConstant of" 27074941.324810438
//    // mercury, pluto and jupiter seem to be weichensteller in astrology
//    double multiplyAndAddFactor1 = 2.0; // error correction through mercury
//    double multiplyAndAddFactor2 = 2.0; // error correction through jupiter
//    double multiplyAndAddFactor3 = 2.0; // error correction through pluto
//    double correctionKonstante1 = 20000000.0;    // error correction through mercury
//    double correctionKonstante2 = 200000000.0;   // error correction through jupiter
//    double correctionKonstante3 = 2000000000.0;  // error correction through pluto
//};
//
//double checky(double energyScales, checkySetup &k) {
//    double t = 0.5 * pi;
//    double r = sin(t) * k.assumedSunConstant * energyScales; // dunno what this value means
//    double mercury = r * k.multiplyAndAddFactor1 - k.correctionKonstante1; // not sure about this correction constant, but should be ok
//    double venus = mercury * k.multiplikationsKonstante1;
//    double earth = pow(venus , k.potenzKonstante1);
//    double mars = pow(earth, k.potenzKonstante2);
//    double jupiter = mars * k.multiplyAndAddFactor2 + k.correctionKonstante2; // seems to be a fehlerkorrektur
//    double saturn = pow(jupiter , k.potenzKonstante3);
//    double uranus = saturn * k.multiplikationsKonstante2;
//    double neptune = pow(uranus, k.potenzKonstante4);
//    double pluto = neptune * k.multiplyAndAddFactor3 - k.correctionKonstante3;
//    double error = (pluto - k.distanceToSun_pluto) * (pluto - k.distanceToSun_pluto);
//    return error;
//}
//
//#include <functional>
//
//double errorReduce(double rangeLow, double rangeHi, double threshold, std::function<double(double v, double threshold)> checky) {
//    double rv = 0;
//    int steps = 10000000;
//    double lastError = 10000000000000000.0;
//    double lowerBound = rangeLow;
//    double upperBound = rangeHi;
//    double calcConstant = 0.0;
//    while (lastError > 0.0001) { // kilometers
//        double k = (lowerBound + upperBound) / 2.0;
//        double error1 = checky((k + upperBound) / 2.0, threshold);
//        double error2 = checky((k + lowerBound) / 2.0, threshold);
//        error1 *= error1; // least squares
//        error2 *= error2; // least squares
//        if (error1 > error2) {
//            upperBound = k;
//            lastError = error1;
//        }
//        else {
//            lowerBound = k;
//            lastError = error2;
//        }
//        double lastErrorB = lastError;
//        rv = k;
//    }
//    return rv;
//}
//
//double reCalcCurrentSunConstant() { // whatever that is
//    checkySetup a1;
//    int steps = 10000000;
//    double lastError = 10000000000000000.0;
//    double lowerBound = 0.900000003;
//    double upperBound = 1.100000001;
//    double calcConstant = 0.0;
//    while (lastError > 0.0001) { // kilometers
//        double k = (lowerBound + upperBound) / 2.0;
//        double error1 = checky((k + upperBound) / 2.0, a1);
//        double error2 = checky((k + lowerBound) / 2.0, a1);
//        if (error1 > error2) {
//            upperBound = k;
//            lastError = error1;
//        }
//        else {
//            lowerBound = k;
//            lastError = error2;
//        }
//        double lastErrorB = lastError;
//        calcConstant = k * a1.assumedSunConstant;
//    }
//    return calcConstant;
//}
//
//double factorial(int v) {
//    double k = 1;
//    for (int i = 1; i <= v; ++i)
//        k *= i;
//    return k;
//}
//
//double keralaSin(double v, int iterations = 25) {
//    double p = v;
//    for (int i = 1; i <= iterations - 1; ++i) {
//        const double k = 2.0 * i + 1;
//        p += pow(v, k) / factorial(k) * ((i & 1) == 0 ? -1 : 1);
//    }
//    return p;
//}
//
//bool printStats = true;
//bool printNumber = true;
//bool printEsoteric = true;
//bool printYesNo = true;
//bool printEmptySection = true;
//FILE* numberOut = nullptr;
//FILE* additionalsOut = nullptr;
//FILE* processOut = nullptr;
//FILE* formulasOut = nullptr;
//FILE* esotericOut = nullptr;
//FILE* advancersOut = nullptr;
//
//#define PRINT(...) printf(__VA_ARGS__)
////#define PRINT(...) fprintf(out,__VA_ARGS__)
//
//#define FILEPRINT(_fileId_,...) if (_fileId_ != nullptr) {fprintf(_fileId_,__VA_ARGS__);}
//#define PRINTNUMBER(...) FILEPRINT(numberOut,__VA_ARGS__)
//#define PRINTADDITIONALS(...) FILEPRINT(additionalsOut,__VA_ARGS__)
//#define PRINTPROCESS(...) FILEPRINT(processOut,__VA_ARGS__)
//#define PRINTFORMULA(...) FILEPRINT(formulasOut,__VA_ARGS__)
//#define PRINTESOTERIC(...) FILEPRINT(esotericOut,__VA_ARGS__)
//#define PRINTADVANCERS(...) FILEPRINT(advancersOut,__VA_ARGS__)
//
//bool addLastSection = true;
//
//int printAllDividersOf(double k, bool onlyEsotericOutput = false, int esotericThreshold = 12, int *sectionsAboveThresh = nullptr, int *dividerCount = nullptr) {
//    if (printStats) {
//        PRINT(">>>>>>>NUMBER:%d\n", (int)k);
//        PRINTNUMBER("%d\n", (int)k);
//        PRINTADDITIONALS("\n[%d]\n", (int)k);
//    }
//    int sectionsAboveThreshold = 0;
//    int maxEsotericLineLength = 0;
//    int z = 1;
//    int lastOne = 0;
//    int il = 0;
//    int sectionLength = 0;
//    int yesCount = 0;
//    int noCount = 0;
//    int yesNumberQuerSum = 0;
//    int noNumberQuerSum = 0;
//    std::string esotericLine, biggestEsotericLine, lastEsotericLine, firstEsotericLine, secondEsotericLine;
//    std::string yesNumber, noNumber;
//    std::string lastClass;
//    int counti = 0;
//    int allSectionsTogether = 0;
//    std::vector<int> dividerProgress;
//    const int ki = k;
//    for (int i = 1; i <= k; ++i) {
//        if ((ki % i) == 0) {
//            allSectionsTogether++; // cycle count for periodof i
//            std::string classs = ((i - lastOne) == il ? "=" : (i - lastOne) > il ? "y" : ".");
//            std::string prefix = "";
//            if (classs == "=") {
//                if (lastClass == "y") {
//                    yesNumber += std::to_string(counti);
//                    yesNumberQuerSum += counti;
//                }
//                if (lastClass == ".") {
//                    noNumber += std::to_string(counti);
//                    noNumberQuerSum += counti;
//                }
//                if (esotericLine.length() > maxEsotericLineLength) {
//                    maxEsotericLineLength = esotericLine.length();
//                    biggestEsotericLine = esotericLine;
//                }
//
//                if (esotericLine.length() >= esotericThreshold) {
//                    if (firstEsotericLine.empty() && esotericLine != "=")
//                        firstEsotericLine = esotericLine;
//                    else
//                    if (secondEsotericLine.empty() && esotericLine != "=")
//                        secondEsotericLine = esotericLine;
//
//                    if (printStats) {
//                        PRINT(">[%d]->%d\n",z+1,sectionsAboveThreshold+1);
//                        PRINT("esoteric:\"%s\" length:%d\n", esotericLine.c_str(), esotericLine.length());
//                        PRINT("esoYes:%d:\"%s\" length:%d\n", yesCount, yesNumber.c_str(), yesNumber.size());
//                        PRINT("esoNo:%d:\"%s\" length:%d\n", noCount, noNumber.c_str(), noNumber.size());
//                        PRINTADDITIONALS(">[%d]->%d\n", z+1, sectionsAboveThreshold+1);
//                        PRINTADDITIONALS("esoteric:\"%s\" length:%d\n", esotericLine.c_str(), esotericLine.length());
//                        PRINTADDITIONALS("esoYes:%d:\"%s\" length:%d\n", yesCount, yesNumber.c_str(), yesNumber.size());
//                        PRINTADDITIONALS("esoNo:%d:\"%s\" length:%d\n", noCount, noNumber.c_str(), noNumber.size());
//                        PRINTADDITIONALS("DividersTillHere:%d\n", allSectionsTogether);
//                    }
//                    dividerProgress.push_back(allSectionsTogether);
//                    sectionsAboveThreshold++;
//                }
//                esotericLine.clear();
//                yesNumberQuerSum = 0;
//                noNumberQuerSum = 0;
//                yesNumber.clear();
//                noNumber.clear();
//                counti = 0;
//                lastClass = "";
//
//                if (sectionLength == 0) {
//                    //if (printStats)
//                    //    PRINT("empty section\n");
//                }
//                else {
//                    if (printStats) {
//                        PRINT("%d:section length:%d secs:%d yes:%d no:%d difference:%d\n", z - 1, allSectionsTogether, yesCount, noCount, yesCount - noCount);
//                    }
//                }
//                sectionLength = 0;
//                yesCount = 0;
//                noCount = 0;
//            }
//            else {
//                if (classs == "y")
//                    yesCount++;
//                else
//                    noCount++;
//                sectionLength++;
//            }
//            if (!onlyEsotericOutput) {
//                if (printStats) {
//                    PRINT("%d\t\t: delta:%d \t\tbigger:%s [%d]\n", i, i - lastOne, classs.c_str(), z);
//                }
//            }
//            else {
//                esotericLine += classs;
//                if (classs != lastClass) {
//                    if (lastClass == "y") {
//                        yesNumber += std::to_string(counti);
//                        yesNumberQuerSum += counti;
//                    }
//                    if (lastClass == ".") {
//                        noNumber += std::to_string(counti);
//                        noNumberQuerSum += counti;
//                    }
//                    counti = 0;
//                    lastClass = classs;
//                }
//                if (classs == "=")
//                    counti = 0;
//                counti++;
//            }
//            ++z;
//            il = i - lastOne;
//            lastOne = i;
//        }
//    }
//    //dividerProgress.push_back(allSectionsTogether);
//
//    if (lastClass == "y") {
//        yesNumber += std::to_string(counti);
//        yesNumberQuerSum += counti;
//    }
//    if (lastClass == ".") {
//        noNumber += std::to_string(counti);
//        noNumberQuerSum += counti;
//    }
//
//    lastEsotericLine = esotericLine;
//
//    if (esotericLine.length() >= esotericThreshold && addLastSection) {
//        if (printStats) {
//            PRINT(">[%d]->%d\n", z + 1, sectionsAboveThreshold + 1);
//            PRINT("esoteric:\"%s\" length:%d\n", esotericLine.c_str(), esotericLine.length());
//            PRINT("esoYes:%d:\"%s\" length:%d\n", yesCount, yesNumber.c_str(), yesNumber.size());
//            PRINT("esoNo:%d:\"%s\" length:%d\n", noCount, noNumber.c_str(), noNumber.size());
//            PRINTADDITIONALS(">[%d]->%d\n", z + 1, sectionsAboveThreshold + 1);
//            PRINTADDITIONALS("esoteric:\"%s\" length:%d\n", esotericLine.c_str(), esotericLine.length());
//            PRINTADDITIONALS("esoYes:%d:\"%s\" length:%d\n", yesCount, yesNumber.c_str(), yesNumber.size());
//            PRINTADDITIONALS("esoNo:%d:\"%s\" length:%d\n", noCount, noNumber.c_str(), noNumber.size());
//            PRINTADDITIONALS("DividersTillHere:%d\n", allSectionsTogether);
//        }
//        sectionsAboveThreshold++;
//    }
//
//    if (sectionsAboveThresh != nullptr)
//        *sectionsAboveThresh = sectionsAboveThreshold;
//    if (dividerCount != nullptr)
//        *dividerCount = allSectionsTogether;
//
//    if (printStats)
//        PRINTESOTERIC("[%09d]\t\tbiggest:%s\t\tlast:%s first:%s second:%s\n", (int)k,biggestEsotericLine.c_str(), lastEsotericLine.c_str(), firstEsotericLine.c_str(), secondEsotericLine.c_str());
//
//    if (printStats)
//        PRINT("%d:section length:%d secs:%d yes:%d no:%d difference:%d\n", z-1, sectionLength, allSectionsTogether, yesCount, noCount, yesCount - noCount);
//
//    if (printStats) {
//        PRINTADDITIONALS(">divideradvances:");
//        PRINT("divideradvances:");
//        for (auto& c : dividerProgress) {
//            PRINTADDITIONALS("%d,", c);
//            PRINT("%d,", c);
//        }
//        PRINTADDITIONALS(" ; count of dividers:%d\n", allSectionsTogether);
//        PRINT(" ; count of dividers:%d\n", allSectionsTogether);
//    }
//
//    return maxEsotericLineLength;
//}
//
//int querSum(int64_t i) {
//    int q = 0;
//    while (i > 0) {
//        q += i % 10;
//        i /= 10;
//    }
//    return q;
//}
//
//void printEsotericProperties(int aha, int sectionThreshold = 12) {
//    int sections = 0;
//    int esotericLineLength = printAllDividersOf(aha, true, sectionThreshold, &sections);
//    PRINT(">>\tnumber:%d maxEsoteric:%d querSum:%d querSumQuersum:%d sections:%d\n", aha, esotericLineLength, querSum(aha), querSum(querSum(aha)), sections);
//    PRINTADDITIONALS(">>\tnumber:%d maxEsoteric:%d querSum:%d querSumQuersum:%d sections:%d\n", aha, esotericLineLength, querSum(aha), querSum(querSum(aha)), sections);
//}
//
//int fide;
//std::string fileNameBase;
//
//#include <conio.h>
//#include <stdlib.h>
//#include <stdio.h>
//
//std::string input(const std::string& textline) {
//    printf("%s\n", textline.c_str());
//    std::string k;
//    char d;
//    while ((d = _getch())!=13) {
//        printf("%c", d);
//        k += d;
//    }
//    printf("\n");
//    return k;
//}
//
//int inputInt(const std::string& textline) {
//    return std::stoi(input(textline));
//}
//
//#include <set>
//#include <map>
//
//FILE *setupFile(const std::string& fileName, int esotericNumber, int sectionThreshold, int singleSection, int countLastSection, int rangeLo, int rangeHi) {
//    FILE *k = fopen(fileName.c_str(), "w");
//    fprintf(k, "esoteric number:%d sectionThreshold:%d singleSection:%d countLastSection:%d rangeLo:%d rangeHi:%d\n", esotericNumber, sectionThreshold, singleSection, countLastSection, rangeLo, rangeHi);
//    fclose(k);
//    return k;
//}
//
//#include <string>
//
//std::string getRatioString(int found, int lastFound) {
//    const int all = found < lastFound ? found : lastFound;
//    for (int i = all - 1; i > 1; --i) {
//        int k = found % i;
//        int p = lastFound % i;
//        if (k == 0 && p == 0) {
//            found /= i;
//            lastFound /= i;
//        }
//    }
//    return std::to_string(lastFound) + "|" + std::to_string(found);
//}
//
//bool isGanzeZahl(double v) {
//    return fabs(v - trunc(v)) < 0.0000000001;
//}
//
//std::string formatRatio(double ratio) {
//    int integerPart = trunc(ratio);
//    double fractionalPart = ratio - integerPart;
//    std::string k;
//    int r = 0;
//    while (integerPart > 0 || r < 4) {
//        std::string p;
//        p = (integerPart % 10) + '0';
//        k = p + k;
//        integerPart /= 10;
//        ++r;
//    }
//    std::string f = std::to_string(fractionalPart);
//    return k + f.substr(1,f.length()-1);
//}
//
//void printAllNumbersWithEsoteric(const int desiredEsotericNumber = 50, int rangeLo = 1, int rangeHi = 10000000) {
//    const bool o = printStats;
//    int esotericNumber, sectionThreshold, notifySectionsWithMoreThan, alsoAddLastSection, singleSection, fromThisNumber, toThisNumber;
//
//    alsoAddLastSection = inputInt("count last section, too(0 or 1):");
//    fromThisNumber = inputInt("from Number:(-1 means default):");
//    toThisNumber = inputInt("to Number:(-1 means default):");
//    if (fromThisNumber > 0) rangeLo = fromThisNumber;
//    if (toThisNumber > 0) rangeHi = toThisNumber;
//    addLastSection = alsoAddLastSection != 0;
//    esotericNumber = inputInt("esoteric to display (e.g. 50, -1 means another mode):");
//    if (esotericNumber < 0) {
//        singleSection = inputInt("only single section (0 or 1):");
//        sectionThreshold = inputInt("section length threshold (e.g. 12):");
//        if (singleSection == 0)
//            notifySectionsWithMoreThan = inputInt("numbers with more sections than(e.g. 3):");
//        else
//            notifySectionsWithMoreThan = inputInt("numbers with exactly (e.g. 4) sections:");
//    }
//    printStats = false;
//    std::set<int> foundSectionCounts;
//    std::map<int, int> countForFoundSectionId;
//    std::string numberFileName = fileNameBase + "_numbers.txt";
//    std::string additionalsFileName = fileNameBase + "_additionals.txt";
//    std::string processFileName = fileNameBase + "_process.txt";
//    std::string formulasFileName = fileNameBase + "_formulas.txt";
//    std::string esotericFileName = fileNameBase + "_esoteric.txt";
//    std::string advancersFileName = fileNameBase + "_advancers.txt";
//    numberOut = setupFile(numberFileName, esotericNumber, sectionThreshold, singleSection, addLastSection, rangeLo, rangeHi);
//    additionalsOut = setupFile(additionalsFileName, esotericNumber, sectionThreshold, singleSection, addLastSection, rangeLo, rangeHi);
//    processOut = setupFile(processFileName, esotericNumber, sectionThreshold, singleSection, addLastSection, rangeLo, rangeHi);
//    formulasOut = setupFile(formulasFileName, esotericNumber, sectionThreshold, singleSection, addLastSection, rangeLo, rangeHi);
//    esotericOut = setupFile(esotericFileName, esotericNumber, sectionThreshold, singleSection, addLastSection, rangeLo, rangeHi);
//    advancersOut = setupFile(advancersFileName, esotericNumber, sectionThreshold, singleSection, addLastSection, rangeLo, rangeHi);
//    double lastFound = 1;
//    double dlastFound = 1;
//    std::map<int, int> lastFoundForSection, dlastFoundForSection;
//    for (int i = rangeLo; i <= rangeHi; ++i) {
//        numberOut = fopen(numberFileName.c_str(), "a+");
//        additionalsOut = fopen(additionalsFileName.c_str(), "a+");
//        processOut = fopen(processFileName.c_str(), "a+");
//        formulasOut = fopen(formulasFileName.c_str(), "a+");
//        esotericOut = fopen(esotericFileName.c_str(), "a+");
//        advancersOut = fopen(advancersFileName.c_str(), "a+");
//        int sections = 0;
//        int dividers = 0;
//        int esotericLineLength = printAllDividersOf(i, true, sectionThreshold, &sections, &dividers);
//        bool valid = false;
//        valid |= esotericNumber >= 0 ? esotericLineLength == esotericNumber : false;
//        valid |= sectionThreshold >= 0 ? sections > notifySectionsWithMoreThan : false;
//        if (singleSection != 0) {
//            valid = sections == notifySectionsWithMoreThan;
//        }
//        if (valid) {
//            foundSectionCounts.insert(sections);
//            countForFoundSectionId[sections]++;
//            printStats = true;
//            printEsotericProperties(i, sectionThreshold);
//            int di = dividers;
//            PRINTFORMULA("[%d]:[dividers:%d] -> *=\n", i, di);
//            {
//                double ratio = (double)lastFound / i;
//                PRINT("%s;%f[%s], %f[%s]\n", (isGanzeZahl(ratio) || isGanzeZahl(1.0 / ratio)) ? "*" : " ", ratio, getRatioString(i, lastFound).c_str(), 1.0 / ratio, getRatioString(lastFound, i).c_str());
//                PRINTFORMULA("%s;%f[%s], %f[%s]", (isGanzeZahl(ratio) || isGanzeZahl(1.0 / ratio)) ? "*" : " ", ratio, getRatioString(i, lastFound).c_str(), 1.0 / ratio, getRatioString(lastFound, i).c_str());
//                double lastFoundS = lastFoundForSection[sections];
//                if (lastFoundS == 0)
//                    lastFoundS = 1;
//                double ratioS = (double)lastFoundS / i;
//                PRINTFORMULA(" %ssection%d:before:%d;%f[%s], %f[%s]\n", (isGanzeZahl(ratioS) || isGanzeZahl(1.0 / ratioS)) ? "*" : " ", sections, (int)lastFoundS, ratioS, getRatioString(i, lastFoundS).c_str(), 1.0 / ratioS, getRatioString(lastFoundS, i).c_str());
//            }
//
//            {
//                std::string prefix = "";
//                if (dlastFoundForSection.find(sections) == dlastFoundForSection.end())
//                    prefix = "initial->";
//                double dratio = (double)dlastFound / di;
//                PRINT("%s;%f[%s], %f[%s]\n", (isGanzeZahl(dratio) || isGanzeZahl(1.0 / dratio)) ? "*" : " ", dratio, getRatioString(di, dlastFound).c_str(), 1.0 / dratio, getRatioString(dlastFound, di).c_str());
//                PRINTFORMULA("%s;%f[%s], %f[%s]", (isGanzeZahl(dratio) || isGanzeZahl(1.0 / dratio)) ? "*" : " ", dratio, getRatioString(di, dlastFound).c_str(), 1.0 / dratio, getRatioString(dlastFound, di).c_str());
//                double dlastFoundS = dlastFoundForSection[sections];
//                if (dlastFoundS == 0)
//                    dlastFoundS = 1;
//                double dratioS = (double)dlastFoundS / di;
//                PRINTFORMULA(" %ssection%d:before:%d;%f[%s], %f[%s]\n", (isGanzeZahl(dratioS) || isGanzeZahl(1.0 / dratioS)) ? "*" : " ", sections, (int)dlastFoundS, dratioS, getRatioString(di, dlastFoundS).c_str(), 1.0 / dratioS, getRatioString(dlastFoundS, di).c_str());
//                PRINTADVANCERS("%s%s%d:%04d:ratio:%s\t\t[%s]=%f\t\t[%s]=%f ; before:%d,after:%d\n", prefix.c_str(), (isGanzeZahl(dratioS) || isGanzeZahl(1.0 / dratioS)) ? "*" : " ", i, sections, formatRatio(1.0 / dratioS),getRatioString(di, dlastFoundS).c_str(), dratioS, getRatioString(dlastFoundS,di).c_str(), 1.0/dratioS, (int)dlastFoundS, di);
//            }
//
//            PRINTPROCESS("[%d]\n", i);
//            PRINT("sectionIds:");
//            PRINTPROCESS("found sectionIds:");
//            for (auto& c : foundSectionCounts) {
//                PRINT("%d,", c);
//                PRINTPROCESS("%d[%d],", c, countForFoundSectionId[c]);
//            }
//            PRINT("\n");
//            PRINTPROCESS("\n");
//            printStats = false;
//            lastFound = i;
//            dlastFound = di;
//            lastFoundForSection[sections] = lastFound;
//            dlastFoundForSection[sections] = dlastFound;
//        }
//        fclose(advancersOut);
//        fclose(formulasOut);
//        fclose(numberOut);
//        fclose(additionalsOut);
//        fclose(processOut);
//        fclose(esotericOut);
//    }
//    printStats = o;
//}
//
//#include <set>
//#include <chrono>
//
//int getCipherCount(double v) {
//    //v = trunc(v);
//    double r = 0;
//    while (v >= 1.0) {
//        v /= 10;
//        r++;
//    }
//    return r;
//}
//
//double clipToNDigitsFraction(double v, int digits = 10, bool rounding = false) {
//    v *= pow(10.0, digits);
//    v = trunc(v);
//    v /= pow(10.0, digits);
//    return v;
//}
//
//int getCipherCountInverse(double v) {
//    v = fmod(v, 1.0);
//    int r = 0;
//    if (v == 0)
//        return 0;
//    while (v < 1) {
//        v *= 10;
//        r++;
//    }
//    return r;
//}
//
//std::string formatRatio(double ratio, int ciphers) {
//    int integerPart = trunc(ratio);
//    double fractionalPart = ratio - integerPart;
//    std::string k;
//    int r = 0;
//    while (integerPart > 0 || r < ciphers) {
//        std::string p;
//        p = (integerPart % 10) + '0';
//        k = p + k;
//        integerPart /= 10;
//        ++r;
//    }
//    std::string f = std::to_string(fractionalPart);
//    return k + f.substr(1, f.length() - 1);
//}
//
//#include <stdint.h>
//
//FILE* checkyx;
//int yugaDNASize = 42+22+72;
//void decipherYuga(const std::string& caption, double baseYuga, int steps = 2000, bool forward = true, int letterCountFuture = 10, int letterCountPast = 10) {
//    forward = true;
//    letterCountFuture = yugaDNASize;
//    letterCountPast = yugaDNASize;
//    getCipherCount(baseYuga);
//    FILEPRINT(checkyx, "%s %s\n", caption.c_str(), forward ? "> shows actual program direction" : "> shows previous steps");
//    printf("%s %s\n", caption.c_str(), forward ? "> shows actual program direction" : "> shows previous steps");
//    double currentYuga = baseYuga;
//    int cipherStepsBefore = 0;
//    for (int i = 0; i < 10000; ++i) {
//        int ciphers = getCipherCount(currentYuga);
//        int prefix = letterCountFuture - ciphers;
//        if (prefix <= 0)
//            break;
//        currentYuga *= 2.5;
//        cipherStepsBefore++;
//    }
//    int64_t lastNachkomma = -1210;
//    //currentYuga /= 2.5;
//    for (int i = 0; i < 10000; ++i) {
//        std::string prefix2 = " ";
//        if (cipherStepsBefore-- == 0) {
//            prefix2 = "*";
//            currentYuga = baseYuga;
//        }
//
//        double c = currentYuga;// clipToNDigitsFraction(currentYuga, letterCountPast);
//        double vorkomma = trunc(c);
//        double k = c - vorkomma;
//        const double round = 0.5;
//        int64_t nachkomma = k * pow(10, letterCountPast) + round;
//        std::string prefix3 = " ";
//        if (nachkomma != 0 && lastNachkomma == 0)
//            prefix3 = "-";
//
//        int ciphers = getCipherCount(vorkomma);
//        int prefix = letterCountFuture - ciphers;
//        if ((!forward) && prefix < 0)
//            break;
//        if (vorkomma == 0 && nachkomma == 0)
//            break;
//        std::string formatString2 = prefix2 + prefix3 + "%03d:";
//        for (int j = 0; j < prefix - (ciphers == 0 ? 1 : 0); ++j)
//            formatString2 += "0";
//        formatString2 += "%." + std::to_string(letterCountPast) + "f\n";
//        FILEPRINT(checkyx, formatString2.c_str(), i, currentYuga);
//        printf(formatString2.c_str(), i, currentYuga);
//        if (forward)
//            currentYuga /= 2.5;
//        else
//            currentYuga *= 2.5;
//        lastNachkomma = nachkomma;
//    }
//}
//
//bool isInteger(const Number& b) {
//    Number c = cleanUpAll(b);
//    return !b.empty() && b[b.size()-1] == '.';
//}
//
//bool decipherYugaInverse = false;
//void decipherYuga2(const std::string& caption, double baseYuga, int steps = 2000, bool forward = true, int letterCountFuture = 10, int letterCountPast = 10) {
//    letterCountFuture = yugaDNASize;
//    letterCountPast = yugaDNASize;
//    Number currentYuga = toNumber(baseYuga);
//    double currentYuga2 = baseYuga;
//    int cipherStepsBefore = 0;
//    FILEPRINT(checkyx, "%s %s\n", caption.c_str(), forward ? "> shows actual program direction" : "> shows previous steps");
//    printf("%s %s\n", caption.c_str(), forward ? "> shows actual program direction" : "> shows previous steps");
//    for (int i = 0; i < 10000; ++i) {
//        int ciphers = dotPos(currentYuga);
//        int prefix = letterCountFuture - ciphers;
//        if (prefix <= 0)
//            break;
//        currentYuga = multiplyBy2Point5(currentYuga);
//        currentYuga2 *= 2.5;
//        cipherStepsBefore++;
//    }
//    std::string prefix3 = "";
//    for (int i = 0; i < 10000; ++i) {
//        std::string prefix2 = " ";
//        if (cipherStepsBefore-- == 0) {
//            prefix2 = "*";
//            currentYuga = toNumber(baseYuga);
//        }
//        currentYuga = cleanUpAll(currentYuga);
//        Number b = concatNumberTo(currentYuga, letterCountFuture, letterCountPast);
//        Number c = cleanUpAll(b);
//        if (c == "0.0")
//            break;
//        prefix3 = " ";
//        if (isInteger(currentYuga)) {
//            Number checkYuga = currentYuga;
//            DNumber d = ditherNumberFromRight(cleanUpAll(checkYuga), 1.0 / 1.5, 1000);
//            if ((!d.second.empty() && d.second[0] == '*') || d.second.empty())
//                prefix3 = "+";
//        }
//        FILEPRINT(checkyx, (prefix3 + prefix2 + "%03d:%s\n").c_str(),i, b.c_str());
//        printf((prefix3 + prefix2 + "%03d:%s\n").c_str(), i, b.c_str());
//        Number yugaBefore = currentYuga;
//        currentYuga = divideBy2Point5(currentYuga);
//    }
//}
////   for (int i = 0; i < 100; ++i) {
////       printNumber2(currentYuga);
////       currentYuga = divideBy2Point5(currentYuga);
////   for (int i = 0; i < 10000; ++i) {
////       std::string prefix2 = " ";
////       if (cipherStepsBefore-- == 0) {
////           prefix2 = "*";
////           currentYuga = baseYuga;
////       }
////   }
////   int debug = 1;
////
//
//void unitTest() {
//    srand(1);
//    for (int k = 0; k < 1000; ++k) {
//        double a1 = rand() / 100.0;
//        double a2 = rand() / 100.0;
//        Number n1 = toNumber(a1);
//        Number n2 = toNumber(a2);
//        //double ra1 = fromNumber(n1);
//        //double ra2 = fromNumber(n2);
//        double v1 = a1 + a2;
//        Number nv1 = add(n1, n2);
//    }
//}
//
//
//int iterations = 0;
//int main(int argc, char**argv)
//{
//    testi();
//    //DNumber k7 = ditherNumberFromRight("707", 1.0/1.5, 1000);
//    double k9 = 1.0 / 1.5;
//    DNumber k7 = ditherNumberFromRight("9999", k9, 1000);
//    double k10 = 1.5 / 1.0;
//    DNumber k8 = ditherNumberFromRight("9999", k10, 1000, true);
//    printDNumber(k7);
//
//    Number r6 = divideByTwo("101.11");
//
//    unitTest();
//    //Number r5 = add("123", "99");
//    //Number d2 = divideByFour("100");
//    //Number d3 = multiplyWithTen("112.030");
//    Number d4 = add("112.030", "1.95");
//    Number d5 = add("111223423.030", "2345234231.95");
//
//    printf("%.50f\n", 273 / 2.5);
//
//    double a = 27074941.324810438;
//    double b = 7391920.7250077892;
//    //double loga = log(a);
//    //double logb = log(b);
//    //double logma = log(-a);
//    //double logmb = log(-b);
//    //double lalb = loga * logb;
//    a = 9;
//    b = 3;
//    double log2a = log2(a);
//    double log2b = log2(b);
//    double log2ma = log2(-a);
//    double log2mb = log2(-b);
//    //double lalb2 = log2a * log2b;
//    //double k9 = lalb2;
//
//    double stepsOfMoonAnglePhaseShiftToReach1 = 1.0 / (fmod(moonAngleToEarthAfterOneSolarYear,360));
//    double k = stepsOfMoonAnglePhaseShiftToReach1 * 360;
//    double solarDaysForTheMoonToReach360 = k * 366.0;
//    double p = (solarDaysForTheMoonToReach360/366.0) * fmod(moonAngleToEarthAfterOneSolarYear, 360);
//    // solarDaysForTheMoonToReach360 means that the moon is ?fixed? again to the same position around the earth almost every 922 days
//    // around 190st day of 2nd solar year
//    double moonSolarYear = trunc(solarDaysForTheMoonToReach360 / 366.0); // not gregorian
//    double moonSolarDay = fmod(solarDaysForTheMoonToReach360, 366.0); // not gregorian
//    // check when the solar day would be a full day not a fractional day now
//    //iterations = 0;
//    // this function has many 0 stellen but only converges to an epsilon of 0.0001 in certain cases e.g. (0.1, 10000.0)
//    const double round = 0.5;
//    double iterations_22 = errorReduce(0.1, 10000.0, 0.0, [=](double v, double thresh) -> double {
//        iterations++;
//        return fabs(solarDaysForTheMoonToReach360 * v - trunc(solarDaysForTheMoonToReach360 * v)-round);
//        });
//    //iterations = 0;
//    // this function has many 0 stellen but only converges to an epsilon of 0.0001 in certain cases e.g. (1.0, 100000.0)
//    double iterations_72 = errorReduce(1.0, 100000.0, 0.0, [=](double v, double thresh) -> double {
//        iterations++;
//        return fabs(solarDaysForTheMoonToReach360 * v - trunc(solarDaysForTheMoonToReach360 * v)-round);
//        });
//    // dunno whatever numbers could be good, could be checked here?
//    // these numbers below represent the complete days the sun would need that the configuration of sun and moon is the same on a direct daystart
//    double check_22 = solarDaysForTheMoonToReach360 * iterations_22;
//    double check_72 = solarDaysForTheMoonToReach360 * iterations_72;
//    const double yugaPerfection4 = solarDaysForTheMoonToReach360 * 1728000;
//    const double yugaPerfection3 = solarDaysForTheMoonToReach360 * (1728000 * 3 / 4);
//    const double yugaPerfection2 = solarDaysForTheMoonToReach360 * (1728000 * 2 / 4);
//    const double yugaPerfection1 = solarDaysForTheMoonToReach360 * (1728000 * 1 / 4);
//    const double completeYugCycle = yugaPerfection1 + yugaPerfection2 + yugaPerfection3 + yugaPerfection4;
//    // period of nearly 42 days and 1563 solar years (it's really bad, because these numbers are very rough estimates)
//    double check22_moonSolarYear = trunc(check_22 / 366.0); // not gregorian
//    double check22_moonSolarDay = fmod(check_22, 366.0); // not gregorian
//    // period of nearly 277 days and 88605 solar year (it's really bad, because these numbers are very rough estimates)
//    double check72_moonSolarYear = trunc(check_72 / 366.0); // not gregorian
//    double check72_moonSolarDay = fmod(check_72, 366.0); // not gregorian
//
//
//    //double checkthis = check_72 / check_22;
//
//    iterations = 0;
//    double iterations_solarYear = errorReduce(0.1, 10000.0, 0.0, [=](double v, double thresh) -> double {
//        iterations++;
//        if (iterations > 20000000)
//            return 0;
//        const double maybe = solarDaysForTheMoonToReach360 / 366.0;
//        return fabs(maybe * v - trunc(maybe * v)-round);
//        });
//    // number represents the number of solar years it would need that the moon is on the same position to earth on a year start like on the absolute ?start? (the angle offsets of the moon,earth,sun at start would be needed for a meaningful number e.g. age)
//    double r5t = solarDaysForTheMoonToReach360 * iterations_solarYear;
//
//    //checkyx = NULL;
//    //decipherYuga("moon constant1:", r5t, 60, true);
//
//    std::vector<double> kommastellentobeshown = {
//        10,
//        20,
//        22,
//        26,
//        28,
//        42,
//        50,
//        55,
//        60,
//        66,
//        72,
//        77,
//        90,
//        100,
//        216,
//        42+22+72,
//        42+72,
//        22+72,
//        42+22
//    };
//    std::string name = "corrected/second/checknr3_";
//
//    for (auto k : kommastellentobeshown) {
//        yugaDNASize = k;
//        std::string fileName = name + std::to_string(yugaDNASize) + ".txt";
//        checkyx = fopen(fileName.c_str(), "w");
//        const int showSteps = 1000;
//        //decipherYuga("real constellations in current sun system (with current data, and guessedperfect correction constants) backward:", assumedSunConstant, showSteps, false);
//        decipherYuga2("real constellations in current sun system (with current data, and guessedperfect correction constants) forward:", assumedSunConstant, showSteps, true);
//        decipherYuga2("moon constant represents when the moon is on the same angle around earth and sun", r5t, showSteps);
//        for (int i = 1; i <= 4; ++i) {
//            //decipherYuga("perfection yuga as defined backward:" + std::to_string(i), 1728000 * i / 4, showSteps, false);
//            decipherYuga2("perfection yuga as defined forward:" + std::to_string(i), 1728000 * i / 4, showSteps, true);
//        }
//        decipherYuga2("3080 just for a test", 3080, showSteps);
//        decipherYuga2("273 just for a test", 273, showSteps);
//        fclose(checkyx);
//    }
//
//
//    double stepsOfDayShiftPerYearToReach1 = 1.0 / fmod(solarDaysForTheMoonToReach360/366.0,1.0);
//    double oneCompleteSunMoonCycleEveryNthSolarYears = 1.0 / stepsOfDayShiftPerYearToReach1;
//    double daysForThat = oneCompleteSunMoonCycleEveryNthSolarYears * 366.0;
//    //double gregorian = daysForThat * 365 / 366;
//
//    double moonSunRatio1 = assumedSunConstant / r5t;
//    double moonSunRatio2 = r5t / assumedSunConstant;
//    // could be 0.01
//    double meaningConstant1 = assumedSunConstant / strange_estimated_NumberOfSunYears_NeededForAlignmentOf_SunMoonAndEarth_OnAYearStart;
//    // could be 2
//    double meaningConstant2 = strange_estimated22_NumberOfSunYears_NeededForAlignmentOf_SunMoonAndEarth_OnADayStart / assumedSunConstant;
//    // could be 0.001
//    double meaningConstant3 = strange_estimated72_NumberOfSunYears_NeededForAlignmentOf_SunMoonAndEarth_OnADayStart / assumedSunConstant;
//
//    srand(std::chrono::high_resolution_clock::now().time_since_epoch().count());
//    fide = rand();
//    fileNameBase = "checkys/esoteric" + std::to_string(fide);
//    checkySetup a1;
//    double standardPi = 3.141592653589793238;
//    double keralaPi = 104348.0 / 33215.0;
//    double keralaPiRatio = standardPi / keralaPi;
//    double pi = standardPi;
//    double earthRotationsPerYear = 366;
//    double perfectionHinduYuga = 1296000;
//    double kaliYuga = 1296000/4;
//    double secondYuga = perfectionHinduYuga - kaliYuga * 1;
//    double thirdYuga = perfectionHinduYuga - kaliYuga * 2;
//    double sunConstant = a1.assumedSunConstant;
//    //----------------------------------
//    //double sunConstant = reCalcCurrentSunConstant();
//    double perfectionHinduYugaEarthRotations = perfectionHinduYuga * earthRotationsPerYear;
//    double sunConstantWithYugaRotations = sunConstant / perfectionHinduYugaEarthRotations;
//
//    //printAllDividersOf(299880, true);
//
//    printAllNumbersWithEsoteric(50);
//
//    //printStats = false;
//    //int esotericId = 0;
//    //int lastEsoteric = 0;
//    //std::vector<std::set<int>> sectionQuerSums(100);
//    //for (int i = 0; i < 10000000; ++i) {
//    //    out = fopen(("checkys/check" + std::to_string(fide)).c_str(), "a+");
//    //    int niceNumber = i;
//    //    int esotericLineLength = 0;
//    //    int sections = 0;
//    //    esotericLineLength = printAllDividersOf(niceNumber, true, 11, &sections);
//    //    sectionQuerSums[sections].insert(querSum(niceNumber));
//    //    if (sections > 3) {
//    //        PRINT(">>>>>>>>>>>>>>>>>>>>>>?WHATSSS DISSS?: (next number)");
//    //        // ?special number?
//    //        printStats = true;
//    //        printEsotericProperties(niceNumber);
//    //        printStats = false;
//    //    }
//    //    if (esotericLineLength > lastEsoteric) {
//    //        printStats = true;
//    //        printEsotericProperties(niceNumber);
//    //        printStats = false;
//    //        lastEsoteric = esotericLineLength;
//    //        esotericId++;
//    //    }
//    //
//    //    if ((i % 100000) == 0) {
//    //        for (int j = 0; j < sectionQuerSums.size(); ++j) {
//    //            if (!sectionQuerSums[j].empty()) {
//    //                PRINT("quersums for sectioncount[%d]:", j);
//    //                for (auto& r : sectionQuerSums[j]) {
//    //                    PRINT(",%d", r);
//    //                }
//    //                PRINT("\n");
//    //            }
//    //            ++i;
//    //        }
//    //    }
//    //    fclose(out);
//    //}
//    //out = fopen(("checkys/check" + std::to_string(fide)).c_str(), "a+");
//    //PRINT("perfection:\n");
//    //printEsotericProperties(perfectionHinduYuga);
//    //PRINT("second:\n");
//    //printEsotericProperties(secondYuga);
//    //PRINT("third:\n");
//    //printEsotericProperties(thirdYuga);
//    //PRINT("kali:\n");
//    //printEsotericProperties(kaliYuga);
//    ////printf("31416:\n");
//    ////printEsotericProperties(31416);
//    //fclose(out);
//
//    double correctNorm = errorReduce(0.5, 1.5, 55, [=](double v, double threshold) -> double {
//        double r = standardPi / (sunConstantWithYugaRotations * v);
//        return (r - threshold);
//        });
//    double correctWhatever = perfectionHinduYuga * correctNorm;
//    int debug = 1;
//}
//
//
///*
//0,
//1,
//2,
//3,
//4,
//5,
//6,
//7,
//8,
//9,
//10,
//11,
//12,
//13,
//14,
//15,
//16,
//17,
//18,
//19,
//20,
//21,
//22,
//23,
//24,
//25,
//26,
//27,
//28,
//29,
//30,
//31,
//32,
//33,
//34,
//35,
//36,
//37,
//38,
//39,
//40,
//41,
//42,
//43,
//44,
//45,
//46,
//47,
//48,
//49,
//50,
//51,
//52,
//53
//; 54 values
//
//3,
//6,
//7,
//8,
//9,
//10,
//11,
//12,
//13,
//14,
//15,
//16,
//17,
//18,
//19,
//20,
//21,
//22,
//23,
//24,
//25,
//26,
//27,
//28,
//29,
//30,
//31,
//32,
//33,
//34,
//35,
//36,
//37,
//38,
//39,
//40,
//41,
//42,
//44,
//45,
//46,
//48,
//51
//; 43 values ; 11 values difference
//
//3,
//9,
//10,
//12,
//15,
//18,
//21,
//23,
//24,
//27,
//30,
//33,
//36,
//39,
//42,
//45
//; 16 values ; 27 values difference
//*/

//const double r3 = calculateConstantForDesiredDNumberSingleDigitPeriod('3', "10", 60, 1);
//const double r7 = calculateConstantForDesiredDNumberSingleDigitPeriod('7', "10", 60, 1);
//const double r5 = calculateConstantForDesiredDNumberInitialStringPeriod(kr, "10", -1, 1);
//{
//    int precision = 5000;
//    const double rb3_1 = calculateConstantForDesiredDNumberSingleDigitPeriod('3', "10", precision, 0, true);
//    const double rb3_1b = calculateConstantForDesiredDNumberSingleDigitPeriod('3', "10", precision, 1, false);
//    precision /= 10;
//    const double rb3_2 = calculateConstantForDesiredDNumberSingleDigitPeriod('3', "10", precision, 0, true);
//    const double rb3_2b = calculateConstantForDesiredDNumberSingleDigitPeriod('3', "10", precision, 1, false);
//    precision /= 10;
//    const double rb3_3 = calculateConstantForDesiredDNumberSingleDigitPeriod('3', "10", precision, 0, true);
//    const double rb3_3b = calculateConstantForDesiredDNumberSingleDigitPeriod('3', "10", precision, 1, false);
//}
