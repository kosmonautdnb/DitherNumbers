#include "string_math1.hpp"
#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <map>

int dotPos(const Number& n) {
    for (int i = 0; i < n.size(); ++i) {
        if (n[i] == '.')
            return i;
    }
    return n.size();
}

Number cleanUp(const Number& n) {
    Number b;
    bool trailing = true;
    for (int i = 0; i < n.size(); ++i) {
        if (n[i] != '0')
            trailing = false;
        if (!trailing)
            b += n[i];
    }
    if (!b.empty() && b[0] == '.')
        b = "0" + b;
    return b;
}

Number cleanUpAll(const Number& n) {
    Number b = cleanUp(n);
    if (dotPos(b) != b.size()) {
        for (int i = b.size() - 1; i >= 0; --i) {
            if (b[i] != '0') {
                b = b.substr(0, i + 1);
                break;
            }
        }
    }
    if (b == "0.")
        b = "0.0";
    return b;
}

Number fixAll(const Number& n) {
    Number b = cleanUpAll(n);
    if (!b.empty() && b[b.size() - 1] == '.')
        b += "0";
    if (dotPos(b) == b.size())
        b += ".0";
    return b;
}

Number concatNumberTo(const Number& n, int vorkomma, int nachkomma) {
    Number b = cleanUpAll(n);
    int dot = dotPos(b);
    int neededBehind = nachkomma - (b.size() - dot) + 1;
    int neededBefore = vorkomma - dot;
    Number k;
    for (int j = 0; j < neededBefore; ++j)
        k += "0";
    k += b;
    for (int j = 0; j < neededBehind; ++j)
        k += "0";
    if (neededBehind < 0) {
        k = k.substr(0, k.size() + neededBehind);
    }
    return k;
}
#include <set>

Number toNumber(double n) {
    int64_t k = (int64_t)floor(n);
    Number r = std::to_string(k);
    Number fractional = std::to_string(n - k);
    if (fractional.size() > 2 && fractional[1] == '.')
        fractional = fractional.substr(2, fractional.size() - 2);
    r += ".";
    r += fractional;
    //printNumber2(r);
    return r;
}

double fromNumber(Number n) {
    n = cleanUpAll(n);
    int64_t vorkomma = 0;
    double nachkomma = 0;
    int i = 0;
    for (; i < n.size(); ++i) {
        char k = n[i];
        if (k == '.')
            break;
        vorkomma *= 10;
        vorkomma += k - '0';
    }
    double divisor = 1.0;
    for (int j = i; j < n.size(); ++j) {
        char k = n[j];
        if (k != '.') {
            nachkomma *= 10;
            nachkomma += k - '0';
            divisor *= 10;
        }
    }
    return vorkomma + (nachkomma / divisor);
}

double cleanDouble(double n) {
    return fromNumber(toNumber(n));
}


Number makeWhole(const Number& n) {
    Number r;
    int prepend = 0;
    for (int i = 0; i < n.size(); ++i) {
        if (n[i] != '.')
            r += n[i];
        else
            prepend = i;
    }
    for (int i = 0; i < prepend; ++i) {
        r = "0" + r;
    }
    return r;
}

int makeSameBase(Number& a, Number& b) {
    int da = a.size() - dotPos(a); // the longer the more digits after comma
    int db = b.size() - dotPos(b);
    int dotLength = da;
    if (da > db) {
        for (int i = db; i < da; ++i) {
            b = b + "0";
        }
    } else {
        for (int i = da; i < db; ++i) {
            a = a + "0";
        }
        dotLength = db;
    }
    if (a.size() > b.size()) {
        for (int i = b.size(); i < a.size(); ++i) {
            b = "0" + b;
        }
    }
    else {
        for (int i = a.size(); i < b.size(); ++i) {
            a = "0" + a;
        }
    }
    return dotLength;
}

Number add(const Number& n1, const Number& n2) {
    Number a1 = n1;
    Number a2 = n2;
    int dotLen = makeSameBase(a1, a2);

    Number result;
    int carryCount = 0;
    int pos1 = a1.size() - 1;
    int pos2 = a2.size() - 1;
    while (pos1 >= 0 || pos2 >= 0 || carryCount != 0) {
        int c1 = 0;
        int c2 = 0;
        bool dot1 = false;
        bool dot2 = false;
        if (pos1 == a1.size() - dotLen) {
            c1 = 0;
            dot1 = true;
        }
        else {
            if (pos1 >= 0)
                c1 = a1[pos1] - '0';
        }

        if (pos2 == a2.size() - dotLen) {
            c2 = 0;
            dot2 = true;
        }
        else {
            if (pos2 >= 0)
                c2 = a2[pos2] - '0';
        }
        if ((!dot1) && (!dot2)) {
            int c = c1 + c2 + carryCount;
            carryCount = c / 10;
            c -= carryCount * 10;
            char d = c + '0';
            result = d + result;
        }
        else {
            result = '.' + result;
        }
        pos1--;
        pos2--;
    }
    return cleanUp(result);
}

Number multiplyWithTen(const Number& number) {
    Number r;
    int l = number.length();
    bool found = false;
    for (int k = 0; k < l; ++k) {
        if ((k - 2 >= 0) && number[k-2] == '.') {
            r += ".";
            found = true;
        }
        if (number[k] != '.')
            r += number[k];
    }
    if (!found)
        r += ".0";
    return fixAll(r);
}

Number divideByTen(const Number& number) {
    Number result;
    for (int i = 0; i < number.size(); ++i) {
        if (i + 1 < number.size() && number[i + 1] == '.') {
            if (result.empty())
                result += "0";
            result += ".";
        }
        if (number[i] != '.')
            result += number[i];
    }
    return fixAll(result);
}

Number multiplyWithTwo(const Number& number) {
    return add(number,number);
}

Number multiplyWithThree(const Number& number) {
    return add(number, multiplyWithTwo(number));
}

Number multiplyWithFour(const Number& number) {
    return add(multiplyWithTwo(number), multiplyWithTwo(number));
}

Number multiplyWithFive(const Number& number) {
    return add(number, multiplyWithFour(number));
}

Number multiplyWithSix(const Number& number) {
    return multiplyWithTwo(multiplyWithThree(number));
}

Number multiplyWithSeven(const Number& number) {
    return add(multiplyWithThree(number), multiplyWithFour(number));
}

Number multiplyWithEight(const Number& number) {
    return multiplyWithTwo(multiplyWithFour(number));
}

Number multiplyWithNine(const Number& number) {
    return add(multiplyWithFive(number), multiplyWithFour(number));
}

Number divideByTwo(const Number& number) {
    Number result = number;
    for (int i = number.size() - 1; i >= 0; --i) {
        if (number[i] != '.') {
            result[i] = ((number[i] - '0') / 2) + '0';
        }
    }
    std::string numberk = number;
    std::string point5 = "0.0";
    makeSameBase(point5, result);
    point5 += "5";
    makeSameBase(point5, result);
    for (int i = numberk.size() - 1; i >= 0; --i) {
        if (numberk[i] != '.') {
            bool addPoint5 = ((numberk[i] - '0') % 2);
            if (addPoint5) {
                result = add(result, point5);
            }
            point5 = multiplyWithTen(point5);
        }
    }
    return cleanUp(result);
}

Number divideByLessThan2Test(const Number& number, double divider) {
    Number result = number;
    for (int i = number.size() - 1; i >= 0; --i) {
        if (number[i] != '.') {
            result[i] = (trunc((number[i] - '0') / divider)) + '0';
        }
    }
    std::string numberk = number;
    std::string point5 = "0.0";
    makeSameBase(point5, result);
    Number r = toNumber(1.0 / divider);
    point5 += r.substr(2, r.size() - 2); // todo: place propper "fractional" function here
    makeSameBase(point5, result);
    Number val = point5;
    Number muls[10];
    muls[0] = "0.0";
    muls[1] = val;
    muls[2] = multiplyWithTwo(val);
    muls[3] = multiplyWithThree(val);
    muls[4] = multiplyWithFour(val);
    muls[5] = multiplyWithFive(val);
    muls[6] = multiplyWithSix(val);
    muls[7] = multiplyWithSeven(val);
    muls[8] = multiplyWithEight(val);
    muls[9] = multiplyWithNine(val);

    
    for (int i = numberk.size() - 1; i >= 0; --i) {
        if (numberk[i] != '.') {
            double r7 = fmod((numberk[i] - '0'), divider);
            bool addPoint5 = r7 > 0.1;
            if (addPoint5) {
                result = add(result, point5);
            }
            point5 = multiplyWithTen(point5);
        }
    }
    return cleanUp(result);
}


Number divideByFour(const Number& number) {
    return divideByTwo(divideByTwo(number));
}


Number multiplyBy2Point5(const Number& number) {
    Number oneHalve = divideByTwo(number);
    return add(multiplyWithTwo(number),oneHalve);
}

Number divideBy2Point5(const Number& number) {
    Number result1 = multiplyWithFour(number);
    Number result3 = divideByTen(result1);
    return cleanUp(result3);
}

Number divideBy1Point5(const Number& number) {
    // 10 / 1.5 = 10 * 0.666 = 6.66666
    Number result1 = multiplyWithThree(number);
    Number result3 = divideByTen(result1);
    return cleanUp(result3);
}

int maxI = 0;
int countSameDigits(Number& a, int digit) { // < 0 returns max same digits
    int ret = 0;
    if (digit != -1) {
        for (int i = 0; i < a.size(); ++i) {
            ret += a[i] == digit ? 1 : 0;
        }
    }
    else {
        std::map<int, int> digitCounts;
        for (int i = 0; i < a.size(); ++i) {
            digitCounts[a[i]]++;
        }
        int maxCount = 0;
        for (auto k : digitCounts) {
            if (k.second > maxCount) {
                maxCount = k.second;
                maxI = k.first - '0';
            }
        }
        return maxCount;
    }
    return ret;
}

void printNumber2(const Number& n) {
    printf("%s\n", n.c_str());
}

static std::string intensity(const std::string& n) {
    const char digitsToIntensities[] = {
        ' ',
        '.',
        ',',
        ';',
        'i',
        'I',
        'O',
        'Q',
        '0',
        'B',
        '&',
        '$',
    };
    std::string r;
    for (int i = 0; i < n.size(); ++i) {
        char c = n[i];
        if (c >= '0' && c <= '9')
            c = digitsToIntensities[c - '0'];
        else if (c == '.') c = '|';
        r += c;
    }
    return r;
}

void printNumber3(const Number& n) {
    printf("%s\n", intensity(n).c_str());
}

void printDNumber(const DNumber& n) {
    printf("%s:%s\n", n.second.c_str(), n.first.c_str());
}

Number multiplyWith0to9(const Number &n, int number) {
    switch (number) {
        case 0:return toNumber(0);
        case 1:return n;
        case 2:return multiplyWithTwo(n);
        case 3:return multiplyWithThree(n);
        case 4:return multiplyWithFour(n);
        case 5:return multiplyWithFive(n);
        case 6:return multiplyWithSix(n);
        case 7:return multiplyWithSeven(n);
        case 8:return multiplyWithEight(n);
        case 9:return multiplyWithNine(n);
    }
    return toNumber(-1);
}

// the binary way but as tenary :)
Number multiply(Number a, Number b) {
    Number result;
    Number adder = fixAll(b);
    Number number = fixAll(a);
    makeSameBase(adder, number);
    int dp = number.length() - dotPos(number);
    number.erase(number.length() - dp,1);
    adder.erase(adder.length() - dp,1);
    dp--;
    adder = fixAll(adder);
    for (int i = 0; i < dp*2; ++i) {
        adder = divideByTen(adder);
    }
    for (int i = 0; i < number.size(); ++i) {
        int letterHere = number[number.size() - 1 - i];
        if (letterHere != '.') {
            int multiplyBy = letterHere - '0';
            result = add(result, multiplyWith0to9(adder, multiplyBy));
        }
        adder = multiplyWithTen(adder);
    }
    double res = fromNumber(a) * fromNumber(b);
    fixAll(result);
    return result;
}
