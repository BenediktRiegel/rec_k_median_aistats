#include "fixedDouble.h"
#include <math.h>

using namespace std;


pair<fixedDouble, fixedDouble> getToSameP(fixedDouble a, fixedDouble b) {
    int diffp = a.precision - b.precision;
    if(diffp < 0) {
        diffp *= -1;
        a.value *= pow(10, diffp);
        a.precision = b.precision;
    }
    else {
        b.value *= pow(10, diffp);
        b.precision = a.precision;
    }
    return pair<fixedDouble, fixedDouble> (a, b);
}


fixedDouble reducePrecision(fixedDouble a) {
    while(a.value % 10 == 0 && a.precision != 0) {
        a.value = a.value / 10;
        --a.precision;
    }
    return a;
}


fixedDouble fixedDouble::operator+(fixedDouble b) const {
    long av = this->value;
    int ap = this->precision;
    long bv = b.value;
    int bp = b.precision;
    if (bp < ap) {
        av = b.value;
        ap = b.precision;
        bv = this->value;
        bp = this->precision;
    }
    int difp = bp - ap;
    av = av * pow(10, difp);
    //now ap = bp, but it's not needed
    fixedDouble newDouble;
    newDouble.value = av + bv;
    newDouble.precision = bp;
    return reducePrecision(newDouble);
}


bool fixedDouble::operator<(fixedDouble b) const {
    pair<fixedDouble, fixedDouble> sameP = getToSameP(*this, b);
    return sameP.first.value < sameP.second.value;
}

bool fixedDouble::operator>(fixedDouble b) const {
    pair<fixedDouble, fixedDouble> sameP = getToSameP(*this, b);
    return sameP.first.value > sameP.second.value;
}

fixedDouble fixedDouble::operator/(fixedDouble b) {
    // a / b = (av / ap) / (bv / bp)
    // if ap == bp, then a / b = av / bv = c = cv / cp
    // => cv == k * av and vp == k * cp
    //wrong, since cp = 10^y with y in Z

    //instead get to same p
    fixedDouble a = *this;
    fixedDouble result;
    if (b.precision > this->precision) {
        pair<fixedDouble, fixedDouble> sameP = getToSameP(*this, b);
        b = sameP.second;
        a = sameP.first;
    }
    result.value = (a.value * pow(10, 2)) / b.value;
    result.precision = 2 + a.precision - b.precision;
    return reducePrecision(result);
}

fixedDouble fixedDouble::operator*(fixedDouble b) {
    fixedDouble result;
    result.value = (*this).value * b.value;
    result.precision = (*this).precision + b.precision;
    return reducePrecision(result);
}

bool fixedDouble::operator==(double b) const {
    return this->operator==(stofixedd(to_string(b)));
}

fixedDouble::fixedDouble(string value) {
    string myvalue;
    this->precision = 0;
    for (int i = 0; i < value.size(); ++i) {
        if (value.at(i) == '.') {
            this->precision = value.size() - 1 - i;
        } else {
            myvalue.push_back(value.at(i));
        }
    }
    this->value = stol(myvalue);
    *this = reducePrecision(*this);
}

fixedDouble stofixedd(const string& inString) {
    return fixedDouble(inString);
}

string to_string(fixedDouble a) {
    return "v: " + to_string(a.value) + ", p: " + to_string(a.precision)
            + " = " + to_string(a.getDouble());
}

string stringValue(fixedDouble a) {
    string result;
    string value = to_string(a.value);
    if (a.precision == 0) {
        return value;
    }
    int numbeforpoint = value.size() - a.precision;
    if (numbeforpoint <= 0) {
        result = "0.";
        numbeforpoint *= -1;
        for (int i = 0; i < numbeforpoint; ++i) {
            result.push_back('0');
        }
        result += value;
    }
    else {
        for (int i = 0; i < numbeforpoint; ++i) {
            result.push_back(value[i]);
        }
        result.push_back('.');
        for (int i = numbeforpoint; i < value.size(); ++i) {
            result.push_back(value[i]);
        }
    }
    return result;
}
