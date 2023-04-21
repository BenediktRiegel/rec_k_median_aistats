#ifndef BAAlgo_FIXEDDOUBLE_H
#define BAAlgo_FIXEDDOUBLE_H

#include <math.h>
#include <string>
using namespace std;


class fixedDouble {
public:
    long value;
    int precision;

    fixedDouble(){

    }

    fixedDouble(long myvalue, int myprecision){
        this->value = myvalue;
        this->precision = myprecision;
    }

    fixedDouble(string value);

    fixedDouble(double value) {
        *this = fixedDouble(to_string(value));
    }

    double getDouble() const {
        return (double) (value / pow(10, precision));
    }


    /*
    fixedDouble operator + (fixedDouble b) {
        int c = this->value + b.value;
        return fixedDouble (c, 1);
    }*/

    fixedDouble operator+(fixedDouble b) const;

    void operator+=(fixedDouble b) {
        (*this) = (*this) + b;
    }

    fixedDouble operator-(fixedDouble b) const {
        b.value *= -1;
        return (*this) + b;
    }

    void operator-=(fixedDouble b) {
        b.value *= -1;
        (*this) = (*this) + b;
    }

    fixedDouble operator/(fixedDouble b);

    fixedDouble operator*(fixedDouble b);

    bool operator==(fixedDouble b) const {
        return this->getDouble() == b.getDouble();
    }

    bool operator==(double b) const;

    bool operator==(const int b) const {
        return this->getDouble() == b;
    }


    bool operator!=(fixedDouble b) const {
        return !((*this) == b);
    }

    bool operator!=(double b) const {
        return !((*this) == b);
    }


    bool operator<(fixedDouble b) const;

    bool operator>(fixedDouble b) const;

    bool operator>=(fixedDouble b) const{
        return (*this).operator>(b) || (*this).operator==(b);
    }

    bool operator<=(fixedDouble b) const{
        return (*this).operator<(b) || (*this).operator==(b);
    }
};

fixedDouble stofixedd(const string& inString);

string to_string(fixedDouble a);

string stringValue(fixedDouble a);

#endif //BAAlgo_FIXEDDOUBLE_H
