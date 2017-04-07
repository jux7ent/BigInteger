#include <iostream>
#include <cstring>
#include <stdexcept>

using namespace std;

typedef long long val_type;

class BigIntegerDivisionByZero : std::logic_error {
public:
    BigIntegerDivisionByZero();
};

BigIntegerDivisionByZero::BigIntegerDivisionByZero() : std::logic_error("Division by zero") {}

class BigIntegerOverflow : std::runtime_error {
public:
    BigIntegerOverflow();
};

BigIntegerOverflow::BigIntegerOverflow() : std::runtime_error("Overflow") {}


class BigInteger {

    class array {

    private:

        val_type *value_;
        size_t size_;
        size_t capacity_;

    public:

        array();

        array(val_type val);

        array(const array &copy);

        ~array();

        void push_back(val_type a);

        const val_type &operator[](size_t index) const;

        val_type &operator[](size_t index);

        array &operator=(const array &val);

        size_t size() const;

        bool empty() const;

        void insert(size_t pos, val_type val);

        void reserve(size_t newCapacity);

        void pop_back();
    };

private:

    array digits_;
    bool negative_;
    static const char base_len_ = 3;
    static const long base_ = 1e+3;
    static const size_t max_ = 20000 / base_len_ + 10;

    bool isNegative() const;

    bool isZero() const;

    void normalize();

    void add(const BigInteger &a, const bool substract = 0, const bool invert = 0);

    void divide(const BigInteger &a, const bool modulo = 0);

public:

    BigInteger();

    BigInteger(const BigInteger &copy);

    BigInteger(long long integer);

    explicit BigInteger(const char *string);

    size_t getSize() const;

    char getBaseLength() const;

    long long int getDigit(size_t pos) const;

    friend std::istream &operator>>(std::istream &is, BigInteger &a);

    const BigInteger operator-() const;

    BigInteger &operator++();

    BigInteger &operator--();

    const BigInteger operator++(int);

    const BigInteger operator--(int);

    BigInteger &operator+=(const BigInteger &a);

    BigInteger &operator-=(const BigInteger &a);

    BigInteger &operator*=(const BigInteger &a);

    BigInteger &operator=(const BigInteger &a);

    BigInteger &operator/=(const BigInteger &a);

    BigInteger &operator%=(const BigInteger &a);

    int compare(const BigInteger &a) const;

    BigInteger sqrt();

    BigInteger gcd(const BigInteger &q);

};

std::ostream &operator<<(std::ostream &os, const BigInteger &a);

const BigInteger operator+(const BigInteger &a, const BigInteger &b);

const BigInteger operator-(const BigInteger &a, const BigInteger &b);

const BigInteger operator*(const BigInteger &a, const BigInteger &b);

const BigInteger operator/(const BigInteger &a, const BigInteger &b);

const BigInteger operator%(const BigInteger &a, const BigInteger &b);

bool operator>(const BigInteger &a, const BigInteger &b);

bool operator<(const BigInteger &a, const BigInteger &b);

bool operator>=(const BigInteger &a, const BigInteger &b);

bool operator<=(const BigInteger &a, const BigInteger &b);

bool operator==(const BigInteger &a, const BigInteger &b);

bool operator!=(const BigInteger &a, const BigInteger &b);

BigInteger::array::array() : size_(0), capacity_(1) {
    value_ = new val_type[1];
}

BigInteger::array::array(val_type val) : size_(1), capacity_(1) {
    value_ = new val_type[1];
    value_[0] = val;
}

BigInteger::array::array(const array &copy) : size_(0), capacity_(1) {
    value_ = new val_type[1];
    *this = copy;
}

BigInteger::array::~array() {
    delete[] value_;
}

void BigInteger::array::reserve(size_t newCapacity) {
    if (newCapacity > capacity_) {
        val_type *temp = new val_type[newCapacity];

        for (size_t i = 0; i < size_; ++i) {
            temp[i] = value_[i];
        }

        delete[] value_;

        value_ = temp;
        capacity_ = newCapacity;
    }
}

void BigInteger::array::push_back(val_type a) {
    if (size_ == capacity_) {
        reserve(capacity_ * 2);
    }

    ++size_;

    if (size_ > BigInteger::max_)
        throw BigIntegerOverflow();

    value_[size_ - 1] = a;
}

size_t BigInteger::array::size() const {
    return size_;
}

bool BigInteger::array::empty() const {
    return size_ == 0;
}

const val_type &BigInteger::array::operator[](size_t index) const {
    return value_[index];
}

val_type &BigInteger::array::operator[](size_t index) {
    return value_[index];
}

void BigInteger::array::insert(size_t pos, val_type val) {
    if (size_ == capacity_) {
        reserve(capacity_ * 2);
    }


    for (size_t i = size_; i > pos; --i) {
        value_[i] = value_[i - 1];
    }
    value_[pos] = val;

    ++size_;
}

void BigInteger::array::pop_back() {
    if (size_ > 0) {
        --size_;
    }
}

BigInteger::array &BigInteger::array::operator=(const array &val) {
    size_ = 0;

    for (size_t i = 0; i < val.size(); ++i) {
        push_back(val[i]);
    }

    return *this;
}


BigInteger::BigInteger() : digits_({0}), negative_(0) {}

BigInteger::BigInteger(const BigInteger &copy) : digits_(copy.digits_), negative_(copy.negative_) {}

BigInteger::BigInteger(long long integer) : digits_({integer}), negative_(integer < 0) {
    if (digits_[0] < 0) {
        digits_[0] *= -1;
    }
    for (size_t i = 0; digits_[i] >= base_; ++i) {
        digits_.push_back(digits_[i] / base_);
        digits_[i] %= base_;
    }
}

BigInteger::BigInteger(const char *string) {
    negative_ = (string[0] == '-');

    long long digit = 0, multiplier = 1;
    for (size_t i = strlen(string) - 1; i >= negative_; --i) {
        digit += multiplier * (string[i] - '0');
        multiplier *= 10;
        if (multiplier >= base_) {
            digits_.push_back(digit);
            digit = 0;
            multiplier = 1;
        }
    }

    if (digit || !getSize())
        digits_.push_back(digit);
}

size_t BigInteger::getSize() const {
    return digits_.size();
}

char BigInteger::getBaseLength() const {
    return base_len_;
}

long long int BigInteger::getDigit(size_t pos) const {
    return pos < digits_.size() ? digits_[pos] : 0;
}

bool BigInteger::isNegative() const {
    return negative_;
}

bool BigInteger::isZero() const {
    return getSize() == 1 && digits_[0] == 0;
}


void BigInteger::normalize() {
    while (getSize() > 1 && digits_[getSize() - 1] == 0)
        digits_.pop_back();
}

void BigInteger::add(const BigInteger &a, const bool substract, const bool invert) {
    for (size_t i = 0; i < a.getSize(); ++i) {
        if (i == getSize())
            digits_.push_back(0);
        digits_[i] += a.digits_[i] * (2 * !substract - 1);
        digits_[i] *= (2 * !invert - 1);
    }

    for (size_t i = 0; i < getSize(); ++i) {
        if (digits_[i] < 0) {
            --digits_[i + 1];
            digits_[i] += base_;
        } else if (digits_[i] >= base_) {
            if (i + 1 == getSize())
                digits_.push_back(0);
            ++digits_[i + 1];
            digits_[i] -= base_;
        }
    }

    negative_ = negative_ != invert;

    normalize();
}

const BigInteger BigInteger::operator-() const {
    BigInteger a = *this;
    a.negative_ = !isZero() && !negative_;
    return a;
}

BigInteger &BigInteger::operator++() {
    *this += 1;
    return *this;
}

BigInteger &BigInteger::operator--() {
    *this += 1;
    return *this;
}

const BigInteger BigInteger::operator++(int) {
    BigInteger a;
    *this += 1;
    return *this;
}

const BigInteger BigInteger::operator--(int) {
    *this += 1;
    return *this;
}

BigInteger &BigInteger::operator+=(const BigInteger &a) {
    if (negative_ == a.negative_)
        add(a, 0, 0);
    else if ((*this >= -a) != negative_)
        add(a, 1, 0);
    else
        add(a, 1, 1);
    return *this;
}

BigInteger &BigInteger::operator-=(const BigInteger &a) {
    return *this += (-a);
}

BigInteger &BigInteger::operator*=(const BigInteger &a) {
    BigInteger c;
    for (size_t i = 0; i < getSize(); ++i) {
        long long carry = 0;
        for (size_t j = 0; j < a.getSize() || carry; ++j) {
            if (c.getSize() <= i + j)
                c.digits_.push_back(0);

            long long temp = c.digits_[i + j] + digits_[i] * a.getDigit(j) + carry;

            c.digits_[i + j] = temp % base_;
            carry = temp / base_;
        }
    }

    c.normalize();

    digits_ = c.digits_;
    negative_ = !isZero() && (negative_ != a.negative_);

    return *this;
}

BigInteger &BigInteger::operator/=(const BigInteger &a) {
    if (a.isNegative())
        divide(-a);
    else
        divide(a);
    negative_ = !isZero() && (negative_ != a.negative_);
    return *this;
}

BigInteger &BigInteger::operator%=(const BigInteger &a) {
    if (a.isNegative())
        divide(-a, 1);
    else
        divide(a, 1);
    return *this;
}

void BigInteger::divide(const BigInteger &divisor, const bool modulo) {
    if (divisor == 0) {
        throw BigIntegerDivisionByZero();
    }

    BigInteger quotient, remainder;

    for (size_t i = digits_.size() - 1; i <= digits_.size() - 1; --i) {
        remainder.digits_.insert(0, digits_[i]);
        remainder.normalize();

        quotient.digits_.insert(0, 0);
        quotient.normalize();

        while (divisor <= remainder) {
            if (divisor.getSize() == 1) {
                long long count = (remainder.getSize() != 1 ? base_ : remainder.getDigit(0)) / divisor.getDigit(0);
                remainder -= count * divisor.getDigit(0);
                quotient += count;
            } else {
                ++quotient;
                remainder -= divisor;
            }
        }
    }

    if (!modulo)
        digits_ = quotient.digits_;
    else
        digits_ = remainder.digits_;

    normalize();
}

BigInteger &BigInteger::operator=(const BigInteger &a) {
    digits_ = a.digits_;
    negative_ = a.negative_;
    return *this;
}

std::ostream &operator<<(std::ostream &os, const BigInteger &a) {
    if (a < 0)
        os << '-';
    os.fill('0');
    for (size_t i = a.getSize(); i > 0; --i) {
        if (i != a.getSize())
            os.width(a.getBaseLength());
        os << a.getDigit(i - 1);
    }
    return os;
}

std::istream &operator>>(std::istream &is, BigInteger &a) {
    char input[10010];
    is >> input;
    a = BigInteger(input);
    return is;
}

const BigInteger operator+(const BigInteger &a, const BigInteger &b) {
    BigInteger c = a;
    c += b;
    return c;
}

const BigInteger operator-(const BigInteger &a, const BigInteger &b) {
    BigInteger c = a;
    c -= b;
    return c;
}

const BigInteger operator*(const BigInteger &a, const BigInteger &b) {
    BigInteger c = a;
    c *= b;
    return c;
}

const BigInteger operator/(const BigInteger &a, const BigInteger &b) {
    BigInteger c = a;
    c /= b;
    return c;
}

const BigInteger operator%(const BigInteger &a, const BigInteger &b) {
    BigInteger c = a;
    c %= b;
    return c;
}

int BigInteger::compare(const BigInteger &a) const {
    if (negative_ != a.negative_)
        return a.negative_ ? 1 : -1;

    bool invert = negative_;

    if (getSize() != a.getSize())
        return (getSize() > a.getSize()) != invert ? 1 : -1;

    for (size_t i = getSize() - 1; i < getSize(); --i)
        if (digits_[i] != a.digits_[i])
            return (digits_[i] > a.digits_[i]) != invert ? 1 : -1;

    return 0;
}

BigInteger BigInteger::sqrt() {
    BigInteger left = 1, right = *this;
    BigInteger middle;

    if (*this != 1) {
        while (left + 1 < right) {
            middle = left + ((right - left) / 2);

            if (middle * middle <= *this) {
                left = middle;
            } else {
                right = middle;
            }
        }
    }

    return left;
}

BigInteger BigInteger::gcd(const BigInteger &a) {
    BigInteger q = *this, w = a;

    while (q != 0 && w != 0) {
        if (q > w) {
            q = q % w;
        }
        else {
            w = w % q;
        }
    }

    return q + w;
}

bool operator>(const BigInteger &a, const BigInteger &b) {
    return a.compare(b) == 1;
}

bool operator<(const BigInteger &a, const BigInteger &b) {
    return a.compare(b) == -1;
}

bool operator>=(const BigInteger &a, const BigInteger &b) {
    return a.compare(b) >= 0;
}

bool operator<=(const BigInteger &a, const BigInteger &b) {
    return a.compare(b) <= 0;
}

bool operator==(const BigInteger &a, const BigInteger &b) {
    return a.compare(b) == 0;
}

bool operator!=(const BigInteger &a, const BigInteger &b) {
    return a.compare(b) != 0;
}

int main() {}