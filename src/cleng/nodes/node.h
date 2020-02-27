#pragma once

#include "../matrix.h"

class Node {
public:
    virtual Point point() const = 0;

    virtual void shift(const Point &shift) = 0;

    virtual void shift(const Matrix<Real> &matrix) = 0;

    virtual void pushSystemPoints(std::map<int, Point> &pointsById) const = 0;

    virtual bool operator<(const Node &other) const {return point() < other.point();}

    virtual std::string to_string() const = 0;

    virtual bool inSubBoxRange(const Point &subBoxRange) const = 0;

    virtual Point _returnSystemPoint() const = 0;

    virtual bool isIdInside(const int &ID) const = 0;

    virtual ~Node() = default;
};
