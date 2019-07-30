#pragma once

#include "ecf/ECF_base.h"
#include <ecf/floatingpoint/FloatingPoint.h>


class MyFloatingPoint : public FloatingPoint::FloatingPoint
{
public:
    MyFloatingPoint() {
        name_ = "MyFPGenotype";
    }


    MyFloatingPoint* copy() {
        MyFloatingPoint *newObject = new MyFloatingPoint(*this);
        return newObject;
    }


    std::vector<MutationOpP> getMutationOp();
    void registerParameters(StateP state);
    bool initialize(StateP state);
    
    const std::vector<double>& getMinVec() const { return _minVec; }
    const std::vector<double>& getMaxVec() const { return _maxVec; }

private:
    std::vector<double> _minVec;
    std::vector<double> _maxVec;

};
typedef boost::shared_ptr<MyFloatingPoint> MyFloatingPointP;
