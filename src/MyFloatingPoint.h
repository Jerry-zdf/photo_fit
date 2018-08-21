#ifndef MYFLOATINGPOINT_H
#define MYFLOATINGPOINT_H

#include "ecf/ECF_base.h"
#include <ecf/floatingpoint/FloatingPoint.h>
#include "GaussMutOp.h"


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
    std::vector<double> getMinVec() { return minVec; }
    std::vector<double> getMaxVec() { return maxVec; }

private:
    std::vector<double> minVec;
    std::vector<double> maxVec;

};
typedef boost::shared_ptr<MyFloatingPoint> MyFloatingPointP;

#endif // MYFLOATINGPOINT_H
