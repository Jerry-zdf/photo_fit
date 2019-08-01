#pragma once

#include <vector>

#include "ecf/ECF_base.h"

#include <ecf/floatingpoint/FloatingPoint.h>

class MyFloatingPoint : public FloatingPoint::FloatingPoint {
   public:
    MyFloatingPoint();

    MyFloatingPoint* copy();

    std::vector<MutationOpP> getMutationOp() override;
    void registerParameters(StateP state) override;
    bool initialize(StateP state) override;

    const std::vector<double>& getMinVec() const { return _minVec; }
    const std::vector<double>& getMaxVec() const { return _maxVec; }

   private:
    std::vector<double> _minVec;
    std::vector<double> _maxVec;
};
typedef boost::shared_ptr<MyFloatingPoint> MyFloatingPointP;
