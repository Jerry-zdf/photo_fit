#include "EvalOp.h"
#include "MyFloatingPoint.h"

#include <Eigen/Dense>
#include <cassert>


FitnessP EvalOp::evaluate(IndividualP individual) {
    FitnessP fitness (new FitnessMin);
    MyFloatingPoint* gen = static_cast<MyFloatingPoint*> (individual->getGenotype().get());

    Eigen::VectorXd fvec(functor.values());
    Eigen::VectorXd xvec(functor.inputs());

    assert(xvec.size() == gen->realValue.size());
    for(uint i = 0; i < xvec.size(); ++i)
        xvec[i] = gen->realValue[i];

    functor(xvec, fvec);
    double val = fvec.dot(fvec);
    std::isnan(val) ? fitness->setValue(1.e+20) : fitness->setValue(val);
    return fitness;
}
