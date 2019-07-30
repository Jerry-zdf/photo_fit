#include "EvalOp.h"
#include "MyFloatingPoint.h"

#include <Eigen/Dense>
#include <limits>

FitnessP EvalOp::evaluate(IndividualP individual) {
    FitnessP fitness (new FitnessMin);
    MyFloatingPoint* gen = static_cast<MyFloatingPoint*> (individual->getGenotype().get());

    Eigen::VectorXd fvec(_functor.values());
    Eigen::VectorXd xvec(_functor.inputs());

    assert(xvec.size() == gen->realValue.size());
    for(uint i = 0; i < xvec.size(); ++i)
        xvec[i] = gen->realValue[i];

    _functor(xvec, fvec);
    const double val = fvec.dot(fvec);
    std::isnan(val) ? fitness->setValue(std::numeric_limits<double>::max()) : fitness->setValue(val);
    return fitness;
}
