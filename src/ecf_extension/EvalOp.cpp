#include "EvalOp.h"

#include <limits>

#include <Eigen/Dense>

#include "MyFloatingPoint.h"

EvalOp::EvalOp(const Gaussian_fit& functor) : _functor(functor) {}

FitnessP EvalOp::evaluate(IndividualP individual) {
    FitnessP fitness(new FitnessMin);
    auto gen = static_cast<MyFloatingPoint*>(individual->getGenotype().get());

    assert(_functor.inputs() == gen->realValue.size());

    Eigen::VectorXd xvec = Eigen::Map<Eigen::VectorXd>(gen->realValue.data(), gen->realValue.size());

    Eigen::VectorXd fvec;
    _functor(xvec, fvec);
    const double val = fvec.dot(fvec);
    std::isnan(val) ? fitness->setValue(std::numeric_limits<double>::max()) : fitness->setValue(val);
    return fitness;
}
