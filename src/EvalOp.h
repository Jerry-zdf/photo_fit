#ifndef EVALOP_H
#define EVALOP_H

#include <ecf/ECF.h>
#include "GaussianFitFunctor.h"


class EvalOp : public EvaluateOp
{
public:
    EvalOp(const GaussianFit &functor_) : functor(functor_) {}
    FitnessP evaluate(IndividualP individual);

private:
    GaussianFit functor;
};
typedef boost::shared_ptr<EvalOp> EvalOpP;

#endif // EVALOP_H
