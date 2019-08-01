#pragma once

#include <ecf/ECF.h>

#include "gaussian_functor.h"

class EvalOp : public EvaluateOp {
   public:
    EvalOp(const Gaussian_fit &functor);
    FitnessP evaluate(IndividualP individual) override;

   private:
    const Gaussian_fit& _functor;
};
typedef boost::shared_ptr<EvalOp> EvalOpP;
