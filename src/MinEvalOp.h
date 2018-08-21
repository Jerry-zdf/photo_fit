/** A simple evaluation oerator class for function minimization problem
 */

#ifndef MINEVALOP_H
#define MINEVALOP_H

#include <ecf/ECF.h>
#include "ArrayGenotype.h"

class MinEvalOp : public EvaluateOp
{
public:
    FitnessP evaluate(IndividualP individual);
};

#endif // MINEVALOP_H
