#ifndef ARRAYCRXOP_H
#define ARRAYCRXOP_H

#include <ecf/ECF_base.h>
#include <ecf/Crossover.h>

class ArrayCrxOp : public CrossoverOp
{
public:
    bool mate(GenotypeP gen1, GenotypeP gen2, GenotypeP child);
};

typedef boost::shared_ptr <ArrayCrxOp> ArrayCrxOpP;

#endif // ARRAYCRXOP_H
