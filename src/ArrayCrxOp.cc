#include <array>

#include "ArrayCrxOp.h"
#include "ArrayGenotype.h"


bool ArrayCrxOp::mate(GenotypeP gen1, GenotypeP gen2, GenotypeP child) {
    ArrayGenotype* parent1   = static_cast <ArrayGenotype*> (gen1.get());
    ArrayGenotype* parent2   = static_cast <ArrayGenotype*> (gen2.get());
    ArrayGenotype* offspring = static_cast <ArrayGenotype*> (child.get());

    for (uint i = 0; i < parent1->genome.size(); ++i ) {
        switch ( state_->getRandomizer()->getRandomInteger(0, 1) ) {
        case (0): offspring->genome[i] = parent1->genome[i];
            break;
        case (1): offspring->genome[i] = parent2->genome[i];
            break;
        }
    }
    return true;
}
