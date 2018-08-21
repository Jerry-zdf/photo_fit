/** Asimple program for finding a minimum of f(x_1, ... , x_k) = cos(x_1) + ... + cos(x_k)
 * over the domain [0, 2Pi]^k.
 * Representation is set as a floating point vector.
 * Set the parameters file for details of genotype.
 */

#include "MinEvalOp.h"

#include <cmath>

FitnessP MinEvalOp::evaluate(IndividualP individual)
{
    FitnessP fitness (new FitnessMin);
    ArrayGenotype<4>* gen = static_cast <ArrayGenotype<4>*> (individual->getGenotype().get());
    double val = 0.0;

    for (uint i = 0; i < gen->genome.size(); i++)
        val += std::cos( gen->genome[i] );

    fitness->setValue(val);
    return fitness;
}
