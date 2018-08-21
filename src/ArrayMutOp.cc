#include <random>

#include "ArrayMutOp.h"
#include "ArrayGenotype.h"


void VectorMutOp::registerParameters( StateP state )
{
    myGenotype_->registerParameter( state, "mut.vector.stdeviation",
                                    static_cast<voidP> (new double(0.1)), ECF::DOUBLE );

    myGenotype_->registerParameter( state, "mut.vector.allrate",
                                    static_cast<voidP> (new double(0.3)), ECF::DOUBLE );

    myGenotype_->registerParameter( state, "mut.vector.indrate",
                                    static_cast<voidP> (new double(0.3)), ECF::DOUBLE );
}


bool VectorMutOp::initialize( StateP state )
{
    voidP sptr = myGenotype_->getParameterValue( state, "mut.vector.indrate" );
    probability_ = *static_cast <double*> (sptr.get());

    sptr = myGenotype_->getParameterValue( state, "mut.vector.allrate" );
    allRate = *static_cast <double*> (sptr.get());

    sptr = myGenotype_->getParameterValue( state, "mut.vector.stdeviation" );
    standardDviation = *static_cast <double*> (sptr.get());

    engine.seed( static_cast <uint> (time(NULL)) );

    return true;
}


bool VectorMutOp::mutate( GenotypeP gene )
{
    VectorGenotype* GT = static_cast <VectorGenotype*> (gene.get());

    std::normal_distribution <double> Distribution( 0.0, standardDviation );
    double offset;

    for (uint i = 0; i < (uint) GT->genome.size(); ++i)
    {
        if ( state_->getRandomizer()->getRandomDouble() < allRate)
        {
            for (uint j = 0; j < 3; ++j)
            {
                offset = Distribution(engine);
                GT->genome[i][j] += offset;

                if( GT->genome[i][j] > GT->getMaxVector()[j] )
                    GT->genome[i][j] = GT->getMaxVector()[j];

                if( GT->genome[i][j] < GT->getMinVector()[j] )
                    GT->genome[i][j] = GT->getMinVector()[j];
            }
        }
    }

    return true;
}
