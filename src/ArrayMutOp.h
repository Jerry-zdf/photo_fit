#ifndef ARRAYMUTOP_H
#define ARRAYMUTOP_H

#include <ecf/ECF_base.h>
#include <ecf/Mutation.h>

#include <array>
#include <string>
#include <sstream>

#include "RandomNumbers.h"


template<uint size>
class ArrayMutOp: public MutationOp
{
public:
    bool mutate( GenotypeP gene ) {
        ArrayGenotype* GT = static_cast <ArrayGenotype*> (gene.get());

        double offset;
        for(uint i = 0; i < size; ++i) {
            if ( state_->getRandomizer()->getRandomDouble() < allRate[i]) {
                std::normal_distribution <double> dist( 0.0, stdDev[i] );
                offset = dist( getRandomEngine() );
                GT->genome[i] += offset;

                GT->genome[i] > GT->getMaxArray()[i] ? GT->genome[i] = GT->getMaxArray()[i];
                GT->genome[i] < GT->getMinArray()[i] ? GT->genome[i] = GT->getMinArray()[i];
            }
        }
        return true;
    }


    bool initialize( StateP state ) {
        voidP sptr = myGenotype_->getParameterValue( state, "mut.array.indrate" );
        probability_ = *static_cast <double*> (sptr.get());

        std::stringstream ss;
        sptr = myGenotype_->getParameterValue( state, "mut.array.allrate" );
        allRate = *static_cast <double*> (sptr.get());

        ss << *static_cast <std::string*> (sptr.get());
        for (uint i = 0; i < size; ++i) ss >> allRate[i];
        ss.clear();

        sptr = myGenotype_->getParameterValue( state, "mut.array.stdeviation" );
        standardDviation = *static_cast <double*> (sptr.get());

        ss << *static_cast <std::string*> (sptr.get());
        for (uint i = 0; i < size; ++i) ss >> stdDev[i];
        ss.clear();

        return true;
    }


    void registerParameters( StateP state ) {
        std::string stdDev_array("0.1");
        std::string allRate_array("0.3");

        for (uint k = 1; k < size; ++k) {
            stdDev_array  += stdDev_array;
            allRaet_array += allRate_array;
        }
        myGenotype_->registerParameter( state, "mut.array.stdeviation",
                                        static_cast<voidP> (&stdDev_array), ECF::STRING);

        myGenotype_->registerParameter( state, "mut.array.allrate",
                                        static_cast<voidP> (&allRate_array), ECF::STRING );

        myGenotype_->registerParameter( state, "mut.array.indrate",
                                        static_cast<voidP> (new double(0.3)), ECF::DOUBLE );
    }


    std::array<double, size> stdDev;
    std::array<double, size> allRate;
};

template<uint size>
using ArrayMutOpP = boost::shared_ptr<ArrayMutOp<size>>;

#endif // ARRAYMUTOP_H
