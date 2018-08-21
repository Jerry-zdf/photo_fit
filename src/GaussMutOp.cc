#include "RandomNumbers.h"
#include "GaussMutOp.h"
#include "MyFloatingPoint.h"

#include <string>
#include <sstream>



void GaussMutOp::registerParameters( StateP state ) {
    myGenotype_->registerParameter( state, "mut.gauss.stdev",
                                    static_cast <voidP> (new std::string("0.2")), ECF::STRING );
    myGenotype_->registerParameter( state, "mut.gauss.vecrate",
                                    static_cast <voidP> (new std::string("0.3")), ECF::STRING );
    myGenotype_->registerParameter( state, "mut.gauss.indrate",
                                    static_cast <voidP> (new double(0.3)), ECF::DOUBLE );
}


bool GaussMutOp::initialize( StateP state ) {
    voidP sptr = myGenotype_->getParameterValue( state, "mut.gauss.indrate" );
    probability_ = *static_cast<double*> (sptr.get());

    sptr = myGenotype_->getParameterValue( state, "size" );
    uint gen_size = *static_cast<uint*> (sptr.get());

    vecRate.resize(gen_size);
    stdDev.resize(gen_size);

    std::stringstream ss;
    sptr = myGenotype_->getParameterValue( state, "mut.gauss.vecrate" );
    ss << *static_cast<std::string*> (sptr.get());

    uint i = 0;
    while(ss) {
        ss >> vecRate[i];
        if(++i == gen_size) break;
    }
    ss.clear();

    sptr = myGenotype_->getParameterValue( state, "mut.gauss.stdev" );
    ss << *static_cast<std::string*> (sptr.get());
    i = 0;
    while(ss) {
        ss >> stdDev[i];
        if(++i == gen_size) break;
    }
    ss.clear();
    return true;
}


bool GaussMutOp::mutate( GenotypeP gene ) {
    MyFloatingPoint* FP = static_cast<MyFloatingPoint*> (gene.get());
    std::normal_distribution <double> dist;
    double offset;

    for (uint i = 0; i < (uint) FP->realValue.size(); i++) {
        if ( state_->getRandomizer()->getRandomDouble() < vecRate[i]) {
            dist = std::normal_distribution <double>(0.0, stdDev[i]);
            offset = dist(getRandomEngine());
            FP->realValue[i] += offset;
            if(FP->realValue[i] > FP->getMaxVec()[i]) FP->realValue[i] = FP->getMaxVec()[i];
            if(FP->realValue[i] < FP->getMinVec()[i]) FP->realValue[i] = FP->getMinVec()[i];
        }
    }
    return true;
}

