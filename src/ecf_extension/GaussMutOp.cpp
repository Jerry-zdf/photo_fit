#include "GaussMutOp.h"

#include <random>
#include <sstream>
#include <string>

#include "MyFloatingPoint.h"
#include "utils.h"

void GaussMutOp::registerParameters(StateP state) {
    myGenotype_->registerParameter(state, "mut.gauss.stdev",
                                   static_cast<voidP>(new std::string("0.2")), ECF::STRING);
    myGenotype_->registerParameter(state, "mut.gauss.vecrate",
                                   static_cast<voidP>(new std::string("0.3")), ECF::STRING);
    myGenotype_->registerParameter(state, "mut.gauss.indrate",
                                   static_cast<voidP>(new double(0.3)), ECF::DOUBLE);
}

bool GaussMutOp::initialize(StateP state) {
    voidP sptr   = myGenotype_->getParameterValue(state, "mut.gauss.indrate");
    probability_ = *static_cast<double*>(sptr.get());

    sptr          = myGenotype_->getParameterValue(state, "size");
    auto gen_size = *static_cast<uint*>(sptr.get());

    _vecRate.resize(gen_size);
    _stdDev.resize(gen_size);

    sptr = myGenotype_->getParameterValue(state, "mut.gauss.vecrate");
    std::stringstream ss_rate(*static_cast<std::string*>(sptr.get()));
    sptr = myGenotype_->getParameterValue(state, "mut.gauss.stdev");
    std::stringstream ss_dev(*static_cast<std::string*>(sptr.get()));

    for (uint i = 0; i < gen_size; ++i) {
        ss_rate >> _vecRate[i];
        ss_dev >> _stdDev[i];
    }
    return true;
}

bool GaussMutOp::mutate(GenotypeP gene) {
    auto FP = static_cast<MyFloatingPoint*>(gene.get());

    for (uint i = 0; i < (uint)FP->realValue.size(); ++i) {
        if (state_->getRandomizer()->getRandomDouble() < _vecRate[i]) {
            std::normal_distribution<double> dist(0.0, _stdDev[i]);
            FP->realValue[i] += dist(get_random_engine());
            if (FP->realValue[i] > FP->getMaxVec()[i])
                FP->realValue[i] = FP->getMaxVec()[i];
            if (FP->realValue[i] < FP->getMinVec()[i])
                FP->realValue[i] = FP->getMinVec()[i];
        }
    }
    return true;
}
