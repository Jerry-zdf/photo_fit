#include "MyFloatingPoint.h"

#include <vector>
#include <string>
#include <sstream>


std::vector<MutationOpP> MyFloatingPoint::getMutationOp() {
    std::vector<MutationOpP> mut;
    mut.push_back(static_cast<MutationOpP> (new GaussMutOp));
    return mut;
}


void MyFloatingPoint::registerParameters(StateP state) {
    registerParameter(state, "size",
                      static_cast <voidP> (new uint(1)), ECF::UINT);

    registerParameter(state, "min.vec",
                      static_cast <voidP> (new std::string("0.0")), ECF::STRING);

    registerParameter(state, "max.vec",
                      static_cast <voidP> (new std::string("10.0")), ECF::STRING);
}


bool MyFloatingPoint::initialize(StateP state) {
    if(!isParameterDefined(state, "min.vec") ||
            !isParameterDefined(state, "max.vec") ||
            !isParameterDefined(state, "size")) {
        ECF_LOG_ERROR(state, "Error: required parameters for MyFloatingPoint genotype not defined (min.vec, max.vec, size)!");
        throw("");
    }

    voidP sptr = getParameterValue(state, "size");
    nDimension_ = *static_cast<uint*> (sptr.get());

    if(nDimension_ < 1) {
        ECF_LOG_ERROR(state, "Error: 'size' must be > 0 for FloatingPoint genotype!");
        throw("");
    }
    realValue.resize(nDimension_);
    minVec.resize(nDimension_);
    maxVec.resize(nDimension_);

    std::stringstream ss;

    sptr = getParameterValue(state, "min.vec");
    ss << *static_cast<std::string*> (sptr.get());

    for (uint i = 0; i < nDimension_; ++i) ss >> minVec[i];
    ss.clear();

    sptr = getParameterValue(state, "max.vec");
    ss << *static_cast <std::string*> (sptr.get());

    for (uint i = 0; i < nDimension_; ++i) ss >> maxVec[i];
    ss.clear();

    for (uint i = 0; i < nDimension_; ++i){
            realValue[i] = ( minVec[i] + (maxVec[i] - minVec[i]) * state->getRandomizer()->getRandomDouble() );
    }
    return true;
}
