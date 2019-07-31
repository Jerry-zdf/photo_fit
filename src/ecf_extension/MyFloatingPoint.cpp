#include "MyFloatingPoint.h"

#include <sstream>
#include <string>
#include <vector>

#include "GaussMutOp.h"

MyFloatingPoint::MyFloatingPoint() {
    name_ = "MyFPGenotype";
}

MyFloatingPoint* MyFloatingPoint::copy() {
    MyFloatingPoint* newObject = new MyFloatingPoint(*this);
    return newObject;
}

std::vector<MutationOpP> MyFloatingPoint::getMutationOp() {
    std::vector<MutationOpP> mut;
    mut.push_back(static_cast<MutationOpP>(new GaussMutOp));
    return mut;
}

void MyFloatingPoint::registerParameters(StateP state) {
    registerParameter(state, "size",
                      static_cast<voidP>(new uint(1)), ECF::UINT);

    registerParameter(state, "min.vec",
                      static_cast<voidP>(new std::string("0.0")), ECF::STRING);

    registerParameter(state, "max.vec",
                      static_cast<voidP>(new std::string("10.0")), ECF::STRING);
}

bool MyFloatingPoint::initialize(StateP state) {
    if (!isParameterDefined(state, "min.vec") ||
        !isParameterDefined(state, "max.vec") ||
        !isParameterDefined(state, "size")) {
        ECF_LOG_ERROR(state, "Error: Parameters required for MyFloatingPoint genotype not defined (min.vec, max.vec, size)!");
        throw std::runtime_error("Error: Parameters required for MyFloatingPoint genotype not defined (min.vec, max.vec, size)!");
    }

    voidP sptr  = getParameterValue(state, "size");
    nDimension_ = *static_cast<uint*>(sptr.get());

    if (nDimension_ < 1) {
        ECF_LOG_ERROR(state, "Error: 'size' must be > 0 for FloatingPoint genotype!");
        throw std::runtime_error("Error: 'size' must be > 0 for FloatingPoint genotype!");
    }
    realValue.resize(nDimension_);
    _minVec.resize(nDimension_);
    _maxVec.resize(nDimension_);

    sptr = getParameterValue(state, "min.vec");
    std::stringstream ss_min(*static_cast<std::string*>(sptr.get()));
    sptr = getParameterValue(state, "max.vec");
    std::stringstream ss_max(*static_cast<std::string*>(sptr.get()));

    for (uint i = 0; i < nDimension_; ++i) {
        ss_min >> _minVec[i];
        ss_max >> _maxVec[i];
        realValue[i] = (_minVec[i] + (_maxVec[i] - _minVec[i]) * state->getRandomizer()->getRandomDouble());
    }
    return true;
}
