#include "ArrayGenotype.h"

#include <string>
#include <sstream>


ArrayGenotype* ArrayGenotype::copy() {
    ArrayGenotype *newObject = new ArrayGenotype(*this);
    return newObject;
}


std::vector<MutationOpP> ArrayGenotype::getMutationOp()
{
    std::vector<MutationOpP> mut;
    mut.push_back(static_cast<MutationOpP> (new ArrayMutOp));
    return mut;
}

std::vector<CrossoverOpP> ArrayGenotype::getCrossoverOp()
{
    std::vector<CrossoverOpP> crx;
    crx.push_back(static_cast<CrossoverOpP> (new ArrayCrxOp));
    return crx;
}


void ArrayGenotype::registerParameters(StateP state)
{
    registerParameter(state, "size",
                      static_cast <voidP> (new uint(1)), ECF::UINT);

    registerParameter(state, "min.vector",
                      static_cast <voidP> (new std::string("0.0 0.0 0.0 0.0")), ECF::STRING);

    registerParameter(state, "max.vector",
                      static_cast <voidP> (new std::string("10.0 10.0 1.0 1.0")), ECF::STRING);
}


bool ArrayGenotype::initialize(StateP state)
{
    if(!isParameterDefined(state, "size")) {
        state->getLogger()->log(1, "Error: Genotype size not defined");
        throw("");
    }

    voidP sizeP = getParameterValue(state, "size");
    uint size = *static_cast <uint*> (sizeP.get());
    genome.resize(size);

    std::stringstream ss;

    voidP sptr = getParameterValue(state, "min.vector");
    ss << *static_cast <std::string*> (sptr.get());

    for (uint i = 0; i < 3; ++i)
        ss >> minArray[i];

    ss.clear();

    sptr = getParameterValue(state, "max.vector");
    ss << *static_cast <std::string*> (sptr.get());

    for (uint i = 0; i < 3; ++i)
        ss >> maxArray[i];

    ss.clear();

    for(uint i = 0; i < size; ++i)
        genome[i] = Array <3>::generateRandomArray( minArray, maxArray );

    return true;
}


void ArrayGenotype::write(XMLNode &xArrayGenotype)
{
    xArrayGenotype = XMLNode::createXMLTopNode("ArrayGenotype");
    std::stringstream ss;
    ss << genome.size();
    xArrayGenotype.addAttribute("size", ss.str().c_str());

    ss.str("");
    for(uint i = 0; i < genome.size(); ++i)
        ss << "\t" << genome[i];

    xArrayGenotype.addText(ss.str().c_str());
}
