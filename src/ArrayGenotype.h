#ifndef ARRAYGENOTYPE_H
#define ARRAYGENOTYPE_H

#include <array>
#include <vector>
#include <string>
#include <sstream>
#include <iomanip>

#include <ecf/ECF_base.h>
#include <ecf/Genotype.h>
#include "ArrayMutOp.h"
#include "ArrayCrxOp.h"


template<uint size>
class ArrayGenotype : public Genotype
{
public:
    static_assert(size > 0, "Size should be greater than 0.");

    ArrayGenotype() {
        name_ = "ArrayGenotype";
    }


    ArrayGenotype<size>* copy(){
        ArrayGenotype<size> *newObject = new ArrayGenotype<size> (*this);
        return newObject;
    }


    std::vector<MutationOpP> getMutationOp() {
        std::vector<MutationOpP> mut;
        mut.push_back(static_cast<MutationOpP> (new ArrayMutOp<size>));
        return mut;
    }


    std::vector<CrossoverOpP> getCrossoverOp() {
        std::vector<CrossoverOpP> crx;
        crx.push_back(static_cast<CrossoverOpP> (new ArrayCrxOp<size>));
        return crx;
    }


    void registerParameters(StateP state) {
        registerParameter(state, "size",
                          static_cast <voidP> (new uint(size)), ECF::UINT);
        std::string min_array("0.0");
        std::string max_array("10.0");

        for (uint k = 1; k < size; ++k) {
            min_array += min_array;
            max_array += max_array;
        }
        registerParameter(state, "min.array",
                          static_cast <voidP> (&min_array), ECF::STRING);

        registerParameter(state, "max.array",
                          static_cast <voidP> (&max_array), ECF::STRING);
    }


    bool initialize(StateP state) {

        std::stringstream ss;
        voidP sptr = getParameterValue(state, "min.array");
        ss << *static_cast <std::string*> (sptr.get());

        for (uint i = 0; i < size; ++i) ss >> minArray[i];

        ss.clear();

        sptr = getParameterValue(state, "max.array");
        ss << *static_cast <std::string*> (sptr.get());

        for (uint i = 0; i < size; ++i) ss >> maxArray[i];

        ss.clear();

        genome = generateRandomArray<size>( minArray, maxArray );
        return true;
    }

    void write(XMLNode &xArrayGenotype) {
        xArrayGenotype = XMLNode::createXMLTopNode("ArrayGenotype");
        std::stringstream ss;
        ss.str("");
        for(uint i = 0; i < size; ++i)
            ss << std::setw(5) << genome[i];

        ss << std::setw(5);
        xArrayGenotype.addText(ss.str().c_str());
    }

    void read(XMLNode& xArrayGenotype);
    /*   {
        XMLCSTR entry = xArrayGenotype.getText();
        std::stringstream ss;
        ss << entry;

        for(uint i = 0; i < genome.size(); ++i)
            ss >> genome[i];
    }
*/

    inline std::array<double, size>& getMinArray() {
        return minArray;
    }


    inline std::array<double, size>& getMaxArray() {
        return maxArray;
    }


    std::array<double, size> genome;

protected:

    std::array<double, size> minArray;
    std::array<double, size> maxArray;
};

template<uint size>
using ArrayGenotypeP = boost::shared_ptr<ArrayGenotype<size>>;

#endif // ARRAYGENOTYPE_H
