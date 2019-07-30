#pragma once

#include <vector>
#include <ecf/ECF_base.h>
#include <ecf/Mutation.h>


class GaussMutOp: public MutationOp
{
public:
    bool mutate( GenotypeP gene );
    bool initialize( StateP state );
    void registerParameters( StateP state );

private:
    std::vector<double> _stdDev;
    std::vector<double> _vecRate;

};
typedef boost::shared_ptr< GaussMutOp > GaussMutOpP;
