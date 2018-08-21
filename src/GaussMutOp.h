#ifndef GAUSSMUTOP_H
#define GAUSSMUTOP_H

#include <vector>
#include <ecf/ECF_base.h>
#include <ecf/Mutation.h>


class GaussMutOp: public MutationOp
{
public:
    bool mutate( GenotypeP gene );
    bool initialize( StateP state );
    void registerParameters( StateP state );

    std::vector<double> stdDev;
    std::vector<double> vecRate;

};
typedef boost::shared_ptr< GaussMutOp > GaussMutOpP;

#endif // GAUSSMUTOP_H
