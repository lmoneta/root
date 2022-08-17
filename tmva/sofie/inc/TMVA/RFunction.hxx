#ifndef TMVA_SOFIE_RFUNCTION
#define TMVA_SOFIE_RFUNCTION


#include "TMVA/ROperator.hxx"

namespace TMVA{
namespace Experimental{
namespace SOFIE{

class RModel;

class RFunction: public ROperator{
    enum class FunctionType{
        UPDATION=0, AGGREGATE=2
    };
    std::unique_ptr<RModel> function_block;
};


}//SOFIE
}//Experimental
}//TMVA


#endif //TMVA_SOFIE_RFUNCTION