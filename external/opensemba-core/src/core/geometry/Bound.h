
#include "Box.h"
#include "element/Element.h"

namespace semba {
namespace geometry {

template< class It >
BoxR3 getBound(It first, It last) {
        
    if (first == last) {
        return BoxR3().setInfinity();
    }

    BoxR3 bound;
        
    while (first != last) {
        if ((*first)->template is<ElemR>()) {
            bound << (*first)->template castTo<ElemR>()->getBound();
        }

        if ((*first)->template is<ElemI>()) {
            BoxI3 boxI = (*first)->template castTo<ElemI>()->getBound();
            math::CVecI3 minP = boxI.getMin();
            math::CVecI3 maxP = boxI.getMax();
            using math::CVecR3;
            using namespace math::Constants;
            bound << BoxR3(CVecR3(minP(x), minP(y), minP(z)), CVecR3(maxP(x), maxP(y), maxP(z)));
        }
            
        first++;
    }

    return bound;
}

} /* namespace Geometry */
} /* namespace SEMBA */
