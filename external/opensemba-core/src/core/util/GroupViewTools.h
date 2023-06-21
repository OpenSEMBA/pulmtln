#pragma once

#include <algorithm>
#include <vector>
#include <functional>

namespace semba {
namespace util {
    
class View {
public:
    static inline geometry::BoxR3 getBound(const geometry::ElemView& view)
    {
        if (view.size() == 0) {
            return geometry::BoxR3().setInfinity();
        }

        geometry::BoxR3 bound;

        for (const auto& elem : view) {
            if (!elem->is<geometry::ElemR>()) {
                continue;
            }

            bound << elem->castTo<geometry::ElemR>()->getBound();
        }

        for (auto const& elem : view) {
            if (!elem->is<geometry::ElemI>()) {
                continue;
            }

            geometry::BoxI3 boxI = elem->castTo<geometry::ElemI>()->getBound();
            math::CVecI3 minP = boxI.getMin();
            math::CVecI3 maxP = boxI.getMax();
            using math::CVecR3;
            using namespace math::Constants;
            bound << geometry::BoxR3(CVecR3(minP(x), minP(y), minP(z)), CVecR3(maxP(x), maxP(y), maxP(z)));
        }

        return bound;
    }

    static inline std::vector<const geometry::ElemR*> filterView(
        const std::vector<const geometry::ElemR*>& view,
        const std::function<bool(const geometry::ElemR*)>& map
    ) {
        std::vector<const geometry::ElemR*> res;

        std::copy_if(
            view.begin(),
            view.end(),
            std::back_inserter(res),
            map
        );

        return res;
    }

    // TODO: To update Structured mesh and consider using Integer inside and communicate with Real outside
    static inline std::vector<const geometry::ElemR*> castToReal(
        const geometry::ElemView& view
    ) {

        for (const auto& elem : view) {
            if (!elem->is<geometry::ElemR>()) {
                throw std::logic_error("View contains elements of type different than Geometry::ElemR");
            }
        }

        std::vector<const geometry::ElemR*> res(view.size());
        std::transform(
            view.cbegin(),
            view.cend(),
            res.begin(),
            [](const geometry::Elem* elem) { return elem->castTo<geometry::ElemR>(); }
        );

        return res;
    }
};

}
}