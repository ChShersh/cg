#pragma once

namespace cg {

    struct point3d {
        double x, y, z;

        point3d(double x, double y, double z)
            : x(x), y(y), z(z)
        {}

        point3d(point3d const & o)
            : x(o.x), y(o.y), z(o.z)
        {}

        point3d()
            : x(0), y(0), z(0)
        {}
   };
   
}
