#pragma once

#include <boost/numeric/interval.hpp>
#include <gmpxx.h>

#include <cg/primitives/point3d.h>

namespace cg {

   enum dead_sign
   {
      NEG_DEAD = -1,
      ZERO_DEAD = 0,
      POS_DEAD = 1
   };

   // 2D determinant sigh
   inline dead_sign orientation_2d(double a, double b, double c, double d)
   {
       /*
        * | a b |
        * | c d |
        */

       double l = a * d;
       double r = b * c;
       double dres = l - r;
       double eps = (fabs(l) + fabs(r)) * 4 * std::numeric_limits<double>::epsilon();

       if (dres > eps) {
           return POS_DEAD;
       }

       if (dres < -eps) {
           return NEG_DEAD;
       }

       typedef boost::numeric::interval_lib::unprotect<boost::numeric::interval<double>>::type interval;
       boost::numeric::interval<double>::traits_type::rounding _;

       interval res = interval(a) * interval(d) - interval(b) * interval(c);
       if (res.lower() > 0) {
           return POS_DEAD;
       }

       if (res.upper() < 0) {
           return NEG_DEAD;
       }

       mpq_class mres = mpq_class(a) * mpq_class(d) - mpq_class(b) * mpq_class(c);
       int cres = cmp(mres, 0);
       if (cres > 0) {
           return POS_DEAD;
       }

       if (cres < 0) {
           return NEG_DEAD;
       }

       return ZERO_DEAD;
   }

   // 3D orientation
   inline dead_sign orientation_3d(const point3d & a, const point3d & b, const point3d & p)
   {
       double l = a.x * (b.y * p.z - b.z * p.y);
       double m = a.y * (b.z * p.x - b.x * p.z);
       double r = a.z * (b.x * p.y - b.y * p.x);
       double dres = l + m + r;
       double eps = (fabs(l) + fabs(r) + fabs(m)) * 16 * std::numeric_limits<double>::epsilon();

       if (dres > eps) {
           return POS_DEAD;
       }

       if (dres < -eps) {
           return NEG_DEAD;
       }

       typedef boost::numeric::interval_lib::unprotect<boost::numeric::interval<double>>::type interval;

       boost::numeric::interval<double>::traits_type::rounding _;

       interval res = interval(a.x) * (interval(b.y) * interval(p.z) - interval(b.z) * interval(p.y))
                    + interval(a.y) * (interval(b.z) * interval(p.x) - interval(b.x) * interval(p.z))
                    + interval(a.z) * (interval(b.x) * interval(p.y) - interval(b.y) * interval(p.x));

       if (res.lower() > 0) {
           return POS_DEAD;
       }

       if (res.upper() < 0) {
           return NEG_DEAD;
       }

       mpq_class mres = mpq_class(a.x) * (mpq_class(b.y) * mpq_class(p.z) - mpq_class(b.z) * mpq_class(p.y))
                      + mpq_class(a.y) * (mpq_class(b.z) * mpq_class(p.x) - mpq_class(b.x) * mpq_class(p.z))
                      + mpq_class(a.z) * (mpq_class(b.x) * mpq_class(p.y) - mpq_class(b.y) * mpq_class(p.x));

       int cres = cmp(mres, 0);

       if (cres > 0) {
           return POS_DEAD;
       }

       if (cres < 0) {
           return NEG_DEAD;
       }

       return ZERO_DEAD;
   }

}
