#pragma once

#include <cg/primitives/line.h>

#include <array>
#include <memory>

namespace cg {
    struct line_cross
    {
        std::shared_ptr<line> l1, l2;

        line_cross() {}

        line_cross(const std::shared_ptr<line> & line1, const std::shared_ptr<line> & line2)
            : l1(line1), l2(line2)
        {}
    };

    struct triangle_k
    {
       triangle_k() {}

       triangle_k(line_cross const & a, line_cross const & b, line_cross const & c)
          : pts_( {{a, b, c}} ) {}

       triangle_k(const triangle_k & other)
           : pts_({{other[0], other[1], other[2]}}) {}

       line_cross &         operator [] (size_t id)       { return pts_[id]; }
       line_cross const &   operator [] (size_t id) const { return pts_[id]; }

       void reverse() {
          std::swap(pts_[0], pts_[2]);
       }

    private:
       std::array<line_cross, 3> pts_;
    };

    bool triangle_contains_point(const triangle_k & t, const point_2 & p)
    {
        bool inside = true;
        for (int i = 0; i < 3 && inside; i++) {
            inside &= (CG_RIGHT != point_segment_orientation(*t[i].l1, *t[i].l2, *t[(i + 1) % 3].l1, *t[(i + 1) % 3].l2, p));
        }
        return inside;
    }

    bool triangle_contains_star_point(const triangle_k & t, const line & l1, const line & l2)
    {
        bool inside = true;
        for (int i = 0; i < 3 && inside; i++) {
            inside &= (CG_LEFT == precise_turn_predicate(*t[i].l1, *t[i].l2, *t[(i + 1) % 3].l1, *t[(i + 1) % 3].l2, l1, l2));
        }
        return inside;
    }

    bool triangle_contains_convex_point(const triangle_k & t, const line & l1, const line & l2)
    {
        bool inside = true;
        for (int i = 0; i < 3 && inside; i++) {
            inside &= (CG_RIGHT != precise_turn_predicate(*t[i].l1, *t[i].l2, *t[(i + 1) % 3].l1, *t[(i + 1) % 3].l2, l1, l2));
        }
        return inside;
    }

    bool triangle_intersection(const triangle_k & t1, const triangle_k & t2)
    {
        bool has_intersection = false;
        for (int i = 0; i < 3 && !has_intersection; i++) {
            has_intersection = triangle_contains_convex_point(t1, *t2[i].l1, *t2[i].l2);
        }

        if (!has_intersection) {
            for (int i = 0; i < 3 && !has_intersection; i++) {
                has_intersection = triangle_contains_convex_point(t2, *t1[i].l1, *t1[i].l2);
            }
        }

        if (!has_intersection) {
            for (int i = 0; i < 3 && !has_intersection; i++) {
                for (int j = 0; j < 3 && !has_intersection; j++) {
                    std::shared_ptr<line> t1l11 = t1[i].l1, t1l12 = t1[i].l2, t1l21 = t1[(i + 1) % 3].l1, t1l22 = t1[(i + 1) % 3].l2;
                    std::shared_ptr<line> t2l11 = t2[j].l1, t2l12 = t2[j].l2, t2l21 = t2[(j + 1) % 3].l1, t2l22 = t2[(j + 1) % 3].l2;

                    orientation_t turn1 = precise_turn_predicate(*t1l11, *t1l12, *t1l21, *t1l22, *t2l11, *t2l12);
                    orientation_t turn2 = precise_turn_predicate(*t1l11, *t1l12, *t1l21, *t1l22, *t2l21, *t2l22);

                    if (turn1 == turn2 && turn1 == CG_COLLINEAR) {
                        std::shared_ptr<line> xminl1 = t1l11, xminl2 = t1l12, xmaxl1 = t1l21, xmaxl2 = t1l22;
                        std::shared_ptr<line> xmins1 = t2l11, xmins2 = t2l12, xmaxs1 = t2l21, xmaxs2 = t2l22;

                        if (x_dif(*xminl1, *xminl2, *xmaxl1, *xmaxl2) > 0) {
                            std::swap(xminl1, xmaxl1);
                            std::swap(xminl2, xmaxl2);
                        }
                        if (x_dif(*xmins1, *xmins2, *xmaxs1, *xmaxs2) > 0) {
                            std::swap(xmins1, xmaxs1);
                            std::swap(xmins2, xmaxs2);
                        }

                        //if (ab[0] == ab[1] && ab[0] == CG_COLLINEAR)
                        //   return (min(a) <= b[0] && max(a) >= b[0])
                        //      || (min(a) <= b[1] && max(a) >= b[1])
                        //      || (min(b) <= a[0] && max(b) >= a[0])
                        //      || (min(b) <= a[1] && max(b) >= a[1]);

                        bool bound1 = x_dif(*t2l11, *t2l12, *xminl1, *xminl2) >= 0 && x_dif(*t2l11, *t2l12, *xmaxl1, *xmaxl2) <= 0;
                        bool bound2 = x_dif(*t2l21, *t2l22, *xminl1, *xminl2) >= 0 && x_dif(*t2l21, *t2l22, *xmaxl1, *xmaxl2) <= 0;
                        bool bound3 = x_dif(*t1l11, *t1l12, *xmins1, *xmins2) >= 0 && x_dif(*t1l11, *t1l12, *xmaxs1, *xmaxs2) <= 0;
                        bool bound4 = x_dif(*t1l21, *t1l22, *xmins1, *xmins2) >= 0 && x_dif(*t1l21, *t1l22, *xmaxs1, *xmaxs2) <= 0;

                        has_intersection = bound1 || bound2 || bound3 || bound4;
                    } else if (turn1 != turn2) {
                        orientation_t turn3 = precise_turn_predicate(*t2l11, *t2l12, *t2l21, *t2l22, *t1l11, *t1l12);
                        orientation_t turn4 = precise_turn_predicate(*t2l11, *t2l12, *t2l21, *t2l22, *t1l21, *t1l22);

                        has_intersection = turn3 != turn4;
                    }
                }
            }
        }

        return has_intersection;
    }
}
