#include <cg/operations/orientation.h>
#include <cg/operations/orientation_3d.h>

namespace cg {
    struct line
    {
        // line: a * x + b * y + c = 0
        double a, b, c;

        line() {}

        line(const line & l)
            : a(l.a), b(l.b), c(l.c)
        {}

        line(double a, double b, double c)
            : a(a), b(b), c(c)
        {}

        void inverse_vector()
        {
            a *= -1;
            b *= -1;
            c *= -1;
        }
    };

    inline bool is_direct_vector_right(const line & l)
    {
        return l.b < 0 || (l.b == 0 && l.a > 0);
    }

    inline bool is_normal_vector_up(const line & l)
    {
        return !is_direct_vector_right(l);
    }

    bool ray_line_intersection(const line & cross_line, const line & edge_line,
                               const line & sl1, const line & sl2)
    {
        line l = line(cross_line);
        if (!is_normal_vector_up(l)) {
            l.inverse_vector();
        }

        auto lp = point3d{l.a, l.b, l.c};
        auto v1 = point3d{sl1.a, sl1.b, sl1.c};
        auto v2 = point3d{sl2.a, sl2.b, sl2.c};

        int vpos = orientation_3d(lp, v1, v2);
        int vdet = orientation_2d(sl1.a, sl1.b, sl2.a, sl2.b);

        int res = vpos * vdet;

        auto s = point_2{-edge_line.b, edge_line.a};
        auto p = point_2{-l.b, l.a}; // always look left here
        if (is_direct_vector_right(edge_line)) {
            p.x *= -1;
            p.y *= -1;
        }

        int orient = orientation_2d(s.x, s.y, p.x, p.y);

        if (orient == CG_COLLINEAR) return false;

        if (res > 0) { //
            return (!is_direct_vector_right(edge_line) ? orient == CG_RIGHT : orient == CG_LEFT);
        } else {
            return (!is_direct_vector_right(edge_line) ? orient == CG_LEFT : orient == CG_RIGHT);
        }
    }

    int line_point_sign(const line & l, const line & sl1, const line & sl2)
    {
        auto lp = point3d{l.a, l.b, l.c};
        auto v1 = point3d{sl1.a, sl1.b, sl1.c};
        auto v2 = point3d{sl2.a, sl2.b, sl2.c};

        int vpos = orientation_3d(lp, v1, v2);
        int vdet = orientation_2d(sl1.a, sl1.b, sl2.a, sl2.b);

        return vpos * vdet;
    }

    bool segment_line_intersection(const line & l,
                                   const line & sl1, const line & sl2,
                                   const line & dl1, const line & dl2)
    {
        return line_point_sign(l, sl1, sl2) != line_point_sign(l, dl1, dl2);
    }

    template<typename Scalar>
    point_2t<Scalar> intersection_point(const line & l1, const line & l2)
    {
        Scalar det = Scalar(l1.a) * Scalar(l2.b) - Scalar(l2.a) * Scalar(l1.b);
        Scalar x = -(Scalar(l1.c) * Scalar(l2.b) - Scalar(l2.c) * Scalar(l1.b)) / det;
        Scalar y = -(Scalar(l1.a) * Scalar(l2.c) - Scalar(l2.a) * Scalar(l1.c)) / det;

        return point_2t<Scalar>(x, y);
    }

    dead_sign line_position(const line & un_line, const point_2 & p)
    {
        line l = line(un_line);
        if (!is_normal_vector_up(l)) {
            l.inverse_vector();
        }

        int res = line_point_sign(l, line{1, 0, -p.x}, line{0, 1, -p.y});
        switch (res) {
            case  1: return POS_DEAD; break;
            case  0: return ZERO_DEAD; break;
            case -1: return NEG_DEAD; break;
        }

        return ZERO_DEAD;
    }

    orientation_t precise_turn_predicate(const line & l1, const line & l2,
                                         const line & s1, const line & s2,
                                         const line & t1, const line & t2)
    {
        double a11 = l1.a, b11 = l1.b, c11 = l1.c;
        double a12 = l2.a, b12 = l2.b, c12 = l2.c;
        double a21 = s1.a, b21 = s1.b, c21 = s1.c;
        double a22 = s2.a, b22 = s2.b, c22 = s2.c;
        double a31 = t1.a, b31 = t1.b, c31 = t1.c;
        double a32 = t2.a, b32 = t2.b, c32 = t2.c;

        int detl = orientation_2d(a11, b11, a12, b12);
        int dets = orientation_2d(a21, b21, a22, b22);
        int dett = orientation_2d(a31, b31, a32, b32);

        int det_sign = detl * detl * dets * dett;

        double p1x = -c11 * b12 + b11 * c12;
        double p1y = -a11 * c12 + c11 * a12;
        double p2x = -c21 * b22 + b21 * c22;
        double p2y = -a21 * c22 + c21 * a22;
        double p3x = -c31 * b32 + b31 * c32;
        double p3y = -a31 * c32 + c31 * a32;

        double det1 = a11 * b12 - b11 * a12;
        double det2 = a21 * b22 - b21 * a22;
        double det3 = a31 * b32 - b31 * a32;

        double x1 = p2x * det1 - p1x * det2;
        double x2 = p3y * det1 - p1y * det3;
        double x3 = p2y * det1 - p1y * det2;
        double x4 = p3x * det1 - p1x * det3;
        double res_d = x1 * x2 - x3 * x4;
        double eps = (fabs(x1 * x2) + fabs(x3 * x4)) * 45 * std::numeric_limits<double>::epsilon();

        if (res_d > eps)
            return (det_sign > 0 ? CG_LEFT : CG_RIGHT);

        if (res_d < -eps)
            return (det_sign > 0 ? CG_RIGHT : CG_LEFT);

        typedef boost::numeric::interval_lib::unprotect<boost::numeric::interval<double> >::type interval;
        boost::numeric::interval<double>::traits_type::rounding _;

        interval p1xi = interval(-c11) * interval(b12) + interval(b11) * interval(c12);
        interval p1yi = interval(-a11) * interval(c12) + interval(c11) * interval(a12);
        interval p2xi = interval(-c21) * interval(b22) + interval(b21) * interval(c22);
        interval p2yi = interval(-a21) * interval(c22) + interval(c21) * interval(a22);
        interval p3xi = interval(-c31) * interval(b32) + interval(b31) * interval(c32);
        interval p3yi = interval(-a31) * interval(c32) + interval(c31) * interval(a32);

        interval det1i = interval(a11) * interval(b12) - interval(b11) * interval(a12);
        interval det2i = interval(a21) * interval(b22) - interval(b21) * interval(a22);
        interval det3i = interval(a31) * interval(b32) - interval(b31) * interval(a32);

        interval x1i = p2xi * det1i - p1xi * det2i;
        interval x2i = p3yi * det1i - p1yi * det3i;
        interval x3i = p2yi * det1i - p1yi * det2i;
        interval x4i = p3xi * det1i - p1xi * det3i;
        interval res_i = x1i * x2i - x3i * x4i;

        if (res_i.lower() > 0)
            return (det_sign > 0 ? CG_LEFT : CG_RIGHT);

        if (res_i.upper() < 0)
            return (det_sign > 0 ? CG_RIGHT : CG_LEFT);

        if (res_i.upper() == res_i.lower())
            return CG_COLLINEAR;

        mpq_class p1xr = mpq_class(-c11) * mpq_class(b12) + mpq_class(b11) * mpq_class(c12);
        mpq_class p1yr = mpq_class(-a11) * mpq_class(c12) + mpq_class(c11) * mpq_class(a12);
        mpq_class p2xr = mpq_class(-c21) * mpq_class(b22) + mpq_class(b21) * mpq_class(c22);
        mpq_class p2yr = mpq_class(-a21) * mpq_class(c22) + mpq_class(c21) * mpq_class(a22);
        mpq_class p3xr = mpq_class(-c31) * mpq_class(b32) + mpq_class(b31) * mpq_class(c32);
        mpq_class p3yr = mpq_class(-a31) * mpq_class(c32) + mpq_class(c31) * mpq_class(a32);

        mpq_class det1r = mpq_class(a11) * mpq_class(b12) - mpq_class(b11) * mpq_class(a12);
        mpq_class det2r = mpq_class(a21) * mpq_class(b22) - mpq_class(b21) * mpq_class(a22);
        mpq_class det3r = mpq_class(a31) * mpq_class(b32) - mpq_class(b31) * mpq_class(a32);

        mpq_class x1r = p2xr * det1r - p1xr * det2r;
        mpq_class x2r = p3yr * det1r - p1yr * det3r;
        mpq_class x3r = p2yr * det1r - p1yr * det2r;
        mpq_class x4r = p3xr * det1r - p1xr * det3r;
        mpq_class res_r = x1r * x2r - x3r * x4r;

        int cres = cmp(res_r, 0);

        if (cres > 0)
           return (det_sign > 0 ? CG_LEFT : CG_RIGHT);

        if (cres < 0)
           return (det_sign > 0 ? CG_RIGHT : CG_LEFT);

        return CG_COLLINEAR;
    }

    orientation_t point_segment_orientation(const line & sl1, const line & sl2,
                                            const line & dl1, const line & dl2,
                                            const point_2 & c)
    {
        return precise_turn_predicate(sl1, sl2, dl1, dl2, line{1, 0, -c.x}, line{0, 1, -c.y});
    }

    dead_sign x_dif(const line & l1, const line & l2,
                    const line & s1, const line & s2)
    {
        int det_ls = orientation_2d(l1.a, l1.b, l2.a, l2.b) * orientation_2d(s1.a, s1.b, s2.a, s2.b);

        double det1 = -l1.c * l2.b + l1.b * l2.c;
        double det2 = s1.a * s2.b - s1.b * s2.a;
        double det3 = -s1.c * s2.b + s1.b * s2.c;
        double det4 = l1.a * l2.b - l1.b * l2.a;

        double res_d = det1 * det2 - det3 * det4;
        double eps = (fabs(det1 * det2) + fabs(det3 * det4)) * 18 * std::numeric_limits<double>::epsilon();

        if (res_d > eps)
            return (det_ls > 0 ? POS_DEAD : NEG_DEAD);

        if (res_d < -eps)
            return (det_ls > 0 ? NEG_DEAD : POS_DEAD);

        typedef boost::numeric::interval_lib::unprotect<boost::numeric::interval<double> >::type interval;
        boost::numeric::interval<double>::traits_type::rounding _;

        interval res_i = (interval(-l1.c) * interval(l2.b) + interval(l1.b) * interval(l2.c))
                       * (interval(s1.a) * interval(s2.b) - interval(s1.b) * interval(s2.a))
                       - (interval(-s1.c) * interval(s2.b) + interval(s1.b) * interval(s2.c))
                       * (interval(l1.a) * interval(l2.b) - interval(l1.b) * interval(l2.a));

        if (res_i.lower() > 0)
            return (det_ls > 0 ? POS_DEAD : NEG_DEAD);

        if (res_i.upper() < 0)
            return (det_ls > 0 ? NEG_DEAD : POS_DEAD);

        mpq_class res_r = (mpq_class(-l1.c) * mpq_class(l2.b) + mpq_class(l1.b) * mpq_class(l2.c))
                        * (mpq_class(s1.a) * mpq_class(s2.b) - mpq_class(s1.b) * mpq_class(s2.a))
                        - (mpq_class(-s1.c) * mpq_class(s2.b) + mpq_class(s1.b) * mpq_class(s2.c))
                        * (mpq_class(l1.a) * mpq_class(l2.b) - mpq_class(l1.b) * mpq_class(l2.a));

        int cres = cmp(res_r, 0);

        if (cres > 0)
           return (det_ls > 0 ? POS_DEAD : NEG_DEAD);

        if (cres < 0)
           return (det_ls > 0 ? NEG_DEAD : POS_DEAD);

        return ZERO_DEAD;
    }

    /*dead_sign y_dif(const line & l1, const line & l2,
                    const line & s1, const line & s2)
    {
        int det_ls = orientation_2d(l1.a, l1.b, l2.a, l2.b) * orientation_2d(s1.a, s1.b, s2.a, s2.b);

        double det1 = -l1.a * l2.c + l1.c * l2.a;
        double det2 = s1.a * s2.b - s1.b * s2.a;
        double det3 = -s1.a * s2.c + s1.c * s2.a;
        double det4 = l1.a * l2.b - l1.b * l2.a;

        double res_d = det1 * det2 - det3 * det4;
        double eps = (fabs(det1 * det2) + fabs(det3 * det4)) * 18 * std::numeric_limits<double>::epsilon();

        if (res_d > eps)
            return (det_ls > 0 ? POS_DEAD : NEG_DEAD);

        if (res_d < -eps)
            return (det_ls > 0 ? NEG_DEAD : POS_DEAD);

        typedef boost::numeric::interval_lib::unprotect<boost::numeric::interval<double> >::type interval;
        boost::numeric::interval<double>::traits_type::rounding _;

        interval res_i = (interval(-l1.a) * interval(l2.c) + interval(l1.c) * interval(l2.a))
                       * (interval(s1.a) * interval(s2.b) - interval(s1.b) * interval(s2.a))
                       - (interval(-s1.a) * interval(s2.c) + interval(s1.c) * interval(s2.a))
                       * (interval(l1.a) * interval(l2.b) - interval(l1.b) * interval(l2.a));

        if (res_i.lower() > 0)
            return (det_ls > 0 ? POS_DEAD : NEG_DEAD);

        if (res_i.upper() < 0)
            return (det_ls > 0 ? NEG_DEAD : POS_DEAD);

        mpq_class res_r = (mpq_class(-l1.a) * mpq_class(l2.c) + mpq_class(l1.c) * mpq_class(l2.a))
                        * (mpq_class(s1.a) * mpq_class(s2.b) - mpq_class(s1.b) * mpq_class(s2.a))
                        - (mpq_class(-s1.a) * mpq_class(s2.c) + mpq_class(s1.c) * mpq_class(s2.a))
                        * (mpq_class(l1.a) * mpq_class(l2.b) - mpq_class(l1.b) * mpq_class(l2.a));

        int cres = cmp(res_r, 0);

        if (cres > 0)
           return (det_ls > 0 ? POS_DEAD : NEG_DEAD);

        if (cres < 0)
           return (det_ls > 0 ? NEG_DEAD : POS_DEAD);

        return ZERO_DEAD;
    }*/
}
