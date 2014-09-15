#pragma once

#include <cg/primitives/point.h>

#include <random>
#include <string>

namespace util
{

    template<class ValueType, class Distribution, class RandomDevice>
    class random_generator
    {
        RandomDevice generator;
        Distribution d;
    public:
        random_generator(ValueType min, ValueType max) : d(min, max)
        {}

        random_generator() : random_generator(std::numeric_limits<ValueType>::min(),
                                                  std::numeric_limits<ValueType>::max())
        {}

        ValueType operator() ()
        {
            return d(generator);
        }

        void reset(ValueType min, ValueType max)
        {
            d = Distribution(min, max);
        }

        random_generator& operator>> (ValueType& rhs)
        {
            rhs = (*this)();
            return *this;
        }
    };

    template<class Int, class RandomDevice = std::random_device>
    using uniform_random_int = random_generator<Int,
          std::uniform_int_distribution<Int>, RandomDevice>;

    template<class Real, class RandomDevice = std::random_device>
    using uniform_random_real = random_generator<Real,
          std::uniform_real_distribution<Real>, RandomDevice>;

    inline std::string randomString(int length)
    {
        uniform_random_int<char> rand('a', 'z');
        std::string result;
        result.reserve(length);

        for (int i = 0; i < length; ++i)
        {
            result += rand();
        }

        return result;
    }

    inline cg::point_2 random_point(double minx, double maxx, double miny, double maxy)
    {
        uniform_random_real<double> xdomain(minx, maxx);
        uniform_random_real<double> ydomain(miny, maxy);

        cg::point_2 res;
        xdomain >> res.x;
        ydomain >> res.y;

        return res;
    }

    inline std::vector<cg::point_2> uniform_points(size_t count)
    {
        std::vector<cg::point_2> res(count);

        for (size_t l = 0; l != count; ++l)
        {
            res[l] = random_point(-200, 200, -200, 200);
        }

        return res;
    }
}
