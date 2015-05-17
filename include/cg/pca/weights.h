#pragma once

#include <cg/primitives/point.h>

#include <set>
#include <vector>
#include <cmath>

namespace cg
{
    struct weight
    {
        virtual double operator() (point_2 const & p) const = 0;
    };

    struct unit_weight : weight
    {
        double operator() (point_2 const & p) const
        {
            return 1;
        }
    };

    struct manhattan_weight : weight
    {
        double operator() (point_2 const & p) const
        {
            return std::abs(p.x) + std::abs(p.y);
        }
    };

    struct euclidian_weight : weight
    {
        double operator() (point_2 const & p) const
        {
            return sqrt(p.x * p.x + p.y * p.y); // + 1
        }
    };

    struct mahalanobis_weight : weight
    {
        mahalanobis_weight() {}

        mahalanobis_weight(const std::set<point_2> & points)
        {
            std::vector<double> mean(2);
            std::vector<std::vector<double>> s(2, std::vector<double>(2));

            for (point_2 p : points) {
                mean[0] += p.x;
                mean[1] += p.y;
            }
            mean[0] /= points.size() - 1;
            mean[1] /= points.size() - 1;

            for (int i = 0; i < 2; i++) {
                for (int j = 0; j < 2; j++) {
                    for (point_2 p : points) {
                        std::vector<double> point{p.x, p.y};
                        s[i][j] += (mean[i] - point[i]) * (mean[j] - point[j]);
                    }
                    s[i][j] /= points.size() - 1;
                }
            }

            cov_matrix = s;
            double det = s[0][0] * s[1][1] - s[0][1] * s[1][0];
            cov_matrix[0][0] = s[1][1] / det;
            cov_matrix[0][1] = -s[0][1] / det;
            cov_matrix[1][0] = -s[1][0] / det;
            cov_matrix[1][1] = s[0][0] / det;
        }

        double operator() (point_2 const & p) const
        {
            std::vector<double> point{p.x, p.y};
            std::vector<double> res(2);
            for (int i = 0; i < 2; i++) {
                for (int j = 0; j < 2; j++) {
                    res[i] += point[j] * cov_matrix[i][j];
                }
            }
            return std::sqrt(res[0] * p.x + res[1] * p.y);
        }

    private:
        std::vector<std::vector<double>> cov_matrix;
    };
}
