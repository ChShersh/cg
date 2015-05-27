#pragma once

#include <cg/pca/weights.h>
#include <boost/optional.hpp>

#include <set>
#include <cmath>
#include <memory>
#include <functional>

namespace cg
{
    struct eigen
    {
        point_2 dir;
        double val;
    };

    enum weight_type
    {
        standard_pca = 1,
        manhattan_pca = 2,
        euclidian_pca = 3,
        mahalanobis_pca = 4
    };

    struct weighted_pca
    {
        weighted_pca()
        {}

        void set_weight(const weight_t & weight_func)
        {
            w = weight_func;
        }

        std::pair<point_2, point_2> find_principal_component(const std::set<point_2> & points) const
        {
            double w_norm = 0;
            std::vector<double> mean(2);
            std::vector<std::vector<double>> cov_matrix(2, std::vector<double>(2));

            for (point_2 p : points) {
                w_norm += w(p);
                mean[0] += p.x * w(p);
                mean[1] += p.y * w(p);
            }
            mean[0] /= w_norm;
            mean[1] /= w_norm;

            for (int i = 0; i < 2; i++) {
                for (int j = 0; j < 2; j++) {
                    for (point_2 p : points) {
                        std::vector<double> point{p.x, p.y};
                        cov_matrix[i][j] += w(p) * (mean[i] - point[i]) * (mean[j] - point[j]);
                    }
                    cov_matrix[i][j] /= w_norm;
                }
            }

            std::vector<eigen> eigens = find_eigen(cov_matrix);
            eigen principal_eigen = ((eigens[0].val >= eigens[1].val) ? eigens[0] : eigens[1]);

            return std::make_pair(principal_eigen.dir, point_2{mean[0], mean[1]});
        }

        std::vector<eigen> find_eigen(const std::vector<std::vector<double>> & m) const
        {
            const double tolerance = 0.1e-20;
            double A = m[0][0], B = m[0][1], C = m[1][0], D = m[1][1];
            eigen eig1, eig2;

            if (B * C <= tolerance  ) {
                eig1.val = A; eig1.dir.x = 1; eig1.dir.y = 0;
                eig2.val = D; eig2.dir.x = 0; eig2.dir.y = 1;
                return std::vector<eigen>{eig1, eig2};
            }

            double tr = A + D;
            double det = A * D - B * C;
            double S = std::sqrt((tr / 2) * (tr / 2) - det);
            eig1.val = tr / 2 + S;
            eig2.val = tr / 2 - S;

            double SS = sqrt(std::max((A - D) * (A - D) / 4 + B * C, 0.0));
            if( A - D < 0 ) {
                eig1.dir.x = C;
                eig1.dir.y = - (A - D) / 2 + SS;
                eig2.dir.x = + (A - D) / 2 - SS;
                eig2.dir.y = B;
            } else {
                eig2.dir.x = C;
                eig2.dir.y = - (A - D) / 2 - SS;
                eig1.dir.x = + (A - D) / 2 + SS;
                eig1.dir.y = B;
            }

            double n1 = sqrt(eig1.dir.x * eig1.dir.x + eig1.dir.y * eig1.dir.y);
            eig1.dir.x /= n1;
            eig1.dir.y /= n1;
            double n2 = sqrt(eig2.dir.x * eig2.dir.x + eig2.dir.y * eig2.dir.y);
            eig2.dir.x /= n2;
            eig2.dir.y /= n2;

            return std::vector<eigen>{eig1, eig2};
        }

    private:
        weight_t w;
    };
}
