#ifndef PATH_SMOOTHER_HPP
#define PATH_SMOOTHER_HPP

#include "cubic_spline.hpp"
#include "lbfgs.hpp"

#include <Eigen/Eigen>

#include <cmath>
#include <cfloat>
#include <iostream>
#include <vector>
#include <algorithm>

namespace path_smoother
{

    class PathSmoother
    {
    private:
        cubic_spline::CubicSpline cubSpline;
        const int maxItr_ = 3000;

        int pieceN;
        Eigen::Matrix3Xd diskObstacles;
        double penaltyWeight;
        Eigen::Vector2d headP;
        Eigen::Vector2d tailP;
        Eigen::Matrix2Xd points;
        Eigen::Matrix2Xd gradByPoints;

        lbfgs::lbfgs_parameter_t lbfgs_params;

    private:
        static inline double costFunction(void *ptr,
                                          const Eigen::VectorXd &x,
                                          Eigen::VectorXd &tau_d)
        {
            // TODO
            double cost = INFINITY;
            int n_minus_one = x.size() / 2;

            Eigen::VectorXd x_next = x + tau_d;
            Eigen::Matrix2Xd inPs(Eigen::Map<Eigen::Matrix2Xd>(x_next.data(), 2, n_minus_one));

            cubic_spline::CubicSpline* cubSpline = &((PathSmoother *) ptr)->cubSpline;
            Eigen::Matrix3Xd* diskObstacles = &((PathSmoother *) ptr)->diskObstacles;
            double penaltyWeight = ((PathSmoother *) ptr)->penaltyWeight;
            cubSpline->setInnerPoints(inPs);

            /*
            // f(x) = stretchEnergy + potential
            double energy;
            double potential;
            cubSpline->getStretchEnergy(energy);

            for (size_t j = 0; j < diskObstacles->cols(); j++)
            {
                // for j-th obstable
                Eigen::VectorXd obstacle = diskObstacles->col(j);
                for (size_t k = 0; k < n_minus_one; k++)
                {
                    // for k-th inner point
                    double ox = obstacle[0];
                    double oy = obstacle[1];
                    double r = obstacle[2];
                    double x = inPs(0, k);
                    double y = inPs(1, k);
                    double dist = std::sqrt((ox - x) * (ox - x) + (oy - y) * (ox - y));
                    // std::cout << "dist: " << dist << std::endl;
                    if (dist < r)
                    {
                        potential += penaltyWeight * (r - dist);
                    }
                }
            }

            cost = energy + potential;
            */

            return cost;
        }

    public:
        inline bool setup(const Eigen::Vector2d &initialP,
                          const Eigen::Vector2d &terminalP,
                          const int &pieceNum,
                          const Eigen::Matrix3Xd &diskObs,
                          const double penaWeight)
        {
            pieceN = pieceNum;
            diskObstacles = diskObs;
            penaltyWeight = penaWeight;
            headP = initialP;
            tailP = terminalP;

            cubSpline.setConditions(headP, tailP, pieceN);

            points.resize(2, pieceN - 1);
            gradByPoints.resize(2, pieceN - 1);

            return true;
        }

        inline double optimize(CubicCurve &curve,
                               const Eigen::Matrix2Xd &iniInPs,
                               const double &relCostTol)
        {
            // init curve
            cubSpline.setInnerPoints(iniInPs);
            cubSpline.getCurve(curve);

            double c = 0.5;

            Eigen::VectorXd tau_d(2 * (pieceN - 1));
            Eigen::Matrix2Xd x_2d = iniInPs;
            Eigen::VectorXd x(Eigen::Map<Eigen::VectorXd>(x_2d.data(), iniInPs.cols() * iniInPs.rows()));
            double f_x = costFunction(this, x, tau_d);
            int step = 0;

            for (size_t i = 0; i < maxItr_; i++)
            {
                // Caculate gradient
                cubSpline.getGrad(gradByPoints);
                // std::cout << "gradByPoints: \n"
                //         << gradByPoints << std::endl;
                for (size_t j = 0; j < diskObstacles.cols(); j++)
                {
                    // for j-th obstable
                    Eigen::VectorXd obstacle = diskObstacles.col(j);
                    Eigen::Matrix2Xd grad(2, pieceN - 1);
                    for (size_t k = 0; k < pieceN - 1; k++)
                    {
                        // for k-th inner point
                        double ox = obstacle[0];
                        double oy = obstacle[1];
                        double r = obstacle[2];
                        double x = iniInPs(0, k);
                        double y = iniInPs(1, k);
                        double dist = std::sqrt((ox - x) * (ox - x) + (oy - y) * (ox - y));
                        std::cout << "dist: " << dist << std::endl;
                        if (dist < r)
                        {
                            // potential += penaltyWeight * (r - dist);
                            grad(0, k) = penaltyWeight * (ox - x) / (dist + 1e-6);
                            grad(1, k) = penaltyWeight * (oy - y) / (dist + 1e-6);
                        }
                    }
                    gradByPoints = gradByPoints + grad;
                }

                // Amijo condition
                double tau = 1.0;
                Eigen::VectorXd g(Eigen::Map<Eigen::VectorXd>(gradByPoints.data(), gradByPoints.cols() * gradByPoints.rows()));
                Eigen::VectorXd d = g.normalized();
                // std::cout << "d: " << d << std::endl;
                
                double f_x_prime;
                double impro_factor;
                do {
                    tau_d = - tau * d;
                    // std::cout << "tau_d: " << tau_d << std::endl;

                    f_x_prime = costFunction(this, x, tau_d);
                    double dot_prod = tau_d.transpose() * g;
                    double impro_factor = - c * tau * dot_prod;
                    tau = tau / 2;
                    std::cout << "delta f: " << f_x - f_x_prime << std::endl;
                    // std::cout << "impro_factor: " << impro_factor << std::endl;
                } while(f_x - f_x_prime < impro_factor);

                f_x = f_x_prime;
                step++;
                std::cout << "step: " << step << std::endl;
                std::cout << "cost: " << f_x << std::endl;
                break;
                // if (gradByPoints.norm())
            }
            return f_x;
        }
    };

}

#endif
