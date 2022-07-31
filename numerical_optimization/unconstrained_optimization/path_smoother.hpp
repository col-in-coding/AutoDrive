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
                                          Eigen::VectorXd &g)
        {
           //TODO
            double cost = 0;
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
            
            // int c = 0.5;
            // double cost;

            // for (size_t i = 0; i < maxItr_; i++)
            // {

            cubSpline.getGrad(gradByPoints);
            std::cout << "gradByPoints: \n" << gradByPoints << std::endl;

            // f(x) = stretchEnergy + potential
            double energy;
            double potential;
            cubSpline.getStretchEnergy(energy);
            // std::cout << "energy: " << energy << std::endl;
            for (size_t j = 0; j < diskObstacles.cols(); j++)
            {
                // for j-th obstable
                Eigen::VectorXd obstacle = diskObstacles.col(j);
                Eigen::Matrix2Xd grad(2, pieceN-1);
                for (size_t k = 0; k < pieceN - 1; k++)
                {
                    // for k-th inner point
                    double ox = obstacle[0];
                    double oy = obstacle[1];
                    double r = obstacle[2];
                    double x = iniInPs(0, k);
                    double y = iniInPs(1, k);
                    double dist = std::sqrt((ox - x)*(ox - x) + (oy - y)*(ox - y));
                    std::cout << "dist: " << dist << std::endl;
                    if (dist < r)
                    {
                        potential += penaltyWeight * (r - dist);
                        grad(0, k) = penaltyWeight * (ox - x) / (dist + 1e-6);
                        grad(1, k) = penaltyWeight * (oy - y) / (dist + 1e-6);
                    }
                }
                std::cout << "potential grad: \n" << grad << std::endl;
                gradByPoints = gradByPoints + grad;
            }
            double fx = energy + potential;

            // Amijo condition
            int tau = 1;

            // double cost = costFunction();

            // }


            double minCost = INFINITY;
            return minCost;
        }
    };

}

#endif
