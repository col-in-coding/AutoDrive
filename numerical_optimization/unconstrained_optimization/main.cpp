#include "path_smoother.hpp"
#include <iostream>


int main(int argc, char const *argv[])
{
    Eigen::Vector2d initialP {0, 0};
    Eigen::Vector2d terminalP {3, 0};
    int N {3};
    double relCostTol {1.0e-6};
    Eigen::Matrix2Xd innerPoints(2, N - 1);
    for (int i = 0; i < N - 1; ++i)
    {
        // innerPoints.col(i) = (initialP - terminalP) * (i + 1.0) / N + initialP;
        innerPoints.col(i) = (terminalP - initialP) * (i + 1.0) / N + initialP;
    }

    // FOR TEST
    // Eigen::Vector2d* 

    // (x, y, r)
    Eigen::Matrix3Xd circleObs(3,1);
    circleObs << 0.9, 0, 0.5;

    double penaltyWeight {100};

    path_smoother::PathSmoother pathSmoother;
    pathSmoother.setup(initialP, terminalP, N, circleObs, penaltyWeight);

    CubicCurve curve;
    if (std::isinf(pathSmoother.optimize(curve, innerPoints, relCostTol)))
    {
        std::cout << "Optimization Faild" << std::endl;
    } else {
        auto positions = curve.getPositions();
        std::cout << positions << std::endl;
    }

    return 0;
}
