#ifndef CUBIC_SPLINE_HPP
#define CUBIC_SPLINE_HPP

#include "cubic_curve.hpp"

#include <Eigen/Eigen>

#include <cmath>
#include <vector>
#include <memory>
#include <iostream>

namespace cubic_spline
{

    // The banded system class is used for solving
    // banded linear system Ax=b efficiently.
    // A is an N*N band matrix with lower band width lowerBw
    // and upper band width upperBw.
    // Banded LU factorization has O(N) time complexity.
    class BandedSystem
    {
    public:
        // The size of A, as well as the lower/upper
        // banded width p/q are needed
        inline void create(const int &n, const int &p, const int &q)
        {
            // In case of re-creating before destroying
            destroy();
            N = n;
            lowerBw = p;
            upperBw = q;
            int actualSize = N * (lowerBw + upperBw + 1);
            ptrData = new double[actualSize];
            std::fill_n(ptrData, actualSize, 0.0);
            return;
        }

        inline void destroy()
        {
            if (ptrData != nullptr)
            {
                delete[] ptrData;
                ptrData = nullptr;
            }
            return;
        }

    private:
        int N;
        int lowerBw;
        int upperBw;
        // Compulsory nullptr initialization here
        double *ptrData = nullptr;

    public:
        // Reset the matrix to zero
        inline void reset(void)
        {
            std::fill_n(ptrData, N * (lowerBw + upperBw + 1), 0.0);
            return;
        }

        // The band matrix is stored as suggested in "Matrix Computation"
        inline const double &operator()(const int &i, const int &j) const
        {
            return ptrData[(i - j + upperBw) * N + j];
        }

        inline double &operator()(const int &i, const int &j)
        {
            return ptrData[(i - j + upperBw) * N + j];
        }

        // This function conducts banded LU factorization in place
        // Note that NO PIVOT is applied on the matrix "A" for efficiency!!!
        inline void factorizeLU()
        {
            int iM, jM;
            double cVl;
            for (int k = 0; k <= N - 2; ++k)
            {
                iM = std::min(k + lowerBw, N - 1);
                cVl = operator()(k, k);
                for (int i = k + 1; i <= iM; ++i)
                {
                    if (operator()(i, k) != 0.0)
                    {
                        operator()(i, k) /= cVl;
                    }
                }
                jM = std::min(k + upperBw, N - 1);
                for (int j = k + 1; j <= jM; ++j)
                {
                    cVl = operator()(k, j);
                    if (cVl != 0.0)
                    {
                        for (int i = k + 1; i <= iM; ++i)
                        {
                            if (operator()(i, k) != 0.0)
                            {
                                operator()(i, j) -= operator()(i, k) * cVl;
                            }
                        }
                    }
                }
            }
            return;
        }

        // This function solves Ax=b, then stores x in b
        // The input b is required to be N*m, i.e.,
        // m vectors to be solved.
        template <typename EIGENMAT>
        inline void solve(EIGENMAT &b) const
        {
            int iM;
            for (int j = 0; j <= N - 1; ++j)
            {
                iM = std::min(j + lowerBw, N - 1);
                for (int i = j + 1; i <= iM; ++i)
                {
                    if (operator()(i, j) != 0.0)
                    {
                        b.row(i) -= operator()(i, j) * b.row(j);
                    }
                }
            }
            for (int j = N - 1; j >= 0; --j)
            {
                b.row(j) /= operator()(j, j);
                iM = std::max(0, j - upperBw);
                for (int i = iM; i <= j - 1; ++i)
                {
                    if (operator()(i, j) != 0.0)
                    {
                        b.row(i) -= operator()(i, j) * b.row(j);
                    }
                }
            }
            return;
        }

        // This function solves ATx=b, then stores x in b
        // The input b is required to be N*m, i.e.,
        // m vectors to be solved.
        template <typename EIGENMAT>
        inline void solveAdj(EIGENMAT &b) const
        {
            int iM;
            for (int j = 0; j <= N - 1; ++j)
            {
                b.row(j) /= operator()(j, j);
                iM = std::min(j + upperBw, N - 1);
                for (int i = j + 1; i <= iM; ++i)
                {
                    if (operator()(j, i) != 0.0)
                    {
                        b.row(i) -= operator()(j, i) * b.row(j);
                    }
                }
            }
            for (int j = N - 1; j >= 0; --j)
            {
                iM = std::max(0, j - lowerBw);
                for (int i = iM; i <= j - 1; ++i)
                {
                    if (operator()(j, i) != 0.0)
                    {
                        b.row(i) -= operator()(j, i) * b.row(j);
                    }
                }
            }
        }
    };

    class CubicSpline
    {
    public:
        CubicSpline() = default;
        ~CubicSpline() { A.destroy(); }

    private:
        int N;
        Eigen::Vector2d headP;
        Eigen::Vector2d tailP;
        BandedSystem A;
        Eigen::MatrixX2d b;

        std::unique_ptr<CubicCurve> curve_;
        const int duration_ = 2;
        Eigen::MatrixXd tmpA_;

    public:
        inline void setConditions(const Eigen::Vector2d &headPos,
                                  const Eigen::Vector2d &tailPos,
                                  const int &pieceNum)
        {
            // set x0(s0),x0(s1)...xn(s0),xn(s1) points
            headP = headPos;
            tailP = tailPos;
            N = pieceNum;
            // init tmpA matrix
            tmpA_.resize(N, N);
            for (size_t i = 0; i < N; i++)
            {
                tmpA_.row(i) = Eigen::VectorXd(N);
                if (i != 0)
                {
                    tmpA_(i, i - 1) = 1;
                }
                tmpA_(i, i) = 4;
                if (i != N - 1)
                {
                    tmpA_(i, i + 1) = 1;
                }
            }
            return;
        }

        inline void setInnerPoints(const Eigen::Ref<const Eigen::Matrix2Xd> &inPs)
        {
            // Init curve
            std::vector<double> durs;
            std::vector<Eigen::Matrix<double, 2, 4>> cMats;
            cMats.resize(N);

            Eigen::Vector2d p_s0 = headP;
            Eigen::Vector2d p_s1 = inPs.col(0);
            cMats[0] << 0.0, 0.0, p_s1(0) - headP(0), headP(0),
                        0.0, 0.0, p_s1(1) - headP(1), headP(1);
            durs.emplace_back(duration_);

            for (size_t i = 1; i < N - 1; i++)
            {
                p_s0 = p_s1;
                p_s1 = inPs.col(0);

                cMats[i] << 0.0, 0.0, p_s1(0) - p_s0(0), p_s0(0),
                            0.0, 0.0, p_s1(1) - p_s0(1), p_s0(1);
                durs.emplace_back(duration_);
            }

            p_s0 = p_s1;

            cMats[N-1] << 0.0, 0.0, tailP(0) - p_s0(0), p_s0(0),
                          0.0, 0.0, tailP(1) - p_s0(1), p_s0(1);
            durs.emplace_back(duration_);

            curve_ = std::make_unique<CubicCurve>(durs, cMats);
            return;
        }

        inline void getCurve(CubicCurve &curve) const
        {
            //TODO
            curve = *curve_.get();
            return;
        }

        inline void getStretchEnergy(double &energy) const
        {
            //TODO
            energy = 0;
            for (size_t i = 0; i < curve_->getPieceNum(); i++)
            {
                auto piece = (*curve_)[i];
                Eigen::Matrix<double, 2, 4> coeffMat = piece.getCoeffMat();
                std::cout << coeffMat << std::endl;

                double ci_x = coeffMat(1, 0);
                double ci_y = coeffMat(1, 1);
                double di_x = coeffMat(0, 0);
                double di_y = coeffMat(0, 1);
                energy += 4 * (ci_x*ci_x + ci_y*ci_y) + 12 * (ci_x*di_x + ci_y*di_y) + 12 * (di_x*di_x + di_y*di_y);
            }
            return;
        }

        inline const Eigen::MatrixX2d &getCoeffs(void) const
        {
            return b;
        }

        // inline void getGrad(Eigen::Ref<Eigen::Matrix2Xd> gradByPoints) const
        inline void getGrad(Eigen::Ref<Eigen::Matrix2Xd> gradByPoints)
        {
            //TODO
            Eigen::VectorXd X(2 * (N + 1));
            std::cout << "X: " << X << std::endl;

            // dD/dX = inverse(A) @ dB/dX
            std::cout << "tmpA" << tmpA_ << std::endl;

            /*
            int i = 0;
            // dEi/dX = partialX @ cf
            auto piece = (*curve_)[i];
            Eigen::Matrix<double, 2, 4> coeffMat = piece.getCoeffMat();
            // cf = [8xci, 8yci, 12xci, 12yci, 12xdi, 12ydi, 24xdi, 24ydi]
            Eigen::VectorXd cf;
            cf << 8 * coeffMat(0, 1), 8 * coeffMat(1, 1),
                  12 * coeffMat(0, 1), 12 * coeffMat(1, 1),
                  12 * coeffMat(0, 0), 12 * coeffMat(1, 0),
                  24 * coeffMat(0, 0), 12 * coeffMat(1, 0);
            // partialX = [dxci/dX, dyci/dX, dxdi/dX, dydi/dX, dxci/dX, dyci/dX, dxdi/dX, dydi/dX]
            // dxci/dX = 3d(x_{i+1} - x_i)/dX + 2 dDi/dX + dD_{i+1}/dX
            */




        }
    };
}

#endif
