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
        // derivative from D1 to D_{N-1}
        Eigen::MatrixXd dDdX_;
        // 3 * d(x_{i+1} - x_i)/dX
        Eigen::MatrixXd tmpM_;

    public:
        inline void setConditions(const Eigen::Vector2d &headPos,
                                  const Eigen::Vector2d &tailPos,
                                  const int &pieceNum)
        {
            // set x0(s0),x0(s1)...xn(s0),xn(s1) points
            headP = headPos;
            tailP = tailPos;
            N = pieceNum;

            // dD/dX = inverse(tmpA) @ dB/dX
            Eigen::MatrixXd tmpA;
            Eigen::MatrixXd dBdX;
            tmpA.resize(N-1, N-1);
            dBdX.resize(N-1, N-1);
            // tmpM = 3 * d(x_i - x_{i+1})/dX
            tmpM_.resize(N, N-1);
            for (size_t i = 0; i < N - 1; i++)
            {
                tmpA.row(i) = Eigen::VectorXd(N-1);
                dBdX.row(i) = Eigen::VectorXd(N-1);
                if (i != 0)
                {
                    tmpA(i, i - 1) = 1;
                    dBdX(i, i - 1) = -3;
                    tmpM_(i, i - 1) = 3;
                }
                tmpA(i, i) = 4;
                tmpM_(i, i) = -3;
                if (i != N - 2)
                {
                    tmpA(i, i + 1) = 1;
                    dBdX(i, i + 1) = 3;
                }
            }
            tmpM_(N-1, N-2) = 3;

            // std::cout << "tmpA: \n" << tmpA << std::endl;
            // std::cout << "dBdX: \n" << dBdX << std::endl;

            dDdX_ = tmpA.inverse() * dBdX;
            // std::cout << "dDdX_: \n" << dDdX_ << std::endl;
            // std::cout << "tmpM_: \n" << tmpM_ << std::endl;
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
                p_s1 = inPs.col(i);

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
                // std::cout << "coeffMat: \n" << coeffMat << std::endl;

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
            Eigen::Matrix2Xd grad(2, N-1);
            Eigen::Matrix2Xd dEidX(2, N-1);
            gradByPoints = Eigen::Matrix2Xd::Zero(2, N-1);

            // dEi/dX = cf @ partialX
            for (size_t i = 0; i < N; i++)
            {
                // size_t i = 2;
                // get coeffMat of curve i
                auto piece = (*curve_)[i];
                Eigen::Matrix<double, 2, 4> coeffMat = piece.getCoeffMat();
                // std::cout << "cMat: " << coeffMat << std::endl;
                // cf = [8ci, 12ci, 12di, 24di]
                Eigen::MatrixXd cf_x(1, 4);
                Eigen::MatrixXd cf_y(1, 4);
                cf_x << 8 * coeffMat(0, 1), 12 * coeffMat(0, 1), 12 * coeffMat(0, 0), 24 * coeffMat(0, 0);
                cf_y << 8 * coeffMat(1, 1), 12 * coeffMat(1, 1), 12 * coeffMat(1, 0), 24 * coeffMat(1, 0);
                // std::cout << "cf_x: \n" << cf_x << std::endl;

                // partialX_mat = [dci/dX, ddi/dX, dci/dX, ddi/dX]^T
                // j = i+1
                // dci/dX = 3 d(x_j - x_i)/dX + 2 dDi/dX + dD_j/dX
                Eigen::VectorXd dcidX;
                Eigen::VectorXd dcidX_part1 = -tmpM_.row(i);

                Eigen::VectorXd dDidX(N-1);
                Eigen::VectorXd dDjdX(N-1);

                if (i != 0)
                {
                    dDidX = dDdX_.row(i - 1);
                }
                if (i != N - 1){
                    dDjdX = dDdX_.row(i);
                }
                dcidX = dcidX_part1 - 2 * dDidX - dDjdX;
                // std::cout << "dcidX: \n" << dcidX << std::endl;

                // ddi/dX = 2 d(x_i - x_j)/dX + dDi/dX + dD_j/dX
                Eigen::VectorXd ddidX;
                Eigen::VectorXd ddidX_part1 = tmpM_.row(i);

                ddidX = ddidX_part1 + dDidX + dDjdX;
                // std::cout << "ddidX: \n" << ddidX << std::endl;

                Eigen::MatrixXd partialX_mat(4, N-1);
                partialX_mat.row(0) = dcidX;
                partialX_mat.row(1) = ddidX;
                partialX_mat.row(2) = dcidX;
                partialX_mat.row(3) = ddidX;
                // std::cout << "partialX_mat: \n" << partialX_mat << std::endl;

                grad.row(0) = cf_x * partialX_mat;
                grad.row(1) = cf_y * partialX_mat;
                
                // std::cout << "grad: \n" << grad << std::endl;
                gradByPoints = gradByPoints + grad;
            }
        }
    };
}

#endif
