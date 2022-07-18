import argparse
import numpy as np


class RosenBroke:
    def __init__(self, N) -> None:
        self.N = N
        self.x1 = np.zeros(N)
        self.x2 = np.zeros(N)
    
    def fx(self):
        z = 100 * np.square(np.square(self.x1) - self.x2) + np.square(self.x1 - 1)
        return z.sum()

    def get_x(self):
        return np.stack((self.x1, self.x2), axis=1).reshape(2*self.N)

    def get_gradiant_T(self):
        odd_vec = 100 * (4 * np.power(self.x1, 3) - 4 * self.x2 * self.x1) + 2 * (self.x1 - 1)
        even_vec = 200 * self.x2 - 200 * np.square(self.x1)
        return np.stack((odd_vec, even_vec), axis=1).reshape(2*self.N)

    def get_hessian(self):
        A = 100 * (12 * np.square(self.x1) - 4 * self.x2) + 2
        B = 100 * (-4 * self.x1)
        C = -4 * self.x1
        D = 2 * np.ones(self.N)
        H = np.zeros((2 * self.N, 2 * self.N))
        for i in range(self.N):
            H[2*i, 2*i] = A[i]
            H[2*i + 1, 2*i + 1] = D[i]
            H[2*i, 2*i + 1] = B[i]
            H[2*i + 1, 2*i] = C[i]
        return H

    def set_x(self, x):
        self.x1, self.x2 = x.reshape(self.N, 2).transpose()
        return self


class LinearSearchSGDSolver:
    def __init__(self, funcObj, c=0.5, max_itr=1000) -> None:
        self.funcObj = funcObj
        self.max_itr = max_itr
        self.c = c

    def solve(self):
        for _ in range(self.max_itr):
            tau = 1
            G_T = self.funcObj.get_gradiant_T()
            # print("G_T: ", G_T)
            d_T = - self._norm(G_T)
            # d_T = - G_T
            x = self.funcObj.get_x()
            # print("x: ", x)
            fx = self.funcObj.fx()
            # print("fx: ", fx)
            while tau > 1e-6:
                x_next = x + tau * d_T
                # print("x_next: ", x_next)
                fx_next = self.funcObj.set_x(x_next).fx()
                # print("fx_next: ", fx_next)
                impro_factor = - self.c * tau * np.dot(d_T, G_T)
                # print("if: ", impro_factor)
                if fx - fx_next >= impro_factor:
                    break
                tau = tau / 2
            print("===> ", x)

            if impro_factor < 1e-16:
                break

        return self.funcObj.get_x()

    def _norm(self, x):
        return x / np.sqrt(np.sum(x**2))


def main(args):
    funcObj = RosenBroke(args.N)
    # print(funcObj.get_x())
    # print(funcObj.get_gradiant())
    # print(funcObj.get_hessian())

    solver = LinearSearchSGDSolver(funcObj, max_itr=10000)
    result = solver.solve()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-N', type=int, default=1)
    args = parser.parse_args()
    main(args)