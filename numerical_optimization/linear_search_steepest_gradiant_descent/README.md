# Linear-search Steepest Gradiant Descent

## Rosenbrock Function
In mathematical optimization, the Rosenbrock function is a non-convex function, introduced by Howard H. Rosenbrock in 1960, which is used as a performance test problem for optimization algorithms. It is also known as Rosenbrock's valley or Rosenbrock's banana function.

$$ f(x)=f(x_1,x_2...x_N)=\sum_{i=1}^{N/2}[100(x_{2i-1}^2-x_{2i})^2+(x_{2i-1}-1)^2] $$

<img src="/images/Rosenbrock-contour.svg" width="50%">

The global minimum of the function is $(1,1,...)_{N+1} $

## Gradiant
$$
\triangledown f(x)=
[\frac{\partial f}{\partial x_1},  \frac{\partial f}{\partial x_2}...\frac{\partial f}{\partial x_N}]^T
$$
$$
=\begin{bmatrix} 100(4x_1^3-4x_2x_1)+2(x_1-1) \\
2x_2-2x_1^2 \\
100(4x_3^3-4x_4x_3)+2(x_3-1) \\
2x_4-2x_3^2 \\
...\end{bmatrix}
$$

## Hessian
$$
\triangledown^2 f(x)=\begin{bmatrix}
 \frac{\partial^2 f}{\partial x_1^2} &  \frac{\partial^2 f}{\partial x_1x_2}&...&\frac{\partial^2 f}{\partial x_1x_N}  \\
 \frac{\partial^2 f}{\partial x_2x_1} &  \frac{\partial^2 f}{\partial x_2^2}&...&\frac{\partial^2 f}{\partial x_2x_N} \\
 \frac{\partial^2 f}{\partial x_Nx_1} &  \frac{\partial^2 f}{\partial x_Nx_2}&...&\frac{\partial^2 f}{\partial x_N^2}
\end{bmatrix}
$$
$$
=\begin{bmatrix}
100(12x_1^2-4x_2)+2 & 100(-4x_1) & 0 & 0 &... \\
-4x_1& 2 & 0 & 0 & ... \\
0&0&100(12x_3^2-4x_4)+2 & 100(-4x_3) &... \\
0 & 0 &-4x_3& 2 & ... \\
... & ... & ... & ... & ...
\end{bmatrix}
$$