# 梯度方程求解

待解的式子：

```math
u - a\nabla^2u = f
```
## 弱形式：

用一个试函数 ``\varphi`` 乘到等号两边，然后积分：

```math
\int_{\Omega}\varphi u - a\int_{\Omega}\varphi \nabla^2 u = \int_{\Omega}\varphi f
```
其中，

```math
\int_{\Omega}\varphi \nabla^2 u = \int_{\partial\Omega}\varphi \mathbf{n}\cdot \nabla u
 - \int_{\Omega}\nabla \varphi\cdot \nabla u
```
边界条件不管是 ``u|_{\partial\Omega} = 0`` 还是
``\nabla u|_{\partial\Omega} = 0`` ,上式在边界上的积分都为0。
此时弱形式化为：
```math
\int_{\Omega}\varphi u + a\int_{\Omega}\nabla \varphi\cdot \nabla u = \int_{\Omega}\varphi f
```
## 离散

离散成网格，在每个单元中，节点上的u值，记为(以三角形网格为例)

```math
u^e =
\begin{bmatrix}
u_1\\
u_2\\
u_3
\end{bmatrix}
```
单元内任一点A的u可用形函数 ``N_i`` (i = 1,2,3) 近似. 记 N
为某点A的形函数的值：

```math
N =
\begin{bmatrix}
N_1(A)\\
N_2(A)\\
N_3(A)
\end{bmatrix}
```
则点A处，u的近似值为：

```math
\hat{u_h} \approx N^T * u^e
```
``\nabla u`` 的近似值：

```math
\nabla u \approx \nabla N^T * u^e=
\begin{bmatrix}
\frac{\partial}{\partial x} \\
\frac{\partial}{\partial y}
\end{bmatrix}
\begin{bmatrix}
N_1&N_2&N_3
\end{bmatrix}
\begin{bmatrix}
u_1\\
u_2\\
u_3
\end{bmatrix}
=
\begin{bmatrix}
\frac{\partial N_1}{\partial x} & \frac{\partial N_2}{\partial x} & \frac{\partial N_3}{\partial x} \\
\frac{\partial N_1}{\partial y} & \frac{\partial N_2}{\partial y} & \frac{\partial N_3}{\partial y}
\end{bmatrix}
\begin{bmatrix}
u_1\\
u_2\\
u_3
\end{bmatrix}
```
记 ``\nabla N`` =
```math
\begin{bmatrix}
\frac{\partial N_1}{\partial x} &
\frac{\partial N_1}{\partial y}\\
\frac{\partial N_2}{\partial x} &
\frac{\partial N_2}{\partial y}\\
\frac{\partial N_3}{\partial x} &
\frac{\partial N_3}{\partial y}
\end{bmatrix}
```
为 ``B`` ：
弱形式在任一单元内的积分近似为：

```math
\int_{\Omega h}\varphi N^T* u^e + a\int_{\Omega h}\nabla \varphi\cdot B^T* u^e = \int_{\Omega h}\varphi f
```
这个方程有三个未知量，即三个节点的u值，所以试函数 ``\varphi`` 需要取三个。

令他们分别为三个形函数， ``\varphi_i = N_i`` ,方程组变为：

```math
\int_{\Omega h}N\cdot N^T* u^e + a\int_{\Omega h}B\cdot B^T* u^e = \int_{\Omega h}N* f
```
## 数值积分

上式积分还是太难算，用高斯积分法(x表示高斯积分点坐标，nx为单元内积分点个数)：

```math
\sum_{x = 1}^{nx} \omega \left[N(x)\cdot N(x)^T + aB(x)\cdot B(x)^T\right] u^e = \sum_{x = 1}^{nx}\omega N(x) * f(x)
```
上式记为：

```math
\begin{gathered}
Ku^e = F \\
K = \sum_{x = 1}^{nx} \omega (N(x)\cdot N(x)^T + aB(x)\cdot B(x)^T)\\
F = \sum_{x = 1}^{nx}\omega N(x) * f(x)
\end{gathered}
```
每个单元刚度矩阵求出来再组装成整体刚度矩阵
