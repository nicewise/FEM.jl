# 弹性理论

弹性问题的控制方程有三个：

```math
\begin{gathered}
  \nabla\cdot\mathbf{\sigma}+\mathbf{f}=0\\
  \mathbf{\varepsilon} = \frac{1}{2}(\nabla \mathbf{u}+\mathbf{u}\nabla)\\
  \mathbf{\sigma} = \mathbb{D}:\mathbf{\varepsilon}
\end{gathered}
```
此外，还有常见的边界条件：
- 位移边界条件: ``\mathbf{u}=\bar{\mathbf{u}}``
- 应力边界条件: ``\mathbf{\sigma}\cdot \mathbf{n}=\bar{\mathbf{f}}``

计算中，高阶张量可以降阶表示。 三维情况下，

```math
\begin{gathered}
\sigma =
\begin{bmatrix}
\sigma_x & \sigma_y & \sigma_z & \sqrt{2}\tau_{yz} & \sqrt{2}\tau_{zx} & \sqrt{2}\tau_{xy}
\end{bmatrix}^T \\
\varepsilon =
\begin{bmatrix}
\varepsilon_x & \varepsilon_y & \varepsilon_z & \sqrt{2}\varepsilon_{yz} & \sqrt{2}\varepsilon_{zx} & \sqrt{2}\varepsilon_{xy}
\end{bmatrix}^T
\end{gathered}
```
## 平衡微分方程：

```math
\begin{bmatrix}
\frac{\partial}{\partial x} & 0 & 0 & 0 & \frac{1}{\sqrt{2}}\frac{\partial}{\partial z} & \frac{1}{\sqrt{2}}\frac{\partial}{\partial y} \\
0 & \frac{\partial}{\partial y} & 0 & \frac{1}{\sqrt{2}}\frac{\partial}{\partial z} & 0 & \frac{1}{\sqrt{2}}\frac{\partial}{\partial x} \\
0 & 0 & \frac{\partial}{\partial z} & \frac{1}{\sqrt{2}}\frac{\partial}{\partial y} & \frac{1}{\sqrt{2}}\frac{\partial}{\partial x} & 0
\end{bmatrix}
\begin{bmatrix}
\sigma_x\\
\sigma_y\\
\sigma_z\\
\sqrt{2}\tau_{yz}\\
\sqrt{2}\tau_{zx}\\
\sqrt{2}\tau_{xy}
\end{bmatrix}
+
\begin{bmatrix}
f_x\\
f_y\\
f_z
\end{bmatrix} = \mathbf{0}
```
记为： ``L^T\sigma+\mathbf{f}=\mathbf{0}``

## 几何方程

```math
\begin{bmatrix}
\varepsilon_x\\
\varepsilon_y\\
\varepsilon_z\\
\sqrt{2}\varepsilon_{yz}\\
\sqrt{2}\varepsilon_{zx}\\
\sqrt{2}\varepsilon_{xy}
\end{bmatrix}
=
\begin{bmatrix}
\frac{\partial}{\partial x} & 0 & 0 \\
0 & \frac{\partial}{\partial y} & 0 \\
0 & 0 & \frac{\partial}{\partial z} \\
0 & \frac{1}{\sqrt{2}}\frac{\partial}{\partial z} & \frac{1}{\sqrt{2}}\frac{\partial}{\partial y} \\
\frac{1}{\sqrt{2}}\frac{\partial}{\partial z} & 0 & \frac{1}{\sqrt{2}}\frac{\partial}{\partial x} \\
\frac{1}{\sqrt{2}}\frac{\partial}{\partial y} & \frac{1}{\sqrt{2}}\frac{\partial}{\partial x} & 0
\end{bmatrix}
\begin{bmatrix}
u_x\\
u_y\\
u_z
\end{bmatrix}
```
记为： ``\varepsilon = L\mathbf{u}``

## 物理方程

在複雜應力狀態下，應力與應變有如下關係（運用小變形的假定）：

```math
\begin{gathered}
\left.\begin{matrix}
\varepsilon_x=\frac{1}{E}\left [ \sigma _x-\nu \left ( \sigma_y+\sigma_z \right ) \right ]\\
\varepsilon_y=\frac{1}{E}\left [ \sigma _y-\nu \left ( \sigma_z+\sigma_x \right ) \right ]\\
\varepsilon_z=\frac{1}{E}\left [ \sigma _z-\nu \left ( \sigma_x+\sigma_y \right ) \right ]\\
\gamma_{xy}=\frac{\tau_{xy}}{G}=\frac{2\left ( 1+\nu  \right )}{E}\tau_{xy}\\
\gamma_{yz}=\frac{\tau_{yz}}{G}=\frac{2\left ( 1+\nu  \right )}{E}\tau_{yz}\\
\gamma_{zx}=\frac{\tau_{zx}}{G}=\frac{2\left ( 1+\nu  \right )}{E}\tau_{zx}
\end{matrix}\right\} \\
\begin{pmatrix}
\varepsilon_x\\
\varepsilon_y\\
\varepsilon_z\\
\gamma_{xy}\\
\gamma_{yz}\\
\gamma_{zx}
\end{pmatrix}=\frac{1}{E}\begin{pmatrix}
1 & -\nu & -\nu & 0 & 0 & 0\\
-\nu & 1 & -\nu & 0 & 0 & 0\\
-\nu & -\nu & 1 & 0 & 0 & 0\\
0 & 0 & 0 & 2\left ( 1+\nu \right ) & 0 & 0\\
0 & 0 & 0 & 0 & 2\left ( 1+\nu \right ) & 0\\
0 & 0 & 0 & 0 & 0 & 2\left ( 1+\nu \right )
\end{pmatrix}\begin{pmatrix}
\sigma_x\\
\sigma_y\\
\sigma_z\\
\tau_{xy}\\
\tau_{yz}\\
\tau_{zx}
\end{pmatrix} \\
\begin{pmatrix}
\sigma_x\\
\sigma_y\\
\sigma_z\\
\tau_{xy}\\
\tau_{yz}\\
\tau_{zx}
\end{pmatrix}=\frac{E}{(1+\nu)(1-2\nu)}\begin{pmatrix}
1-\nu & \nu & \nu & 0 & 0 & 0\\
\nu & 1-\nu & \nu & 0 & 0 & 0\\
\nu & \nu & 1-\nu & 0 & 0 & 0\\
0 & 0 & 0 & \frac{1-2\nu}{2} & 0 & 0\\
0 & 0 & 0 & 0 & \frac{1-2\nu}{2} & 0\\
0 & 0 & 0 & 0 & 0 & \frac{1-2\nu}{2}
\end{pmatrix}\begin{pmatrix}
\varepsilon_x\\
\varepsilon_y\\
\varepsilon_z\\
\gamma_{xy}\\
\gamma_{yz}\\
\gamma_{zx}
\end{pmatrix}\\
\begin{pmatrix}
\sigma_x\\
\sigma_y\\
\sigma_z\\
\sqrt{2}\tau_{yz}\\
\sqrt{2}\tau_{zx}\\
\sqrt{2}\tau_{xy}
\end{pmatrix}=\frac{E}{(1+\nu)(1-2\nu)}\begin{pmatrix}
1-\nu & \nu & \nu & 0 & 0 & 0\\
\nu & 1-\nu & \nu & 0 & 0 & 0\\
\nu & \nu & 1-\nu & 0 & 0 & 0\\
0 & 0 & 0 & 1-2\nu& 0 & 0\\
0 & 0 & 0 & 0 & 1-2\nu & 0\\
0 & 0 & 0 & 0 & 0 & 1-2\nu
\end{pmatrix}\begin{pmatrix}
\varepsilon_x\\
\varepsilon_y\\
\varepsilon_z\\
\sqrt{2}\varepsilon_{yz}\\
\sqrt{2}\varepsilon_{zx}\\
\sqrt{2}\varepsilon_{xy}
\end{pmatrix}
\end{gathered}
```
```math
\mathbb{D} =
\frac{E}{(1+\nu)(1-2\nu)}\begin{pmatrix}
1-\nu & \nu & \nu & 0 & 0 & 0\\
\nu & 1-\nu & \nu & 0 & 0 & 0\\
\nu & \nu & 1-\nu & 0 & 0 & 0\\
0 & 0 & 0 & 1-2\nu& 0 & 0\\
0 & 0 & 0 & 0 & 1-2\nu & 0\\
0 & 0 & 0 & 0 & 0 & 1-2\nu
\end{pmatrix}
```
记为： ``\sigma = D\varepsilon``

## 按位移求解

弹性问题按位移求解的控制方程为：

```math
L^TDL\mathbf{u}+\mathbf{f}=\mathbf{0}
```
要找到一个u的表达式，在所有点满足这个偏微分方程组大部分情况下是做不到的

## 最小势能原理

应变能密度：
```math
v_{\varepsilon} = \frac{1}{2}\sigma:\varepsilon = \frac{1}{2}[\sigma]^T[\varepsilon]= \frac{1}{2}[\varepsilon]^T[\sigma]
```
系统的应变势能为：

```math
V_{\varepsilon} = \int_V v_{\varepsilon} dV = \frac{1}{2}\int_V \frac{1}{2}[\varepsilon]^T[\sigma] dV
 = \frac{1}{2}\int_V [\mathbf{u}]^TL^T[\mathbf{D}]L[\mathbf{u}] dV
```
系统总势能：

```math
E(\mathbf{u}) = V_{\varepsilon}
- \int_V[\mathbf{f}]^T[\mathbf{u}] dV
- \int_S[\bar{\mathbf{f}}]^T[\mathbf{u}] dS
```
推导有时间再补。 反正要找这么一个位移场，使得系统总势能是极小值：

```math
\frac{\partial E(u)}{\partial u} = \int_V L^T[\mathbf{D}]L[\mathbf{u}] dV
- \int_V[\mathbf{f}] dV
- \int_S[\bar{\mathbf{f}}] dS
= \mathbf{0}
```
这么一个 ``\mathbf{u}`` 还是很难找。

# 有限元

将物体所占区域划分为有限个网格，系统总势能这个积分式子变为在所有单元上的求和，
在单个单元里，点 (x, y, z) 位移用形函数 ``N_i(x, y, z)``
来近似，以四面体单元为例:

```math
\begin{bmatrix}
u_x\\
u_y\\
u_z
\end{bmatrix}=
\begin{bmatrix}
N_1&0&0&N_2&0&0&N_3&0&0&N_4&0&0\\
0&N_1&0&0&N_2&0&0&N_3&0&0&N_4&0\\
0&0&N_1&0&0&N_2&0&0&N_3&0&0&N_4
\end{bmatrix}
\begin{bmatrix}
u_{1x}\\
u_{1y}\\
u_{1z}\\
u_{2x}\\
u_{2y}\\
u_{2z}\\
u_{3x}\\
u_{3y}\\
u_{3z}\\
u_{4x}\\
u_{4y}\\
u_{4z}
\end{bmatrix}
```
上式记为： ``[\mathbf{u}] = [N][u^e]``

```math
\begin{bmatrix}
f_x\\
f_y\\
f_z
\end{bmatrix}=
\begin{bmatrix}
N_1&0&0&N_2&0&0&N_3&0&0&N_4&0&0\\
0&N_1&0&0&N_2&0&0&N_3&0&0&N_4&0\\
0&0&N_1&0&0&N_2&0&0&N_3&0&0&N_4
\end{bmatrix}
\begin{bmatrix}
f_{1x}\\
f_{1y}\\
f_{1z}\\
f_{2x}\\
f_{2y}\\
f_{2z}\\
f_{3x}\\
f_{3y}\\
f_{3z}\\
f_{4x}\\
f_{4y}\\
f_{4z}
\end{bmatrix}
```
上式记为： ``[\mathbf{f}] = [N][f^e]``

系统总势能泛函极值表达式：

```math
\sum_H\int_H [N]^TL^T[\mathbf{D}]L[N] dV[u^e]
- \int_H[N][f^e] dV
- \int_SN[\bar{f^e}] dS
= \mathbf{0}
```
记 ``L[N]`` 为``B``：

```math
\mathbf{B} =
\begin{bmatrix}
  N_{1,x} & 0 & 0 & N_{2,x} & 0 & 0 & N_{3,x} & 0 & 0 & N_{4,x} & 0 & 0 \\
  0 & N_{1,y} & 0 & 0 & N_{2,y} & 0 & 0 & N_{3,y} & 0 & 0 & N_{4,y} & 0 \\
  0 & 0 & N_{1,z} & 0 & 0 & N_{2,z} & 0 & 0 & N_{3,z} & 0 & 0 & N_{4,z} \\
  0 & \frac{N_{1,z}}{\sqrt{2}} & \frac{N_{1,y}}{\sqrt{2}} & 0 & \frac{N_{2,z}}{\sqrt{2}} & \frac{N_{2,y}}{\sqrt{2}} & 0 & \frac{N_{3,z}}{\sqrt{2}} & \frac{N_{3,y}}{\sqrt{2}} & 0 & \frac{N_{4,z}}{\sqrt{2}} & \frac{N_{4,y}}{\sqrt{2}} \\
  \frac{N_{1,z}}{\sqrt{2}} & 0 & \frac{N_{1,x}}{\sqrt{2}} & \frac{N_{2,z}}{\sqrt{2}} & 0 & \frac{N_{2,x}}{\sqrt{2}} & \frac{N_{3,z}}{\sqrt{2}} & 0 & \frac{N_{3,x}}{\sqrt{2}} & \frac{N_{4,z}}{\sqrt{2}} & 0 & \frac{N_{4,x}}{\sqrt{2}} \\
  \frac{N_{1,y}}{\sqrt{2}} & \frac{N_{1,x}}{\sqrt{2}} & 0 & \frac{N_{2,y}}{\sqrt{2}} & \frac{N_{2,x}}{\sqrt{2}} & 0 & \frac{N_{3,y}}{\sqrt{2}} & \frac{N_{3,x}}{\sqrt{2}} & 0 & \frac{N_{4,y}}{\sqrt{2}} & \frac{N_{4,x}}{\sqrt{2}} & 0
\end{bmatrix}
```
泛函极值表达式写作：

```math
\sum_H\int_H B^T[\mathbf{D}]B dV[u^e]
- \int_H[N][f^e] dV
- \int_SN[\bar{f^e}] dS
= \mathbf{0}
```
## 高斯积分

在每个单元上，积分还是太难求，用数值积分：

```math
\int_H B^T[\mathbf{D}]B dV[u^e]
= \sum \omega B^T[\mathbf{D}]B [u^e]
```
# 二维问题

```math
\begin{gathered}
\sigma =
\begin{bmatrix}
\sigma_x & \sigma_y & \sqrt{2}\tau_{xy}
\end{bmatrix}^T \\
\varepsilon =
\begin{bmatrix}
\varepsilon_x & \varepsilon_y & \sqrt{2}\varepsilon_{xy}
\end{bmatrix}^T
\end{gathered}
```
## 平衡微分方程：

```math
\begin{bmatrix}
\frac{\partial}{\partial x} & 0 & \frac{1}{\sqrt{2}}\frac{\partial}{\partial y} \\
0 & \frac{\partial}{\partial y} & \frac{1}{\sqrt{2}}\frac{\partial}{\partial x}
\end{bmatrix}
\begin{bmatrix}
\sigma_x\\
\sigma_y\\
\sqrt{2}\tau_{xy}
\end{bmatrix}
+
\begin{bmatrix}
f_x\\
f_y
\end{bmatrix} = \mathbf{0}
```
记为： ``L^T\sigma+\mathbf{f}=\mathbf{0}``

## 几何方程

```math
\begin{bmatrix}
\varepsilon_x\\
\varepsilon_y\\
\sqrt{2}\varepsilon_{xy}
\end{bmatrix}
=
\begin{bmatrix}
\frac{\partial}{\partial x} & 0 \\
0 & \frac{\partial}{\partial y} \\
\frac{1}{\sqrt{2}}\frac{\partial}{\partial y} & \frac{1}{\sqrt{2}}\frac{\partial}{\partial x}
\end{bmatrix}
\begin{bmatrix}
u_x\\
u_y
\end{bmatrix}
```
记为： ``\varepsilon = L\mathbf{u}``

## 物理方程

### 平面应力

```math
0 = \sigma_x = \frac{E}{(1+\nu)(1-2\nu)}(\nu\varepsilon_x + \nu\varepsilon_y+(1-\nu)\varepsilon_z)
```
```math
\varepsilon_z = -\frac{\nu}{1-\nu}(\varepsilon_x + \varepsilon_y)
```
```math
\begin{pmatrix}
\sigma_x\\
\sigma_y\\
\sqrt{2}\tau_{xy}
\end{pmatrix}
=\frac{E}{(1+\nu)(1-2\nu)}\begin{pmatrix}
1-\nu-\frac{\nu^2}{1-\nu} & \nu-\frac{\nu^2}{1-\nu} & 0\\
\nu-\frac{\nu^2}{1-\nu} & 1-\nu-\frac{\nu^2}{1-\nu} & 0\\
0 & 0 & 1-2\nu
\end{pmatrix}\begin{pmatrix}
\varepsilon_x\\
\varepsilon_y\\
\sqrt{2}\varepsilon_{xy}
\end{pmatrix}
```
```math
\begin{pmatrix}
\sigma_x\\
\sigma_y\\
\sqrt{2}\tau_{xy}
\end{pmatrix}
=\frac{E}{1-\nu^2}\begin{pmatrix}
1 & \nu & 0\\
\nu & 1 & 0\\
0 & 0 & 1-\nu
\end{pmatrix}\begin{pmatrix}
\varepsilon_x\\
\varepsilon_y\\
\sqrt{2}\varepsilon_{xy}
\end{pmatrix}
```
### 平面应变

```math
\sigma_z = \frac{E\nu}{(1+\nu)(1-2\nu)}(\varepsilon_x + \varepsilon_y)
```
```math
\begin{pmatrix}
\sigma_x\\
\sigma_y\\
\sqrt{2}\tau_{xy}
\end{pmatrix}
=\frac{E}{(1+\nu)(1-2\nu)}\begin{pmatrix}
1-\nu & \nu & 0\\
\nu & 1-\nu & 0\\
0 & 0 & 1-2\nu
\end{pmatrix}\begin{pmatrix}
\varepsilon_x\\
\varepsilon_y\\
\sqrt{2}\varepsilon_{xy}
\end{pmatrix}
```
