@doc raw"""
形函数``N``对全局坐标系的导数``B``如何求解以一阶三角形等参单元为例：
```math
B =
\begin{bmatrix}
\frac{\partial N_1}{\partial x} & \frac{\partial N_1}{\partial y} \\
\frac{\partial N_2}{\partial x} & \frac{\partial N_2}{\partial y} \\
\frac{\partial N_3}{\partial x} & \frac{\partial N_3}{\partial y}
\end{bmatrix}
 =
\begin{bmatrix}
\frac{\partial N_1}{\partial \xi} & \frac{\partial N_1}{\partial \eta}\\
\frac{\partial N_2}{\partial \xi} & \frac{\partial N_2}{\partial \eta}\\
\frac{\partial N_3}{\partial \xi} & \frac{\partial N_3}{\partial \eta}
\end{bmatrix}
\begin{bmatrix}
\frac{\partial \xi}{\partial x} & \frac{\partial \xi}{\partial y}\\
\frac{\partial \eta}{\partial x} & \frac{\partial \eta}{\partial y}
\end{bmatrix}
```
这里面，形函数对局部坐标系的导数很好求，
但局部坐标系对全局坐标系的导数不好求。

反过来写：

```math
\begin{bmatrix}
\frac{\partial N_1}{\partial \xi} & \frac{\partial N_1}{\partial \eta}\\
\frac{\partial N_2}{\partial \xi} & \frac{\partial N_2}{\partial \eta}\\
\frac{\partial N_3}{\partial \xi} & \frac{\partial N_3}{\partial \eta}
\end{bmatrix}
=
\begin{bmatrix}
\frac{\partial N_1}{\partial x} & \frac{\partial N_1}{\partial y} \\
\frac{\partial N_2}{\partial x} & \frac{\partial N_2}{\partial y} \\
\frac{\partial N_3}{\partial x} & \frac{\partial N_3}{\partial y}
\end{bmatrix}
\begin{bmatrix}
\frac{\partial x}{\partial \xi} & \frac{\partial x}{\partial \eta}\\
\frac{\partial y}{\partial \xi} & \frac{\partial y}{\partial \eta}
\end{bmatrix}
```
这个式子里面全局坐标系对局部坐标系的导数就好求了：

由于：

```math
\begin{gathered}
x = N_1 * x_1 + N_2 * x_2 + N_3 * x_3 \\
y = N_1 * y_1 + N_2 * y_2 + N_3 * y_3
\end{gathered}
```
即：

```math
\begin{bmatrix}
x\\
y
\end{bmatrix}
=
\begin{bmatrix}
x_1 & x_2 & x_3 \\
y_1 & y_2 & y_3
\end{bmatrix}
\begin{bmatrix}
N_1\\
N_2\\
N_3
\end{bmatrix}
```
所以：

```math
\begin{bmatrix}
\frac{\partial x}{\partial \xi} & \frac{\partial x}{\partial \eta}\\
\frac{\partial y}{\partial \xi} & \frac{\partial y}{\partial \eta}
\end{bmatrix}
=
\begin{bmatrix}
x_1 & x_2 & x_3 \\
y_1 & y_2 & y_3
\end{bmatrix}
\begin{bmatrix}
\frac{\partial N_1}{\partial \xi} & \frac{\partial N_1}{\partial \eta}\\
\frac{\partial N_2}{\partial \xi} & \frac{\partial N_2}{\partial \eta}\\
\frac{\partial N_3}{\partial \xi} & \frac{\partial N_3}{\partial \eta}
\end{bmatrix}
```
上式等号左边项即为局部坐标系到全局坐标系的转换矩阵，
或者叫 Jacobi 矩阵，记为 J。

```math
\begin{bmatrix}
\frac{\partial N_1}{\partial x} & \frac{\partial N_1}{\partial y} \\
\frac{\partial N_2}{\partial x} & \frac{\partial N_2}{\partial y} \\
\frac{\partial N_3}{\partial x} & \frac{\partial N_3}{\partial y}
\end{bmatrix}
 =
\begin{bmatrix}
\frac{\partial N_1}{\partial \xi} & \frac{\partial N_1}{\partial \eta}\\
\frac{\partial N_2}{\partial \xi} & \frac{\partial N_2}{\partial \eta}\\
\frac{\partial N_3}{\partial \xi} & \frac{\partial N_3}{\partial \eta}
\end{bmatrix}
\cdot
\text{inv(J)}
```
在参考坐标系下的高斯积分，其结果要乘上Jacobi 矩阵的行列式。
"""
B_gen(C::Type, ref_coords::Vector{Float64}, weight::Float64, coords::Matrix{Float64}) =
    B_gen(N_r_gen(C, ref_coords), weight, coords)
@inline function B_gen(N_r::Matrix{Float64}, weight::Float64, coords::Matrix{Float64})
    J = coords * N_r #N_r --- the derivative of N to the reference coordinates
    N_g = N_r / J #N_g --- the derivative of N to the global coordinates
    return N_g, weight * det(J) #JxW --- the product of the determinant of the Jacobian matrix and the weight of quadrature point
end

@doc raw"""
弹性理论的小变形假设下的一个点的应变是这么求的：
```math
\begin{equation*}
\begin{bmatrix}
\varepsilon_x\\
\varepsilon_y\\
\sqrt{2}\varepsilon_{xy}
\end{bmatrix}=
\begin{bmatrix}
\frac{\partial}{\partial x} & 0\\
0 & \frac{\partial}{\partial y}\\
\frac{1}{\sqrt{2}}\frac{\partial}{\partial y} & \frac{1}{\sqrt{2}}\frac{\partial}{\partial x}
\end{bmatrix}
\begin{bmatrix}
u_x\\
u_y
\end{bmatrix}
\end{equation*}
```
用形函数来近似位移：

```math
\begin{equation*}
\begin{bmatrix}
u_x\\
u_y
\end{bmatrix}=
\begin{bmatrix}
N_1&0&N_2&0&N_3&0\\
0&N_1&0&N_2&0&N_3
\end{bmatrix}
\begin{bmatrix}
u_{1x}\\
u_{1y}\\
u_{2x}\\
u_{2y}\\
u_{3x}\\
u_{3y}
\end{bmatrix}
\end{equation*}
```
带入前面的式子里：

```math
\begin{equation*}
\begin{bmatrix}
\varepsilon_x\\
\varepsilon_y\\
\varepsilon_{xy}
\end{bmatrix}=
\begin{bmatrix}
\frac{\partial N_1}{\partial x} & 0 &
\frac{\partial N_2}{\partial x} & 0 &
\frac{\partial N_3}{\partial x} & 0\\
0 & \frac{\partial N_1}{\partial y} &
0 & \frac{\partial N_2}{\partial y} &
0 & \frac{\partial N_3}{\partial y}\\
\frac{1}{\sqrt{2}}\frac{\partial N_1}{\partial y} &
\frac{1}{\sqrt{2}}\frac{\partial N_1}{\partial x} &
\frac{1}{\sqrt{2}}\frac{\partial N_2}{\partial y} &
\frac{1}{\sqrt{2}}\frac{\partial N_2}{\partial x} &
\frac{1}{\sqrt{2}}\frac{\partial N_3}{\partial y} &
\frac{1}{\sqrt{2}}\frac{\partial N_3}{\partial x} &
\end{bmatrix}
\begin{bmatrix}
u_{1x}\\
u_{1y}\\
u_{2x}\\
u_{2y}\\
u_{3x}\\
u_{3y}
\end{bmatrix}
\end{equation*}
```
上式等号右边第一项是弹性力学问题求解中需要用到的``B``
"""
function elastic_B_gen(C::Type, D::Type, ref_coords::Vector{Float64}, weight::Float64, coords::Matrix{Float64})
    N_g, JxW = B_gen(C, ref_coords, weight, coords)
    B = elastic_B_gen(D, N_g, size(coords, 2))
    return B, JxW
end
function elastic_B_gen(D::Type, N_g::Matrix{Float64}, node_number)
    dim = D <: Dim2 ? 2 : 3
    dof = dim * node_number
    B = dim == 2 ? zeros(3, dof) : zeros(6, dof)
    if dim == 2
        for i in axes(N_g, 1)
            B[1, 2i-1] = B[3,   2i] = N_g[i, 1]
            B[2,   2i] = B[3, 2i-1] = N_g[i, 2]
        end
    else
        for i in axes(N_g, 1)
            B[1, 3i-2] = B[5,   3i] = B[6, 3i-1] = N_g[i, 1]
            B[2, 3i-1] = B[4,   3i] = B[6, 3i-2] = N_g[i, 2]
            B[3,   3i] = B[4, 3i-1] = B[5, 3i-2] = N_g[i, 3]
        end
    end
    B[dim+1:end, :] ./= sqrt(2)
    return B
end



@doc raw"""
```
 ^η
 |
 3 (0,1)
 |'\
 |  '\
 |    '\
 |      '\
 |        '\
 1----------2-->ξ
(0,0)      (1,0)
```
三个结点的面积坐标分别为
```math
1（1,0,0)\quad 2(0,1,0)\quad 3(0,0,1)
```
由上图可以看到，我这里的局部坐标系是这样定义的：
```math
  \xi=L_2,\quad \eta=L_3,\quad 1-\xi-\eta=L_1
```
形函数即面积坐标：

```math
\begin{aligned}
  &N_1 = L_1 = 1-\xi-\eta \\
  &N_2 = L_2 = \xi \\
  &N_3 = L_3 = \eta
\end{aligned}
```

形函数对局部坐标``\xi`` 、``\eta`` 的导数为：
```math
N_{\xi\eta} = \begin{bmatrix}
N_{1,\xi}&N_{1,\eta}\\
N_{2,\xi}&N_{2,\eta}\\
N_{3,\xi}&N_{3,\eta}
\end{bmatrix} = \begin{bmatrix}
-1&-1\\1&0\\0&1
\end{bmatrix}
```
其中``N_{1,\xi}`` 表示``\frac{\partial N_1}{\partial\xi}`` , etc.
"""
N_gen(::Type{t1}, ref_coords::Vector{Float64}) = ref_coords
N_r_gen(::Type{t1}, ref_coords::Vector{Float64}) = [-1 -1.0
                                                     1  0
                                                     0  1]
@doc raw"""
```
  ^η
  |
1 3
  |'\
  |  '\
  6    '5
  |      '\
  |        '\
  1----4-----2-->ξ
  0          1
```
各个点的面积坐标``(L_1, L_2, L_3)``：\
角点：
```math
  1(1, 0,0)\quad 2(0,1,0)\quad 3(0,0,1)
```
边中点：
```math
  4(0.5,0.5,0)\quad 5(0,0.5,0.5)\quad 6(0.5,0,0.5)
```
和一阶的时候一样，局部坐标是这样布置的：
```math
  \xi=L_2,\quad \eta=L_3,\quad 1-\xi-\eta=L_1
```
故形函数为：

```math
\begin{aligned}
N_1 =& (2L_1-1)L_1\\
N_2 =& (2L_2-1)L_2\\
N_3 =& (2L_3-1)L_3\\
N_4 =& 4L_1L_2\\
N_5 =& 4L_2L_3\\
N_6 =& 4L_3L_1
\end{aligned}
```
```math
N_r = \frac{\partial N}{\partial L}\frac{\partial L}{\partial r} =
\begin{bmatrix}
1-4L_1 & 1-4L_1\\
4L_2-1 & 0\\
0 & 4L_3-1\\
4L_1-4L_2& -4L_2\\
4L3 & 4L2\\
-4L3 & 4L1-4L3
\end{bmatrix}
```
"""
function N_gen(::Type{t2}, ref_coords::Vector{Float64})
    L1, L2, L3 = ref_coords
    [(2L1 - 1) * L1
     (2L2 - 1) * L2
     (2L3 - 1) * L3
     4L1*L2
     4L2*L3
     4L3*L1]
end
@inline function N_r_gen(::Type{t2}, ref_coords::Vector{Float64})
    L1, L2, L3 = ref_coords
    [ 1-4L1     1-4L1
      4L2-1     0
      0         4L3-1
      4L1-4L2  -4L2
      4L3       4L2
     -4L3       4L1-4L3]
end

@doc raw"""
```
                      v
                    .
                  ,/
                 /
               3
             ,/|'\
           ,/  |  '\
         ,/    '.   '\
       ,/       |     '\
     ,/         |       '\
    1-----------'.--------2 --> u
     '\.         |       ,/
        '\.       |    ,/
           '\.    '. ,/
              '\.  |/
                 ' 4
                    '\.
                       ' w
```
这里局部坐标是这么布置的：（u、v、w相互正交)\
u = L2, v = L3, w = L4, 1-u-v-w = L1

各个角点的体积坐标为：\
1(1 0 0 0), 2(0, 1, 0, 0), 3(0, 0, 1, 0), 4(0, 0, 0, 1)

形函数为：
N1 = L1 = 1-u-v-w\
N2 = L2 = u\
N3 = L3 = v\
N4 = L4 = w\
N_local = [-1 -1 -1\
           1 0 0\
           0 1 0\
           1 0 1]\
可以看到这里又和积分点坐标无关了，
T1单元是常应变的单元，合理。
"""
N_gen(::Type{T1}, ref_coords::Vector{Float64}) = ref_coords
N_r_gen(::Type{T1}, ref_coords::Vector{Float64}) = [-1 -1 -1.0
                                                     1  0  0
                                                     0  1  0
                                                     0  0  1]

@doc raw"""
```
       ^η
       |1
  4----------3
  |    |     |
  |    |     |
--|----+-----|--->
-1|    |     |1 ξ
  |    |     |
  1----------2
       |-1
```
四个角点在局部坐标系中的坐标``(\xi, \eta)``为：
```math
  1(-1, -1)\quad  2(1, -1)\quad  3(1, 1)\quad  4(-1, 1)
```
形函数N为：
```math
  N_i = (1+\xi_i\xi)(1+\eta_i\eta)/4
```
thus:

```math
  N_{\xi\eta} = \frac{1}{4} \begin{bmatrix}
\xi_1(1+\eta_1\eta) & \eta_1(1+\xi_1\xi)\\
\xi_2(1+\eta_2\eta) & \eta_2(1+\xi_2\xi)\\
\xi_3(1+\eta_3\eta) & \eta_3(1+\xi_3\xi)\\
\xi_4(1+\eta_4\eta) & \eta_4(1+\xi_4\xi)
\end{bmatrix} = \frac{1}{4} \begin{bmatrix}
\eta-1 & \xi-1\\
1-\eta & -1-\xi\\
1+\eta & 1+\xi\\
-1-\eta & 1-\xi
\end{bmatrix}
```
"""
function N_gen(::Type{q1}, ref_coords::Vector{Float64})
    ξ, η = ref_coords
    [(1-ξ)*(1-η)
     (1+ξ)*(1-η)
     (1+ξ)*(1+η)
     (1-ξ)*(1+η)] / 4
end
@inline function N_r_gen(::Type{q1}, ref_coords::Vector{Float64})
    ξ, η = ref_coords
    [ η-1  ξ-1
      1-η -1-ξ
      1+η  1+ξ
     -1-η  1-ξ] / 4
end

@doc raw"""
```
           η
    4----------3
    |\     ^   |\
    | \    |   | \
    |  \   |   |  \
    |   8----------7
    |   |  +-- |-- | -> ξ
    1---+---\--2   |
     \  |    \  \  |
      \ |     \  \ |
       \|      ζ  \|
        5----------6
```
其实这个局部坐标怎么布置不重要。这里三个坐标系的顺序是：ξ-η-ζ\
重要的是结点的排布，我这里的结点排布是和gmsh一致的

各个角点在局部坐标系的坐标为：\
1(-1,-1,-1) 2(1,-1,-1) 3(1,1,-1) 4(-1,1,-1)\
5(-1,-1, 1) 6(1,-1, 1) 7(1,1, 1) 8(-1,1, 1)

形函数为：\
Ni = (1+ξi*ξ)(1+ηi*η)(1+ζi*ζ)/8
"""
function N_gen(::Type{H1}, ref_coords::Vector{Float64})
    ξ, η, ζ = ref_coords
    [(1-ξ)*(1-η)*(1-ζ)
     (1+ξ)*(1-η)*(1-ζ)
     (1+ξ)*(1+η)*(1-ζ)
     (1-ξ)*(1+η)*(1-ζ)
     (1-ξ)*(1-η)*(1+ζ)
     (1+ξ)*(1-η)*(1+ζ)
     (1+ξ)*(1+η)*(1+ζ)
     (1-ξ)*(1+η)*(1+ζ)] / 8
end
@inline function N_r_gen(::Type{H1}, ref_coords::Vector{Float64})
    ξ, η, ζ = ref_coords
    [-(1-η)*(1-ζ) -(1-ξ)*(1-ζ) -(1-ξ)*(1-η)
      (1-η)*(1-ζ) -(1+ξ)*(1-ζ) -(1+ξ)*(1-η)
      (1+η)*(1-ζ)  (1+ξ)*(1-ζ) -(1+ξ)*(1+η)
     -(1+η)*(1-ζ)  (1-ξ)*(1-ζ) -(1-ξ)*(1+η)
     -(1-η)*(1+ζ) -(1-ξ)*(1+ζ)  (1-ξ)*(1-η)
      (1-η)*(1+ζ) -(1+ξ)*(1+ζ)  (1+ξ)*(1-η)
      (1+η)*(1+ζ)  (1+ξ)*(1+ζ)  (1+ξ)*(1+η)
     -(1+η)*(1+ζ)  (1-ξ)*(1+ζ)  (1-ξ)*(1+η)] / 8
end
