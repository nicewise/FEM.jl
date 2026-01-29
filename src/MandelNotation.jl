# MandelNotation.jl
# encoding: utf-8

@doc raw"""
两向量的并矢
```math
a_ib_j = \begin{bmatrix}
a_1b_1&a_1b_2&a_1b_3\\
a_2b_1&a_2b_2&a_2b_3\\
a_3b_1&a_3b_2&a_3b_3
\end{bmatrix}
```
如果这个二阶张量是对称的话，则可以降阶为：
```math
a_ib_j = \begin{bmatrix}
a_1b_1\\
a_2b_2\\
a_3b_3\\
\sqrt{2}a_2b_3\\
\sqrt{2}a_3b_1\\
\sqrt{2}a_1b_2
\end{bmatrix}
```

由两个向量构件二阶对称张量：
```math
a\otimes^sb = \frac{1}{2}(a\otimes b+b\otimes a)
=\frac{1}{2}\begin{bmatrix}
a_1b_1+a_1b_1&a_1b_2+a_2b_1&a_1b_3+a_3b_1\\
a_2b_1+a_1b_2&a_2b_2+a_2b_2&a_2b_3+a_3b_2\\
a_3b_1+a_1b_3&a_3b_2+a_2b_3&a_3b_3+a_3b_3
\end{bmatrix}\overset{\text{降阶}}{=}
\frac{1}{2}\begin{bmatrix}
a_1b_1+a_1b_1\\
a_2b_2+a_2b_2\\
a_3b_3+a_3b_3\\
\sqrt{2}(a_2b_3+a_3b_2)\\
\sqrt{2}(a_3b_1+a_1b_3)\\
\sqrt{2}(a_1b_2+a_2b_1)
\end{bmatrix}
```
"""
function ⊗(a::SVector{3}, b::SVector{3}; symmetric::Bool = false)
    if symmetric
        return (a ⊗ b + b ⊗ a) / 2
    else
        return @SVector [
            a[1] * b[1],
            a[2] * b[2],
            a[3] * b[3],
            a[2] * b[3] * √2,
            a[3] * b[1] * √2,
            a[1] * b[2] * √2
        ]
    end
end

function ⊗(A::SVector{6}, B::SVector{6}; symmetric::Bool = false)
    if symmetric
        c = [1, 1, 1, √2, √2, √2]
        a = A ./ c
        b = B ./ c
        o = [
            a[1]*b[1] a[6]*b[6] a[5]*b[5] a[6]*b[5] a[5]*b[1] a[1]*b[6]
            a[6]*b[6] a[2]*b[2] a[4]*b[4] a[2]*b[4] a[4]*b[6] a[6]*b[2]
            a[5]*b[5] a[4]*b[4] a[3]*b[3] a[4]*b[3] a[3]*b[5] a[5]*b[4]
            a[6]*b[5] a[2]*b[4] a[4]*b[3] a[2]*b[3] a[4]*b[5] a[6]*b[4]
            a[5]*b[1] a[4]*b[6] a[3]*b[5] a[4]*b[5] a[3]*b[1] a[5]*b[6]
            a[1]*b[6] a[6]*b[2] a[5]*b[4] a[6]*b[4] a[5]*b[6] a[1]*b[2]
        ] + [
            a[1]*b[1] a[6]*b[6] a[5]*b[5] a[5]*b[6] a[1]*b[5] a[6]*b[1]
            a[6]*b[6] a[2]*b[2] a[4]*b[4] a[4]*b[2] a[6]*b[4] a[2]*b[6]
            a[5]*b[5] a[4]*b[4] a[3]*b[3] a[3]*b[4] a[5]*b[3] a[4]*b[5]
            a[6]*b[5] a[2]*b[4] a[4]*b[3] a[4]*b[4] a[6]*b[3] a[2]*b[5]
            a[5]*b[1] a[4]*b[6] a[3]*b[5] a[3]*b[6] a[5]*b[5] a[4]*b[1]
            a[1]*b[6] a[6]*b[2] a[5]*b[4] a[5]*b[2] a[1]*b[4] a[6]*b[6]
        ]
        o ./= 2
        o[:, 4:6] .*= √2
        o[4:6, :] .*= √2
        return SMatrix{6, 6}(o)
    else
        return A * B'
    end
end
⋆(A::SVector{3}, B::SVector{3}) = ⊗(A, B, symmetric = true) # ⊗s
⋆(A::SVector{6}, B::SVector{6}) = ⊗(A, B, symmetric = true) # ⊗s

struct MandelNotation
    δ::SVector{6}
    I::SMatrix{6, 6}
    J::SMatrix{6, 6, Float64}
    K::SMatrix{6, 6, Float64}

    function MandelNotation()
        δ = @SVector [1, 1, 1, 0, 0, 0]
        # I = δ ⊗s δ, 但我不想让它变成浮点数矩阵
        I = @SMatrix [ 1  0  0  0  0  0
                       0  1  0  0  0  0
                       0  0  1  0  0  0
                       0  0  0  1  0  0
                       0  0  0  0  1  0
                       0  0  0  0  0  1 ]
        J = δ ⊗ δ / 3
        K = I - J
        return new(δ, I, J, K)
    end
end

@doc raw"""
- Mandel.δ : 二阶单位张量
- Mandel.I : 对称四阶单位张量
- Mandel.J : 体积四阶单位张量
- Mandel.K : 偏差四阶单位张量

二阶对称张量``A``的降阶表示为：
```math
A = \begin{bmatrix}
A_{11}\\
A_{22}\\
A_{33}\\
\sqrt{2}A_{23}\\
\sqrt{2}A_{31}\\
\sqrt{2}A_{12}\\
\end{bmatrix}
```
四阶张量U的降阶表示为：
```math
U = \begin{bmatrix}
U_{1111}&U_{1122}&U_{1133}&\sqrt{2}U_{1123}&\sqrt{2}U_{1131}&\sqrt{2}U_{1112}\\
U_{2211}&U_{2222}&U_{2233}&\sqrt{2}U_{2223}&\sqrt{2}U_{2231}&\sqrt{2}U_{2212}\\
U_{3311}&U_{3322}&U_{3333}&\sqrt{2}U_{3323}&\sqrt{2}U_{3331}&\sqrt{2}U_{3312}\\
\sqrt{2}U_{2311}&\sqrt{2}U_{2322}&\sqrt{2}U_{2333}& 2U_{2323}& 2U_{2331}& 2U_{2312}\\
\sqrt{2}U_{3111}&\sqrt{2}U_{3122}&\sqrt{2}U_{3133}& 2U_{3123}& 2U_{3131}& 2U_{3112}\\
\sqrt{2}U_{1211}&\sqrt{2}U_{1222}&\sqrt{2}U_{1233}& 2U_{1223}& 2U_{1231}& 2U_{1212}
\end{bmatrix}
```
```math
\delta = \begin{bmatrix}
1\\1\\1\\0\\0\\0
\end{bmatrix}
```
```math
I = \delta \otimes^s \delta
  = \frac{1}{2} (\delta_{ik} \delta_{jl} + \delta_{il} \delta_{jk}) e_i\otimes e_j\otimes e_k\otimes e_l
  = \begin{bmatrix}
   1& 0& 0& 0& 0& 0\\
   0& 1& 0& 0& 0& 0\\
   0& 0& 1& 0& 0& 0\\
   0& 0& 0& 1& 0& 0\\
   0& 0& 0& 0& 1& 0\\
   0& 0& 0& 0& 0& 1
\end{bmatrix}
```
```math
J = \frac{1}{3} \delta \otimes \delta
  = \frac{1}{3} \delta_{ij} \delta_{kl}
  = \frac{1}{3}\begin{bmatrix}
   1& 1& 1& 0& 0& 0\\
   1& 1& 1& 0& 0& 0\\
   1& 1& 1& 0& 0& 0\\
   0& 0& 0& 0& 0& 0\\
   0& 0& 0& 0& 0& 0\\
   0& 0& 0& 0& 0& 0
\end{bmatrix}
```
```math
K = I - J
  = \frac{1}{3}\begin{bmatrix}
    2& -1& -1& 0& 0& 0\\
   -1&  2& -1& 0& 0& 0\\
   -1& -1&  2& 0& 0& 0\\
    0&  0&  0& 3& 0& 0\\
    0&  0&  0& 0& 3& 0\\
    0&  0&  0& 0& 0& 3
\end{bmatrix}
```
```jldoctest
julia> Mandel.J
6×6 StaticArraysCore.SMatrix{6, 6, Float64, 36} with indices SOneTo(6)×SOneTo(6):
 0.333333  0.333333  0.333333  0.0  0.0  0.0
 0.333333  0.333333  0.333333  0.0  0.0  0.0
 0.333333  0.333333  0.333333  0.0  0.0  0.0
 0.0       0.0       0.0       0.0  0.0  0.0
 0.0       0.0       0.0       0.0  0.0  0.0
 0.0       0.0       0.0       0.0  0.0  0.0

julia> Mandel.K
6×6 StaticArraysCore.SMatrix{6, 6, Float64, 36} with indices SOneTo(6)×SOneTo(6):
  0.666667  -0.333333  -0.333333  0.0  0.0  0.0
 -0.333333   0.666667  -0.333333  0.0  0.0  0.0
 -0.333333  -0.333333   0.666667  0.0  0.0  0.0
  0.0        0.0        0.0       1.0  0.0  0.0
  0.0        0.0        0.0       0.0  1.0  0.0
  0.0        0.0        0.0       0.0  0.0  1.0
```
"""
const Mandel = MandelNotation()
