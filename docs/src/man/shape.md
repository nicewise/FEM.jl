# 形函数 N 以及对坐标的导数
这里说的是等参单元的形函数。
```@docs
B_gen
```
弹性力学问题中需要的并不是上边所描述的``B``，而是：
```@docs
elastic_B_gen
```
## 平面单元
### 一阶三角形单元
```@docs
N_gen(::Type{t1}, ::Vector{Float64})
```
### 二阶三角形单元
```@docs
N_gen(::Type{t2}, ::Vector{Float64})
```
### 一阶四边形单元
```@docs
N_gen(::Type{q1}, ::Vector{Float64})
```
## 空间单元
### 一阶四面体单元
```@docs
N_gen(::Type{T1}, ::Vector{Float64})
```
### 一阶六面体单元
```@docs
N_gen(::Type{H1}, ::Vector{Float64})
```
