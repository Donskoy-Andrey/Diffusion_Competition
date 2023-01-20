# Competitive task on supercomputer practice

## Formulation of the diffusion problem
The solution of the non-stationary diffusion equation with unknown U is considered.
```math
$U = U(x, y, z, t):$
```

```math
$∂U/∂t - div (D · grad U) = f(x, y, z, t), (x,y,z)$
``` 

принадлежит $Ω = [0; 1]^3$, t на отрезке $[0, T]$

```math
$U(x, y, z, t) = g(x, y, z)$ 
```
на границе области $∂Ω$

```math
$U(x, y, z, 0) = 0$ 
```
в начальный момент времени

Конечный момент времени зададим как $T = 1$.

В задаче будем использовать диагональный тензор $D:$

```math
D=\begin{pmatrix}
d_x & 0 & 0 \\
0 & d_y & 0 \\
0 & 0 & d_z 
\end{pmatrix}
```

где $d_x = 0.25, d_y = 0.15, d_z = 0.1$

Также зададим
```math
$g(x, y, z) = 0$
$f(x, y, z) = (d_x+d_y+d_z)·π²·sin(πx)sin(πy)sin(πz)$
```

У решаемого уравнения существует аналитическое решение:
```math
$U_analityc = sin(πx)sin(πy)sin(πz)·(1 - exp(-(d_x+d_y+d_z)·π²·t))$
```

# Дискретизация
Построим параллелепипедную дискретизацию нашей области.

А именно, выберем числа $Nx, Ny, Nz > 1$ - количество узлов, которые будут укладываться вдоль оси $Ox, Oy и Oz$ соотвественно.
Тогда определим шаг сетки $Δx = \frac{1}{(Nx-1)}, Δy = \frac{1}{(Ny-1)}, Δz = \frac{1}{(Nz-1)}.$
Также определим шаг по времени $Δt$.
Обозначим через $V_{ijk}$ узел сетки с координатами $x_i = i·Δx, y_j = j·Δy, z_k = k·Δz.$

Будем описывать дискретную функцию $[U]^h$ в момент времени $nΔt$ её степенями свободы, которые расположим в узлах сетки и степень свободы в узле V_ijk обозначим как

$U_{ijk}^n, 0 ⩽ i ⩽ Nx-1, 0 ⩽ j ⩽ Ny-1, 0 ⩽ k ⩽ Nz-1.$

Дискретизуем наше уравнение по пространству методом конечных разностей, а по времени явной схемой Эйлера.

Заметим, что для такой дискретизации шаг по времени должен удовлетворять условию Куранта:

$Δt < \frac{0.5}{ \frac{d_x } {(Δx)^2} + \frac{d_y}{(Δy)^2} + \frac{d_z}{(Δz)^2} }$

Введём дискретные операторы вторых пространственных производные:

$Lx_{ijk} U^n = \frac{ U_{(i-1)jk}^n - 2·U_{ijk}^n + U_{(i+1)jk}^n }{Δx^2}$

$Ly_{ijk} U^n = \frac{ U_{i(j-1)k}^n - 2·U_{ijk}^n + U_{i(j+1)k}^n }{Δy)^2}$

$Lz_{ijk} U^n = \frac{ U_{ij(k-1)}^n - 2·U_{ijk}^n + U_{ij(k+1)}^n }{Δz^2}$

Также пусть

$f_{ijk}^n = f(i·Δx, j·Δy, k·Δz, nΔt)$ и $g_{ijk} = g(i·Δx, j·Δy, k·Δz)$

Тогда имеем следующую численную схему:
$U_{ijk}^{(n+1)} = U_{ijk}^n + Δt·(f_{ijk}^n + d_x·Lx_{ijk} U^n + d_y·Ly_{ijk} U^n + d_z·Lz_{ijk} U^n), 

если

1 ⩽ i ⩽ Nx-2, 1 ⩽ j ⩽ Ny-2, 1 ⩽ k ⩽ Nz-2$

$U_{ijk}^{(n+1)} = g_{ijk}$, если $(i%(Nx-1)) · (j%(Ny-1)) · (k%(Nz-1)) = 0$,
где $x%y$ - операция взятия остатка от деления $x$ на $y$.

## Start code

```bash
pip install -r requirements.txt
mkdir build
cd build
cmake ..
make
./main
```

## Analytical solution
<p align="center">
  <img src="data/images/analytical.png">
</p>