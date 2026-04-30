# 第四章：不动点定理及其应用

在泛函分析、微分方程、积分方程以及其他各类方程的理论中，解的存在性、唯一性以及近似解的收敛性等都是相当重要的课题。为了讨论这些方程解的存在性及唯一性，我们常常将它们转化成求某一映射的**不动点问题**。

!!! info "定义 4.1 (不动点)"

    设 $T$ 是定义在空间 $X$ 上的映射。如果存在点 $x \in X$，满足：

    $$
    T(x) = x
    $$

    则称点 $x$ 为映射 $T$ 的**不动点**。

---

## 1. 基本的不动点定理：压缩映射原理

### 1.1 压缩映射与不动点定理

!!! success "定理 4.1 (压缩映射原理)"

    设 $X$ 是**完备的距离空间**，距离为 $\rho$。$T$ 是由 $X$ 到其自身的映射，且对任意的 $x, y \in X$，不等式：

    $$
    \rho(Tx, Ty) \le \theta \rho(x, y)
    $$

    成立，其中 $\theta$ 是满足 $0 \le \theta < 1$ 的常数。
    
    那么，$T$ 在 $X$ 中存在**唯一**的不动点，即存在唯一的 $\overline{x}$ 使得 $T\overline{x} = \overline{x}$，且该不动点可以用迭代法求得。

满足定理条件的映射 $T$ 称为**压缩映射**。

??? proof "压缩映射原理的证明"

    **1. 构造迭代序列并证明其为基本列 (Cauchy 列)：**

    在 $X$ 中任意取定一点 $x_0$，我们通过迭代构造点列：

    $$
    x_1 = Tx_0, \quad x_2 = Tx_1, \quad \dots, \quad x_{n+1} = Tx_n, \quad \dots
    $$

    我们首先计算相邻两项之间的距离。由压缩映射的定义可知：

    $$
    \rho(x_1, x_2) = \rho(Tx_0, Tx_1) \le \theta \rho(x_0, x_1)
    $$

    $$
    \rho(x_2, x_3) = \rho(Tx_1, Tx_2) \le \theta \rho(x_1, x_2) \le \theta^2 \rho(x_0, x_1)
    $$

    应用数学归纳法，可以得到一般性的递推关系：

    $$
    \rho(x_{n-1}, x_n) \le \theta^{n-1} \rho(x_0, x_1)
    $$

    进而：

    $$
    \rho(x_n, x_{n+1}) = \rho(Tx_{n-1}, Tx_n) \le \theta \rho(x_{n-1}, x_n) \le \theta^n \rho(x_0, x_1)
    $$

    对于任意的自然数 $p \ge 1$，利用三角不等式以及上述距离放缩，我们有：

    $$
    \rho(x_n, x_{n+p}) \le \rho(x_n, x_{n+1}) + \rho(x_{n+1}, x_{n+2}) + \dots + \rho(x_{n+p-1}, x_{n+p})
    $$

    $$
    \le (\theta^n + \theta^{n+1} + \dots + \theta^{n+p-1}) \rho(x_0, x_1)
    $$

    利用等比数列求和公式，并注意到 $\theta < 1$：

    $$
    \le \frac{\theta^n}{1 - \theta} \rho(x_0, x_1)
    $$

    由于 $0 \le \theta < 1$，当 $n \rightarrow \infty$ 时，$\theta^n \rightarrow 0$。因此对于任给的 $\epsilon > 0$，当 $n$ 充分大时，对一切 $p \ge 1$ 均有 $\rho(x_n, x_{n+p}) < \epsilon$。
    
    这说明 $\{x_n\}$ 是一个**基本列**。由于空间 $X$ 是**完备**的，$\{x_n\}$ 必然收敛于 $X$ 中的某一点 $\overline{x} \in X$。

    **2. 证明极限点是不动点：**

    由于压缩映射满足 Lipschitz 条件（常数 $\theta < 1$），它必然是一个**连续映射**。
    
    在递推关系 $x_{n+1} = Tx_n$ 的两端同时令 $n \rightarrow \infty$，利用连续性可得：

    $$
    \lim_{n \rightarrow \infty} x_{n+1} = \lim_{n \rightarrow \infty} Tx_n = T(\lim_{n \rightarrow \infty} x_n)
    $$

    即：

    $$
    \overline{x} = T\overline{x}
    $$

    这证明了 $\overline{x}$ 是映射 $T$ 的一个不动点。

    **3. 证明不动点的唯一性：**

    假设存在另外一个不动点 $\overline{y} \in X$ 满足 $\overline{y} = T\overline{y}$。我们计算 $\overline{x}$ 与 $\overline{y}$ 的距离：

    $$
    \rho(\overline{x}, \overline{y}) = \rho(T\overline{x}, T\overline{y}) \le \theta \rho(\overline{x}, \overline{y})
    $$

    即 $(1 - \theta)\rho(\overline{x}, \overline{y}) \le 0$。
    
    由于 $\theta < 1$ 且距离非负，必定有 $\rho(\overline{x}, \overline{y}) = 0$，也就是 $\overline{x} = \overline{y}$。证明完毕。 $\square$

### 1.2 迭代序列的误差估计

由证明过程可以看出，为了获得不动点，我们可以从 $X$ 中的任一点 $x_0$ 出发建立迭代序列 $x_n = T^n x_0$。

在大多数情况下，精确的不动点不易求得，因此我们往往用 $x_n$ 作为其近似值。根据证明过程中的放缩，我们可以得到 $x_n$ 与真实不动点 $\overline{x}$ 之间的误差估计（令 $p \rightarrow \infty$）：

$$
\rho(x_n, \overline{x}) \le \frac{\theta^n}{1 - \theta} \rho(x_0, x_1)
$$

---

## 2. 压缩映射原理的应用

### 2.1 微分方程初值问题解的存在性与唯一性

考察一阶常微分方程初值问题（Picard 存在唯一性定理）：

$$
\frac{dy}{dx} = f(x, y), \quad y|_{x=x_0} = y_0
$$

假设 $f(x, y)$ 在一个区域上连续，并且关于 $y$ 满足 **Lipschitz 条件**：

$$
|f(x, y) - f(x, y')| \le K|y - y'|
$$

其中 $K > 0$ 为常数。

由微积分理论可知，该微分方程初值问题等价于下面的 Volterra 型积分方程：

$$
y(x) = y_0 + \int_{x_0}^x f(t, y(t)) dt
$$

??? proof "利用压缩映射原理论证（点击展开）"

    取充分小的 $\delta > 0$，使得 $K\delta < 1$。
    
    我们在连续函数空间 $C[x_0 - \delta, x_0 + \delta]$ 上定义映射 $T$：

    $$
    (Ty)(x) = y_0 + \int_{x_0}^x f(t, y(t)) dt, \quad x \in [x_0 - \delta, x_0 + \delta]
    $$

    对于任意两个连续函数 $y_1(x), y_2(x) \in C[x_0 - \delta, x_0 + \delta]$，计算它们经过映射 $T$ 后的距离：

    $$
    \rho(Ty_1, Ty_2) = \max_{x \in [x_0 - \delta, x_0 + \delta]} \left| \int_{x_0}^x [f(t, y_1(t)) - f(t, y_2(t))] dt \right|
    $$

    利用 Lipschitz 条件放缩：

    $$
    \le \max_{x} \left| \int_{x_0}^x K |y_1(t) - y_2(t)| dt \right|
    $$

    $$
    \le K \delta \max_{t \in [x_0 - \delta, x_0 + \delta]} |y_1(t) - y_2(t)| = K\delta \rho(y_1, y_2)
    $$

    由于我们选取了 $K\delta < 1$，故 $T$ 是一个**压缩映射**。由压缩映射原理可知，在空间 $C[x_0 - \delta, x_0 + \delta]$ 中存在唯一的连续函数 $y_0(x)$ 满足积分方程，即初值问题在该区间内有唯一解。
    
    随后可将此解进一步向外延拓，最终可将解延拓到整个实轴上。 $\square$

### 2.2 线性积分方程解的存在性与唯一性

设有 Fredholm 第二类线性积分方程：

$$
x(t) = f(t) + \lambda \int_a^b K(t, s)x(s) ds
$$

其中已知函数 $f \in L^2[a, b]$，$\lambda$ 为参数。核函数 $K(t, s)$ 是定义在矩形区域 $a \le t \le b$, $a \le s \le b$ 上的可测函数，且满足平方可积条件：

$$
\int_a^b \int_a^b |K(t, s)|^2 dt ds < \infty
$$

??? proof "证明：当 $|\lambda|$ 充分小时方程有唯一解"

    我们在完备的距离空间 $L^2[a, b]$ 上定义映射 $T$：

    $$
    (Tx)(t) = f(t) + \lambda \int_a^b K(t, s)x(s) ds
    $$

    首先由 Cauchy-Schwarz 不等式易证 $Tx \in L^2[a, b]$，即 $T$ 是由 $L^2[a, b]$ 到其自身的映射：

    $$
    \int_a^b \left| \int_a^b K(t,s)x(s)ds \right|^2 dt \le \int_a^b \left[ \int_a^b |K(t,s)|^2 ds \int_a^b x(s)^2 ds \right] dt < \infty
    $$

    接着，取 $|\lambda|$ 充分小，使得：

    $$
    \theta = |\lambda| \left[ \int_a^b \int_a^b |K(t, s)|^2 ds dt \right]^{\frac{1}{2}} < 1
    $$

    计算任意 $y_1, y_2 \in L^2[a,b]$ 映射后的距离（由 Cauchy-Schwarz 不等式）：

    $$
    \rho(Ty_1, Ty_2) = |\lambda| \left[ \int_a^b \left| \int_a^b K(t, s)(y_1(s) - y_2(s)) ds \right|^2 dt \right]^{\frac{1}{2}}
    $$

    $$
    \le |\lambda| \left[ \int_a^b \int_a^b |K(t, s)|^2 ds dt \right]^{\frac{1}{2}} \left[ \int_a^b |y_1(s) - y_2(s)|^2 ds \right]^{\frac{1}{2}}
    $$

    $$
    = \theta \rho(y_1, y_2)
    $$

    由于 $\theta < 1$，故 $T$ 为压缩映射。由压缩映射原理可知，当参数 $\lambda$ 的模充分小时，方程在 $L^2[a, b]$ 中存在唯一的解 $y$。 $\square$

---

## 3. 推广的压缩映射原理

有时，映射 $T$ 本身并不满足严格的压缩映射条件，但它的某次迭代映射满足该条件（例如后文将介绍的 Volterra 积分方程）。此时需要对压缩映射原理加以推广。

记 $T^2x = T(Tx)$，一般地，对任何自然数 $n$，归纳定义 $T^nx = T(T^{n-1}x)$。

!!! success "定理 4.2 (推广的压缩映射原理)"

    设 $X$ 是完备的距离空间，距离为 $\rho$。$T$ 是由 $X$ 到其自身的映射。如果存在一个自然数 $n_0$，使得对任意的 $x, y \in X$，不等式：

    $$
    \rho(T^{n_0}x, T^{n_0}y) \le \theta \rho(x, y)
    $$

    成立，其中 $\theta$ 是满足 $0 \le \theta < 1$ 的常数。
    
    那么，$T$ 在 $X$ 中存在**唯一**的不动点。

??? proof "推广的压缩映射原理的证明"

    **1. 存在性：**
    
    由已知条件可知，映射 $T^{n_0}$ 满足标准的压缩映射原理。由于 $X$ 完备，故映射 $T^{n_0}$ 在 $X$ 中存在**唯一的不动点**，记为 $\overline{x}$。即：

    $$
    T^{n_0}\overline{x} = \overline{x}
    $$

    我们需要证明这个 $\overline{x}$ 也是原映射 $T$ 的不动点。我们在上式两边同时作用算子 $T$：

    $$
    T(T^{n_0}\overline{x}) = T\overline{x}
    $$

    利用算子结合律：

    $$
    T^{n_0}(T\overline{x}) = T^{n_0+1}(\overline{x}) = T(T^{n_0}\overline{x}) = T\overline{x}
    $$

    这说明 $T\overline{x}$ 这个元素，经过 $T^{n_0}$ 映射后等于自身。即 $T\overline{x}$ 也是映射 $T^{n_0}$ 的不动点。
    
    但是，根据前面的结论，$T^{n_0}$ 的不动点是**唯一的**，所以必定有：

    $$
    T\overline{x} = \overline{x}
    $$

    这证明了 $\overline{x}$ 确实也是映射 $T$ 的不动点。

    **2. 唯一性：**
    
    假设 $T$ 存在另一个不动点 $\overline{y}$，满足 $T\overline{y} = \overline{y}$。那么进行多次迭代自然也有：

    $$
    T^{n_0}\overline{y} = T^{n_0-1}(T\overline{y}) = T^{n_0-1}\overline{y} = \dots = \overline{y}
    $$

    说明 $\overline{y}$ 也是 $T^{n_0}$ 的不动点。由于 $T^{n_0}$ 不动点的唯一性，必然得出 $\overline{y} = \overline{x}$。 $\square$

---

## 4. 推广定理的应用：Volterra 积分方程解的存在性

对于 Volterra 积分方程，其积分上限是变量 $t$，而不是固定的常数 $b$。

!!! success "定理 4.3"

    设核函数 $K(t, s)$ 是定义在三角形区域 $a \le t \le b$, $a \le s \le t$ 上的连续函数。
    
    则 Volterra 积分方程：

    $$
    x(t) = f(t) + \lambda \int_a^t K(t, s)x(s) ds
    $$

    对任何给定的 $f \in C[a, b]$ 以及**任何常数 $\lambda \ne 0$**，在空间 $C[a, b]$ 中存在唯一的解 $x$。

??? proof "证明：利用推广的压缩映射原理（点击展开）"

    我们在完备空间 $C[a, b]$ 上定义映射 $T$：

    $$
    (Tx)(t) = f(t) + \lambda \int_a^t K(t, s)x(s) ds
    $$

    记核函数的最大值为 $M = \max_{a \le t \le b, a \le s \le t} |K(t, s)|$。

    对于任意的 $y_1, y_2 \in C[a, b]$，计算经过一次映射后的偏差：

    $$
    |(Ty_1)(t) - (Ty_2)(t)| = |\lambda| \left| \int_a^t K(t, s)(y_1(s) - y_2(s)) ds \right|
    $$

    $$
    \le |\lambda| M \int_a^t |y_1(s) - y_2(s)| ds
    $$

    $$
    \le |\lambda| M (t - a) \max_{s \in [a, b]} |y_1(s) - y_2(s)| = |\lambda| M (t - a) \rho(y_1, y_2)
    $$

    注意这里得到的是带有变量 $t$ 的上界。当 $\lambda$ 较大或者区间较长时，$|\lambda| M (t - a)$ 可能不小于 1，所以 $T$ 并不一定是压缩映射。
    
    但我们可以继续计算第二次、第三次迭代，利用数学归纳法证明如下不等式：

    $$
    |(T^n y_1)(t) - (T^n y_2)(t)| \le \frac{(|\lambda| M (t - a))^n}{n!} \rho(y_1, y_2)
    $$

    *归纳推导过程：假设对于 $n$ 成立，对于 $n+1$ 步：*

    $$
    |(T^{n+1} y_1)(t) - (T^{n+1} y_2)(t)| = |\lambda| \left| \int_a^t K(t, s) (T^n y_1(s) - T^n y_2(s)) ds \right|
    $$

    $$
    \le |\lambda| M \int_a^t \frac{(|\lambda| M (s - a))^n}{n!} \rho(y_1, y_2) ds
    $$

    由于 $\int_a^t (s - a)^n ds = \frac{(t - a)^{n+1}}{n+1}$，代入上式即得：

    $$
    \le \frac{(|\lambda| M (t - a))^{n+1}}{(n+1)!} \rho(y_1, y_2)
    $$

    *归纳法证明完毕。*

    由于对于任意的 $t \in [a, b]$，都有 $(t-a) \le (b-a)$，我们在两侧对 $t$ 取上确界，可以得到距离：

    $$
    \rho(T^n y_1, T^n y_2) \le \frac{(|\lambda| M (b - a))^n}{n!} \rho(y_1, y_2)
    $$

    由微积分知识可知，对于任意固定的常数 $C = |\lambda| M (b - a)$，随着 $n \rightarrow \infty$，有：

    $$
    \lim_{n \rightarrow \infty} \frac{C^n}{n!} = 0
    $$

    因此，无论 $\lambda$ 有多大，总是可以选取足够大的 $n_0$，使得：

    $$
    \theta = \frac{(|\lambda| M (b - a))^{n_0}}{n_0!} < 1
    $$

    这意味着映射 $T^{n_0}$ 是一个严格的**压缩映射**。根据推广的压缩映射原理，映射 $T$ 在 $C[a, b]$ 中存在唯一的不动点，即原 Volterra 积分方程对任何参数 $\lambda$ 都存在唯一解。 $\square$