# 第一章：度量空间

在数学分析中，我们熟悉了实数系 $\mathbb{R}$ 和欧氏空间 $\mathbb{R}^n$ 中的极限与连续性，而这些概念的本质都依赖于“距离”。泛函分析的第一步，就是将“距离”这一直观概念抽象化，推广到包含函数、数列等更复杂元素的集合中去。本章将介绍距离空间的基本公理，并引入几个泛函分析中极为核心的具体空间与经典不等式。

---

## 1. 距离空间的定义

!!! info "定义 1.1 (距离空间 Metric Space)"

    设 $X$ 为一非空集合。如果对于 $X$ 中任给的两个元素 $x, y$，均有一个确定的实数与它们对应，记为 $\rho(x, y)$，且满足以下三个条件（距离公理）：

    **1. 非负性与正定性：**
    
    \[
    \rho(x, y) \ge 0
    \]
    
    且 $\rho(x, y) = 0$ 的充分必要条件是 $x = y$。

    **2. 对称性：**
    
    \[
    \rho(x, y) = \rho(y, x)
    \]

    **3. 三角不等式：**
    
    \[
    \rho(x, y) \le \rho(x, z) + \rho(z, y), \quad \forall z \in X
    \]

    则称 $\rho$ 是集合 $X$ 上的一个**距离 (Metric)**，称 $X$ 是以 $\rho$ 为距离的**距离空间**，记为 $(X, \rho)$。距离空间中的元素通常称为**点**。

**子空间与离散空间：**

* **子空间**：若 $A$ 为 $X$ 的非空子集，则 $A$ 按照相同的距离 $\rho$ 自身也构成一个距离空间，称为 $X$ 的子空间。
* **离散距离空间**：在任何非空集合 $X$ 上，都可以定义一种平凡的距离：当 $x = y$ 时 $\rho(x,x)=0$；当 $x \ne y$ 时 $\rho(x,y)=1$。这构成了离散距离空间。

---

## 2. 泛函分析中的经典空间实例

### 2.1 欧氏空间 $\mathbb{R}^n$ 与 $\mathbb{C}^n$

对于 $n$ 维欧氏空间 $\mathbb{R}^n$，对于任意两点 $x = (x_1, \dots, x_n)$ 和 $y = (y_1, \dots, y_n)$，其标准距离定义为：

\[
\rho(x, y) = \left( \sum_{j=1}^n (x_j - y_j)^2 \right)^{\frac{1}{2}}
\]

其三角不等式的证明直接依赖于经典代数中的 **Cauchy 不等式**。
复空间 $\mathbb{C}^n$ 也可以定义类似的欧氏距离，或者定义最大值距离 $\rho_\infty(x, y) = \max_j |x_j - y_j|$。

### 2.2 连续函数空间 $C[a,b]$

由闭区间 $[a,b]$ 上所有实（或复）连续函数构成的集合。对于任意 $x, y \in C[a,b]$，定义距离为最大一致偏差：

\[
\rho(x, y) = \max_{t \in [a,b]} |x(t) - y(t)|
\]

??? proof "证明：$C[a,b]$ 满足三角不等式"

    设 $x, y, z \in C[a,b]$。对于任意固定的 $t \in [a,b]$，根据实数绝对值的三角不等式：

    \[
    |x(t) - y(t)| \le |x(t) - z(t)| + |z(t) - y(t)|
    \]

    将右侧各项放大为其在整个区间上的最大值：

    \[
    |x(t) - y(t)| \le \max_{t \in [a,b]} |x(t) - z(t)| + \max_{t \in [a,b]} |z(t) - y(t)| = \rho(x, z) + \rho(z, y)
    \]

    由于上式右端是一个不依赖于 $t$ 的常数，且对所有 $t \in [a,b]$ 都成立，因此对其左端取最大值依然成立：

    \[
    \max_{t \in [a,b]} |x(t) - y(t)| \le \rho(x, z) + \rho(z, y)
    \]

    即 $\rho(x, y) \le \rho(x, z) + \rho(z, y)$。 $\square$

### 2.3 数列空间 $l^p$ 与 $l^\infty$

* **空间 $l^p$ ($1 \le p < \infty$)**：由所有满足 $\sum_{j=1}^\infty |x_j|^p < \infty$ 的数列 $x = \{x_1, x_2, \dots\}$ 构成。其距离定义为：

    \[
    \rho(x, y) = \left( \sum_{j=1}^\infty |x_j - y_j|^p \right)^{\frac{1}{p}}
    \]

* **空间 $l^\infty$**：由一切**有界**数列构成，距离定义为上确界：

    \[
    \rho(x, y) = \sup_{j} |x_j - y_j|
    \]

为了证明 $l^p$ 空间满足三角不等式，我们必须引入泛函分析中最基石的两个不等式：Hölder 不等式与 Minkowski 不等式。

---

## 3. 经典不等式及其证明 (Hölder & Minkowski)

!!! success "定理 3.1 (Hölder 不等式)"

    设 $p > 1, q > 1$ 且满足共轭指数关系 $\frac{1}{p} + \frac{1}{q} = 1$。
    对于任意数列 $a = \{a_1, a_2, \dots\} \in l^p$ 和 $b = \{b_1, b_2, \dots\} \in l^q$，有：

    \[
    \sum_{j=1}^\infty |a_j b_j| \le \left( \sum_{j=1}^\infty |a_j|^p \right)^{\frac{1}{p}} \cdot \left( \sum_{j=1}^\infty |b_j|^q \right)^{\frac{1}{q}}
    \]

??? proof "Hölder 不等式的证明（点击展开）"

    **Step 1: 引入 Young 不等式**
    
    考虑函数 $f(x) = x^\alpha$ (其中 $0 < \alpha < 1, x \ge 0$)。由于 $f''(x) = \alpha(\alpha-1)x^{\alpha-2} < 0$，该函数是上凸的。它在点 $(1,1)$ 处的切线方程为 $y = \alpha(x-1) + 1 = \alpha x + \beta$ （其中 $\beta = 1 - \alpha$）。
    由上凸性知函数图像位于切线下方，故：

    \[
    x^\alpha \le \alpha x + \beta
    \]

    令 $x = \frac{u}{v} \; (u, v > 0)$，代入上式同乘 $v$，得到 **Young 不等式**：

    \[
    u^\alpha v^\beta \le \alpha u + \beta v \quad (\alpha + \beta = 1)
    \]

    **Step 2: 构造单位化变量**
    
    令 $\alpha = \frac{1}{p}, \beta = \frac{1}{q}$。对于任意给定的 $j$，取：

    \[
    u = \frac{|a_j|^p}{\sum_{j=1}^\infty |a_j|^p}, \quad v = \frac{|b_j|^q}{\sum_{j=1}^\infty |b_j|^q}
    \]

    **Step 3: 应用 Young 不等式并求和**
    
    代入 Young 不等式，得到：

    \[
    \frac{|a_j b_j|}{ \left( \sum |a_j|^p \right)^{\frac{1}{p}} \left( \sum |b_j|^q \right)^{\frac{1}{q}} } \le \frac{1}{p} \frac{|a_j|^p}{\sum |a_j|^p} + \frac{1}{q} \frac{|b_j|^q}{\sum |b_j|^q}
    \]

    对上式两边关于 $j$ 从 $1$ 到 $\infty$求和：

    \[
    \frac{\sum_{j=1}^\infty |a_j b_j|}{ \left( \sum |a_j|^p \right)^{\frac{1}{p}} \left( \sum |b_j|^q \right)^{\frac{1}{q}} } \le \frac{1}{p} \cdot 1 + \frac{1}{q} \cdot 1 = 1
    \]

    移项即得证 Hölder 不等式。 $\square$

!!! success "定理 3.2 (Minkowski 不等式)"

    对于 $p \ge 1$，任取 $l^p$ 中的两个元素 $a$ 和 $b$，则 $a+b \in l^p$，且满足：

    \[
    \left( \sum_{j=1}^\infty |a_j + b_j|^p \right)^{\frac{1}{p}} \le \left( \sum_{j=1}^\infty |a_j|^p \right)^{\frac{1}{p}} + \left( \sum_{j=1}^\infty |b_j|^p \right)^{\frac{1}{p}}
    \]

    *注：这本质上就是 $l^p$ 空间满足三角不等式的核心证明。*

??? proof "Minkowski 不等式的证明（点击展开）"

    **Step 1: 证明封闭性 ($a+b \in l^p$)**
    
    利用初等不等式 $|a_j + b_j|^p \le (2 \max\{|a_j|, |b_j|\})^p \le 2^p (|a_j|^p + |b_j|^p)$，对其求和可知只要 $a, b \in l^p$，和序列必定也绝对收敛，即 $a+b \in l^p$。

    **Step 2: 拆分与 Hölder 不等式**
    
    当 $p > 1$ 时（$p=1$ 时由绝对值三角不等式直接成立），我们利用共轭指数关系 $(p-1)q = p$ 进行恒等变形：

    \[
    \sum_{j=1}^\infty |a_j + b_j|^p = \sum_{j=1}^\infty |a_j + b_j| \cdot |a_j + b_j|^{p-1} \le \sum_{j=1}^\infty (|a_j| + |b_j|) |a_j + b_j|^{p-1}
    \]

    展开为两项后，分别对这两项应用刚才证明的 **Hölder 不等式**：

    \[
    \sum_{j=1}^\infty |a_j| |a_j + b_j|^{p-1} \le \left( \sum_{j=1}^\infty |a_j|^p \right)^{\frac{1}{p}} \left( \sum_{j=1}^\infty |a_j + b_j|^{(p-1)q} \right)^{\frac{1}{q}}
    \]

    因为 $(p-1)q = p$，所以右边第二个因子就是 $\left( \sum |a_j + b_j|^p \right)^{1/q}$。同理对 $|b_j|$ 项处理。

    **Step 3: 合并同类项**
    
    将两项相加，提取公因式：

    \[
    \sum_{j=1}^\infty |a_j + b_j|^p \le \left[ \left( \sum_{j=1}^\infty |a_j|^p \right)^{\frac{1}{p}} + \left( \sum_{j=1}^\infty |b_j|^p \right)^{\frac{1}{p}} \right] \left( \sum_{j=1}^\infty |a_j + b_j|^p \right)^{\frac{1}{q}}
    \]

    两边同时除以 $\left( \sum |a_j + b_j|^p \right)^{\frac{1}{q}}$，注意到 $1 - \frac{1}{q} = \frac{1}{p}$，即得：

    \[
    \left( \sum_{j=1}^\infty |a_j + b_j|^p \right)^{\frac{1}{p}} \le \left( \sum_{j=1}^\infty |a_j|^p \right)^{\frac{1}{p}} + \left( \sum_{j=1}^\infty |b_j|^p \right)^{\frac{1}{p}}
    \]
    
    证明完毕。 $\square$

*(注：对于可测函数空间 $L^p(F)$，将上述的求和号 $\sum$ 替换为在可测集 $F$ 上的勒贝格积分 $\int_F$，即可得到相应的积分形式的 Hölder 与 Minkowski 不等式。)*

---

## 4. 距离空间中的收敛性

在建立起距离的度量标准后，极限操作得以被自然地定义。

!!! info "定义 4.1 (点列的收敛)"

    设 $\{x_n\}$ 为距离空间 $(X, \rho)$ 中的一个点列。如果存在 $X$ 中的点 $y$，使得当 $n \rightarrow \infty$ 时，有：

    \[
    \rho(x_n, y) \rightarrow 0
    \]

    则称点列 $\{x_n\}$ 收敛到 $y$，记为 $x_n \rightarrow y$ 或 $\lim_{n \rightarrow \infty} x_n = y$。$y$ 称为该点列的极限。

### 4.1 收敛点列的基本性质

对于距离空间中的收敛点列 $\{x_n\}$，它天然继承了实数轴上极限的若干优良性质：

1. **极限唯一性**：如果 $x_n \rightarrow y$ 且 $x_n \rightarrow y'$，则必定有 $y = y'$。
2. **有界性**：对于任意给定的点 $y_0 \in X$，实数列 $\rho(x_n, y_0)$ 是有界的。
3. **子列继承性**：$\{x_n\}$ 的任一子列也必然收敛，且收敛于同一个极限。反之，若任一子列均收敛，则原点列必定收敛。

??? proof "证明：极限的唯一性与有界性"

    **1. 唯一性证明：**
    
    利用三角不等式：

    \[
    0 \le \rho(y, y') \le \rho(y, x_n) + \rho(x_n, y') = \rho(x_n, y) + \rho(x_n, y')
    \]

    由于 $x_n \rightarrow y$ 且 $x_n \rightarrow y'$，当 $n \rightarrow \infty$ 时，右端趋于 0。由夹逼定理得 $\rho(y, y') = 0$，根据距离公理第一条，得 $y = y'$。

    **2. 有界性证明：**
    
    设 $x_n \rightarrow y$，即 $\rho(x_n, y) \rightarrow 0$。由收敛定义，必存在 $M > 0$ 使得 $\rho(x_n, y) < M$ 对全体 $n$ 恒成立。
    任取 $y_0 \in X$，由三角不等式：

    \[
    \rho(x_n, y_0) \le \rho(x_n, y) + \rho(y, y_0) < M + \rho(y, y_0)
    \]

    右侧是一个不依赖于 $n$ 的常数，故实数列 $\rho(x_n, y_0)$ 有界。 $\square$

### 4.2 不同空间中的“收敛”等价性辨析

泛函分析的精妙之处在于：**空间的元素相同，但赋予的距离不同，导致的“收敛”概念可能完全不同。**

* **在有限维欧氏空间 $\mathbb{R}^n$ 中：** 点列按欧氏距离收敛的 **充分必要条件** 是它的每个坐标（分量）分别在实数域中收敛。此时距离空间的收敛退化为逐点收敛。

* **在连续函数空间 $C[a,b]$ 中：**
  * 若赋予最大值距离 $\rho(x, y) = \max_{t \in [a,b]} |x(t) - y(t)|$，点列收敛等价于函数列在 $[a,b]$ 上的 **一致收敛 (Uniform Convergence)**。
  * **反例：** 若在 $C[a,b]$ 中赋予均方误差距离 $\rho_1(x, y) = \left( \int_a^b |x(t) - y(t)|^2 dt \right)^{1/2}$。考虑 $[0,1]$ 上的函数列 $x_n(t) = t^n$。
    显然 $\rho_1(x_n, 0) \rightarrow 0$（按积分距离收敛于 $0$），但作为函数列它在 $[0,1]$ 上并不一致收敛于 $0$（在 $t=1$ 处始终为 $1$）。

这说明，距离的具体定义形式深刻影响着空间内部的拓扑结构和极限行为，这也是后续我们要研究拓扑性质与完备性的重要动机。