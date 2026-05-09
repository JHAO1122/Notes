# 5.1 $L^p$ 空间、不等式与完备性

在建立了完备的 Lebesgue 积分理论之后，我们将探讨由可积函数构成的函数空间——$L^p$ 空间。它们不仅是实分析的结晶，更是泛函分析的核心研究对象。本节我们将重点讨论 $L^p$ 空间中几个极为重要的不等式及其结构包含关系。

---

## 1. Chebyshev 不等式

在探讨 $L^p$ 空间之前，我们先引入测度论中最基础的尾概率控制不等式——Chebyshev 不等式（又称 Markov-Chebyshev 不等式）。它建立了函数的 $L^p$ 积分与其水平集测度之间的联系。

!!! info "定理 5.1.1 (Chebyshev 不等式)"

    设 $f \in L^p(X, \mu)$（其中 $0 < p < \infty$），对于任意 $\alpha > 0$，定义水平集 $E_\alpha = \{x \in X : |f(x)| > \alpha\}$。那么有：

    \[
    \mu(\{x \in X : |f(x)| > \alpha\}) \le \left( \frac{\|f\|_p}{\alpha} \right)^p = \frac{1}{\alpha^p} \int_X |f|^p \, d\mu
    \]

??? proof "Chebyshev 不等式的证明"

    利用特征函数和积分的单调性，证明非常直接。
    
    显然，在集合 $E_\alpha$ 上，我们有 $|f(x)| > \alpha$，从而 $|f(x)|^p > \alpha^p$。
    
    于是我们可以放大积分：

    \[
    \int_X |f|^p \, d\mu \ge \int_{E_\alpha} |f|^p \, d\mu \ge \int_{E_\alpha} \alpha^p \, d\mu = \alpha^p \mu(E_\alpha)
    \]

    两边同除以 $\alpha^p$ 即可得到结论。这说明，如果一个函数在 $L^p$ 空间中，那么它取大值的区域的测度必然衰减得非常快。

---

## 2. $L^p$ 空间的包含与插值关系

不同指数 $p$ 对应的 $L^p$ 空间之间存在着深刻的联系。它们之间的包含关系很大程度上取决于底层测度空间的有限性。

### 2.1 有限测度下的包含关系

如果测度空间是有限的，低阶可积性由高阶可积性保证。

!!! success "定理 5.1.2 (有限测度的 $L^p$ 包含)"

    如果 $\mu(A) < \infty$，且 $0 < p < q \le \infty$，那么：

    \[
    L^q(A) \subset L^p(A)
    \]

    并且存在常数 $C$ 使得 $\|f\|_p \le C \|f\|_q$。

    *(提示：这可以直接通过 Hölder 不等式，将 $\int |f|^p \cdot 1 \, d\mu$ 中的常数 1 视为共轭函数来证明。)*

### 2.2 空间的代数分解与插值

当空间测度无限时（如 $\mathbb{R}^n$），$L^p$ 和 $L^q$ 之间没有绝对的包含关系。但在 $0 < p < q < r \le \infty$ 的条件下，存在以下优美的代数关系。

!!! success "定理 5.1.3 ($L^q$ 函数的分解)"

    对于 $0 < p < q < r \le \infty$，任何处于中间空间 $L^q$ 的函数都可以被分解为属于 $L^p$ 和 $L^r$ 的两部分，即：

    \[
    L^q \subset L^p + L^r
    \]

??? proof "分解定理的证明"

    任取 $f \in L^q$。我们通过高度阈值 1 来截断函数。
    
    令 $E = \{x \in X : |f(x)| \ge 1\}$。定义分解 $f = g + h$，其中：

    * 尖峰部分：$g = f \cdot \chi_E$

    * 平缓部分：$h = f \cdot (1 - \chi_E)$

    **证明 $g \in L^p$**：
    由于在 $E$ 上 $|f| \ge 1$，且 $p < q$，所以 $|f|^p \le |f|^q$。

    \[
    \int_X |g|^p \, d\mu = \int_E |f|^p \, d\mu \le \int_E |f|^q \, d\mu \le \int_X |f|^q \, d\mu < \infty
    \]

    所以 $g \in L^p$。

    **证明 $h \in L^r$**：
    由于在 $E^c$ 上 $|f| < 1$，且 $q < r$，所以 $|f|^r \le |f|^q$。

    \[
    \int_X |h|^r \, d\mu = \int_{E^c} |f|^r \, d\mu \le \int_{E^c} |f|^q \, d\mu \le \int_X |f|^q \, d\mu < \infty
    \]

    所以 $h \in L^r$。由此得证 $f = g + h \in L^p + L^r$。

!!! success "定理 5.1.4 ($L^p$ 空间的插值性质 / 绝对凸性)"

    对于 $0 < p < q < r \le \infty$，两个极端空间的交集必然包含中间空间：

    \[
    L^p \cap L^r \subset L^q
    \]

    更定量地，存在 $\lambda \in (0, 1)$ 满足 $\frac{1}{q} = \frac{\lambda}{p} + \frac{1-\lambda}{r}$（即 $\frac{\lambda q}{p} + \frac{(1-\lambda)q}{r} = 1$），此时对任意 $f \in L^p \cap L^r$，可以通过 Hölder 不等式得到**对数凸性界 (Log-convexity bound)**：

    \[
    \|f\|_q \le \|f\|_p^\lambda \|f\|_r^{1-\lambda}
    \]

---

## 3. 积分算子与 Schur 检验

泛函分析中，我们经常研究由核函数 (Kernel) 定义的积分算子。Schur 检验给出了这类算子在 $L^p$ 空间上有界的一个非常实用的充分条件。

!!! info "定义 5.1.2 (积分算子)"

    设 $(X, \mathcal{M}, \mu)$ 和 $(Y, \mathcal{N}, \nu)$ 为有限或 $\sigma$-有限测度空间。$k: X \times Y \to \mathbb{C}$ 为乘积可测函数。定义积分算子 $T$ 为：

    \[
    Tf(x) = \int_Y k(x, y) f(y) \, d\nu(y)
    \]

!!! success "定理 5.1.5 (Schur 检验 / 广义 Young 不等式)"

    假设存在常数 $C > 0$，使得核函数 $k$ 满足以下两个一致行/列积分为界的条件：

    * 对于几乎所有的 $y \in Y$：

    \[
    \int_X |k(x, y)| \, d\mu(x) \le C
    \]

    * 对于几乎所有的 $x \in X$：

    \[
    \int_Y |k(x, y)| \, d\nu(y) \le C
    \]

    那么对于任意 $1 \le p \le \infty$，算子 $T$ 是从 $L^p(\nu)$ 到 $L^p(\mu)$ 的有界线性算子，且算子范数满足：

    \[
    \|Tf\|_{L^p(\mu)} \le C \|f\|_{L^p(\nu)}
    \]

??? proof "Schur 检验的证明核心 (基于 Hölder 不等式)"

    我们利用 Hölder 不等式的巧妙拆分（针对 $1 < p < \infty$ 的情形，共轭指数 $p'$ 满足 $1/p + 1/p' = 1$）。

    将积分中的 $|k(x, y)|$ 拆分为 $|k(x,y)|^{1/p'} \cdot |k(x,y)|^{1/p}$：

    \[
    |Tf(x)| \le \int_Y |k(x, y)|^{1/p'} \cdot \left(|k(x, y)|^{1/p} |f(y)|\right) d\nu(y)
    \]

    关于测度 $\nu$，对这两项应用 Hölder 不等式：

    \[
    |Tf(x)|^p \le \left( \int_Y |k(x, y)| \, d\nu(y) \right)^{p/p'} \left( \int_Y |k(x, y)| |f(y)|^p \, d\nu(y) \right)
    \]

    利用第二个条件 $\int_Y |k(x, y)| d\nu(y) \le C$，代入并化简（注意到 $p/p' = p-1$）：

    \[
    |Tf(x)|^p \le C^{p-1} \int_Y |k(x, y)| |f(y)|^p \, d\nu(y)
    \]

    现在对两边关于 $x$ 在 $X$ 上积分，并使用 Tonelli 定理交换积分顺序：

    \[
    \int_X |Tf(x)|^p \, d\mu(x) \le C^{p-1} \int_Y |f(y)|^p \left( \int_X |k(x, y)| \, d\mu(x) \right) d\nu(y)
    \]

    利用第一个条件 $\int_X |k(x, y)| d\mu(x) \le C$：

    \[
    \|Tf\|_{L^p(\mu)}^p \le C^{p-1} \cdot C \int_Y |f(y)|^p \, d\nu(y) = C^p \|f\|_{L^p(\nu)}^p
    \]

    两边开 $p$ 次方即得证 $\|Tf\|_p \le C \|f\|_p$。

---

## 4. 广义 Minkowski 积分不等式

经典的 Minkowski 不等式确立了 $L^p$ 空间满足三角不等式 ($\|f + g\|_p \le \|f\|_p + \|g\|_p$)，从而使其成为一个赋范向量空间。

而在乘积测度空间中，这一不等式可以推广到对连续参数求积分的形式，即“积分的 $L^p$ 范数小于等于 $L^p$ 范数的积分”。

!!! success "定理 5.1.6 (广义 Minkowski 积分不等式)"

    设 $f(x, y)$ 是 $X \times Y$ 上的非负乘积可测函数，$1 \le p < \infty$。那么有：

    \[
    \left[ \int_X \left( \int_Y f(x, y) \, d\nu(y) \right)^p d\mu(x) \right]^{1/p} \le \int_Y \left( \int_X f(x, y)^p \, d\mu(x) \right)^{1/p} d\nu(y)
    \]

    用算子范数的语言更直观地表达为：假设 $y \mapsto \|f(\cdot, y)\|_{L^p(\mu)}$ 在 $Y$ 上是可积的，则函数 $F(x) = \int_Y f(x, y) d\nu(y)$ 在 $L^p(\mu)$ 中，并且满足：

    \[
    \left\| \int_Y f(\cdot, y) \, d\nu(y) \right\|_{L^p(\mu)} \le \int_Y \|f(\cdot, y)\|_{L^p(\mu)} \, d\nu(y)
    \]

这个深刻的不等式在处理偏微分方程、算子半群理论以及卷积算子的范数估计时，是一个不可或缺的工具。