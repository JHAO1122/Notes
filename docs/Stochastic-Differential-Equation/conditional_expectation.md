# 第一章：条件期望与鞅基础

在开始随机微分方程前，我们需要抛弃初等概率论中条件期望的定义方法，用测度论的语言严格定义条件期望。这是我们后续探索鞅理论、随机分析的基础。

## 1. 条件期望的测度论定义

在古典概率论中，定义事件 $A$ 在事件 $B$ 发生的条件下的条件概率为 $P(A|B) \triangleq \frac{P(AB)}{P(B)}.$

相应的条件期望为 $E[X|B] \triangleq \frac{1}{P(B)}\int_B X dP.$

但当给定的是一个连续随机变量 $Y$ 时，事件 $\{Y=y\}$ 的概率为 0，古典定义失效。

此时我们需要借助 Radon-Nikodym 定理在 $\sigma$-代数层面重新定义。

!!! abstract "定义：条件期望 (Conditional Expectation)"

    设 $X$ 为概率空间 $(\Omega, \mathcal{F}, P)$ 上的可积随机变量（即 $E[|X|] < \infty$），$\mathcal{G}$ 是 $\mathcal{F}$ 的一个子 $\sigma$-代数（常记为 $\mathcal{G} = \sigma(Y)$）。
    我们称随机变量 $Z$ 为 $X$ 给定 $\mathcal{G}$ 下的**条件期望**，记作 $E[X|\mathcal{G}] = Z$，如果它满足以下两个条件：

    **1. 可测性**：$Z$ 是 $\mathcal{G}$-可测的随机变量。

    **2. 局部平均相等（积分匹配）**：对于任意的 $A \in \mathcal{G}$，都有：

    $$
    \int_A X dP = \int_A Z dP
    $$

    由 Radon-Nikodym 定理保证，这样的 $Z$ 存在且在 almost sure (a.s.) 的意义下唯一。

## 2. 核心性质

设 $\mathcal{G}, \mathcal{H}$ 为子 $\sigma$-代数。

!!! info "基本性质"

    **1. 线性性**：

    $$
    E[aX + bY | \mathcal{G}] = aE[X|\mathcal{G}] + bE[Y|\mathcal{G}] \quad a.s.
    $$

    **2. 保期望性**：

    $$
    E[E[X|\mathcal{G}]] = E[X]
    $$

    > （只需在定义中取 $A = \Omega$ 即可证明）。

    **3. 已知即常数（提出已知因子）**：

    如果 $X$ 是 $\mathcal{G}$-可测的，那么：

    $$
    E[XY|\mathcal{G}] = X E[Y|\mathcal{G}] \quad a.s.
    $$

    **4. 独立则无关**：

    如果 $X$ 与 $\mathcal{G}$ 独立，那么：

    $$
    E[X|\mathcal{G}] = E[X] \quad a.s.
    $$

    **5. 塔牌性质 (Tower Property)**：

    若 $\mathcal{H} \subset \mathcal{G} \subset \mathcal{F}$，则：

    $$
    E\big[ E[X|\mathcal{G}] \big| \mathcal{H} \big] = E[X|\mathcal{H}] \quad a.s.
    $$

    > *(直观理解：小信息集说了算。)*

!!! tip "几何解释：$L^2$ 空间上的正交投影"

    条件期望有一个极其优美的几何直观。

    考虑平方可积空间 $L^2(\Omega, \mathcal{F}, P)$，它是一个 Hilbert 空间，内积定义为：

    $$
    (X,Y) = E[XY]
    $$

    由于 $\mathcal{G} \subset \mathcal{F}$，子空间 $L^2(\Omega, \mathcal{G}, P)$ 也是一个闭子空间。

    在这个视角下，**条件期望 $E[X|\mathcal{G}]$ 实际上就是元素 $X$ 在子空间 $L^2(\Omega, \mathcal{G}, P)$ 上的正交投影**。

    它使得均方误差 $E[(X - Z)^2]$ 在所有 $\mathcal{G}$-可测的随机变量 $Z$ 中达到最小：

    $$
    \min_{Z \in \mathcal{G}} E[(X - Z)^2] = E\big[(X - E[X|\mathcal{G}])^2\big]
    $$

## 3. Jensen 不等式

在高级概率论中，Jensen 不等式是我们处理凸函数和极限定理的利器。

!!! success "定理：条件 Jensen 不等式"

    设 $\Phi: \mathbb{R} \rightarrow \mathbb{R}$ 为凸函数，且 $E[|\Phi(X)|] < \infty$，则有：

    $$
    \Phi(E[X|\mathcal{G}]) \le E[\Phi(X)|\mathcal{G}] \quad a.s.
    $$

    **证明思路 (基于简单函数逼近)**：
    利用凸函数的性质 $\Phi(\sum a_i b_i) \le \sum b_i \Phi(a_i)$（其中 $\sum b_i = 1$）。由于条件期望可以看作是一种加权平均，我们可以用指示函数构造简单函数来逼近，最终通过极限推广至一般情况。这是实变函数/测度论中常用的证明方式，指示函数到简单函数再到可测函数。

## 4. 鞅 (Martingale) 与 Doob 不等式

为了研究随时间演化的随机现象（如布朗运动），我们引入信息流和鞅的概念。

!!! abstract "定义：流 (Filtration) 与 鞅 (Martingale)"

    **1. 信息流 (Filtration)**：一族单调递增的 $\sigma$-代数 $\{\mathcal{F}_t\}_{t \ge 0}$，即当 $s \le t$ 时，$\mathcal{F}_s \subset \mathcal{F}_t$。它代表了随着时间 $t$ 累积的信息历史。

    **2. 鞅**：设随机过程 $\{X_t\}$ 适应于流 $\mathcal{F}_t$，且 $E[|X_t|] < \infty$。如果对于任意 $s \le t$，都有：

    $$
    E[X_t | \mathcal{F}_s] = X_s \quad a.s.
    $$

    则称 $\{X_t\}$ 为一个**鞅**。（代表公平游戏，未来期望等于现在）。

    **3. 上下鞅**：若 $E[X_t | \mathcal{F}_s] \ge X_s$，称为**下鞅 (Submartingale)**；若 $E[X_t | \mathcal{F}_s] \le X_s$，称为**上鞅 (Supermartingale)**。

在随机分析中，我们常常需要控制整个过程在其路径上的最大值。Doob 不等式成为解决此类问题的重要工具。

!!! info "定理：Doob 极大值不等式 (Doob's Maximal Inequality)"

    **(1) 次线性界**：设 $\{X_n\}_{n=1}^N$ 为非负下鞅，则对于任意 $\lambda > 0$：

    $$
    P\left(\max_{1 \le k \le n} X_k \ge \lambda\right) \le \frac{1}{\lambda} E[X_n I_{\{\max X_k \ge \lambda\}}] \le \frac{1}{\lambda} E[X_n]
    $$

    > *(注意它与 Chebyshev 不等式的形式相似，但控制的是整个路径的极值)*

    **(2) $L^p$ 极大不等式**：设 $p > 1, \frac{1}{p} + \frac{1}{q} = 1$。如果 $\{X_n\}$ 是非负下鞅且 $X_n \in L^p$，则：

    $$
    E\left[\max_{1 \le k \le n} X_k^p\right] \le \left(\frac{p}{p-1}\right)^p E[X_n^p]
    $$

    **$L^p$ 边界推导关键步骤补全**：

    证明的核心技巧在于利用 Fubini 定理交换积分顺序，并结合 Hölder 不等式。
    设 $X^* = \max_{1 \le k \le n} X_k$，利用恒等式：

    $$
    E[(X^*)^p] = \int_0^\infty p \lambda^{p-1} P(X^* \ge \lambda) d\lambda
    $$

    代入次线性界后，使用 Fubini 定理将对 $\lambda$ 的积分交换到期望内部：

    $$
    \le \int_{\Omega} X_n \left( \int_0^{X^*} p \lambda^{p-2} d\lambda \right) dP = \frac{p}{p-1} E[X_n (X^*)^{p-1}]
    $$

    最后，利用 Hölder 不等式 $\int f g \le (\int f^p)^{1/p} (\int g^q)^{1/q}$ 剥离 $X_n$ 和 $(X^*)^{p-1}$ 即可证得结论。