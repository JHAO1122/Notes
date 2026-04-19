# 第五章：弱相依数据（一）(Weakly Dependent Data I)

在经典的极限理论中，我们通常假设数据是独立同分布的。然而在实际应用（如金融日度、周度、年度数据等等间隔采样的时间序列）中，数据往往存在序列相关性。本章我们将探讨平稳时间序列、刻画数据相依性的各种**混合系数 (Mixing Coefficients)**，以及在弱相依条件下非常重要的一系列协方差不等式。

设 $Z_1, \dots, Z_n \in \mathbb{R}^d$ 为等间隔采样的时间序列数据，其中 $d$ 为观测向量的维度。

---

## 1. 平稳性与时间序列模型

为了对相依数据进行统计推断，我们通常需要数据满足某种平稳性条件。

!!! info "平稳性的定义 (Stationarity)"

    * **严格平稳 (Strictly Stationary)**：对于任意整数 $l$ 和 $m$，随机向量序列 $(Z_{i_1}, \dots, Z_{i_m})^T$ 与平移后的序列 $(Z_{i_1+l}, \dots, Z_{i_m+l})^T$ 具有完全相同的联合分布（即**强平移不变性**）。
    * **弱平稳 / 二阶平稳 (Weak Stationary / Second Order Stationary)**：过程的一阶和二阶矩仅依赖于时间间隔而不随时间改变，即：
      
    \[
    E(Z_i) = E(Z_{i+1}), \quad Var(Z_i) = Var(Z_{i+l}), \quad Cov(Z_i, Z_j) = Cov(Z_{i+l}, Z_{j+l})
    \]

    （即**弱平移不变性**）。

    *注：在实际操作中，使非平稳时间序列平稳化的常见方法包括：差分 (Difference)、取平方根或对数转换等。*

### 1.1 线性时间序列模型

!!! info "定义 4.1：自回归滑动平均模型 ARMA(p, q)"

    序列 $\{Z_t\}_{t \in \mathbb{Z}}$ 被称为 **ARMA(p, q)** 模型，如果 $\{Z_t\}$ 是弱平稳的，并且对于任意 $t$ 满足：
    
    \[
    Z_t - \theta_1 Z_{t-1} - \dots - \theta_p Z_{t-p} = \epsilon_t - \eta_1 \epsilon_{t-1} - \dots - \eta_q \epsilon_{t-q}
    \]

    其中 $\{\epsilon_t\}$ 是一个独立的白噪声过程 (Independent White Noise)，记为 $WN(0, \sigma^2)$。

设 $\theta(z)$ 和 $\eta(z)$ 分别为 $p$ 阶和 $q$ 阶多项式，定义为：

\[
\theta(z) = 1 - \theta_1 z - \dots - \theta_p z^p, \quad \eta(z) = 1 - \eta_1 z - \dots - \eta_q z^q
\]

引入后移算子 (Backward shift operator) $B$，使得 $B^j Z_t = Z_{t-j}$。

!!! info "定义 4.2：因果性 (Causality)"

    一个 ARMA(p, q) 过程 $\{Z_t\}$ 被称为是**因果的 (causal)**，如果存在绝对可和的系数序列 $\{\psi_j\}_{j=0}^\infty$ （即 $\sum_{j=0}^\infty |\psi_j| < \infty$），使得对于所有 $t \in \mathbb{Z}$，有：
    
    \[
    Z_t = \sum_{j=0}^\infty \psi_j \epsilon_{t-j}
    \]

    *注：因果的 ARMA(p, q) 模型实际上就是一个 $MA(\infty)$ 过程。*

!!! success "定理 4.3 (ARMA 过程的因果条件)"

    设 $\{Z_t\}$ 为 ARMA(p, q) 过程，表示为 $\theta(B)Z_t = \eta(B)\epsilon_t$。假设多项式 $\theta(z)$ 和 $\eta(z)$ 没有公共根，那么 $\{Z_t\}$ 是因果的**当且仅当**：
    
    \[
    \theta(z) \ne 0, \quad \forall z \in \mathbb{C}, \ |z| \le 1
    \]

    即特征多项式的根都严格落在单位圆外。此时系数 $\{\psi_j\}$ 由展开式 $\psi(z) = \eta(z) / \theta(z)$ 决定。

如果允许求和下标扩展到负无穷，我们就可以定义一般化的**线性过程 (Linear Process)**：

\[
Z_t = \sum_{j=-\infty}^\infty \psi_j \epsilon_{t-j}, \quad \epsilon_t \sim i.i.d.\ F(0, \sigma^2)
\]

*(在一定条件下，ARMA 过程属于线性过程的一种特例。)*

### 1.2 ARCH(p) 模型：条件异方差

**自回归条件异方差模型 (Auto-Regression Conditionally Heterogeneity, ARCH)** 是金融时间序列中极其重要的非线性模型：

\[
Z_t = m(\vec{Z}_{t,p}) + \sigma(\vec{Z}_{t,p})\epsilon_t
\]

其中 $\vec{Z}_{t,p} = (Z_{t-1}, \dots, Z_{t-p})^T$。
它在两个层面上推广了经典的 AR(p) 模型（$Z_t = \theta_0 + \theta_1 Z_{t-1} + \dots + \theta_p Z_{t-p} + \epsilon_t$）：
1. 将线性条件均值推广为**非线性条件均值函数** $m(\cdot)$。
2. 将常数条件方差推广为**状态依赖的条件方差函数** $\sigma(\cdot)$。

*注：ARCH 模型及一般的非线性时间序列模型未必总是平稳的，但在一定参数约束下，它们可以保证是“渐近平稳的” (Asymptotically stationary) —— 即在经过一段时间的“预热 (pre-burning)”后趋于平稳。*

---

## 2. 混合系数：刻画相依性的测度

对于弱相依序列，随着时间间隔的拉长，过去的事件对未来的影响应当逐渐衰减为零。我们用混合系数 (Mixing Coefficients) 来精确度量这种衰减速度。

设 $\{Z_1, \dots, Z_t, \dots\}$ 是一个严格平稳过程。对于正整数 $l \le m$，定义 $\mathcal{F}_l^m$ 为由随机变量 $\{Z_i\}_{i=l}^m$ 生成的 $\sigma$-代数。对于时间间隔 $k \ge 1$，定义以下混合系数：

1. **$\alpha$-混合系数 (或强混合, strong mixing)**:
   
    \[
    \alpha(k) = \sup_{B \in \mathcal{F}_{-\infty}^t, C \in \mathcal{F}_{t+k}^\infty} |P(B \cap C) - P(B)P(C)|
    \]

2. **$\beta$-混合系数 (绝对正则, absolute regularity)**:
   
    \[
    \beta(k) = E \left[ \sup_{C \in \mathcal{F}_{t+k}^\infty} |P(C) - P(C | \mathcal{F}_{-\infty}^t)| \right]
    \]

3. **$\phi$-混合系数 (一致混合, uniform mixing)**:
   
    \[
    \phi(k) = \sup_{B \in \mathcal{F}_{-\infty}^t, C \in \mathcal{F}_{t+k}^\infty} |P(C) - P(C | B)|
    \]

4. **$\rho$-混合系数 (最大相关, maximal correlation)**:
   
    \[
    \rho(k) = \sup_{X \in L^2(\mathcal{F}_{-\infty}^t), Y \in L^2(\mathcal{F}_{t+k}^\infty)} |Corr(X, Y)| = \sup_{X, Y} \left| \frac{Cov(X, Y)}{\sqrt{Var(X)Var(Y)}} \right|
    \]
   
   *(其中 $L^2(\mathcal{F})$ 表示在 $\sigma$-代数 $\mathcal{F}$ 上可测且二阶矩有限 ($EX^2 < \infty$) 的全体随机变量集合。)*

### 2.1 混合系数的关系与性质

各种混合系数之间存在如下极其重要的严格不等式关系：

\[
2\alpha(k) \le \beta(k) \le \phi(k), \quad 4\alpha(k) \le \rho(k) \le 2\phi^{1/2}(k)
\]

**定义：**
如果 $\lim_{k \to \infty} \alpha(k) = 0$，则称该过程 $\{Z_t\}$ 是 **$\alpha$-混合的**。同理可定义 $\phi$-混合、$\rho$-混合等。这些混合性质度量了过去 ($\mathcal{F}_{-\infty}^t$) 与未来 ($\mathcal{F}_{t+k}^\infty$) 之间的相依性。当时间间隔 $k \to \infty$ 时，混合系数趋于零，意味着系统展现出**渐近独立性**。

由上述不等式可知，蕴含关系如下：

\[
\phi\text{-混合} \implies \beta\text{-混合}, \quad \rho\text{-混合} \implies \alpha\text{-混合}
\]

*注：$\alpha$-混合在以上所有混合条件中是**最弱的条件**（最容易被满足），但颇具讽刺意味的是，它在历史上被命名为“强混合 (Strong Mixing)”。*

### 2.2 线性过程和马尔可夫过程何时是混合的？

* **线性过程**：对于因果的线性过程 $Z_t = \sum \psi_j \epsilon_{t-j}$，Gorodetskii 证明了在特定条件下其具备 $\alpha$-混合性质。Pham 和 Tran 进一步证明，如果当 $j \to \infty$ 时系数呈指数衰减 $\psi_j = O(r^j)$ (其中 $0 < r < 1$)，则该过程是**几何 $\alpha$-混合 (Geometric $\alpha$-mixing)** 的，即存在常数 $C$ 和 $\rho \in [0, 1)$，使得 $\alpha(k) \le C \rho^k$。
* **马尔可夫过程 / ARCH(p) 模型**：对于马尔可夫形式的过程 $Y_i = m(X_i) + \sigma(X_i)\epsilon_i$，Masry 和 Tjosheim 给出了该过程遍历且满足几何 $\alpha$-混合的具体条件。

---

## 3. 核心不等式：混合序列的协方差界

为了研究相依变量的极限理论（如极限定理中的方差收敛），我们需要对相隔 $k$ 个时间步的随机变量的协方差进行控制。以下四个引理是弱相依数据理论的基石。

### 3.1 Billingsley 不等式 (有界变量的协方差界)

!!! success "引理 4.4 (Billingsley 不等式)"

    假设 $\{Z_i\}$ 是 $\alpha$-混合的（**不要求平稳**），且随机变量 $X \in \mathcal{F}_{-\infty}^t$， $Y \in \mathcal{F}_{t+k}^\infty$ 是一致有界的，即 $|X| \le C_1$ 且 $|Y| \le C_2$。那么它们的协方差有如下界限：
    
    \[
    |Cov(X, Y)| \le 4C_1 C_2 \alpha(k)
    \]

??? proof "Billingsley 不等式的证明（点击展开）"

    考虑协方差的定义，我们可以利用条件期望塔性质展开：

    \[
    |Cov(X, Y)| = |E(XY) - E(X)E(Y)| = |E[X \{ E(Y | \mathcal{F}_{-\infty}^t) - EY \}]|
    \]

    由于 $|X| \le C_1$ 是 $\mathcal{F}_{-\infty}^t$ 可测的，我们有：

    \[
    |Cov(X, Y)| \le C_1 E|E(Y | \mathcal{F}_{-\infty}^t) - EY| = C_1 E[\xi (E(Y | \mathcal{F}_{-\infty}^t) - EY)]
    \]

    这里我们引入符号函数 $\xi = \text{sgn}(E(Y | \mathcal{F}_{-\infty}^t) - EY)$。由于条件期望本身是 $\mathcal{F}_{-\infty}^t$ 可测的，因此指示方向的变量 $\xi$ 也是 $\mathcal{F}_{-\infty}^t$ 可测的。

    将其还原为无条件期望的形式：

    \[
    |Cov(X, Y)| \le C_1 |E(\xi Y) - E\xi EY| = C_1 |Cov(\xi, Y)|  \quad \dots (1)
    \]

    对于 $|Cov(\xi, Y)|$，我们可以对 $Y$ 应用完全相同的技巧。定义 $\eta = \text{sgn}(E(\xi | \mathcal{F}_{t+k}^\infty) - E\xi) \in \mathcal{F}_{t+k}^\infty$，由于 $|Y| \le C_2$，可得：

    \[
    |Cov(\xi, Y)| \le C_2 |E(\xi \eta) - E\xi E\eta| \quad \dots (2)
    \]

    现在我们需要界定 $|E(\xi \eta) - E\xi E\eta|$。注意 $\xi, \eta$ 只能取值 $\{1, -1\}$。
    定义集合 $A = \{\xi = 1\}$，则 $A^c = \{\xi = -1\}$；同理 $B = \{\eta = 1\}$ 和 $B^c = \{\eta = -1\}$。
    显然 $A, A^c \in \mathcal{F}_{-\infty}^t$ 且 $B, B^c \in \mathcal{F}_{t+k}^\infty$。

    展开协方差：

    \[
    \begin{aligned}
    |E(\xi \eta) - E\xi E\eta| &= |P(AB) + P(A^c B^c) - P(A^c B) - P(A B^c) - (P(A) - P(A^c))(P(B) - P(B^c))| \\
    &= |4 (P(A \cap B) - P(A)P(B))| \\
    &\le 4 \alpha(k) \quad \dots (3)
    \end{aligned}
    \]

    *(第二步通过代入 $P(A^c) = 1 - P(A)$ 和 $P(B^c) = 1 - P(B)$ 即可消项化简得到。最后一步直接应用了 $\alpha$-混合的定义。)*

    结合 (1), (2), (3) 式，即证明了 Billingsley 不等式：

    \[
    |Cov(X, Y)| \le 4C_1 C_2 \alpha(k)
    \]

    $\square$

---

### 3.2 截断技术与 $L^p$ 范数界

!!! success "引理 4.5 ($X$ 有高阶矩，$Y$ 有界的情形)"

    假设 $\{Z_i\}$ 是 $\alpha$-混合的（不要求平稳），$X \in \mathcal{F}_{-\infty}^t$ 且存在 $p > 1$ 使得 $E|X|^p < \infty$；同时 $Y \in \mathcal{F}_{t+k}^\infty$ 且有界 $|Y| \le C$。则：
    
    \[
    |Cov(X, Y)| \le 6C \|X\|_p \alpha^{1/q}(k)
    \]

    其中 $q$ 是 $p$ 的共轭数（满足 $\frac{1}{p} + \frac{1}{q} = 1$），且 $\|X\|_p = (E|X|^p)^{1/p}$。

??? proof "引理 4.5 的证明：截断法（点击展开）"

    为了利用引理 4.4（变量有界），我们对无界的 $X$ 引入截断常数 $M > 0$。
    令 $X_M = X \mathbb{I}(|X| \le M)$ 以及尾部 $X'_M = X - X_M = X \mathbb{I}(|X| > M)$。
    由双线性性质和三角不等式：

    \[
    |Cov(X, Y)| = |Cov(X_M, Y) + Cov(X'_M, Y)| \le |Cov(X_M, Y)| + |Cov(X'_M, Y)|
    \]

    对于有界的第一项，直接应用引理 4.4（此时 $X_M$ 受到 $M$ 约束，$Y$ 受到 $C$ 约束）：

    \[
    |Cov(X_M, Y)| \le 4CM \alpha(k)
    \]

    对于尾部的第二项，首先我们来控制其期望。对于 $X'_M$，由于积分区域是 $|X| > M$，因此有比值 $|x|/M > 1$：

    \[
    E|X'_M| = \int_{|x| > M} |x| dF(x) \le \int_{|x| > M} |x| \left(\frac{|x|}{M}\right)^{p-1} dF(x)
    \]

    化简后得到尾部绝对期望的界：

    \[
    E|X'_M| \le M^{-p+1} E|X|^p \quad \dots (4)
    \]

    据此，我们可以控制尾部的协方差：

    \[
    \begin{aligned}
    |Cov(X'_M, Y)| &= |E(X'_M Y) - EX'_M EY| \\
    &\le E|X'_M Y| + C E|X'_M| \\
    &\le 2C E|X'_M| \\
    &\le 2C M^{-p+1} E|X|^p \quad \dots (5)
    \end{aligned}
    \]

    综合两项，得到总协方差的上界为：$4CM\alpha(k) + 2C M^{-p+1} E|X|^p$。
    为了最小化这个上界，我们巧妙地选取截断点 $M$ 为：

    \[
    M = \|X\|_p \{ \alpha(k) \}^{-1/p}
    \]

    代入后即可得到：

    \[
    |Cov(X, Y)| \le 4C \|X\|_p \alpha(k)^{1 - 1/p} + 2C \|X\|_p^p \left(\|X\|_p \alpha(k)^{-1/p}\right)^{1-p}
    \]

    由于 $1 - 1/p = 1/q$，上式化简为：

    \[
    |Cov(X, Y)| \le 4C \|X\|_p \alpha^{1/q}(k) + 2C \|X\|_p \alpha^{1/q}(k) = 6C \|X\|_p \alpha^{1/q}(k)
    \]

    得证。 $\square$

---

### 3.3 Rio 不等式与分位数函数

!!! success "引理 4.6 (Rio's Inequality)"

    设 $X$ 和 $Y$ 为两个可积的实值随机变量，满足极限 $lim_{c \to \infty} E\{|X| \mathbb{I}(|X| > c)\} = 0$。定义 $Q_X(u) = \inf\{t: P(|X| > t) \le u\}$ 为 $|X|$ 的**上分位数函数 (Upper quantile function)**。
    
    如果 $Q_X Q_Y$ 在区间 $(0, 1)$ 上是可积的，则有：
    
    \[
    |Cov(X, Y)| \le 2 \int_0^{2\alpha} Q_X(u) Q_Y(u) du
    \]
    
    其中 $\alpha = \alpha(\sigma(X), \sigma(Y))$ 是由 $X$ 和 $Y$ 生成的 $\sigma$-代数之间的 $\alpha$-混合系数。

*(该不等式的证明涉及非参数统计和随机过程深度知识，详见 Bosq, D. (1998)。它是推导后续更强的不等式的关键基石。)*

---

### 3.4 Davydov 不等式 (一般化的 $L^q-L^r$ 协方差界)

Davydov 不等式是最常用的协方差不等式，它解除了变量必须有界的限制，仅要求它们具有适当的高阶矩。

!!! success "引理 4.7 (Davydov Inequality)"

    设 $X$ 和 $Y$ 是两个实值随机变量，满足 $X \in L^q(\mathcal{F}_{-\infty}^t)$， $Y \in L^r(\mathcal{F}_{t+k}^\infty)$。其中 $q > 1$, $r > 1$，且存在 $p$ 满足 Hölder 形式的关系：
    
    \[
    \frac{1}{q} + \frac{1}{r} = 1 - \frac{1}{p} \quad (\text{其中 } p > 1)
    \]
    
    那么有：
    
    \[
    |Cov(X, Y)| \le 2p (2\alpha(k))^{1/p} \|X\|_q \|Y\|_r
    \]

??? proof "Davydov 不等式的证明（点击展开）"

    证明的核心是利用前面的 Rio 不等式（引理 4.6），按参数的有穷性分三种情况讨论。

    **(i) 假设 $q$ 和 $r$ 都是有限的：**
    根据马尔可夫不等式 (Markov's Inequality)，对于任意 $u \in (0, 1]$，我们有：
    
    \[
    P\left( |X| > \frac{\|X\|_q}{u^{1/q}} \right) \le \frac{E|X|^q}{(\|X\|_q / u^{1/q})^q} = u
    \]
    
    根据上分位数函数的定义，这直接推导出对于所有的 $0 < u \le 1$：
    
    \[
    Q_X(u) \le \frac{\|X\|_q}{u^{1/q}}
    \]
    
    对称地，对于 $Y$ 也有 $Q_Y(u) \le \frac{\|Y\|_r}{u^{1/r}}$。
    将这两项代入 Rio 不等式的积分中：
    
    \[
    \begin{aligned}
    |Cov(X,Y)| &\le 2 \int_0^{2\alpha(k)} \frac{\|X\|_q}{u^{1/q}} \frac{\|Y\|_r}{u^{1/r}} du \\
    &= 2 \|X\|_q \|Y\|_r \int_0^{2\alpha(k)} u^{-\frac{1}{q} - \frac{1}{r}} du
    \end{aligned}
    \]
    
    由于 $\frac{1}{q} + \frac{1}{r} = 1 - \frac{1}{p}$，积分内部的指数部分为 $u^{\frac{1}{p} - 1}$。对该幂函数进行积分：
    
    \[
    \int_0^{2\alpha(k)} u^{\frac{1}{p} - 1} du = \left[ p \cdot u^{1/p} \right]_0^{2\alpha(k)} = p (2\alpha(k))^{1/p}
    \]
    
    代回原式即得：
    
    \[
    |Cov(X,Y)| \le 2p (2\alpha(k))^{1/p} \|X\|_q \|Y\|_r
    \]

    **(ii) 假设 $r = +\infty$，且 $q$ 是有限的：**
    此时依据参数关系，$\frac{1}{q} + \frac{1}{p} = 1$。
    对于 $L^\infty$ 空间中的 $Y$，其上分位数函数存在一个天然的硬上界 $Q_Y(u) \le Q_Y(0) = \|Y\|_\infty = \sup |Y|$。
    *(因为使得 $P(|Y| > t) = 0$ 的最小 $t$ 恰好就是 $\|Y\|_\infty$)*。
    
    再次应用 Rio 不等式：
    
    \[
    \begin{aligned}
    |Cov(X,Y)| &\le 2 \int_0^{2\alpha(k)} \frac{\|X\|_q}{u^{1/q}} \|Y\|_\infty du \\
    &= 2 \|X\|_q \|Y\|_\infty \int_0^{2\alpha(k)} u^{\frac{1}{p} - 1} du \\
    &= 2p (2\alpha(k))^{1/p} \|X\|_q \|Y\|_\infty
    \end{aligned}
    \]

    *(注：这与前面推导的引理 4.5 结果十分相似，但常数因子略有不同。)*

    **(iii) 假设 $r = +\infty$ 且 $q = +\infty$：**
    此时两个变量都是本质有界的，代入参数关系必然有 $1/p = 1 \implies p = 1$。
    此时我们对两个分位数函数都采用硬界 $Q_X(u) \le \|X\|_\infty$ 且 $Q_Y(u) \le \|Y\|_\infty$。
    代入 Rio 不等式：
    
    \[
    |Cov(X,Y)| \le 2 \int_0^{2\alpha(k)} \|X\|_\infty \|Y\|_\infty du = 4 \|X\|_\infty \|Y\|_\infty \alpha(k)
    \]
    
    这恰好完全退化为 Billingsley 不等式（引理 4.4）的形式。
    
    综上所有情况，引理得证。 $\square$