# 第六章：解的存在唯一性定理

在常微分方程 (ODE) 的理论中，我们通过 Picard 迭代法结合 Lipschitz 条件来证明初值问题 $y' = f(x,y), y(x_0)=y_0$ 解的存在唯一性。

对于随机微分方程 (SDE)，我们面临着包含布朗运动积分的复杂情形：

$$
dX_t = b(X_t, t)dt + \sigma(X_t, t)dW_t, \quad X_0 = x_0
$$

为了保证该方程有唯一强解，我们需要对漂移系数 $b$ 和扩散系数 $\sigma$ 施加类似的限制。本节将详细梳理这一核心定理及其严谨证明。

---

## 1. 存在唯一性定理 (Existence and Uniqueness Theorem)

!!! abstract "定理：SDE 强解的存在唯一性"

    设 $T > 0$，函数 $b(x,t): \mathbb{R}^n \times [0,T] \to \mathbb{R}^n$ 和 $\sigma(x,t): \mathbb{R}^n \times [0,T] \to \mathbb{R}^{n \times m}$ 满足以下两个基本条件：

    **1. 全局利普希茨条件 (Global Lipschitz Condition)**：
    存在常数 $L > 0$，使得对于任意的 $x, y \in \mathbb{R}^n$ 和 $t \in [0,T]$，有：

    $$
    |b(x,t) - b(y,t)| + |\sigma(x,t) - \sigma(y,t)| \le L|x - y|
    $$

    **2. 线性增长条件 (Linear Growth Condition)**：
    存在常数 $L > 0$，使得对于任意的 $x \in \mathbb{R}^n$ 和 $t \in [0,T]$，有：

    $$
    |b(x,t)|^2 + |\sigma(x,t)|^2 \le L^2(1 + |x|^2)
    $$

    **初始条件**：
    设 $\xi$ 是一个与布朗运动 $W_t$ 独立的随机变量，且满足二阶矩有限 $E|\xi|^2 < \infty$。

    **结论**：
    那么，初值问题 $X_t = \xi + \int_0^t b(X_s, s)ds + \int_0^t \sigma(X_s, s)dW_s$ 在区间 $[0,T]$ 上**存在唯一的平方可积连续强解** $X_t$，且满足：

    $$
    E\left[ \sup_{0 \le t \le T} |X_t|^2 \right] < \infty
    $$

---

## 2. 定理的证明 (Proof of the Theorem)

证明分为两大部分：**唯一性**（利用 Gronwall 不等式）和**存在性**（利用 Picard 迭代和 Borel-Cantelli 引理）。

### 2.1 唯一性的证明

??? proof "唯一性证明：Gronwall 不等式的应用（点击展开）"

    假设存在两个解 $X_t$ 和 $\tilde{X}_t$ 满足相同的初始条件 $X_0 = \tilde{X}_0 = \xi$。

    令误差过程为 $Y_t = X_t - \tilde{X}_t$，则 $Y_0 = 0$。根据方程的积分形式，我们有：

    $$
    Y_t = \int_0^t [b(X_s, s) - b(\tilde{X}_s, s)]ds + \int_0^t [\sigma(X_s, s) - \sigma(\tilde{X}_s, s)]dW_s
    $$

    利用基本不等式 $(a+b)^2 \le 2a^2 + 2b^2$，两边取绝对值平方并求期望：

    $$
    E|Y_t|^2 \le 2E \left| \int_0^t [b(X_s, s) - b(\tilde{X}_s, s)]ds \right|^2 + 2E \left| \int_0^t [\sigma(X_s, s) - \sigma(\tilde{X}_s, s)]dW_s \right|^2
    $$

    对于第一项（黎曼积分），利用 Cauchy-Schwarz 不等式：

    $$
    E \left| \int_0^t [b(X_s, s) - b(\tilde{X}_s, s)]ds \right|^2 \le t E \int_0^t |b(X_s, s) - b(\tilde{X}_s, s)|^2 ds
    $$

    对于第二项（Itô 积分），利用 Itô 等距同构 (Itô Isometry)：

    $$
    E \left| \int_0^t [\sigma(X_s, s) - \sigma(\tilde{X}_s, s)]dW_s \right|^2 = E \int_0^t |\sigma(X_s, s) - \sigma(\tilde{X}_s, s)|^2 ds
    $$

    将 Lipschitz 条件 $|b(X) - b(\tilde{X})|^2 \le L^2 |X - \tilde{X}|^2 = L^2 |Y_s|^2$ 以及 $|\sigma(X) - \sigma(\tilde{X})|^2 \le L^2 |Y_s|^2$ 代入上式：

    $$
    E|Y_t|^2 \le 2t L^2 \int_0^t E|Y_s|^2 ds + 2L^2 \int_0^t E|Y_s|^2 ds \le 2L^2(T+1) \int_0^t E|Y_s|^2 ds
    $$

    令常数 $C = 2L^2(T+1)$，我们得到：

    $$
    E|Y_t|^2 \le C \int_0^t E|Y_s|^2 ds
    $$

    由于 $E|Y_t|^2$ 是非负的，**由 Gronwall 不等式直接可得**，对于所有的 $t \in [0,T]$：

    $$
    E|Y_t|^2 = 0
    $$

    这意味着对于任意给定的 $t$，几乎必然有 $X_t = \tilde{X}_t$。结合样本路径的连续性，得出这两个解在整个区间 $[0,T]$ 上是不可区分的 (indistinguishable)。唯一性得证。 $\square$

---

### 2.2 存在性的证明 (Picard 迭代法)

存在性的核心思想是构造逼近序列。这一部分的推导逻辑链条较长，需要步步为营。

??? proof "存在性证明：Picard 迭代法（点击展开）"

    **步骤 1：构造 Picard 迭代序列**

    定义初始近似为 $X^{(0)}_t = \xi$。
    
    对于 $k \ge 0$，递归地定义：

    $$
    X^{(k+1)}_t = \xi + \int_0^t b(X^{(k)}_s, s)ds + \int_0^t \sigma(X^{(k)}_s, s)dW_s
    $$

    **步骤 2：证明迭代序列的均方有界性**

    我们需要确保序列中每一个元素的二阶矩 $E|X^{(k)}_t|^2$ 是有限的。这可以通过数学归纳法和线性增长条件来证明：

    $$
    E|X^{(k)}_t|^2 \le C(1 + E|\xi|^2) e^{Ct}
    $$
    
    （此处的推导与唯一性证明中的 Cauchy-Schwarz 和 Itô 等距类似，利用了线性增长条件，笔记中缩写为 $E|X^k_t|^2 \le Ce^{ct}$，此处略去重复步骤）。

    **步骤 3：估计相邻迭代项差异的衰减速度**

    定义相邻两项的均方误差为 $d_k(t) = E|X^{(k+1)}_t - X^{(k)}_t|^2$。
    
    利用与唯一性证明完全相同的放缩技巧（Cauchy-Schwarz + Itô Isometry + Lipschitz 条件），可以得到递推关系：

    $$
    d_k(t) \le 2L^2(T+1) \int_0^t d_{k-1}(s) ds = M \int_0^t d_{k-1}(s) ds
    $$

    其中 $M = 2L^2(T+1)$。

    对于 $k=0$ 时：
    
    $$
    d_0(t) = E|X^{(1)}_t - X^{(0)}_t|^2 \le M t (1 + E|\xi|^2) = C_0 t
    $$

    不断向后递推积分：

    $$
    d_k(t) \le M \int_0^t d_{k-1}(s) ds \le M^k C_0 \int_0^t \frac{s^{k-1}}{(k-1)!} ds = C_0 \frac{M^k t^k}{k!}
    $$

    这说明误差的期望以极快的速度（阶乘级别）衰减。

    **步骤 4：利用 Doob 极大值不等式控制路径上确界**

    为了证明一致收敛，我们需要考察路径在整个区间上的最大误差。利用 Doob 鞅不等式 (Doob's Martingale Inequality) 和 Burkholder-Davis-Gundy (BDG) 不等式处理 Itô 积分项的最大值：

    $$
    E\left[ \sup_{0 \le s \le T} |X^{(k+1)}_s - X^{(k)}_s|^2 \right] \le C_1 \int_0^T E|X^{(k)}_s - X^{(k-1)}_s|^2 ds \le C_2 \frac{(MT)^k}{k!}
    $$

    根据马尔可夫不等式 (Markov's Inequality)：

    $$
    P\left( \sup_{0 \le t \le T} |X^{(k+1)}_t - X^{(k)}_t| > 2^{-k} \right) \le \frac{E[\sup_{0 \le t \le T} |X^{(k+1)}_t - X^{(k)}_t|^2]}{(2^{-k})^2} \le 2^{2k} C_2 \frac{(MT)^k}{k!}
    $$

    **步骤 5：Borel-Cantelli 引理与一致收敛**

    由于阶乘的增长速度远快于指数，级数是收敛的：

    $$
    \sum_{k=0}^{\infty} 2^{2k} C_2 \frac{(MT)^k}{k!} < \infty
    $$

    由 **Borel-Cantelli 引理** 可知，事件 $\{\sup_{0 \le t \le T} |X^{(k+1)}_t - X^{(k)}_t| > 2^{-k}\}$ 发生无限次的概率为 0。

    这意味着，对于几乎所有的样本路径 $\omega$，存在一个 $K(\omega)$，当 $k > K(\omega)$ 时，$\sup_{0 \le t \le T} |X^{(k+1)}_t - X^{(k)}_t| \le 2^{-k}$。
    
    由于 $\sum 2^{-k}$ 收敛，这保证了级数 $\sum (X^{(k+1)}_t - X^{(k)}_t)$ 绝对且一致收敛。因此，极限过程存在：

    $$
    X_t = X^{(0)}_t + \sum_{k=0}^{\infty} (X^{(k+1)}_t - X^{(k)}_t) = \lim_{k \to \infty} X^{(k)}_t
    $$

    且该收敛在 $t \in [0,T]$ 上是一致的，保证了极限过程 $X_t$ 是连续的。

    **步骤 6：相容性 (Consistency) - 验证极限是解**

    最后，需要验证这个几乎必然一致收敛的极限 $X_t$ 确实满足原方程。
    
    在 Picard 迭代公式中对两边取极限 $k \to \infty$。由于 $b$ 和 $\sigma$ 满足 Lipschitz 连续，结合控制收敛定理 (Dominated Convergence Theorem) 和 Itô 积分的等距同构性质，极限可以穿透积分号：

    $$
    \int_0^t \sigma(X^{(k)}_s) dW_s \xrightarrow{L^2} \int_0^t \sigma(X_s) dW_s
    $$

    由此验证了 $X_t = \xi + \int_0^t b(X_s, s)ds + \int_0^t \sigma(X_s, s)dW_s$。存在性得证。 $\square$

---

## 3. 解的高阶矩与初值依赖性 (Properties of the Solution)

在笔记的最后部分，提到了解的高阶矩界限以及对初始值的连续依赖性。这是 SDE 理论中非常重要的推论。

!!! info "性质：高阶矩与依赖性"

    **1. 高阶矩的界限**
    
    如果初始条件 $\xi$ 满足 $E|\xi|^{2p} < \infty$ ($p \ge 1$)，那么在全局 Lipschitz 和线性增长条件下，强解 $X_t$ 满足：

    $$
    E\left[ \sup_{0 \le t \le T} |X_t|^{2p} \right] \le C(1 + E|\xi|^{2p}) e^{CT}
    $$
    
    *(证明同样利用 Itô 公式、BDG 不等式和 Gronwall 不等式)。*

    **2. 对初值的连续依赖性**
    
    设 $X^x_t$ 和 $X^y_t$ 分别是初值为 $X_0 = x$ 和 $X_0 = y$ 的强解，那么它们之间的差异受初值差异的控制：

    $$
    E\left[ \sup_{0 \le t \le T} |X^x_t - X^y_t|^2 \right] \le C |x - y|^2 e^{CT}
    $$
    
    这说明 SDE 的解对其初始状态是 Lipschitz 连续依赖的，这在研究马尔可夫半群和偏微分方程 (PDEs) 的联系时具有关键作用。