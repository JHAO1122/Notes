---
tags:
  - 高等概率论
  - 讲座讲义
  - 大偏差原理
  - LDP
---

# 🎙️ 高等概率论讨论班：从大数定律到大偏差原理

> **"The Strong Law of Large Numbers tells us where we are going; the Large Deviation Principle tells us the cost of getting lost."**

!!! abstract "讲座综述 (Talk Abstract)"
    * **主题**：随机序列的极限行为：从 SLLN 到 Cramér 大偏差定理
    * **时长**：90 分钟 (Full Session)
    * **逻辑架构**：
        1. **收敛基础**：回顾尾代数、0-1 律及级数收敛判定（Kolmogorov 三级数定理）。
        2. **收敛速率**：利用 Kronecker 引理探讨 SLLN 的收敛阶（Marcinkiewicz-Zygmund 定理）。
        3. **大偏差 (LDP) 核心**：研究偏离均值的指数级小概率事件，推导 Cramér 定理与率函数 $I(a)$。
    * **核心难点**：勒让德变换的直观理解、测度变换（Cramér Transform）在下界证明中的应用。

---

## 1. 尾 $\sigma$-代数与 0-1 律

在研究随机变量序列的渐近行为时，我们往往只关心那些“不受有限个随机变量影响”的事件。这类事件的集合构成了一个极其重要的代数结构——尾 $\sigma$-代数。

!!! info "定义 1.1 (尾 $\sigma$-代数 Tail $\sigma$-field)"

    设 $\{X_n\}_{n \ge 1}$ 是一列随机变量。定义 $\mathcal{F}_n' = \sigma(X_n, X_{n+1}, \dots)$ 为从第 $n$ 个变量开始生成的 $\sigma$-代数。
    定义**尾 $\sigma$-代数** $\mathcal{T}$ 为：

    \[
    \mathcal{T} = \bigcap_{n=1}^\infty \mathcal{F}_n'
    \]

    如果一个事件 $A \in \mathcal{T}$，那么对于任意的 $n$，事件 $A$ 的发生与否完全由 $\{X_n, X_{n+1}, \dots\}$ 决定，而与前 $n-1$ 个变量无关。我们称 $\mathcal{T}$ 中的事件为**尾事件**。

**经典推论：Borel-Cantelli 引理 I**

如果 $\sum_{n=1}^\infty P(A_n) < \infty$，那么 $P(A_n \text{ i.o.}) = 0$。（其中 $\text{i.o.}$ 表示 infinitely often，即发生无限次）。
极限上确界 $\limsup A_n = \{A_n \text{ i.o.}\}$ 也是一个典型的尾事件。

### 1.1 Kolmogorov 0-1 律

!!! success "定理 1.1 (Kolmogorov's 0-1 Law)"

    如果序列 $X_1, X_2, \dots$ 是**相互独立**的随机变量，且 $A \in \mathcal{T}$ 是一个尾事件，那么：

    \[
    P(A) = 0 \quad \text{或} \quad P(A) = 1
    \]

??? proof "Kolmogorov 0-1 律的证明（点击展开）"

    核心思想：证明事件 $A$ 与其自身是独立的。

    设 $\mathcal{F}_n = \sigma(X_1, \dots, X_n)$，$\mathcal{F}_{n+1}' = \sigma(X_{n+1}, X_{n+2}, \dots)$。
    因为序列 $\{X_n\}$ 是相互独立的，所以由不相交的变量集合生成的 $\sigma$-代数 $\mathcal{F}_n$ 和 $\mathcal{F}_{n+1}'$ 也是相互独立的。

    对于任意的尾事件 $A \in \mathcal{T}$，由于 $\mathcal{T} \subset \mathcal{F}_{n+1}'$，事件 $A$ 必定属于 $\mathcal{F}_{n+1}'$。
    因此，$A$ 独立于 $\mathcal{F}_n$。由于这一结论对所有的 $n$ 都成立，我们得出 $A$ 独立于所有前有限个变量生成的代数的并集 $\bigcup_{n=1}^\infty \mathcal{F}_n$。

    根据测度论中的 $\pi-\lambda$ 定理（或单调类定理），既然 $A$ 独立于并集代数，那么 $A$ 必然独立于其生成的 $\sigma$-代数：

    \[
    \mathcal{F}_\infty = \sigma(X_1, X_2, \dots)
    \]

    然而，尾 $\sigma$-代数 $\mathcal{T}$ 本身就是 $\mathcal{F}_\infty$ 的一个子 $\sigma$-代数，所以 $A \in \mathcal{F}_\infty$。
    综上所述，$A$ 与 $\mathcal{F}_\infty$ 独立，同时也属于 $\mathcal{F}_\infty$。这意味着 $A$ 必须与**自身**独立：

    \[
    P(A \cap A) = P(A)P(A) \implies P(A) = [P(A)]^2
    \]

    解此方程，只能得到 $P(A) = 0$ 或 $P(A) = 1$。 $\square$

### 1.2 可置换事件与 Hewitt-Savage 0-1 律

除了舍弃前有限项的尾事件，还有一类事件对有限个元素的排列顺序不敏感。

!!! info "定义 1.2 (有限置换与可置换事件)"

    * **有限置换 (Finite Permutation)**：一个映射 $\pi: \mathbb{N} \rightarrow \mathbb{N}$，如果它是一一对应的，且只有有限个 $i$ 满足 $\pi(i) \ne i$，则称 $\pi$ 为有限置换。
    * **可置换事件 (Permutable Event)**：如果对于任意的有限置换 $\pi$，事件 $A$ 的原像 $\pi^{-1}(A) := \{\omega : \pi(\omega) \in A\}$ 总是等于 $A$，则称 $A$ 为可置换事件。
    
    所有可置换事件构成了**可置换 $\sigma$-代数 (Exchangeable $\sigma$-field)**，记为 $\mathcal{E}$。显然有 $\mathcal{T} \subset \mathcal{E}$。

!!! success "定理 1.2 (Hewitt-Savage 0-1 Law)"

    如果 $X_1, X_2, \dots$ 是**独立同分布 (i.i.d.)** 的，且 $A \in \mathcal{E}$ 是一个可置换事件，那么：

    \[
    P(A) \in \{0, 1\}
    \]

??? proof "Hewitt-Savage 0-1 律的证明（点击展开）"

    基本思路同 Kolmogorov 0-1 律：证明 $P(A) = [P(A)]^2$。

    对于可置换事件 $A \in \mathcal{E} \subset \sigma(X_1, X_2, \dots)$，根据测度逼近定理，对于任给的 $\epsilon > 0$，必定存在一个依赖于前 $n$ 个变量的“柱集”事件 $A_n \in \sigma(X_1, \dots, X_n)$，使得对称差的概率极小：

    \[
    P(A \Delta A_n) \rightarrow 0 \quad (n \rightarrow \infty)
    \]

    这同样意味着 $P(A_n) \rightarrow P(A)$。

    现在，构造一个特定的有限置换 $\pi_n$，它将前 $n$ 个坐标与接下来的 $n$ 个坐标对调：
    $\pi_n(1, \dots, n, n+1, \dots, 2n) = (n+1, \dots, 2n, 1, \dots, n)$。

    设 $A_n' = \pi_n(A_n)$。由于 $A_n$ 只依赖于 $X_1, \dots, X_n$，那么 $A_n'$ 就只依赖于 $X_{n+1}, \dots, X_{2n}$。
    因为序列是 i.i.d. 的，所以 $A_n$ 和 $A_n'$ 是**相互独立且同分布**的。因此：

    \[
    P(A_n \cap A_n') = P(A_n) P(A_n') \rightarrow P(A) \cdot P(A) = P(A)^2
    \]

    另一方面，由于 $A$ 是可置换事件，$\pi_n(A) = A$。因此对于对称差运算，置换不改变其概率：

    \[
    P(A_n' \Delta A) = P(\pi_n(A_n) \Delta \pi_n(A)) = P(\pi_n(A_n \Delta A)) = P(A_n \Delta A) \rightarrow 0
    \]

    既然 $A_n$ 和 $A_n'$ 在概率上都逼近于同一个事件 $A$，那么它们的交集 $A_n \cap A_n'$ 在概率上也必然逼近于 $A$ 本身（即 $P(A_n \cap A_n') \rightarrow P(A)$）。

    结合上下两式，得到：

    \[
    P(A) = P(A)^2 \implies P(A) \in \{0, 1\}
    \]
    
    证明完毕。 $\square$

---

## 2. 随机级数的收敛定理

为了研究 $\sum X_n$ 的收敛性，我们需要一个强大的不等式工具来控制局部波动的最大值。

!!! success "引理 2.1 (Kolmogorov 极大值不等式)"

    假设 $X_1, \dots, X_n$ 相互独立，均值为 0，且方差存在。记部分和 $S_k = \sum_{i=1}^k X_i$。对于任意 $x > 0$：

    \[
    P\left( \max_{1 \le k \le n} |S_k| \ge x \right) \le \frac{Var(S_n)}{x^2}
    \]

    *(注：对比 Chebyshev 不等式 $P(|S_n| \ge x) \le x^{-2} Var(S_n)$，极大值不等式给出了更强的一致界。)*

### 2.1 Kolmogorov 级数收敛定理

!!! success "定理 2.1 (Kolmogorov's Convergence Theorem)"

    假设 $\{X_n\}$ 是相互独立的随机变量序列，且 $E(X_n) = 0$。如果方差级数收敛：

    \[
    \sum_{n=1}^\infty Var(X_n) < \infty
    \]

    那么随机级数 $\sum_{n=1}^\infty X_n$ **几乎处处 (a.s.)** 收敛。

??? proof "Kolmogorov 收敛定理的证明（点击展开）"

    令部分和 $S_N = \sum_{n=1}^N X_n$。我们要证明序列 $\{S_N\}$ 在 $\mathbb{R}$ 中是一个 Cauchy 列 a.s.。

    应用 Kolmogorov 极大值不等式考察区间 $(M, N]$ 上的波动：

    \[
    P\left( \max_{M < m \le N} |S_m - S_M| > \epsilon \right) \le \frac{1}{\epsilon^2} Var(S_N - S_M) = \frac{1}{\epsilon^2} \sum_{n=M+1}^N Var(X_n)
    \]

    由于 $\sum Var(X_n) < \infty$，当 $M, N \rightarrow \infty$ 时，级数余项趋于 0。令 $N \rightarrow \infty$，由连续性：

    \[
    P\left( \sup_{m > M} |S_m - S_M| > \epsilon \right) \le \frac{1}{\epsilon^2} \sum_{n=M+1}^\infty Var(X_n) \xrightarrow{M \rightarrow \infty} 0
    \]

    这意味着对于任意 $\epsilon > 0$，尾部最大波动的概率趋于 0。这等价于 $\{S_N\}$ 是一个 Cauchy 列的概率为 1，故级数几乎处处收敛。 $\square$

### 2.2 Kolmogorov 三级数定理

并非所有随机变量都具有方差或期望，此时我们需要运用**截断方法 (Truncation Method)**。

!!! success "定理 2.2 (Kolmogorov's Three-Series Theorem)"

    设 $X_1, X_2, \dots$ 是相互独立的随机变量，任取常数 $A > 0$。定义截断变量：
    
    \[
    Y_n = X_n \mathbb{I}_{\{|X_n| \le A\}}
    \]
    
    则随机级数 $\sum_{n=1}^\infty X_n$ 几乎处处收敛的**充分必要条件**是以下三个级数同时收敛：

    (i) $\sum_{n=1}^\infty P(|X_n| > A) < \infty$

    (ii) $\sum_{n=1}^\infty E(Y_n)$ 收敛

    (iii) $\sum_{n=1}^\infty Var(Y_n) < \infty$

??? proof "三级数定理（充分性）的证明（点击展开）"

    已知条件 (i), (ii), (iii) 收敛，我们要证明原级数收敛。

    根据条件 (iii)，截断变量的方差级数收敛：$\sum Var(Y_n) < \infty$。
    显然 $Y_n - E(Y_n)$ 是相互独立且均值为 0 的随机变量。
    根据 **Kolmogorov 级数收敛定理 (定理 2.1)**，级数 $\sum_{n=1}^\infty (Y_n - E(Y_n))$ 几乎处处收敛。

    再结合条件 (ii)，确定性级数 $\sum_{n=1}^\infty E(Y_n)$ 是收敛的。
    将这两者相加，即可得出：
    
    \[
    \sum_{n=1}^\infty Y_n \quad \text{几乎处处收敛。}
    \]

    接下来，我们要把 $Y_n$ 的收敛性过渡回 $X_n$。
    根据条件 (i)，$\sum_{n=1}^\infty P(|X_n| > A) < \infty$。
    由于 $P(|X_n| > A) = P(X_n \ne Y_n)$，根据 **Borel-Cantelli 引理 I**：

    \[
    P(X_n \ne Y_n \text{ i.o.}) = 0
    \]

    这意味着，在概率为 1 的样本空间中，至多只有有限个 $n$ 满足 $X_n \ne Y_n$。换句话说，对于充分大的 $n$，$X_n$ 与 $Y_n$ 最终是完全一致的。
    因此，级数 $\sum X_n$ 与 $\sum Y_n$ 的敛散性必定完全相同。
    既然 $\sum Y_n$ a.s. 收敛，那么 $\sum X_n$ 也必然 a.s. 收敛。充分性证毕。 $\square$

---

## 3. 强大数定律 (Strong Law of Large Numbers)

在讨论大数定律前，我们需要一个分析学中的基础引理，它能将“级数收敛”转化为“加权平均收敛”。

!!! info "引理 3.1 (Kronecker's Lemma)"

    假设 $\{a_n\}$ 是严格递增且趋于无穷的正数列，即 $0 < a_1 \le a_2 \dots \uparrow \infty$。
    如果实数级数 $\sum_{n=1}^\infty \frac{x_n}{a_n}$ 收敛，那么必然有：

    \[
    \lim_{n \rightarrow \infty} \frac{1}{a_n} \sum_{m=1}^n x_m = 0
    \]

### 3.1 基于方差的强大数定律

!!! success "定理 3.1 (SLLN under finite condition)"

    设 $X_1, X_2, \dots$ 是 **i.i.d.** 随机变量，且 $E(X_1) = \mu, E(X_1^2) < \infty$。
    记部分和 $S_n = \sum_{i=1}^n X_i$，那么当 $n \rightarrow \infty$ 时：

    \[
    \frac{S_n}{n} \xrightarrow{a.s.} \mu
    \]

??? proof "强大数定律的截断证明法（点击展开）"

    **Step 1: 构造截断变量并验证方差级数**
    定义随着 $k$ 动态截断的变量：

    \[
    Y_k = X_k \mathbb{I}_{\{|X_k| \le k\}}
    \]
    
    定义中心化变量 $Z_k = Y_k - E(Y_k)$。显然 $E(Z_k) = 0$ 且相互独立。
    考察其方差级数：

    \[
    \sum_{k=1}^\infty \frac{Var(Z_k)}{k^2} \le \sum_{k=1}^\infty \frac{E(Y_k^2)}{k^2} = \sum_{k=1}^\infty \frac{1}{k^2} E[X_k^2 \mathbb{I}_{\{|X_k| \le k\}}]
    \]

    由于 $\{X_k\}$ 同分布，上式等于：

    \[
    \sum_{k=1}^\infty \frac{1}{k^2} E[X_1^2 \mathbb{I}_{\{|X_1| \le k\}}] < \infty \quad \text{（基于 } E(X_1^2) < \infty \text{ 的假设）}
    \]

    **Step 2: 应用级数收敛与 Kronecker 引理**
    根据 Kolmogorov 级数收敛定理，级数 $\sum_{k=1}^\infty \frac{Z_k}{k}$ a.s. 收敛。
    应用 **Kronecker 引理**（取 $a_n = n$）：

    \[
    \frac{1}{n} \sum_{k=1}^n Z_k = \frac{1}{n} \sum_{k=1}^n (Y_k - E(Y_k)) \xrightarrow{a.s.} 0
    \]

    **Step 3: 期望渐近性与 B-C 引理过渡**
    由于由控制收敛定理 (DCT)，$E(Y_k) \rightarrow E(X_1) = \mu$。由 Césaro 均值性质，$\frac{1}{n} \sum_{k=1}^n E(Y_k) \rightarrow \mu$。
    因此：

    \[
    \frac{1}{n} \sum_{k=1}^n Y_k \xrightarrow{a.s.} \mu
    \]

    最后，由 Chebyshev 不等式和 B-C 引理：
    
    \[
    \sum_{k=1}^\infty P(X_k \ne Y_k) = \sum_{k=1}^\infty P(|X_1| > k) \le E(|X_1|) \le (E(X_1^2))^{1/2} < \infty
    \]

    故 $P(X_k \ne Y_k \text{ i.o.}) = 0$。这意味着在极限意义下，用 $Y_k$ 替代 $X_k$ 不改变平均值的收敛性，因此 $\frac{S_n}{n} \xrightarrow{a.s.} \mu$。 $\square$

*(注：Khinchin 弱化了方差存在的前提，证明了仅需一阶矩 $E|X_1| < \infty$ 即可保证强大数定律成立，这涉及更精细的积分放缩控制。)*

---

## 4. 广义收敛速率与无穷期望

当我们研究超出 $\frac{1}{n}$ 缩放比例的收敛时，往往需要调节矩存在的阶数条件。

### 4.1 方差存在时的收敛速率

!!! success "定理 4.1"

    设 $X_1, X_2, \dots$ 为 i.i.d. 且 $E(X_i) = 0, E(X_i^2) = \sigma^2 < \infty$。
    对于任意给定的 $\delta > 0$，有：

    \[
    \frac{S_n}{n^{1/2}(\log n)^{1/2 + \delta}} \xrightarrow{a.s.} 0
    \]

??? proof "证明简述（点击展开）"

    令归一化系数 $a_n = n^{1/2}(\log n)^{1/2 + \delta}$。
    考察级数 $\sum_{n=1}^\infty \frac{Var(X_n)}{a_n^2} = \sigma^2 \sum_{n=1}^\infty \frac{1}{n (\log n)^{1+2\delta}}$。
    由于积分 $\int \frac{1}{x (\log x)^{1+p}} dx < \infty$ ($p>0$)，该级数收敛。
    根据定理 2.1，$\sum_{n=1}^\infty \frac{X_n}{a_n}$ 几乎处处收敛。再应用 Kronecker 引理，直接得到结论。 $\square$

### 4.2 Marcinkiewicz-Zygmund 强大数定律

如果随机变量的二阶矩不存在，但存在介于一阶和二阶之间的 $p$ 阶矩，会有什么样的收敛速率？

!!! success "定理 4.2"

    设 $X_1, X_2, \dots$ 为 i.i.d. 且 $E(X_i) = 0, E|X_i|^p < \infty$ 其中 $1 < p < 2$。那么：

    \[
    \frac{S_n}{n^{1/p}} \xrightarrow{a.s.} 0
    \]

??? proof "证明：精细截断法（点击展开）"

    定义截断变量 $Y_m = X_m \mathbb{I}_{\{|X_m| \le m^{1/p}\}}$。
    
    **(1) 过渡的合法性：**
    
    \[
    \sum_{m=1}^\infty P(X_m \ne Y_m) = \sum_{m=1}^\infty P(|X_1| > m^{1/p}) = \sum_{m=1}^\infty P(|X_1|^p > m) \le E(|X_1|^p) < \infty
    \]
    
    由 B-C 引理 I，我们只需要证明截断后的序列 $\sum Y_m / m^{1/p} \rightarrow 0$ 即可。

    **(2) 截断方差级数的收敛：**
    中心化 $Z_m = Y_m - E(Y_m)$。我们要证明 $\sum_{m=1}^\infty \frac{E(Y_m^2)}{m^{2/p}} < \infty$。

    \[
    \sum_{m=1}^\infty \frac{E(Y_m^2)}{m^{2/p}} = \sum_{m=1}^\infty m^{-2/p} \int_0^{m^{1/p}} y^2 dP(|X_1| \le y)
    \]
    
    交换求和与积分顺序（Fubini 定理）：
    
    \[
    = \int_0^\infty y^2 \left( \sum_{m \ge y^p} m^{-2/p} \right) dP
    \]
    
    注意尾和估计 $\sum_{m \ge k} m^{-2/p} \le C \cdot k^{1 - 2/p}$。代入得：
    
    \[
    \le C \int y^2 (y^p)^{1 - 2/p} dP = C \int y^2 y^{p-2} dP = C E(|X_1|^p) < \infty
    \]

    **(3) 运用 Kronecker 引理：**
    由于 $\sum Var(Z_m) / (m^{1/p})^2 < \infty$，级数 $\sum \frac{Z_m}{m^{1/p}}$ 收敛 a.s.。
    由 Kronecker 引理，$\frac{1}{n^{1/p}} \sum_{m=1}^n Z_m \rightarrow 0$ a.s.

    **(4) 期望漂移的处理：**
    由于 $E(X_m) = 0$，有 $E(Y_m) = -E[X_1 \mathbb{I}_{\{|X_1| > m^{1/p}\}}]$。
    通过类似的积分放缩，可以证明 $\frac{1}{n^{1/p}} \sum_{m=1}^n E(Y_m) \rightarrow 0$。
    结合 (3) 即得 $\frac{T_n}{n^{1/p}} \rightarrow 0$，最终推得 $\frac{S_n}{n^{1/p}} \xrightarrow{a.s.} 0$。 $\square$

### 4.3 期望无穷大的发散性

如果连一阶期望都不存在（即期望为无穷大），大数定律还能给出一个界限吗？

!!! success "定理 4.3 (Infinite Mean Divergence)"

    设 $X_1, X_2, \dots$ 为 i.i.d. 且 $E|X_i| = \infty$。
    令 $\{a_n\}$ 为任意一个正数列，且满足 $a_n / n$ 是递增的。那么极限上确界必然是极端的：

    \[
    \limsup_{n \rightarrow \infty} \frac{|S_n|}{a_n} = 0 \quad \text{或} \quad \infty \quad \text{a.s.}
    \]

??? proof "证明：运用 Borel-Cantelli 引理 II（点击展开）"

    对于任意常数 $k > 0$，由于 $a_n / n$ 递增且 $a_n > 0$（可以假设 $a_n$ 至少是线性增长），加上 $E|X_1| = \infty$，必定有积分发散：
    
    \[
    \sum_{n=1}^\infty P(|X_n| > k a_n) = \infty
    \]
    
    由于序列是相互独立的，根据 **Borel-Cantelli 引理 II**，事件 $\{|X_n| \ge k a_n\}$ 会无限次发生（i.o.）。
    这意味着：

    \[
    \limsup_{n \rightarrow \infty} \frac{|X_n|}{a_n} = \infty \quad a.s.
    \]

    现在考察部分和与单项的关系：$X_n = S_n - S_{n-1}$。由三角不等式：
    
    \[
    \frac{|X_n|}{a_n} \le \frac{|S_n|}{a_n} + \frac{|S_{n-1}|}{a_n} \le \frac{|S_n|}{a_n} + \frac{|S_{n-1}|}{a_{n-1}} \quad (\text{因为 } a_{n-1} \le a_n)
    \]

    如果 $\limsup \frac{|S_n|}{a_n} = M < \infty$，那么不等式右边就会被有界化为 $2M$，这与左边上确界为 $\infty$ 的事实相矛盾！
    因此，只能得出：
    
    \[
    \limsup_{n \rightarrow \infty} \frac{|S_n|}{a_n} = \infty \quad a.s.
    \]
    
    （如果 $a_n$ 增长得太快，比如超指数增长使得级数收敛，那么极限就可能坍缩为 0。所以结果只能是这两种极端情况）。 $\square$

## 5. 大偏差原理 (Large Deviation Principle)

强大数定律告诉我们 $\frac{S_n}{n}$ 会依概率 1 收敛到 $\mu$。但在统计物理或保险精算中，我们更关心：**当 $a > \mu$ 时，$P(\frac{S_n}{n} \ge a)$ 趋于 0 的速度到底有多快？**

### 5.1 动差生成函数与率函数

在研究指数级收敛前，我们需要刻画分布的“尾部性质”。

!!! info "定义 5.1 (矩母函数与对数母函数)"

    设 $X, X_1, X_2, \dots$ 为 i.i.d. 随机变量：

    * **矩母函数 (MGF)**：$\phi(\theta) = E[e^{\theta X}]$

    * **对数矩母函数**：$\Lambda(\theta) = \log \phi(\theta)$

!!! info "定义 5.2 (率函数 / 勒让德变换)"

    定义 $\Lambda(\theta)$ 的**勒让德变换 (Legendre Transform)** $I(a)$ 为：

    \[
    I(a) = \sup_{\theta \in \mathbb{R}} \{ \theta a - \Lambda(\theta) \}
    \]

    在凸分析中，$I(a)$ 也称为对偶函数。它具有下半连续性和凸性，且在 $a = E[X]$ 处取得最小值 0。

---

### 5.2 Cramér 定理 (Cramér's Theorem)

这是大偏差理论的奠基性定理，揭示了样本均值偏离期望的概率呈指数级衰减。

!!! success "定理 5.1 (Cramér's Theorem in $\mathbb{R}$)"

    设 $\{X_n\}$ 为 i.i.d. 序列，满足对所有 $\theta \in \mathbb{R}$ 都有 $\phi(\theta) < \infty$。记 $S_n = \sum_{i=1}^n X_i$。
    则对于任意 $a > E[X_1]$，有：

    \[
    \lim_{n \rightarrow \infty} \frac{1}{n} \log P\left( \frac{S_n}{n} \ge a \right) = -I(a)
    \]

    这意味着当 $n$ 很大时：$P(S_n \ge na) \approx e^{-n I(a)}$。

---

### 5.3 Cramér 定理的证明

证明分为两部分：上界（利用 Chernoff 估计）和下界（利用测度变换）。

#### 5.3.1 上界部分的证明 (The Upper Bound)

??? proof "证明：$\limsup_{n \rightarrow \infty} \frac{1}{n} \log P(S_n \ge na) \le -I(a)$（点击展开）"

    对于任意 $\theta > 0$，利用 Markov 不等式（Chernoff Bound）：

    \[
    P(S_n \ge na) = P(e^{\theta S_n} \ge e^{\theta na}) \le \frac{E[e^{\theta S_n}]}{e^{n \theta a}}
    \]

    由于 $\{X_n\}$ 是 i.i.d. 的，$E[e^{\theta S_n}] = (E[e^{\theta X_1}])^n = [\phi(\theta)]^n$。代入上式得：

    \[
    P(S_n \ge na) \le \frac{[\phi(\theta)]^n}{e^{n \theta a}} = \exp\left( -n [\theta a - \log \phi(\theta)] \right) = e^{-n (\theta a - \Lambda(\theta))}
    \]

    由于上述不等式对所有 $\theta > 0$ 都成立，我们取上确界：

    \[
    P(S_n \ge na) \le \exp\left( -n \sup_{\theta > 0} \{\theta a - \Lambda(\theta)\} \right)
    \]

    两边取对数并除以 $n$：

    \[
    \frac{1}{n} \log P(S_n \ge na) \le - \sup_{\theta > 0} \{\theta a - \Lambda(\theta)\}
    \]

    当 $a > E[X_1]$ 时，可以证明 $\sup_{\theta > 0}$ 与全局 $\sup_{\theta \in \mathbb{R}}$ 相等（因为 $\theta \le 0$ 时 $\theta a - \Lambda(\theta)$ 在 $a > E[X]$ 处不可能是最大值）。
    因此：
    
    \[
    \limsup_{n \rightarrow \infty} \frac{1}{n} \log P(S_n \ge na) \le -I(a)
    \]

    $\square$

#### 5.3.2 下界部分的证明 (The Lower Bound)

这是讲义中最精彩的部分，采用了**测度变换 (Change of Measure)**，其核心思想是：构造一个新的概率测度，使得原来的“稀有事件”在新的测度下变成“高频事件（大数定律成立的地方）”。

??? proof "证明：$\liminf_{n \rightarrow \infty} \frac{1}{n} \log P(S_n \ge na) \ge -I(a)$（点击展开）"

    **Step 1: 构造新测度 (Cramér Transform)**
    设 $F$ 为 $X$ 的原分布函数。对于选定的参数 $\lambda$，定义新的分布函数 $F_\lambda$：

    \[
    dF_\lambda(x) = \frac{e^{\lambda x}}{\phi(\lambda)} dF(x)
    \]

    在新测度 $P_\lambda$ 下，随机变量 $X$ 的期望为：
    
    \[
    E_\lambda[X] = \int x dF_\lambda(x) = \frac{\phi'(\lambda)}{\phi(\lambda)} = \Lambda'(\lambda)
    \]

    **Step 2: 选取最优参数 $\lambda_a$**
    我们选取特定的 $\lambda_a$，使得在新测度下 $E_{\lambda_a}[X] = a$。
    此时根据凸分析理论，刚好有 $I(a) = \lambda_a a - \Lambda(\lambda_a)$。

    **Step 3: 利用新测度下的 SLLN 估计原概率**
    记 $P_{\lambda_a}^n$ 为 $n$ 个独立同分布于 $F_{\lambda_a}$ 的随机变量的联合分布。则有 Radon-Nikodym 导数：

    \[
    \frac{dP^n}{dP_{\lambda_a}^n} = \prod_{i=1}^n \frac{\phi(\lambda_a)}{e^{\lambda_a X_i}} = [\phi(\lambda_a)]^n e^{-\lambda_a S_n}
    \]

    考虑小邻域 $(a, a+\epsilon)$ 上的概率：

    \[
    P(S_n \in [na, n(a+\epsilon)]) = \int_{na \le S_n \le n(a+\epsilon)} [\phi(\lambda_a)]^n e^{-\lambda_a S_n} dP_{\lambda_a}^n
    \]

    由于在该积分区域内 $S_n \le n(a+\epsilon)$，且 $e^{-\lambda_a S_n} \ge e^{-\lambda_a n(a+\epsilon)}$，得到：

    \[
    P \ge [\phi(\lambda_a)]^n e^{-\lambda_a n(a+\epsilon)} P_{\lambda_a}^n( na \le S_n \le n(a+\epsilon) )
    \]

    注意到在新测度 $P_{\lambda_a}$ 下，$E_{\lambda_a}[S_n/n] = a$。根据 **大数定律 (LLN)**：

    \[
    P_{\lambda_a}^n( a \le \frac{S_n}{n} \le a+\epsilon ) \xrightarrow{n \rightarrow \infty} \frac{1}{2} \text{ (或由 CLT 趋于常数)}
    \]

    两边取对数除以 $n$：

    \[
    \liminf \frac{1}{n} \log P \ge \log \phi(\lambda_a) - \lambda_a(a+\epsilon) + \liminf \frac{1}{n} \log \left(\frac{1}{2}\right)
    \]

    \[
    = - (\lambda_a a - \Lambda(\lambda_a)) - \lambda_a \epsilon = -I(a) - \lambda_a \epsilon
    \]

    令 $\epsilon \rightarrow 0$，即得下界 $-I(a)$。 $\square$

---

## 6. LDP 的直观理解与应用

大偏差原理在统计推断中扮演着极其重要的角色。

### 6.1 常见分布的率函数 $I(a)$

| 分布类型 | 参数 | 率函数 $I(a)$ | 物理/统计意义 |
| :--- | :--- | :--- | :--- |
| **正态分布** | $N(\mu, \sigma^2)$ | $\frac{(a-\mu)^2}{2\sigma^2}$ | 概率随欧氏距离平方呈指数衰减 |
| **伯努利分布** | $p$ | $a \log \frac{a}{p} + (1-a) \log \frac{1-a}{1-p}$ | 即 $a$ 与 $p$ 之间的 **KL 散度** |
| **泊松分布** | $\lambda$ | $a \log \frac{a}{\lambda} - a + \lambda$ | 稀有事件计数的偏差代价 |

### 6.2 在统计推断中的启示

1.  **极大似然估计 (MLE)**：在渐近意义下，MLE 的一致性可以用 LDP 来刻画其错误率收敛的速度。
2.  **假设检验 (Sanov 定理)**：第一类错误和第二类错误的概率衰减速度通常由 LDP 决定。
3.  **大样本下的“小概率”决策**：在大数据背景下（$n$ 很大），虽然样本均值趋于期望，但如果我们必须面对极端的 $a$，LDP 给出了风险的具体量化指标。

---

> *"The essence of LDP is the change of measure: it turns a miracle into a mundane reality."*