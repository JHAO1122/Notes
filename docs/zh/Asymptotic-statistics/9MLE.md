# 第十章：极大似然估计（二）

在上一章中，我们探讨了极大似然估计 (MLE) 的相合性与渐近正态性。本章我们将从著名的 Wilks 定理出发，深入研究似然比检验的渐近性质。接着，我们会讨论如何利用 Bartlett 修正来提高似然比检验的精度，并引入参数不可识别时的反例。最后，我们将跨越参数模型的限制，介绍现代非参数统计的核心工具——**经验似然 (Empirical Likelihood)**，并展示其优美的数学结构与 Wilks 定理的非参数版本。

---

## 1. 似然比的 Wilks 定理 (Wilks' Theorem)

考虑复合假设检验问题 $H_0: \theta \in \Theta_0$ VS $H_1: \theta \in \Theta_1$，我们需要明确原假设空间 $\Theta_0$ 的数学结构。

### 1.1 $\Theta_0$ 的结构设定 (Condition S1)

假设存在一个具有内点的集合 $A \subset \mathbb{R}^r$（即 $r$ 维集合），以及定义在 $\mathbb{R}^r$ 上的 $p$ 个函数 $g_i$，使得 $\theta_i = g_i(\phi_1, ..., \phi_r)$ 对 $i=1,...,p$ 成立。
这为 $\Theta_0$ 和 $A$ 之间建立了 $1-1$ 的对应关系。并且每个 $g_i$ 在 $A$ 的任何内点处都是二次连续可导的。在 $H_0$ 下的真实参数 $\theta_0 \in \Theta_0$ 对应于 $A$ 中的一个内点 $\phi_0 \in A$。

这一条件意味着，原假设下的分布族可以等价地参数化为：

\[
\{ \tilde{f}(x, \phi) = f(x, g_1(\phi_1, ..., \phi_r), ..., g_p(\phi_1, ..., \phi_r)) : \phi = (\phi_1, ..., \phi_r) \in A \}
\]

**等价表述：**
$\Theta_0$ 也可以由 $p-r$ 个约束条件等价定义：

\[
R_i(\theta) = 0, \quad 1 \le i \le p-r
\]

* **简单假设的例子**：对于 $H_0: \theta = \theta_0$，有 $\Theta_0 = \{\theta_0\}$，约束函数可取为 $R_i(\theta) = \theta_i - \theta_{i0}$ ($1 \le i \le p$)。
* **复合假设的例子**：若 $p=3$ 且 $H_0: \theta_1 = \theta_{10}$（只限制了第一个参数），则约束条件为 $R_1(\theta) = \theta_1 - \theta_{10}$。

### 1.2 Wilks 定理与三大检验统计量

令 $\hat{\theta}_n$ 为全局极大似然估计，$\hat{\phi}_n$ 为在族 $\{f(x, \phi), \phi \in A\}$ 下的受限 MLE，且对应的受限参数估计为 $\tilde{\theta}_n = g(\hat{\phi}_n)$。

定义对数似然比为：

\[
LR_n = -2 \log \frac{L_n(\tilde{\theta}_n)}{L_n(\hat{\theta}_n)}
\]

!!! success "定理 9.1 (Wilks 定理, Wilks 1938)"

    在条件 S1 以及全局 MLE 渐近正态性的正则条件下，原假设 $H_0$ 成立时有：

    \[
    LR_n \xrightarrow{d} \chi^2_{p-r}
    \]

在原假设下，以下三大检验统计量是**渐近等价**的：

1.  **似然比统计量 (LR statistic)**: $LR_n$

2.  **Wald 检验统计量 (基于 Delta 方法)**:

    \[
    \mathcal{W}_n = R(\hat{\theta}_n)^T \left( \frac{\partial}{\partial \theta} R(\hat{\theta}_n) \left( -\frac{\partial^2 \log L_n(\hat{\theta}_n)}{\partial \theta \partial \theta^T} \right)^{-1} \frac{\partial}{\partial \theta} R(\hat{\theta}_n)^T \right)^{-1} R(\hat{\theta}_n)
    \]

3.  **Rao 得分检验统计量 (Lagrangian 乘子法)**:

    \[
    R_n = \frac{\partial \log L_n(\tilde{\theta}_n)}{\partial \theta^T} \left( -\frac{\partial^2 \log L_n(\tilde{\theta}_n)}{\partial \theta \partial \theta^T} \right)^{-1} \frac{\partial \log L_n(\tilde{\theta}_n)}{\partial \theta}
    \]


??? proof "Wilks 定理的证明概要（点击展开）"

    在 $H_0$ 下，通过对数似然函数的二阶泰勒展开，可以得到：

    \[
    2(\log L_n(\hat{\theta}_n) - \log L_n(\theta_0)) = \frac{1}{n} \left(\frac{\partial \log L_n(\theta_0)}{\partial \theta}\right)^T I_1(\theta_0)^{-1} \left(\frac{\partial \log L_n(\theta_0)}{\partial \theta}\right) + o_p(1)
    \]

    类似地，使用 $H_0$ 下的参数化 $\phi$，受限对数似然的展开为：

    \[
    2(\log L_n(g(\hat{\phi}_n)) - \log L_n(g(\phi_0))) = \frac{1}{n} \left(\frac{\partial \log L_n(g(\phi_0))}{\partial \phi}\right)^T \tilde{I}_1(\phi_0)^{-1} \left(\frac{\partial \log L_n(g(\phi_0))}{\partial \phi}\right) + o_p(1)
    \]

    根据链式法则：

    \[
    \frac{\partial \log L_n(g(\phi_0))}{\partial \phi} = \frac{\partial g(\phi_0)^T}{\partial \phi} \frac{\partial \log L_n(\theta_0)}{\partial \theta}
    \]

    结合上述结果，似然比统计量 $LR_n$ 可以写为：

    \[
    LR_n = \frac{1}{n} \left( \frac{\partial \log L_n(\theta_0)}{\partial \theta} \right)^T \left( I_1(\theta_0)^{-1} - \frac{\partial g(\phi_0)}{\partial \phi^T} \tilde{I}_1(\phi_0)^{-1} \frac{\partial g(\phi_0)^T}{\partial \phi} \right) \left( \frac{\partial \log L_n(\theta_0)}{\partial \theta} \right) + o_p(1)
    \]

    由中心极限定理 (CLT) 可知，得分函数服从：

    \[
    \frac{1}{\sqrt{n}} I_1(\theta_0)^{-1/2} \frac{\partial \log L_n(\theta_0)}{\partial \theta} \xrightarrow{d} Z \sim N(0, I_{p \times p})
    \]

    记 $A = I_1(\theta_0)$，$C = \tilde{I}_1(\phi_0)$，$D = \frac{\partial g(\phi_0)^T}{\partial \phi}$，以及 $B = A^{-1} - D^T C^{-1} D$。可以验证 $C = D A D^T$。
    因此 $LR_n \Rightarrow Z^T A^{1/2} B A^{1/2} Z$。

    我们需要证明矩阵 $A^{1/2} B A^{1/2}$ 是**幂等矩阵 (Idempotent)**：

    \[
    (A^{1/2} B A^{1/2})^2 = A^{1/2} B A B A^{1/2}
    \]

    \[
    = A^{1/2} (A^{-1} - D^T C^{-1} D) A (A^{-1} - D^T C^{-1} D) A^{1/2}
    \]
    
    \[
    = (I_p - A^{1/2} D^T C^{-1} D A^{1/2})(I_p - A^{1/2} D^T C^{-1} D A^{1/2})
    \]
    
    \[
    = I_p - 2 A^{1/2} D^T C^{-1} D A^{1/2} + A^{1/2} D^T C^{-1} D A D^T C^{-1} D A^{1/2}
    \]

    由于 $D A D^T = C$，最后一项化简为 $A^{1/2} D^T C^{-1} C C^{-1} D A^{1/2} = A^{1/2} D^T C^{-1} D A^{1/2}$。代回得到：

    \[
    = I_p - A^{1/2} D^T C^{-1} D A^{1/2} = A^{1/2} B A^{1/2}
    \]

    由于该矩阵是幂等的，其二次型服从卡方分布，自由度等于其迹 (Trace)，即 $Tr(A^{1/2} B A^{1/2}) = p - r$。因此 $LR_n \xrightarrow{d} \chi^2_{p-r}$。

---

## 2. 讨厌参数、置信区域与反例

### 2.1 讨厌参数情形 (Nuisance Parameter Case)

假设 $\theta = (\beta, \gamma)$，其中 $\beta$ 是我们感兴趣的参数，而 $\gamma$ 是**讨厌参数 (nuisance parameter)**。
测试 $H_0: \beta = \beta_0$ VS $H_1: \beta \ne \beta_0$。

似然比统计量为：

\[
LR_n(\beta_0) = 2 \{ l_n(\hat{\beta}_n, \hat{\gamma}_n) - l_n(\beta_0, \tilde{\gamma}_n) \}
\]

其中 $\hat{\theta}_n = (\hat{\beta}_n, \hat{\gamma}_n)$ 是无约束的全局 MLE，而 $\tilde{\gamma}_n = \arg\sup_{\gamma} l_n(\beta_0, \gamma)$ 是受限 MLE。
假设参数空间满足 $\Theta = \Theta_\beta \oplus \Theta_\gamma$，且 $dim(\Theta_\gamma) = r$。根据 Wilks 定理，有：

\[
LR_n(\beta_0) \xrightarrow{d} \chi^2_{p-r}
\]

假设检验的“局部线性” (local linearity) 结构对于卡方逼近是至关重要的。

### 2.2 似然比置信区域

利用似然比统计量，我们可以构建 $\beta$ 的 $(1-\alpha)$ 级别的置信区域：

\[
I_{1-\alpha}(\beta) = \{ \beta | LR_n(\beta) \le \chi^2_{p-r, 1-\alpha} \}
\]

这种置信区域具有**自然的形状和方向**，完全由数据和模型自动决定，无需主观人为选择边界（这与使用 Wald 检验构造的对称椭球不同）。虽然“画出”这个置信区域可能需要计算或图形工具，但现代计算能力已让其变得十分普及。

### 2.3 Wilks 定理失效的反例

Wilks 定理依赖于严格的正则性条件。如果不满足，会出现异常行为。

**反例**：设 $X_1, \dots, X_n \stackrel{i.i.d}{\sim} \text{Uniform}(0, \theta)$，$\theta > 0$。
这里，密度函数的支撑集边界依赖于参数 $\theta$（违反了定理的基本正则性条件）。

欲检验 $H_0: \theta = \theta_0$ VS $H_1: \theta \ne \theta_0$。
无约束的 MLE 为最大次序统计量 $\hat{\theta}_n = X_{(n)}$。
似然比统计量为：

\[
LR_n = -2 \log \frac{\theta_0^{-n}}{X_{(n)}^{-n}} = -2n \log \frac{X_{(n)}}{\theta_0}
\]

可以证明，此时 $LR_n$ 的渐近分布是 $\chi^2_2$ 的某种变形或指数分布相关，而**非**标准的 $\chi^2_1$。

---

## 3. 似然比的 Bartlett 修正 (Bartlett Correction)

Bartlett (1937) 发现了一种可以提升似然比检验渐近分布逼近精度的修正方法。

在 $H_0$ 下，$LR_n \rightarrow \chi^2_{p-r}$。假设 $LR_n$ 的均值存在如下渐近展开：

\[
E(LR_n) = (p-r) + \frac{b_1}{n} + o(n^{-1}) = (p-r)\left\{ 1 + \frac{b}{n} + o(n^{-1}) \right\}
\]

这里 $b$ 或 $(1 + b/n)$ 被称为 **Bartlett 因子**。设 $\hat{b}_n$ 是 $b$ 的相合估计量（$\hat{b}_n \rightarrow b$），通过 Slutsky 定理，修正后的统计量：

\[
LR_{n, bc} := \frac{LR_n}{1 + b/n} \xrightarrow{d} \chi^2_{p-r}
\]

\[
\hat{LR}_{n, bc} := \frac{LR_n}{1 + \hat{b}/n} \xrightarrow{d} \chi^2_{p-r}
\]

**Bartlett 猜想与 Edgeworth 展开**：
Wilks 定理仅仅保证了 $P(LR_n \le x) = P(\chi^2_{p-r} \le x) + o(1)$。
如果假设数据尾部具有足够的平滑性且满足更多矩条件 (Cramér 条件)，通过 Edgeworth 展开可以证明（Barndorff-Nielsen and Cox, 1984）：

\[
P(LR_n(\theta_0) \le x) = P(\chi^2_{p-r} \le x) - \beta x g_p(x) n^{-1} + O(n^{-2})
\]

其中 $g_p$ 是 $\chi^2_{p-r}$ 的概率密度函数。
而经过 Bartlett 修正后，收敛速度大幅提高：

\[
P\left(\frac{LR_n}{1 + \beta n^{-1}} \le x\right) = P(\chi^2_{p-r} \le x) + O(n^{-2})
\]

\[
P\left(\frac{LR_n}{1 + \hat{\beta} n^{-1}} \le x\right) = P(\chi^2_{p-r} \le x) + O(n^{-3/2})
\]

**Bartlett 修正的意义**：
利用修正后的统计量构造置信区域 $I_{\alpha, bc}$，其覆盖误差 (coverage error) $P(\theta \in I_{\alpha, bc}) - (1-\alpha) = O(n^{-3/2})$，比未修正的 $O(n^{-1})$ 具有更小的高阶误差。它为似然比检验提供了极为精确的近似显著性水平。

---

## 4. 从参数到非参数：经验似然 (Empirical Likelihood)

在标准参数 MLE 中，我们假设数据来自某个已知分布族 $F_\theta$。然而，如果模型错误设定 (mis-specified)，MLE 可能会严重有偏且不相合。

**问题**：能否不假设具体的参数分布族来构造似然函数？
**直觉**：让数据“为自己发声”。经验累积分布函数 (ECDF)：

\[
F_n(x) = \frac{1}{n} \sum_{i=1}^n I(X_i \le x)
\]

是真实 CDF $F(x)$ 的一个自然的、相合的非参数估计。能否定义一个使得 $F_n$ 达到最大化的似然函数？

!!! info "定义 9.1 (分布的非参数似然 Nonparametric Likelihood)"

    设 $X_1, ..., X_n$ 为来自完全未知分布 $F$ 的 i.i.d. 样本。给定数据，分布 $F$ 的非参数似然定义为：

    \[
    L(F) = \prod_{i=1}^n dF(X_i)
    \]

    其中 $dF(X_i) = P_F(X = X_i)$ 是分布 $F$ 赋予观测点 $X_i$ 的概率质量。

为了最大化 $L(F)$，分布 $F$ 必须将其所有的概率质量严格放置在已观测到的数据点上。
令 $p_i = dF(X_i)$。我们需要满足 $p_i \ge 0$ 且 $\sum_{i=1}^n p_i = 1$。似然函数简化为 $L(p) = \prod_{i=1}^n p_i$。

!!! success "定理 9.2 (ECDF 最大化非参数似然)"

    经验分布 $F_n(x)$ 在所有可能的分布中最大化了非参数似然 $L(F)$。

??? proof "证明：利用算术-几何均值不等式 (AM-GM Inequality)"

    我们要在约束 $p_i \ge 0$ 且 $\sum_{i=1}^n p_i = 1$ 下最大化 $\prod_{i=1}^n p_i$。
    根据 AM-GM 不等式：

    \[
    \left( \prod_{i=1}^n p_i \right)^{1/n} \le \frac{1}{n} \sum_{i=1}^n p_i = \frac{1}{n}
    \]

    当且仅当 $p_1 = p_2 = \dots = p_n = \frac{1}{n}$ 时，等号成立。
    这意味着在没有任何约束下，非参数似然在赋予每个观测点 $1/n$ 的质量时达到最大，这正好恢复了经验分布 $F_n$。此时最大似然值为 $(1/n)^n$。

!!! info "定义 9.2 (经验似然比 Empirical Likelihood Ratio)"

    由 Art Owen (1988, 1990) 提出，候选分布 $F$ 的经验似然比定义为其似然值与无约束最大似然值的比值：

    \[
    R(F) = \frac{L(F)}{L(F_n)} = \frac{\prod_{i=1}^n p_i}{\prod_{i=1}^n (1/n)} = \prod_{i=1}^n (n p_i)
    \]

    注意：对于所有合法的分布 $F$，都有 $R(F) \le 1$。这个比值衡量了候选分布 $F$ 相较于无约束的经验分布 $F_n$ 的“合理性”。

---

## 5. 均值的经验似然与优化求解

假设我们对参数向量感兴趣，例如均值 $\mu = E_F(X) = \int x dF(x)$。
$\mu$ 的**轮廓经验似然比 (Profile empirical likelihood ratio)** $\mathcal{R}(\mu)$ 定义为在满足参数约束的所有分布 $F$ 中最大化 $R(F)$：

\[
\mathcal{R}(\mu) = \sup \left\{ \prod_{i=1}^n (n p_i) \Bigg| p_i \ge 0, \sum_{i=1}^n p_i = 1, \sum_{i=1}^n p_i X_i = \mu \right\}
\]

这把统计推断转变成了一个**约束优化问题**。我们需要寻找一组在似然意义上最接近 $1/n$ 的权重 $p_i$，同时迫使其期望值等于候选的 $\mu$。

### 拉格朗日乘子法求解 (Lagrange Multipliers)

我们要最大化对数似然比：

\[
\log \mathcal{R}(\mu) = \sum_{i=1}^n \log(n p_i)
\]

受限于 $\sum p_i = 1$ 和 $\sum p_i(X_i - \mu) = 0$。引入拉格朗日乘子 $\gamma$ 和 $n\lambda$（这里用 $n$ 缩放 $\lambda$ 纯粹为了后续代数化简方便），构造拉格朗日函数：

\[
\mathcal{L}(p, \gamma, \lambda) = \sum_{i=1}^n \log(n p_i) - \gamma \left( \sum_{i=1}^n p_i - 1 \right) - n\lambda^T \left( \sum_{i=1}^n p_i(X_i - \mu) \right)
\]

对 $p_i$ 求偏导并令其为 0：

\[
\frac{\partial \mathcal{L}}{\partial p_i} = \frac{1}{p_i} - \gamma - n\lambda^T (X_i - \mu) = 0
\]

将整个方程同乘 $p_i$ 并对所有的 $i$ 求和：

\[
\sum_{i=1}^n 1 - \gamma \sum_{i=1}^n p_i - n\lambda^T \sum_{i=1}^n p_i(X_i - \mu) = 0
\]

利用约束条件 $\sum p_i = 1$ 和 $\sum p_i(X_i - \mu) = 0$，方程极为优美地坍缩为：

\[
n - \gamma(1) - n\lambda^T(0) = 0 \implies \gamma = n
\]

将 $\gamma = n$ 代回一阶条件并解出 $p_i$，我们得到了最优的概率权重：

\[
p_i = \frac{1}{n} \frac{1}{1 + \lambda^T(X_i - \mu)}
\]

拉格朗日乘子 $\lambda$ 则由剩下的矩约束方程隐式决定。代入 $p_i$ 后得到关于 $\lambda$ 的**得分方程**：

\[
\frac{1}{n} \sum_{i=1}^n \frac{X_i - \mu}{1 + \lambda^T(X_i - \mu)} = 0
\]

将最优权重 $p_i$ 代回目标函数，我们得到了**经验对数似然比统计量 (Profile Empirical Log-Likelihood Ratio Statistic)**：

\[
l_E(\mu) = -2 \log \mathcal{R}(\mu) = -2 \sum_{i=1}^n \log(n p_i) = -2 \sum_{i=1}^n \log \left( \frac{1}{1 + \lambda^T(X_i - \mu)} \right) = 2 \sum_{i=1}^n \log \{ 1 + \lambda^T(X_i - \mu) \}
\]

令人惊奇的是，这个函数在渐近行为上与传统的参数似然比统计量极为相似。

---

## 6. 经验似然的 Wilks 定理

!!! success "定理 9.3 (经验似然的 Wilks 定理, Owen 1988, 1990)"

    设 $X_1, ..., X_n$ 为 $\mathbb{R}^d$ 中的 i.i.d. 随机向量，具有真实的均值 $\mu_0$ 以及有限且正定的协方差矩阵 $\Sigma$。
    那么，当 $n \rightarrow \infty$ 时：

    \[
    l_E(\mu_0) = -2 \log \mathcal{R}(\mu_0) \xrightarrow{d} \chi^2_d
    \]

**该定理的非凡意义**：我们获得了一个非参数的枢轴量 (pivot)。要为 $\mu$ 构建置信区间，我们只需要依赖于 $\chi^2$ 的分位数，既不需要假设任何参数分布，**也根本不需要显式地去估计协方差矩阵 $\Sigma$**！

??? proof "经验似然 Wilks 定理的证明（以一维 $d=1$ 为例，点击展开）"

    为了代数上的简便，假设维度 $d=1$。令 $Z_i = X_i - \mu_0$。由定义，$E(Z_i) = 0$ 且 $E(Z_i^2) = \sigma^2$。

    决定 $\lambda$ 的方程为 $g(\lambda) = \frac{1}{n} \sum_{i=1}^n \frac{Z_i}{1 + \lambda Z_i} = 0$。应用几何级数展开 $\frac{1}{1+x} = 1 - x + \frac{x^2}{1+x}$：

    \[
    0 = \frac{1}{n} \sum_{i=1}^n Z_i \left( 1 - \lambda Z_i + \frac{\lambda^2 Z_i^2}{1 + \lambda Z_i} \right) = \bar{Z} - \lambda \left( \frac{1}{n} \sum_{i=1}^n Z_i^2 \right) + \frac{1}{n} \sum_{i=1}^n \frac{\lambda^2 Z_i^3}{1 + \lambda Z_i}
    \]

    由于方差有限，有 $\max_i |Z_i| = o_p(n^{1/2})$，误差项可以被控制。解出 $\lambda$ 得到其渐近表示：

    \[
    \lambda = \frac{\bar{Z}}{S_Z^2} + o_p(n^{-1/2}) = O_p(n^{-1/2})
    \]

    其中 $S_Z^2 = \frac{1}{n}\sum Z_i^2 \xrightarrow{p} \sigma^2$。

    接下来，我们使用泰勒级数展开经验对数似然比统计量 $l_E(\mu_0)$：$\log(1+x) = x - x^2/2 + x^3/3 - \dots$
    因为 $\lambda = O_p(n^{-1/2})$ 且 $\max |Z_i| = o_p(n^{1/2})$，所以 $\max |\lambda Z_i| \rightarrow 0$，允许安全的泰勒展开：

    \[
    l_E(\mu_0) = 2 \sum_{i=1}^n \log(1 + \lambda Z_i) = 2 \sum_{i=1}^n \left( \lambda Z_i - \frac{1}{2}\lambda^2 Z_i^2 + \eta_i \right)
    \]

    余项 $2 \sum \eta_i \le C |\lambda|^3 \sum |Z_i|^3 = o_p(1)$。因此：

    \[
    l_E(\mu_0) = 2\lambda \sum_{i=1}^n Z_i - \lambda^2 \sum_{i=1}^n Z_i^2 + o_p(1)
    \]

    将主导项 $\lambda \approx \frac{\bar{Z}}{\sigma^2}$ 和 $\sum Z_i^2 \approx n\sigma^2$ 代入展开式：

    \[
    l_E(\mu_0) \approx 2 \left(\frac{\bar{Z}}{\sigma^2}\right)(n\bar{Z}) - \left(\frac{\bar{Z}}{\sigma^2}\right)^2 (n\sigma^2) = \frac{2n\bar{Z}^2}{\sigma^2} - \frac{n\bar{Z}^2}{\sigma^2} = \frac{n\bar{Z}^2}{\sigma^2}
    \]

    我们通过代数将其成功化简为：

    \[
    l_E(\mu_0) \approx \left( \frac{\sqrt{n}\bar{Z}}{\sigma} \right)^2 = \left( \frac{\sqrt{n}(\bar{X} - \mu_0)}{\sigma} \right)^2
    \]

    由经典的中心极限定理 (CLT) 已知 $\frac{\sqrt{n}(\bar{X} - \mu_0)}{\sigma} \xrightarrow{d} N(0,1)$。
    根据连续映射定理，标准正态分布随机变量的平方收敛于自由度为 1 的卡方分布：

    \[
    -2 \log \mathcal{R}(\mu_0) \xrightarrow{d} \chi^2_1
    \]

### 经验似然的核心优势

1.  **数据驱动的形状 (Data-Driven Shape)**：与 Wald 类型的区间（如 $\bar{X} \pm 1.96 \frac{\sigma}{\sqrt{n}}$）不同，经验似然的置信区域完全由数据决定。它们不需要是对称的，能够完美地捕捉到潜在分布的偏态 (skewness)。

2.  **保域性 (Range Preserving)**：EL 区间严格限制在数据的凸包 (convex hull) 内部。估计结果在逻辑上永远不会违背自然边界（例如，不会出现负的概率值）。

3.  **内在学生化 (Internal Studentization)**：方差的缩放（学生化）在经验似然中是隐含自动完成的。我们无需显式地去插入一个方差 $\Sigma$ 的估计量。

4.  **变换不变性 (Transformation Respecting)**：经验似然区域在平滑的参数重参数化 (reparameterization) 之下保持不变。