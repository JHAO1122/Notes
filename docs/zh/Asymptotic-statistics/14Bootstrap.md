# 第十四章：Bootstrap (自助法)

---

## 1. 从 Wald 推断到 Bootstrap

在前面的章节中，对于许多估计量，我们已经得出了它们满足渐近正态性的结论。例如：

\[
\sqrt{n}(\hat{\theta}_n - \theta_0) \rightarrow N(0, \Sigma)
\]

对于正则的 M-估计量，我们有：

\[
\sqrt{n}(\hat{\theta}_n - \theta_0) \rightarrow N(0, A^{-1}BA^{-1})
\]

对于正则的极大似然估计量 (MLE)，我们有：

\[
\sqrt{n}(\hat{\theta}_n - \theta_0) \rightarrow N(0, I(\theta_0)^{-1})
\]

**提出问题：** 既然我们已经有了渐近正态性，为什么我们仍然需要 Bootstrap 呢？

### 1.1 一阶近似 (First-order approximation)

令 $T_n$ 为标准化或学生化 (studentized) 的统计量：

\[
T_n = \frac{\hat{\theta}_n - \theta_0}{\hat{se}_n}
\]

通常的 Wald 近似为：

\[
P(T_n \le x) \approx \Phi(x)
\]

这被称为**一阶近似 (first-order approximation)**。在许多正则问题中，误差阶数为：

\[
P(T_n \le x) = \Phi(x) + O(n^{-1/2})
\]

* 这种近似成功捕捉了极限分布的中心和方差。

* 但是它通常忽略了偏度 (skewness)、偏差 (bias) 以及曲率效应 (curvature effects)。

* 这些被忽略的效应通常正是 $O(n^{-1/2})$ 阶的。

### 1.2 Edgeworth 展开 (Edgeworth Expansion)

为了获得更高阶的近似，我们可以使用 Edgeworth 展开：

\[
P(T_n \le x) = \Phi(x) + n^{-1/2} p_1(x) \phi(x) + n^{-1} p_2(x) \phi(x) + \dots + O(n^{-k/2})
\]

* 其中 $p_1(x)$ 是一个与总体的**偏度**相关的多项式。

* $p_2(x)$ 是一个与总体的**峰度 (kurtosis)** 相关的多项式。

* 然而，这些多项式 $p_1, p_2$ 依赖于未知的总体矩 (population moments)，在实际应用中很难解析地计算出来。

---

## 2. Bootstrap 的核心思想

由 Efron (1979) 提出的 Bootstrap（重抽样/自助法）提供了一种极其绝妙的解决思路。

**核心思想：** 用经验分布 $\mathbb{P}_n$ 去替代未知的总体分布 $P$。
从 $\mathbb{P}_n$ 中进行抽样，以此来模仿从真实分布 $P$ 中抽样的过程。

最令人惊叹的是，这种极其简单的数据驱动操作，**能够自动地捕获高阶项 (higher-order terms)！**

### 2.1 真实世界与 Bootstrap 世界的对应关系

* **真实世界 (The Real World)**：
  未知的总体分布 $P \implies$ 抽取样本 $X_1, \dots, X_n \implies$ 计算估计量 $\hat{\theta}_n$

* **Bootstrap 世界 (The Bootstrap World)**：
  经验分布 $\mathbb{P}_n \implies$ 抽取 Bootstrap 样本 $X_1^*, \dots, X_n^* \implies$ 计算 Bootstrap 估计量 $\hat{\theta}_n^*$

Bootstrap 推断的基石在于：**在 Bootstrap 世界中 $\hat{\theta}_n^*$ 相对于 $\hat{\theta}_n$ 的分布关系，能够极好地模仿真实世界中 $\hat{\theta}_n$ 相对于 $\theta_0$ 的分布关系。**

### 2.2 Bootstrap 的两大主要用途

* **1. 估计分布和标准误 (Standard Errors)**：当估计量的解析方差形式过于复杂或难以推导时，Bootstrap 提供了一种基于计算的非参数估计方法。

* **2. 改进近似精度 (Refining approximations)**：用于替代一阶的正态近似，自动实现 Edgeworth 修正，从而得到具有更小误差界（如 $O(n^{-1})$ 而不是 $O(n^{-1/2})$）的置信区间。

### 2.3 符号说明 (Notations)

* **原始数据**：$X_1, \dots, X_n \sim P$

* **经验测度**：$\mathbb{P}_n = \frac{1}{n} \sum_{i=1}^n \delta_{X_i}$

* **Bootstrap 样本**：$X_1^*, \dots, X_n^*$ 是从 $\mathbb{P}_n$ 中**有放回地 (with replacement)** 独立抽取的样本。即 $X_i^* \sim \mathbb{P}_n$。

* **Bootstrap 经验测度**：$\mathbb{P}_n^* = \frac{1}{n} \sum_{i=1}^n \delta_{X_i^*}$

* **估计量**：$\hat{\theta}_n = T(\mathbb{P}_n)$

* **Bootstrap 估计量**：$\hat{\theta}_n^* = T(\mathbb{P}_n^*)$

---

## 3. Bootstrap 的相合性 (Consistency)

!!! info "定义 10.1 (Bootstrap 的相合性)"

    如果给定原始数据下，$\sqrt{n}(\hat{\theta}_n^* - \hat{\theta}_n)$ 的条件分布弱收敛于与无条件分布 $\sqrt{n}(\hat{\theta}_n - \theta_0)$ 相同的极限，那么我们称该 Bootstrap 是**相合的 (consistent)**。

具体而言，定义真实的累积分布函数为：

\[
R_n(x, P) = P(\sqrt{n}(\hat{\theta}_n - \theta_0) \le x)
\]

Bootstrap 估计的分布函数（条件分布）为：

\[
R_n(x, \mathbb{P}_n) = P^*(\sqrt{n}(\hat{\theta}_n^* - \hat{\theta}_n) \le x \mid X_1, \dots, X_n)
\]

那么 Bootstrap 相合意味着：

\[
\sup_x |R_n(x, \mathbb{P}_n) - R_n(x, P)| \xrightarrow{P} 0 \quad (\text{或者几乎必然成立 a.s.})
\]

### 3.1 示例：样本均值的 Bootstrap

!!! example "例 10.2 (样本均值 Sample Mean)"

    设 $X_i \sim P$，均值为 $\mu$，方差为 $\sigma^2$。目标参数 $\theta_0 = \mu$，估计量 $\hat{\theta}_n = \bar{X}_n$。
    
    由中心极限定理 (CLT) 可知，在真实世界中：

    \[
    \sqrt{n}(\bar{X}_n - \mu) \xrightarrow{d} N(0, \sigma^2)
    \]

    现在考察 Bootstrap 世界：抽取 $X_1^*, \dots, X_n^* \sim \mathbb{P}_n$。
    
    Bootstrap 估计量为 $\hat{\theta}_n^* = \bar{X}_n^* = \frac{1}{n}\sum_{i=1}^n X_i^*$。

    给定原始数据 $(X_1, \dots, X_n)$，计算条件期望和条件方差：

    * **条件期望**：
    
    \[
    E^*(\bar{X}_n^*) = E^*(X_1^*) = \frac{1}{n} \sum_{i=1}^n X_i = \bar{X}_n
    \]

    * **条件方差**：
    
    \[
    Var^*(\sqrt{n}\bar{X}_n^*) = Var^*(X_1^*) = \frac{1}{n} \sum_{i=1}^n (X_i - \bar{X}_n)^2 = \hat{\sigma}_n^2
    \]

### 3.2 样本均值的 Bootstrap 中心极限定理

由于 $X_1^*, \dots, X_n^*$ 在给定原始数据的条件下是独立同分布 (i.i.d.) 的，我们可以对 Bootstrap 样本应用 Lindeberg-Feller 中心极限定理或 Berry-Esseen 定理。

又因为由大数定律可知样本方差是相合的：$\hat{\sigma}_n^2 \xrightarrow{p} \sigma^2$，所以：

\[
\sup_x \left| P^*(\sqrt{n}(\bar{X}_n^* - \bar{X}_n) \le x) - \Phi(x/\hat{\sigma}_n) \right| \xrightarrow{P} 0
\]

同时，由于 $\hat{\sigma}_n \xrightarrow{p} \sigma$ 以及正态累积分布函数 $\Phi$ 的连续性：

\[
\Phi(x/\hat{\sigma}_n) \xrightarrow{P} \Phi(x/\sigma)
\]

由此得出结论，样本均值的 Bootstrap 是相合的。

**注记：为什么不直接使用 $N(0, \hat{\sigma}_n^2)$ 呢？**
对于简单的样本均值而言，直接使用正态近似确实已经足够。但是，对于复杂的估计量（如样本中位数、一般的 M-估计量等），其极限方差 $\Sigma$ 的解析形式可能非常复杂或难以准确估计。Bootstrap 巧妙地绕过了显式计算极限方差这一步骤。

---

## 4. Bootstrap 失效的例子：最大值

Bootstrap **并不是**在所有情况下都相合的。对于极值或处于参数边界上的统计量，标准的 Bootstrap 会失效。

!!! warning "例 10.3 (Bootstrap 失效：最大值 The Maximum)"

    假设 $X_i \sim U(0, \theta)$ 是来自均匀分布的 i.i.d. 样本。
    
    我们要估计参数 $\theta$，自然的选择是最大值统计量：
    
    \[
    \hat{\theta}_n = X_{(n)} = \max_{1 \le i \le n}(X_i)
    \]

    在真实世界中，最大值具有如下的极值极限分布：

    \[
    n(\theta - X_{(n)}) \xrightarrow{d} \text{Exp}(1/\theta)
    \]

    其中极限分布是一个指数分布 (Exponential distribution)。

### 4.1 最大值在 Bootstrap 世界中的表现

现在让我们看看在 Bootstrap 世界中会发生什么。
Bootstrap 估计量为 $\hat{\theta}_n^* = X_{(n)}^* = \max_{1 \le i \le n}(X_i^*)$。

我们来计算 Bootstrap 最大值恰好等于原始样本最大值的概率：

\[
P^*(X_{(n)}^* = X_{(n)}) = P^*(\text{至少有一个 } X_i^* \text{ 等于 } X_{(n)})
\]

这等于 $1$ 减去所有 Bootstrap 样本都不等于 $X_{(n)}$ 的概率。因为每次抽样抽不到 $X_{(n)}$ 的概率是 $1 - 1/n$：

\[
P^*(X_{(n)}^* = X_{(n)}) = 1 - \left(1 - \frac{1}{n}\right)^n
\]

当 $n \rightarrow \infty$ 时，我们有极其著名的极限：

\[
1 - \left(1 - \frac{1}{n}\right)^n \rightarrow 1 - e^{-1} \approx 0.632
\]

* **真实世界**：由于总体是连续分布，$P(X_{(n)} = \theta) = 0$。

* **Bootstrap 世界**：Bootstrap 分布在点 0 处具有一个概率约为 $0.632$ 的巨大**点质量 (point mass)**！

显然，包含一个巨大点质量的 Bootstrap 分布，绝对不可能收敛到一个连续的指数分布。因此，**Bootstrap 在此失效**。

### 4.2 为什么会失效？

其根本原因在于离散性与连续性的冲突：
经验分布 $\mathbb{P}_n$ 本质上是一个**离散分布 (discrete distribution)**，而真实的总体分布 $P$ 在边界处是**连续分布 (continuous distribution)**。
极值统计量对于经验分布这种局部的离散性质极其敏感。

*(解决这类问题的常用方法是采用 "m-out-of-n bootstrap" 或者子抽样法 Subsampling，即重抽样规模 $m$ 远小于原样本量 $n$。)*

---

## 5. Bootstrap 的 Delta 方法

如果一个统计量 $T_n$ 是 Bootstrap 相合的，那么经过非线性光滑变换后的 $g(T_n)$ 还是相合的吗？答案是肯定的。

!!! success "定理 10.4 (Bootstrap 的 Delta 方法)"

    假设在真实世界中 $\sqrt{n}(T_n - \theta) \xrightarrow{d} T$。
    
    设函数 $g$ 在 $\theta$ 处可微，且其导数为 $g'(\theta)$。
    
    如果 $T_n$ 是 Bootstrap 相合的，那么对于 $g(T_n)$，Bootstrap 也是相合的。即在给定数据的条件下，有：

    \[
    \sqrt{n}(g(T_n^*) - g(T_n)) \xrightarrow{d} g'(\theta) T \quad \text{(在 } P^* \text{ 测度下)}
    \]

??? proof "定理 10.4 的证明概要（点击展开）"

    利用均值定理 (Mean Value Theorem) 或泰勒展开：

    \[
    g(T_n^*) - g(T_n) = g'(\tilde{T}_n^*)(T_n^* - T_n)
    \]

    其中 $\tilde{T}_n^*$ 介于 $T_n^*$ 和 $T_n$ 之间。

    因为 $T_n$ 是 Bootstrap 相合的，在条件测度下 $\sqrt{n}(T_n^* - T_n) = O_{P^*}(1)$，即它在 $P^*$ 概率下有界。
    
    同时由于真实世界中 $T_n \xrightarrow{p} \theta$，这蕴含了 $T_n^* \xrightarrow{p^*} \theta$。

    由于 $\tilde{T}_n^*$ 被夹在两者之间，所以 $\tilde{T}_n^* \xrightarrow{p^*} \theta$。

    因为导函数连续，我们有 $g'(\tilde{T}_n^*) \xrightarrow{p^*} g'(\theta)$。

    最后，对乘积应用 Slutsky 定理，由于一部分依条件概率收敛于常数矩阵，另一部分依条件分布收敛于 $T$，即证：

    \[
    \sqrt{n}(g(T_n^*) - g(T_n)) \xrightarrow{d} g'(\theta) T
    \]

    证明完毕。 $\square$

---

## 6. 推广 Delta 方法与泛函导数

对于参数模型，参数通常仅仅是矩的函数，使用经典的 Delta 方法即可。
但是，在非参数问题中，我们感兴趣的参数往往是整个分布的泛函：

\[
\theta = \phi(P)
\]

对应的插入式 (plug-in) 估计量为：

\[
\hat{\theta}_n = \phi(\mathbb{P}_n)
\]

为了对这类的泛函应用 Delta 方法，我们需要引入**泛函导数 (generalized derivative for functionals)** 的概念。

### 6.1 Gateaux 导数 (Gateaux derivative)

Gateaux 导数类似于多变量微积分中的**方向导数 (directional derivative)**。

!!! info "定义 (Gateaux 导数)"

    如果极限：

    \[
    \lim_{t \downarrow 0} \frac{\phi((1-t)P + tH) - \phi(P)}{t} = \int \phi'_P d(H-P)
    \]

    存在，并且它关于方向 $(H-P)$ 是线性的，则称泛函 $\phi$ 在 $P$ 处沿方向 $H-P$ 是 **Gateaux 可微的**。

### 6.2 Frechet 与 Hadamard 可微性

Gateaux 可微仅仅保证了在单个确定方向上的收敛性。这在泛函分析中往往是不够的，我们需要某种一致收敛性。

* **Frechet 可微性 (Frechet differentiability)**：要求极限在**所有方向上是一致收敛**的 (uniform over all directions)。这通常要求过严。

* **Hadamard 可微性 (Hadamard differentiability)**：也被称为**紧可微 (compact differentiability)**。它要求极限在任何**紧集方向上是一致收敛**的 (uniform over compact sets of directions)。

!!! info "定义 (Hadamard 可微性)"

    设 $\phi: \mathcal{P} \to \mathbb{R}$。如果对于任意满足当 $t \to 0$ 时 $H_t \to H$ 的方向序列，都有：

    \[
    \frac{\phi(P + t H_t) - \phi(P)}{t} \to \phi'_P(H)
    \]

    则称 $\phi$ 在 $P$ 处是 Hadamard 可微的。

**为何统计学偏爱 Hadamard 可微性？**
因为根据 Prohorov 定理，经验过程 (Empirical processes) 的弱收敛性本质上与无穷维空间中的紧集 (compact sets) 紧密相关。Hadamard 可微性正好提供了在应用连续映射定理 (Continuous Mapping Theorem) 和非参数 Delta 方法时所需的恰到好处的一致性要求。


## 7. 渐近线性估计量与 Bootstrap (Asymptotic Linear Estimators)

很多在实际中广泛使用的估计量都可以被表示为（或近似为）独立同分布随机变量的均值形式。这就引出了渐近线性估计量的概念。

!!! info "定义 10.5 (渐近线性估计量 Asymptotically linear estimator)"

    估计量 $\hat{\theta}_n$ 被称为是参数 $\theta(P_0)$ 的一个**渐近线性估计量**，如果它可以被展开为如下形式：

    \[
    \hat{\theta}_n - \theta_0 = \frac{1}{n} \sum_{i=1}^n \dot{\phi}_{P_0}(X_i) + o_P(n^{-1/2})
    \]

    其中，函数 $\dot{\phi}_{P_0}(x)$ 被称为**影响函数 (Influence function)**，并且满足以下条件：

    \[
    E_{P_0} \dot{\phi}_{P_0}(X) = 0, \quad E_{P_0} \|\dot{\phi}_{P_0}(X)\|^2 < \infty
    \]

### 7.1 Z-估计量作为渐近线性估计量

Z-估计量是渐近线性估计量的一个极佳示例。回顾前面关于 Z-估计量渐近正态性的定理，如果估计量满足估计方程 $\Psi_n(\hat{\theta}_n) = 0$，我们在一定的正则条件下有：

\[
\sqrt{n}(\hat{\theta}_n - \theta_0) = \frac{1}{\sqrt{n}} \sum_{i=1}^n (-V_{\theta_0}^{-1} \psi_{\theta_0}(X_i)) + o_P(1)
\]

通过对比定义可以看出，Z-估计量是渐近线性的，并且其对应的**影响函数**为：

\[
\dot{\phi}_{P_0}(x) = -V_{\theta_0}^{-1} \psi_{\theta_0}(x)
\]

---

## 8. 渐近线性估计量的 Bootstrap 相合性

!!! success "定理 10.5 (渐近线性估计量的 Bootstrap 相合性)"

    假设 $\hat{\theta}_n$ 是 $\theta(P_0)$ 的一个渐近线性估计量，其影响函数为 $\dot{\phi}_{P_0}(x)$。设 $\hat{\theta}_n^*$ 是将相同的估计过程应用到 Bootstrap 样本 $\mathbb{P}_n^*$ 上得到的 Bootstrap 估计量。
    
    如果 $E_{P_0} \|\dot{\phi}_{P_0}(X)\|^2 < \infty$，并且在给定样本的 Bootstrap 概率测度 $P^*$ 下有如下的渐近线性展开：

    \[
    \sqrt{n}(\hat{\theta}_n^* - \hat{\theta}_n) = \frac{1}{\sqrt{n}} \sum_{i=1}^n (\dot{\phi}_{P_0}(X_i^*) - \mathbb{P}_n \dot{\phi}_{P_0}) + o_{P^*}(1)
    \]

    那么，该 Bootstrap 是相合的 (consistent)。

??? proof "定理 10.5 的证明（点击展开）"

    **第一步：转换为独立同分布均值序列**

    令 $Y_i = \dot{\phi}_{P_0}(X_i)$。由大数定律和已知条件可知：
    
    \[
    \bar{Y}_n = \mathbb{P}_n \dot{\phi}_{P_0} \xrightarrow{p} E_{P_0} \dot{\phi}_{P_0}(X) = 0
    \]
    
    在 Bootstrap 抽样中，$X_1^*, \dots, X_n^*$ 是从经验分布 $\mathbb{P}_n$ 中有放回抽取的独立同分布样本。定义相应的 Bootstrap 随机变量：
    
    \[
    Y_i^* = \dot{\phi}_{P_0}(X_i^*)
    \]
    
    给定原始数据 $X_1, \dots, X_n$，变量 $Y_1^*, \dots, Y_n^*$ 是独立同分布的，并且它们在 Bootstrap 测度下的期望和方差为：
    
    * **条件期望**：
    
    \[
    E^*(Y_1^*) = \bar{Y}_n = \mathbb{P}_n \dot{\phi}_{P_0}
    \]
    
    * **条件方差**：
    
    \[
    Var^*(Y_1^*) = \frac{1}{n} \sum_{i=1}^n (Y_i - \bar{Y}_n)^2 = \frac{1}{n} \sum_{i=1}^n Y_i^2 - (\bar{Y}_n)^2
    \]

    由于 $E_{P_0} Y_1^2 < \infty$ 且 $E_{P_0} Y_1 = 0$，根据大数定律：
    
    \[
    Var^*(Y_1^*) \xrightarrow{p} E_{P_0} Y_1^2 \quad \text{(当 } n \rightarrow \infty \text{ 时)}
    \]

    **第二步：应用 Lindeberg-Feller 中心极限定理**

    我们需要验证在给定数据下，序列 $\frac{1}{\sqrt{n}} \sum_{i=1}^n (Y_i^* - \bar{Y}_n)$ 满足 Lindeberg 条件。具体地，对于任意的 $\epsilon > 0$，我们需要证明：
    
    \[
    E^* \left[ (Y_1^* - \bar{Y}_n)^2 I(|Y_1^* - \bar{Y}_n| > \epsilon \sqrt{n}) \right] \xrightarrow{p} 0
    \]
    
    将上述期望在 Bootstrap 测度下展开：
    
    \[
    = \frac{1}{n} \sum_{i=1}^n (Y_i - \bar{Y}_n)^2 I(|Y_i - \bar{Y}_n| > \epsilon \sqrt{n})
    \]
    
    利用不等式放缩：$(a-b)^2 \le 2a^2 + 2b^2$，以及当 $|Y_i - \bar{Y}_n| > \epsilon \sqrt{n}$ 时，必有 $|Y_i| > \frac{\epsilon}{2}\sqrt{n}$ 或者 $|\bar{Y}_n| > \frac{\epsilon}{2}\sqrt{n}$，我们可以得到：
    
    \[
    \le \frac{1}{n} \sum_{i=1}^n 2(Y_i^2 + \bar{Y}_n^2) \left[ I\left(|Y_i| > \frac{\epsilon}{2}\sqrt{n}\right) + I\left(|\bar{Y}_n| > \frac{\epsilon}{2}\sqrt{n}\right) \right]
    \]
    
    因为 $E_{P_0} Y_1^2 < \infty$，由控制收敛定理 (Dominated Convergence Theorem)：
    
    \[
    E_{P_0} \left[ Y_1^2 I\left(|Y_1| > \frac{\epsilon}{2}\sqrt{n}\right) \right] \rightarrow 0
    \]
    
    并且我们知道 $\bar{Y}_n \xrightarrow{p} 0$。因此，上述放缩项依概率收敛于 0。这证明了 Lindeberg 条件成立。

    **第三步：得出结论**
    
    由 Lindeberg-Feller 中心极限定理，我们得到：
    
    \[
    \frac{1}{\sqrt{n}} \sum_{i=1}^n (Y_i^* - \bar{Y}_n) \xrightarrow{d} N(0, E_{P_0} Y_1^2) \quad \text{(在 } P^* \text{ 测度下)}
    \]
    
    结合前提条件：
    
    \[
    \sqrt{n}(\hat{\theta}_n^* - \hat{\theta}_n) = \frac{1}{\sqrt{n}} \sum_{i=1}^n (Y_i^* - \bar{Y}_n) + o_{P^*}(1)
    \]
    
    由 Slutsky 定理，Bootstrap 分布的极限与真实世界中 $\sqrt{n}(\hat{\theta}_n - \theta_0)$ 的极限完全一致，因此 Bootstrap 是相合的。 $\square$

---

## 9. Z-估计量的 Bootstrap

对于 Z-估计量，原始估计量 $\hat{\theta}_n$ 满足估计方程：

\[
\Psi_n(\theta) = \mathbb{P}_n \psi_\theta = 0 \quad (\text{或者 } = o_P(n^{-1/2}))
\]

在 Bootstrap 世界中，相应的 Bootstrap Z-估计量 $\hat{\theta}_n^*$ 则需要满足基于经验测度 $\mathbb{P}_n^*$ 的估计方程：

\[
\Psi_n^*(\theta) = \mathbb{P}_n^* \psi_\theta = 0 \quad (\text{或者 } = o_{P^*}(n^{-1/2}))
\]

!!! success "定理 10.6 (Z-估计量的 Bootstrap 相合性)"

    假设标准 Z-估计量渐近正态性的所有常规条件均成立。进一步假设 Bootstrap 估计量 $\hat{\theta}_n^*$ 满足：
    
    \[
    \sqrt{n}\mathbb{P}_n^* \psi_{\hat{\theta}_n^*} = o_{P^*}(1)
    \]
    
    并且在条件测度下依概率相合：$\hat{\theta}_n^* \xrightarrow{P} \theta_0$。
    
    那么，在 Bootstrap 测度下，类似于原始的泰勒展开，我们有如下的渐近线性表示：
    
    \[
    \sqrt{n}(\hat{\theta}_n^* - \hat{\theta}_n) = -V_{\theta_0}^{-1} \sqrt{n}(\mathbb{P}_n^* - \mathbb{P}_n)\psi_{\theta_0} + o_{P^*}(1)
    \]

**关于定理 10.6 的注记：**

* 我们需要证明 $\hat{\theta}_n^*$ 在 $P^*$ 的概率意义下相合于 $\theta_0$。
* 在实际的计算机算法实现中，我们并不一定要精确求解这个通常为非凸 (non-convex) 的优化问题。为了避免计算开销，我们经常只需从 $\hat{\theta}_n$ 出发，执行一次**单步 (one-step)** 的 Newton-Raphson 迭代即可，这样计算出的单步 Bootstrap 估计量已经能够保证渐近相合性。

---

## 10. M-估计量的 Bootstrap

对于 M-估计量，原始估计量 $\hat{\theta}_n$ 最大化了随机准则函数 $\mathbb{P}_n m_\theta$。
同样地，Bootstrap M-估计量 $\hat{\theta}_n^*$ 最大化了基于 Bootstrap 样本的准则函数 $\mathbb{P}_n^* m_\theta$。

在 M-估计量的理论中，我们在最大值点附近将目标函数展开成二次型（类似于在 Z-估计量中对导数进行一阶泰勒展开）。

!!! success "定理 10.7 (M-估计量的 Bootstrap 相合性)"

    假设标准 M-估计量渐近正态性的所有常规条件均成立，并且 $\hat{\theta}_n^*$ 在条件概率测度下依概率相合：$\hat{\theta}_n^* \xrightarrow{P} \theta_0$。
    
    那么，在 Bootstrap 测度下，我们有如下的展开：
    
    \[
    \sqrt{n}(\hat{\theta}_n^* - \hat{\theta}_n) = -V_{\theta_0}^{-1} \sqrt{n}(\mathbb{P}_n^* - \mathbb{P}_n)\dot{m}_{\theta_0} + o_{P^*}(1)
    \]
    
    这说明 M-估计量在 Bootstrap 下也是相合的，并且可以写成渐近线性的形式。

---

## 11. 利用 Bootstrap 估计标准误与偏差

通过 Bootstrap，我们可以基于计算和重抽样直接估计统计量的标准误与偏差，而无需推导其复杂的解析渐近表达式。

### 11.1 估计标准误 (Standard Errors)

计算 Bootstrap 标准误的步骤如下：

1. 从原始样本中独立地、有放回地抽取 $B$ 个 Bootstrap 样本。
2. 对每一个 Bootstrap 样本，计算其对应的估计量，记作 $\hat{\theta}_{n,1}^*, \dots, \hat{\theta}_{n,B}^*$。
3. 计算这 $B$ 个估计量的样本标准差，作为原始估计量的标准误估计 $\widehat{se}_B$：

\[
\widehat{se}_B = \sqrt{ \frac{1}{B-1} \sum_{b=1}^B (\hat{\theta}_{n,b}^* - \bar{\hat{\theta}}_n^*)^2 }
\]

其中 $\bar{\hat{\theta}}_n^* = \frac{1}{B} \sum_{b=1}^B \hat{\theta}_{n,b}^*$。

### 11.2 估计偏差 (Bias)

估计量 $\hat{\theta}_n$ 的偏差定义为：

\[
Bias = E(\hat{\theta}_n) - \theta
\]

通过在 Bootstrap 世界中用经验分布替换总体分布，我们可以用 $\bar{\hat{\theta}}_n^*$ 来替代 $E(\hat{\theta}_n)$，用 $\hat{\theta}_n$ 替代 $\theta$。因此，Bootstrap 偏差估计量为：

\[
\widehat{Bias} = E^*(\hat{\theta}_n^*) - \hat{\theta}_n \approx \bar{\hat{\theta}}_n^* - \hat{\theta}_n
\]

*(注意：在实际应用中，如果 $\widehat{Bias}$ 较小，通常不建议直接用来做偏差修正，因为这可能会增加整体的均方误差 (MSE)；只有当偏差显著大于标准误时才进行修正。)*

---

## 12. Bootstrap 置信区间 (Confidence Intervals)

Bootstrap 常被用来构造置信区间，下面介绍三种最常用的方法。

### 12.1 正态置信区间 (Normal Interval)

最简单的方法是假设 Bootstrap 分布近似于正态分布，从而直接利用估计出的标准误来构造：

\[
C_n = \left[ \hat{\theta}_n - z_{1-\alpha/2} \widehat{se}, \quad \hat{\theta}_n + z_{1-\alpha/2} \widehat{se} \right]
\]

* **优点**：计算简单。
* **缺点**：它并不是严格的（非参数的）区间，因为它强制要求分布满足对称性 (symmetry) 和正态性 (normality) 假设。

---

### 12.2 基本 Bootstrap 置信区间 (Basic Bootstrap Interval / Pivotal Method)

这种方法基于构造一个**枢轴量 (pivotal quantity)**。

令根 (Root) 为 $R_n = \hat{\theta}_n - \theta$。我们用 Bootstrap 对应的根 $R_n^* = \hat{\theta}_n^* - \hat{\theta}_n$ 的分布来近似真实世界中 $R_n$ 的分布。

设 $q_{\alpha/2}^*$ 和 $q_{1-\alpha/2}^*$ 分别是 $R_n^*$ 的 $\alpha/2$ 和 $1-\alpha/2$ 分位数 (quantiles)。由于我们用 Bootstrap 分布来近似真实分布，我们有如下近似等式：

\[
P(q_{\alpha/2}^* \le \hat{\theta}_n - \theta \le q_{1-\alpha/2}^*) \approx 1 - \alpha
\]

对不等式中间的 $\theta$ 进行移项并调整符号，我们得到 $\theta$ 的置信区间：

\[
\left[ \hat{\theta}_n - q_{1-\alpha/2}^*, \quad \hat{\theta}_n - q_{\alpha/2}^* \right]
\]

进一步，如果设 $\theta_p^*$ 为 Bootstrap 估计量 $\hat{\theta}_n^*$ 自身的 $p$-分位数，那么 $q_p^* = \theta_p^* - \hat{\theta}_n$。代入上述公式，基本 Bootstrap 置信区间可以重写为一种“翻转”的形式：

\[
C_n = \left[ 2\hat{\theta}_n - \theta_{1-\alpha/2}^*, \quad 2\hat{\theta}_n - \theta_{\alpha/2}^* \right]
\]

---

### 12.3 百分位法置信区间 (Percentile Interval)

百分位区间是最直观也是在实践中最常被使用的方法。它直接取 Bootstrap 分布 $\hat{\theta}_n^*$ 的下分位数和上分位数作为区间的端点：

\[
C_n = \left[ \theta_{\alpha/2}^*, \quad \theta_{1-\alpha/2}^* \right]
\]

**为什么百分位法是合理的？ (Justification)**

百分位法具有极其优秀的**单调变换不变性 (Transformation invariance)**。
假设存在某个未知的单调递增的变换函数 $g$，使得变换后的估计量精确服从正态分布，即：

\[
g(\hat{\theta}_n) - g(\theta) \sim N(0, c^2)
\]

由于 $g$ 是单调的，分位数的相对位置在经过单调变换后保持不变（即分位数的变换等于变换的分位数）：

\[
g(\theta_\alpha^*) = (g(\theta))_\alpha^*
\]

如果我们在变换后的空间 $g(\theta)$ 上构造置信区间，其区间恰好是 $[ (g(\theta))_{\alpha/2}^*, (g(\theta))_{1-\alpha/2}^* ]$。再通过逆变换 $g^{-1}$ 映射回到原始参数空间 $\theta$ 上的置信区间，恰巧就是：

\[
\left[ g^{-1}((g(\theta))_{\alpha/2}^*), \quad g^{-1}((g(\theta))_{1-\alpha/2}^*) \right] = \left[ \theta_{\alpha/2}^*, \quad \theta_{1-\alpha/2}^* \right]
\]

这意味着，**只要存在某种能使分布正态化的单调变换，即使我们并不知道这个变换 $g$ 具体是什么，百分位置信区间都隐式地执行了这种正态化变换**，从而给出了正确的区间界限。

## 13. 学生化 Bootstrap (Studentized Bootstrap / Bootstrap-t)

在标准的正态近似中，我们通常使用标准化或学生化 (studentized) 的统计量（即枢轴量 Pivot）：

\[
T_n = \frac{\hat{\theta}_n - \theta_0}{\widehat{se}}
\]

在 Bootstrap 世界中，我们构造对应的学生化统计量：

\[
T_n^* = \frac{\hat{\theta}_n^* - \hat{\theta}_n}{\widehat{se}^*}
\]

* 注意，这里的 $\widehat{se}^*$ 是从**当前的 Bootstrap 样本**中计算得到的标准误估计量，而不是原始样本的标准误。

### 13.1 学生化 Bootstrap 置信区间

设 $t_{\alpha/2}^*$ 和 $t_{1-\alpha/2}^*$ 分别是 Bootstrap 统计量 $T_n^*$ 的 $\alpha/2$ 和 $1-\alpha/2$ 分位数。
我们利用 Bootstrap 分布来近似真实分布，即假设事件 $t_{\alpha/2}^* \le T_n \le t_{1-\alpha/2}^*$ 的概率近似为 $1-\alpha$。

将其展开并对 $\theta_0$ 进行求解，我们可以得到学生化 Bootstrap 置信区间：

\[
C_n = \left[ \hat{\theta}_n - t_{1-\alpha/2}^* \widehat{se}, \quad \hat{\theta}_n - t_{\alpha/2}^* \widehat{se} \right]
\]

### 13.2 为什么要学生化？(二阶精确性 Second-order accuracy)

* 基本 Bootstrap 区间 (Basic Interval) 和百分位法区间 (Percentile Interval) 的覆盖误差 (coverage error) 通常是 $O(n^{-1/2})$。这与基于正态近似的区间（一阶精确 First-order accurate）处于相同的量级。

* **学生化 Bootstrap 是二阶精确的 (Second-order accurate)**。它的覆盖误差达到了 $O(n^{-1})$。

**原因：Edgeworth 展开 (Edgeworth Expansion)**

学生化统计量 $T_n$ 是一个**渐近枢轴量 (Asymptotic Pivot)**，它的极限分布 $N(0,1)$ 不依赖于任何未知的总体参数。
当我们在 Bootstrap 中重抽样枢轴量时，Bootstrap 能够自动估计并修正 Edgeworth 展开中阶数为 $n^{-1/2}$ 的**偏度修正项 (skewness correction term)**，从而消除了一阶误差，使得最终的近似误差降至 $O(n^{-1})$。

---

## 14. 双重 Bootstrap (Double Bootstrap / Iterated Bootstrap)

在学生化 Bootstrap 中，我们需要计算 $T_n^* = (\hat{\theta}_n^* - \hat{\theta}_n) / \widehat{se}^*$。然而，在某些复杂的模型中，对于每一个 Bootstrap 样本，标准误 $\widehat{se}^*$ 可能没有解析表达式或极难直接计算。

**解决方案：在 Bootstrap 内部再嵌套一层 Bootstrap！**

**具体算法步骤：**

1. 从原始经验分布 $\mathbb{P}_n$ 中有放回地抽取第一层 Bootstrap 样本 $\mathbb{P}_n^*$。利用该样本计算估计量 $\hat{\theta}_n^*$。

2. 从当前的第一层 Bootstrap 样本 $\mathbb{P}_n^*$ 中，再次有放回地抽取第二层（内部）Bootstrap 样本 $\mathbb{P}_n^{**}$，并计算 $\hat{\theta}_n^{**}$。

3. 重复步骤 2 达 $M$ 次，通过这 $M$ 个 $\hat{\theta}_n^{**}$ 的经验标准差来计算标准误 $\widehat{se}^*$。

4. 计算当前的 Bootstrap 枢轴量 $T_n^* = \frac{\hat{\theta}_n^* - \hat{\theta}_n}{\widehat{se}^*}$。

5. 重复步骤 1 到 4 达 $B$ 次，得到 $T_n^*$ 的完整经验分布，进而求得分位数以构造置信区间。

* **缺点**：计算成本极其高昂。如果外层需要 $B$ 次重抽样，内层需要 $M$ 次重抽样，总共需要进行 $B \times M$ 次模型拟合。

---

## 15. 回归模型中的 Bootstrap (Bootstrap for Regression)

考虑标准的线性或非线性回归模型：

\[
Y_i = x_i^T \beta + \epsilon_i
\]

在回归分析中，如何进行重抽样取决于我们对模型的假设程度。通常有两种主要的方法：

### 15.1 方法一：残差 Bootstrap (Residual Bootstrap / Model-based Bootstrap)

如果我们坚信模型结构的正确性，并假设误差项 $\epsilon_i$ 是独立同分布的 (i.i.d.)，且与协变量 $x_i$ 独立（同方差性 Homoscedasticity）。

**步骤：**

1. 对原始数据拟合模型，得到参数估计 $\hat{\beta}$ 以及对应的残差 $\hat{\epsilon}_i = Y_i - x_i^T \hat{\beta}$。

2. 对残差进行中心化处理（推荐）：$\tilde{\epsilon}_i = \hat{\epsilon}_i - \bar{\hat{\epsilon}}$。

3. 从中心化后的残差集 $\{\tilde{\epsilon}_1, \dots, \tilde{\epsilon}_n\}$ 中有放回地随机抽取 $n$ 个 Bootstrap 残差 $\epsilon_1^*, \dots, \epsilon_n^*$。

4. 保持自变量 $x_i$ 不变，生成新的 Bootstrap 响应变量：

\[
Y_i^* = x_i^T \hat{\beta} + \epsilon_i^*
\]

5. 使用生成的 $(x_i, Y_i^*)$ 重新拟合模型，得到 Bootstrap 估计量 $\hat{\beta}^*$。

* **特点**：充分利用了模型的参数结构，在模型设定正确且误差满足同方差假设时极其高效。

### 15.2 方法二：成对 Bootstrap (Pairs Bootstrap / Nonparametric Bootstrap)

如果我们不确定误差是否满足同方差性（即可能存在异方差 Heteroscedasticity），或者怀疑线性模型的正确性，我们可以采用更加非参数的方法。

**步骤：**

1. 将每一个观测视为一个数据对 $Z_i = (x_i, Y_i)$。

2. 直接从数据集 $\{Z_1, \dots, Z_n\}$ 中有放回地随机抽取 $n$ 个数据对，记为 $Z_1^*, \dots, Z_n^*$（其中 $Z_i^* = (x_i^*, Y_i^*)$）。

3. 对这批新的成对样本重新拟合模型，得到 Bootstrap 估计量 $\hat{\beta}^*$。

* **特点**：对异方差和模型误设具有极强的稳健性 (Robustness)。这是因为重抽样保留了 $x_i$ 和 $Y_i$ 之间的联合分布结构。

---

## 16. Bootstrap 什么时候会失效？(When does Bootstrap fail?)

尽管 Bootstrap 极其强大，但它并不是万能的。在以下几种典型场景中，标准的 Bootstrap 会失效（即无法保持相合性）：

### 16.1 1. 重尾分布 (Heavy tails / Infinite variance)

* 如果总体分布的方差无限（例如 $E[X^2] = \infty$），那么样本均值不满足经典的正态极限分布（它可能依分布收敛于某种 $\alpha$-稳定分布）。

* 在这种情况下，Bootstrap 样本均值的条件分布不会弱收敛到一个确定的非随机测度，它在极限下依然保持随机性。

### 16.2 2. 边界参数 (Boundary parameters)

* 我们在前面讨论过的 **最大值统计量** 便是典型的例子（如 $X_i \sim U(0, \theta)$ 且估计量为 $X_{(n)}$）。

* 当参数位于参数空间的边界时，经验分布 $\mathbb{P}_n$ 的离散性会在极值处引发极大概率的“点质量 (point mass)”，从而无法模拟出真实的连续极限分布。

### 16.3 3. 不可微泛函 (Non-differentiable functionals)

* 如果我们感兴趣的参数是经验测度的一个非光滑或不可微的泛函，Bootstrap 可能失效。

* **示例 1：Hodges 估计量** (超有效估计量)。这类估计量在某一点（如 0 点）由于不可微的“收缩”行为，使得 Bootstrap 无法捕捉其局部的突变性质。

* **示例 2：存在结 (ties) 时的样本中位数**。虽然对于连续分布而言，样本中位数是 Bootstrap 相合的；但如果真实分布在总体中位数附近不具备足够的光滑性（例如分布存在断点或高度离散），Bootstrap 中位数的表现将会非常不稳定。

---

## 17. 解决失效的方法：子抽样与 m-out-of-n Bootstrap

当面临上述使得标准 Bootstrap 失效的情境时，我们通常采用以下两种替代策略：

### 17.1 子抽样 (Subsampling)

* **机制**：从容量为 $n$ 的原始样本中，**无放回地 (without replacement)** 抽取容量为 $m$ 的子样本。

* **条件**：要求子样本容量 $m \rightarrow \infty$ 且 $m/n \rightarrow 0$。

* **优点**：其相合性的条件极其微弱。它只要求原始统计量在真实测度下收敛于某个非退化的极限分布即可，而不要求泛函必须可微或不在边界上。

### 17.2 m-out-of-n Bootstrap

* **机制**：从容量为 $n$ 的原始样本中，**有放回地 (with replacement)** 抽取容量为 $m$ 的 Bootstrap 样本。

* **条件**：同样要求 $m \rightarrow \infty$ 且 $m/n \rightarrow 0$。

* **优点**：通过缩小重抽样的规模 $m$，它能够打破因为重抽样 $n$ 次所带来的过多“结 (ties)”和极值处的聚集问题，从而恢复在边界参数和某些重尾情况下的相合性。

---

## 18. 本章最终总结 (Final Summary)

1. **一阶推断**：渐近正态性理论为我们提供了首要的一阶推断工具 (Asymptotic normality gives first-order inference)。

2. **有限样本误差**：Edgeworth 展开在理论上解释了主要的有限样本误差来源 (Edgeworth expansion explains the leading finite-sample error)。

3. **自动修正**：Bootstrap 的强大之处在于，它通常能够利用数据自动估计出 Edgeworth 展开中的高阶修正项 (Bootstrap often estimates Edgeworth correction terms automatically)。

4. **二阶精确性**：**学生化 Bootstrap (Studentized bootstrap)** 是获得具有二阶精确度的置信区间的主要且推荐的途径 (Studentized bootstrap is the main route to second-order accurate intervals)。

5. **理论一致性**：对于 M-估计量、Z-估计量和极大似然估计量 (MLE)，Bootstrap 的相合性可以通过与证明渐近正态性完全相同的线性化 (linearization) 论证来推导。

6. **局限性**：Bootstrap 是一把极其强大的利器，但它并**不是普遍有效**的 (Bootstrap is powerful, but it is not universally valid)。使用前需注意前提条件。

!!! success "一句话核心总结 (One-line summary)"

    Bootstrap 的精髓在于用“平行宇宙”来模仿真实世界：真实的波动 $\hat{\theta}_n - \theta_0$ 被 Bootstrap 世界的 $\hat{\theta}_n^* - \hat{\theta}_n$ 完美模仿；而那不可见的神明——总体分布 $P$，则被我们手中的数据——经验分布 $\mathbb{P}_n$ 所替代。
