# 第七章：Delta 方法

在前面的章节中，我们研究了随机变量序列（如样本均值）本身的收敛性。本章我们将探讨一个在统计学中更为核心的问题：**在平滑变换下，依分布收敛是否能够被保持？** 这就是著名的 **Delta 方法 (Delta Method)** 及其在各种统计推断（如置信区间、假设检验、方差稳定化）中的广泛应用。

假设我们有一系列参数 $\theta \in \mathbb{R}^k$ 的估计量 $\{T_n\}_{n \ge 1}$（取值于 $\mathbb{R}^k$）。

* 对于我们感兴趣的参数函数 $\phi(\theta)$（其中 $\phi: \mathbb{R}^k \rightarrow \mathbb{R}^m$），根据连续映射定理 (Continuous Mapping Theorem)，如果 $T_n \xrightarrow{p} \theta$ 且 $\phi$ 在 $\theta$ 处连续，那么 $\phi(T_n) \xrightarrow{p} \phi(\theta)$。

* 但在统计学中，一个更有趣且更实用的问题是：如果已知 $\sqrt{n}(T_n - \theta) \Rightarrow T$，那么经过非线性变换后的 $\sqrt{n}(\phi(T_n) - \phi(\theta))$ 是否也会收敛到一个确定的分布？

---

## 1. 向量值函数的导数与一阶 Delta 方法

回顾多变量微积分，如果函数 $\phi(\cdot)$ 在 $\theta$ 处可导，意味着存在一个线性映射（矩阵）$\phi'_{\theta}: \mathbb{R}^k \mapsto \mathbb{R}^m$，使得：

\[
\phi(\theta + h) - \phi(\theta) = \phi'(\theta)h + R(h)
\]

其中残差项满足 $R(h) = o(\|h\|)$ 当 $h \rightarrow 0$。

这个导数映射（雅可比矩阵 Jacobian Matrix）具体形式为：

\[
\phi'_{\theta} = \begin{pmatrix} 
\frac{\partial \phi_1}{\partial \theta_1}(\theta) & \cdots & \frac{\partial \phi_1}{\partial \theta_k}(\theta) \\ 
\vdots & \ddots & \vdots \\ 
\frac{\partial \phi_m}{\partial \theta_1}(\theta) & \cdots & \frac{\partial \phi_m}{\partial \theta_k}(\theta) 
\end{pmatrix}_{m \times k}
\]

*(注：如果 $m=1, k>1$，该导数映射即为函数的**梯度 (Gradient)**。)*

### 1.1 一阶 Delta 方法 (First Order Delta Method)

!!! success "定理 5.1 (一阶 Delta 方法)"

    如果 $\phi$ 在 $\theta$ 处可导，且导数矩阵 $\phi'(\theta) \ne 0$。假设存在一个确定的发散数列 $\{r_n\}$（通常 $r_n = \sqrt{n}$）满足 $r_n \rightarrow \infty$，且 $r_n(T_n - \theta) \Rightarrow T$，那么：

    (i) $r_n(\phi(T_n) - \phi(\theta)) - \phi'(\theta)(r_n(T_n - \theta)) \xrightarrow{p} 0$

    (ii) $r_n(\phi(T_n) - \phi(\theta)) \Rightarrow \phi'(\theta)T$

??? proof "定理 5.1 的证明（点击展开）"

    **证明 (i)：**
    
    已知 $r_n(T_n - \theta) \Rightarrow T$。由于 $r_n \rightarrow \infty$，由随机有界性 (Stochastic Boundedness, $O_p(1)$)，必然有：

    \[
    T_n - \theta \xrightarrow{p} 0
    \]

    利用 $\phi$ 在 $\theta$ 处的泰勒展开（可导性），对于充分小的 $h = T_n - \theta$：

    \[
    \phi(T_n) - \phi(\theta) - \phi'(\theta)(T_n - \theta) = o(\|T_n - \theta\|)
    \]

    两边同乘 $r_n$：

    \[
    r_n [ \phi(T_n) - \phi(\theta) - \phi'(\theta)(T_n - \theta) ] = r_n \cdot o_p(\|T_n - \theta\|)
    \]

    将右侧改写为：

    \[
    o_p(1) \cdot r_n \|T_n - \theta\|
    \]

    因为 $r_n(T_n - \theta) = O_p(1)$，所以 $o_p(1) \cdot O_p(1) = o_p(1)$。这就证明了结论 (i)。

    **证明 (ii)：**
    
    将结论 (i) 移项得到：

    \[
    r_n(\phi(T_n) - \phi(\theta)) = \phi'(\theta) r_n(T_n - \theta) + o_p(1)
    \]

    由于 $r_n(T_n - \theta) \Rightarrow T$，根据连续映射定理，$\phi'(\theta) r_n(T_n - \theta) \Rightarrow \phi'(\theta)T$。
    最后，应用 **Slutsky 定理**（加上一个依概率收敛于 0 的项不改变分布收敛），直接得到：

    \[
    r_n(\phi(T_n) - \phi(\theta)) \Rightarrow \phi'(\theta)T
    \]

    证明完毕。 $\square$

**典型应用：正态 Delta 方法**

如果估计量满足渐近正态性：$\sqrt{n}(T_n - \theta) \xrightarrow{d} N(0, \sigma^2(\theta))$。
对于任意在 $\theta$ 处可导且导数 $g'(\theta) \ne 0$ 的标量函数 $g: \mathbb{R} \rightarrow \mathbb{R}$，有：

\[
\sqrt{n}[g(T_n) - g(\theta)] \xrightarrow{d} N(0, [g'(\theta)]^2 \sigma^2(\theta))
\]

---

## 2. 高阶 Delta 方法 (High Order Delta Method)

一阶 Delta 方法极其依赖于 $\phi'(\theta) \ne 0$。如果遇到 $\phi'(\theta) = 0$ 但 $\phi''(\theta) \ne 0$ 的退化情况，一阶方法就会失效（得到退化的点质量分布）。此时我们需要引入高阶泰勒展开。

展开到二阶项：

\[
\phi(T_n) = \phi(\theta) + \frac{1}{2}\phi''(\theta)(T_n - \theta)^2 + \cdots
\]

同乘 $n$（注意这里是 $n$ 而不是 $\sqrt{n}$，因为平方项的存在）：

\[
n(\phi(T_n) - \phi(\theta)) = \frac{1}{2}\phi''(\theta)[\sqrt{n}(T_n - \theta)]^2 \Rightarrow \frac{1}{2}\phi''(\theta)T^2
\]

!!! success "定理 5.2 (高阶 Delta 方法)"

    假设单变量函数 $\phi$ 在 $\theta$ 处 $m$ 次可导，且满足 $\phi^{(m)}(\theta) \ne 0$ 但前面所有的低阶导数均为零（即 $\phi^{(j)}(\theta) = 0, \forall j < m$）。如果 $r_n(T_n - \theta) \Rightarrow T$，那么：

    \[
    \frac{r_n^m (\phi(T_n) - \phi(\theta))}{\frac{1}{m!} \phi^{(m)}(\theta)} \Rightarrow T^m
    \]

### 2.1 高阶 Delta 方法应用示例

假设 $X_1, \dots, X_n$ 是 i.i.d. 序列，均值为 $\mu$，方差已知为 $\sigma^2$。我们要检验原假设 $H_0: \mu = 0$。
在原假设下，统计量 $n\bar{X}_n^2 / \sigma^2 \rightarrow [N(0,1)]^2 = \chi_1^2$。

现在考虑随机变量 $\cos(\bar{X}_n)$ 的极限行为：

* **如果强行使用一阶 Delta 方法**：由于函数 $g(x) = \cos(x)$ 在 $x=0$ 处的导数 $g'(0) = -\sin(0) = 0$，标准化项 $\sqrt{n}$ 会导致：
  
    \[
    \sqrt{n}(\cos(\bar{X}_n) - 1) \xrightarrow{p} 0
    \]
  
  这没有提供任何有用的分布信息，说明 $\sqrt{n}$ 不是正确的收敛速率。

* **使用二阶 Delta 方法**：因为在 $x=0$ 处，二阶导数 $\cos''(0) = -\cos(0) = -1 \ne 0$。展开得：
  
    \[
    \cos(\bar{X}_n) - \cos(0) = (\bar{X}_n - 0) \cdot 0 + \frac{1}{2}(\bar{X}_n - 0)^2 \cdot (-1) + o_p(\bar{X}_n^2)
    \]
  
  同乘 $-2n$：$-2n(\cos(\bar{X}_n) - 1) = n\bar{X}_n^2 + o_p(1) \Rightarrow \chi_1^2 \cdot \sigma^2$
  
  这给出了正确的非退化极限分布。

---

## 3. 渐近正态性与 Delta 方法的经典应用

### 3.1 样本方差与标准差的极限分布

设 $X_1, \dots, X_n \sim i.i.d. F$，具有有限的 4 阶矩。记总体中心矩 $\alpha_i = E(X_1^i)$，样本矩 $m_{ni} = n^{-1}\sum X_j^i$。
样本方差可以写为两个样本矩的函数：

\[
S_n = n^{-1}\sum_{i=1}^n (X_i - \bar{X})^2 = m_{n2} - m_{n1}^2 = \phi(m_{n1}, m_{n2})
\]

其中非线性变换函数为 $\phi(x_1, x_2) = x_2 - x_1^2$。其梯度向量为：

\[
\phi'(\alpha_1, \alpha_2) = (-2\alpha_1, 1)
\]

由多维中心极限定理 (Multivariate CLT)：

\[
\sqrt{n} \left[ \begin{pmatrix} m_{n1} \\ m_{n2} \end{pmatrix} - \begin{pmatrix} \alpha_1 \\ \alpha_2 \end{pmatrix} \right] \xrightarrow{d} N\left( 0, Var \begin{pmatrix} X_1 \\ X_1^2 \end{pmatrix} \right)
\]

应用多元一阶 Delta 方法，样本方差的极限分布为：

\[
\sqrt{n}(S_n - \sigma^2) \xrightarrow{d} N\left( 0, (-2\alpha_1, 1) Var \begin{pmatrix} X_1 \\ X_1^2 \end{pmatrix} \begin{pmatrix} -2\alpha_1 \\ 1 \end{pmatrix} \right)
\]

通过展开二次型，可以巧妙地化简为中心矩的形式：$E(X_1 - \alpha_1)^4 - [E(X_1 - \alpha_1)^2]^2 = c_4 - c_2^2$（即第四中心矩减去方差的平方）。
因此：

\[
\sqrt{n}(S_n - \sigma^2) \xrightarrow{d} N(0, c_4 - c_2^2)
\]

**推论：无偏样本方差与样本标准差**

对于无偏方差 $S_{n-1} = \frac{n}{n-1}S_n$，由于相差的常数在极限下趋于 1，且差异项 $\sqrt{n}(\frac{n}{n-1} - 1)S_n = o_p(1)$，它具有**相同的极限分布**。

对于**样本标准差** $S_n^{1/2}$，应用单变量 Delta 方法，取 $\phi(x) = \sqrt{x}$，导数为 $\phi'(x) = \frac{1}{2}x^{-1/2}$。代入 $\sigma^2$ 处的值：

\[
\sqrt{n}(S_n^{1/2} - \sigma) \xrightarrow{d} N\left( 0, \frac{c_4 - c_2^2}{4\sigma^2} \right)
\]

### 3.2 更多常见变换的例子

假设基础序列 $X_n$ 满足渐近正态性 (Asymptotically Normal, $AN$)：$X_n \sim AN(\mu, \sigma_n^2)$ 且 $\sigma_n \rightarrow 0$。

* (i) $X_n^2 \sim AN(\mu^2, 4\mu^2 \sigma_n^2)$ （要求 $\mu \ne 0$）

* (ii) $\frac{1}{X_n} \sim AN(\mu^{-1}, \frac{\sigma_n^2}{\mu^4})$ （要求 $\mu \ne 0$）

* (iii) $e^{X_n} \sim AN(e^\mu, e^{2\mu} \sigma_n^2)$ （对于任意 $\mu$）

* (iv) $\log|X_n| \sim AN(\log|\mu|, \mu^{-2} \sigma_n^2)$ （要求 $\mu \ne 0$）。如果 $\mu = 0$ 且 $\sigma_n = 1/\sqrt{n}$，则由连续映射定理，极限分布与 $\log|N(0,1)|$ 有关。

**多维二次型的权重 $\chi^2$ 分布：**

设 $X_1, \dots, X_n \sim i.i.d.\ F$ 于 $\mathbb{R}^p$ 空间，均值为 $\mu$，协方差为 $\Sigma$。考察目标 $\hat{\theta} = \bar{X}^T \bar{X}$。

* **若 $\mu \ne 0$**：应用一阶 Delta 方法，$\phi'(\mu) = 2\mu^T$，有 $\sqrt{n}(\bar{X}^T\bar{X} - \mu^T\mu) \xrightarrow{d} N(0, 4\mu^T \Sigma \mu)$。

* **若 $\mu = 0$**：一阶导数为 0（因为 $\mu^T \Sigma \mu = 0$）。此时需要使用高阶映射。因为 $\sqrt{n}\bar{X} \xrightarrow{d} N_p(0, \Sigma)$，所以：
  
    \[
    n\bar{X}^T\bar{X} \Rightarrow N_p^T(0, \Sigma) N_p(0, \Sigma) \stackrel{d}{=} Z^T \Sigma^{1/2} \Sigma^{1/2} Z = Z^T \Sigma Z
    \]
  
  其中 $Z \sim N_p(0, I_p)$。通过特征值分解 $\Sigma = U^T \text{diag}(\lambda_1, \dots, \lambda_p) U$，上式等价于线性组合：$\sum_{i=1}^p \lambda_i \chi_{1i}^2$
  
  这是一个**加权 $\chi^2$ 分布 (Weighted $\chi^2$ distribution)**。

---

## 4. 假设检验中的渐近理论

### 4.1 方差的 $\chi^2$ 检验与超额峰度的影响

设 $X_1, \dots, X_n \sim i.i.d.\ F$，$EX_1^4 < \infty$。我们要检验 $H_0: \sigma^2 \le 1$ VS $H_1: \sigma^2 > 1$。
在正态假设下，检验统计量为 $nS_n$，拒绝域为 $nS_n > \chi^2_{n-1, \alpha}$。检验的 size 恰好为 $\alpha$。

然而，如果数据分布 $F$ **不是正态分布**，存在**超额峰度 (Excessive Kurtosis)** $\kappa = \frac{E(X-\mu)^4}{\sigma^4} - 3 \ne 0$ 时，情况会发生根本改变。

已知对于标准正态变量之和构成的卡方分布，当 $n$ 很大时：

\[
\frac{\chi^2_{n-1} - (n-1)}{\sqrt{2(n-1)}} \xrightarrow{d} N(0, 1)
\]

而真实样本方差的渐近分布（由 3.1 节已知）：

\[
\sqrt{n}\left( \frac{S_n}{\sigma^2} - 1 \right) \xrightarrow{d} N(0, \kappa + 2) \ne N(0, 2)
\]

**检验的实际 Size（第一类错误率）：**

利用卡方临界值的渐近展开 $\chi^2_{n-1, \alpha} \approx (n-1) + Z_\alpha \sqrt{2(n-1)}$，当真实方差位于边界 $\sigma^2 = 1$ 时：

\[
P_{\sigma^2=1}(nS_n > \chi^2_{n-1, \alpha}) \approx P\left( \sqrt{n}(S_n - 1) > \frac{Z_\alpha \sqrt{2n}}{\sqrt{n}} \right) \rightarrow P(N(0, \kappa+2) > \sqrt{2}Z_\alpha)
\]

标准化后：

\[
= 1 - \Phi\left( \frac{\sqrt{2} Z_\alpha}{\sqrt{\kappa + 2}} \right)
\]

* **结论**：对于具有厚尾特征 ($\kappa > 0$) 的分布，真实的 Size 会严格大于名义的 $\alpha$。这就解释了为什么在非正态数据下，传统的方差卡方检验会产生过多的假阳性。

### 4.2 多项分布向量与 Pearson $\chi^2$ 统计量

考虑多项分布 $(n_1, \dots, n_K) \sim Multinomial(n; p_1, \dots, p_K)$。
定义标准化频率向量 $X_n = \sqrt{n}(\frac{n_1}{n} - p_1, \dots, \frac{n_K}{n} - p_K)^T \xrightarrow{d} N(0, \Sigma)$。
其中协方差矩阵 $\Sigma$ 元素为 $\sigma_{ii} = p_i(1-p_i)$ 且 $\sigma_{ij} = -p_i p_j$。

拟合优度 (Goodness-of-fit) 的 Pearson $\chi^2$ 统计量可以写为二次型：

\[
T_n = \sum_{i=1}^K \frac{(n_i - np_i)^2}{np_i} = X_n^T C X_n
\]

其中 $C = \text{diag}(p_1^{-1}, \dots, p_K^{-1})$。

由映射定理，极限分布为二次型 $Z^T \Sigma^{1/2} C \Sigma^{1/2} Z$。
可以证明矩阵 $A = \Sigma^{1/2} C \Sigma^{1/2}$ 是一个**幂等矩阵 (Idempotent matrix, $A^2 = A$)**。
幂等矩阵二次型服从卡方分布，自由度为其迹 (Trace)：

\[
\text{tr}(\Sigma^{1/2} C \Sigma^{1/2}) = \text{tr}(C \Sigma) = \sum_{i=1}^K p_i^{-1} p_i(1-p_i) = K - 1
\]

因此，Pearson 统计量 $T_n \Rightarrow \chi^2_{K-1}$。

### 4.3 Wald 检验 (Wald Test)

对于多维假设检验 $H_0: \mu = \mu_0$ VS $H_1: \mu \ne \mu_0$，常用的 Wald 统计量为：

\[
W_n = n(\bar{X} - \mu_0)^T S_n^{-1} (\bar{X} - \mu_0)
\]

由大数定律，样本协方差矩阵 $S_n \xrightarrow{p} \Sigma$，故 $S_n^{-1} \xrightarrow{p} \Sigma^{-1}$。
通过插入法 (Plug-in) 和渐近展开：

\[
W_n = \sqrt{n}(\bar{X} - \mu_0)^T \Sigma^{-1} \sqrt{n}(\bar{X} - \mu_0) + \sqrt{n}(\bar{X} - \mu_0)^T (S_n^{-1} - \Sigma^{-1}) \sqrt{n}(\bar{X} - \mu_0)
\]

由于 $\sqrt{n}(\bar{X} - \mu_0) = O_p(1)$ 且 $S_n^{-1} - \Sigma^{-1} = o_p(1)$，第二项为 $o_p(1)$。
第一项即为标准的多元正态二次型，故：

\[
W_n \xrightarrow{d} \chi^2_p
\]

---

## 5. 方差稳定变换 (Variance Stabilizing Transform, VST)

在使用渐近正态性构造置信区间时：

\[
T_n \pm Z_{1-\alpha/2} \frac{\sigma(\hat{\theta})}{\sqrt{n}}
\]

我们发现区间的宽度会随着未知参数 $\theta$（体现在 $\sigma(\theta)$ 中）的变化而剧烈波动。
**方差稳定变换 (VST)** 的目的是寻找一个平滑变换 $\phi(\cdot)$，使得变换后的极限方差不再依赖于参数 $\theta$：

\[
\sqrt{n}(\phi(T_n) - \phi(\theta)) \xrightarrow{d} N(0, c^2)
\]

其中 $c > 0$ 是一个常数。

由一阶 Delta 方法已知，变换后的方差为 $(\phi'(\theta))^2 \sigma^2(\theta)$。令其等于常数 $c^2$：

\[
\phi'(\theta) \sigma(\theta) = c \implies \phi'(\theta) = \frac{c}{\sigma(\theta)}
\]

对两边积分，我们得到了 VST 的核心构造公式：

\[
\phi(\theta) = \int \frac{d\theta}{\sigma(\theta)}
\]

### 5.1 VST 应用：Tukey's Hanging Rootgram

在非参数核密度估计 (Kernel Density Estimator, KDE) 中：

\[
\hat{f}_{nh}(x) = \frac{1}{nh} \sum_{i=1}^n K\left( \frac{x - X_i}{h} \right)
\]

已知在适当的带宽条件下：

\[
\sqrt{nh}(\hat{f}_{nh}(x) - f(x)) \Rightarrow N(0, f(x))
\]

即原始估计量的方差正比于密度函数本身 $f(x)$。为了稳定方差以进行统一的误差带绘制，我们应用 VST。这里方差项 $\sigma^2(f) = f$。代入构造公式：

\[
\phi(f) = \int \frac{df}{\sqrt{f}} = f^{1/2}
\]

（忽略积分常数和倍数）。因此，我们对密度估计量开平方（"Root-gram"）：

\[
\hat{f}_{nh}^{1/2}(x) \sim AN\left( f^{1/2}(x), \frac{1}{4nh} \right)
\]

此时，渐近方差仅与样本量和带宽有关，完美消除了对密度值 $f(x)$ 的依赖。

---

## 6. 一致可积性与矩的渐近逼近

Delta 方法不仅可以研究分布的收敛，还能用于近似估计量的期望和方差。但这需要一个连接“依分布收敛”与“矩收敛”的桥梁——**一致可积性**。

!!! info "定义 5.3 (渐近一致可积性 Asymptotic Uniformly Integrable, u.i.)"

    序列 $\{Y_n\}_{n \ge 0}$ 被称为渐近一致可积的，如果满足：

    \[
    \lim_{M \rightarrow \infty} \limsup_{n \rightarrow \infty} E[|Y_n| \mathbb{I}_{\{|Y_n| > M\}}] = 0
    \]

一致可积性是确保期望取极限操作合法的关键。

!!! success "定理 5.4"

    设 $f: \mathbb{R}^k \rightarrow \mathbb{R}$ 在集合 $C$ 上处处连续可测。若 $X_n \xrightarrow{d} X$ 且 $X$ 取值于 $C$。那么：
    
    \[
    E[f(X_n)] \rightarrow E[f(X)] \quad \text{当且仅当序列 } f(X_n) \text{ 是渐近 u.i. 的。}
    \]

**矩的泰勒逼近 (Moment Approximation)**

如果我们想利用二阶泰勒展开来近似 $E[\phi(T_n)]$ 和 $Var(\phi(T_n))$：

\[
\phi(T_n) = \phi(\theta) + \phi'(\theta)(T_n - \theta) + \frac{1}{2}\phi''(\theta)(T_n - \theta)^2 + \cdots
\]

取期望和方差后，我们*期望*得到：

* $E[\phi(T_n)] \approx \phi(\theta) + \phi'(\theta)\text{Bias}(T_n) + \frac{1}{2}\phi''(\theta)\text{MSE}(T_n)$
* $Var(\phi(T_n)) \approx [\phi'(\theta)]^T Var(T_n) [\phi'(\theta)]$

**合法的严谨前提**：为了让上述约等号严格成立，我们必须确保残差项的期望收敛。这要求随机序列 $\phi(T_n) - \phi(\theta)$ 必须是**一致可积 (u.i.)** 的。通常，如果基础偏差 $T_n - \theta$ 是一致可积的，并且函数 $\phi$ 满足 Lipschitz 连续条件，那么变换后的序列也是一致可积的，从而确保了矩逼近的渐近有效性。