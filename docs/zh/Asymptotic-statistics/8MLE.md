# 第八章：极大似然估计（一）

极大似然估计 (Maximum Likelihood Estimator, MLE) 是现代统计学中最核心的参数估计方法之一。本章我们将系统性地探讨极大似然估计的渐近理论，包括 Fisher 信息量、Cramér-Rao 下界、极大似然估计的相合性 (Consistency) 以及它的渐近正态性 (Asymptotic Normality)。最后，我们将引入基于似然的假设检验方法，特别是似然比检验 (Likelihood Ratio Test)。

假设 $X = \{X_1, ..., X_n\}$ 是独立同分布 (i.i.d.) 的样本，其分布 $F_\theta$ 属于某参数族 $\mathcal{F} = \{F_\theta : \theta = (\theta_1, ..., \theta_k)^T \in \Theta\}$，且假设分布 $F_\theta$ 具有概率密度函数 $f_\theta(x)$。样本 $X$ 的**似然函数 (Likelihood function)** 定义为：

\[
L(\theta; X) = \prod_{i=1}^n f_\theta(X_i)
\]

**极大似然估计 (MLE)** $\hat{\theta}$ 就是在参数空间中最大化似然函数的点：

\[
\hat{\theta} = \arg \max_{\theta \in \Theta} \log L(\theta; X)
\]

通常，MLE 可以通过求解**得分方程 (Score equations)** 获得：

\[
\frac{\partial \log L(\theta; X)}{\partial \theta_j} \Bigg|_{\theta=\hat{\theta}} = 0, \quad j=1, 2, \dots, k
\]

得分函数的方差对于 MLE 的渐近正态性起着决定性的作用。

---

## 1. Fisher 信息量 (Fisher Information)

!!! info "定义 8.1 (Fisher 信息量正则性 FI Regularity)"

    假设分布族 $\mathcal{P} = \{P_\theta, \theta \in \Theta\}$ 被一个 $\sigma$-有限测度 $\mu$ 支配 (dominated)。如果存在 $\theta$ 的一个开邻域 $\Theta_\theta$，使得以下条件满足，则称 $\mathcal{P}$ 在 $\theta \in \Theta$ 处是 **Fisher 信息量正则的 (FI regular)**：

    (i) $f_\theta(x) := \frac{dP_\theta(x)}{d\mu} > 0$ 对于任意的 $x$ 和 $\theta \in \Theta_\theta$ 成立。

    (ii) 对于任意 $x$，$f_\theta(x)$ 在 $\theta$ 处可导。

    (iii) 积分 $\int f_\theta(x) \mu(dx)$ 可以在积分号下关于 $\theta$ 求导，即：

    \[
    \int \frac{d}{d\theta'} f_{\theta'}(x) \Big|_{\theta'=\theta} \mu(dx) = 0
    \]

!!! info "定义 8.2 (Fisher 信息量)"

    如果模型 $\mathcal{P} = \{P_\theta, \theta \in \Theta\}$ 是 FI 正则的，那么：

    \[
    I_1(\theta) = E_\theta \left[ \left( \frac{d}{d\theta'} \log f_{\theta'}(X) \Big|_{\theta'=\theta} \right)^2 \right]
    \]

    被称为基于单个观测 $X$ 在 $\theta$ 处的 **Fisher 信息量 (Fisher Information)**。

**关于 Fisher 信息量的重要性质：**

1.  **方差表示**：由定义中积分号下可导的性质可知得分函数的期望为 0，即 $E_\theta[\frac{d}{d\theta'} \log f_{\theta'}(X)|_{\theta'=\theta}] = 0$。因此 Fisher 信息量就是得分函数的方差：

    \[
    I_1(\theta) = var\left( \frac{d}{d\theta'} \log f_{\theta'}(X) \Big|_{\theta'=\theta} \right)
    \]

2.  **多元情形**：如果参数空间 $\Theta \subset \mathbb{R}^K$ ($K > 1$)，其中 $\theta = (\theta_1, \dots, \theta_K)^T$，则得分函数是一个 $K$ 维向量 $\nabla_\theta \log f_{\theta'}(x)$。此时 Fisher 信息量是一个矩阵 (FI matrix)：

    \[
    I(\theta) = E_\theta \left[ \left( \nabla_\theta \log f_{\theta'}(X) \right) \left( \nabla_\theta \log f_{\theta'}(X) \right)^T \Big|_{\theta'=\theta} \right]
    \]

3.  **二阶导数形式**：如果 $\mathcal{P}$ 是 FI 正则的，且对任意 $x$，$f_\theta(x)$ 在 $\theta$ 处**二阶可导**，且恒等式 $1 = \int f_\theta(x)\mu(dx)$ 可以在积分号下求二阶导数（即 $\int \frac{d^2}{d\theta'^2} f_{\theta'}(x)|_{\theta'=\theta} \mu(dx) = 0$），则有：

    \[
    I(\theta) = -E_\theta \left[ \frac{d^2}{d\theta'^2} \log f_{\theta'}(X) \Big|_{\theta'=\theta} \right]
    \]

---

## 2. Cramér-Rao 下界与 Bhattacharya 不等式

### 2.1 Cramér-Rao (C-R) 下界

!!! success "定理 8.1 (Cramér-Rao 下界)"

    设 $(\mathbb{X}, \mathcal{X}, \mathcal{P}=\{P_\theta, \theta \in \Theta\})$ 是随机变量 $X$ 的概率空间，$\mathcal{P}$ 被 $\sigma$-有限测度 $\mu$ 绝对连续支配，密度为 $f_\theta(x) = \frac{dP_\theta}{d\mu}$。假设满足以下条件：
    
    (i) $\Theta \subset \mathbb{R}$ 是开集。
    (ii) $A = \text{support}(f_\theta)$ 不依赖于 $\theta$。
    (iii) $\forall \theta \in \Theta$，偏导数 $\frac{df_\theta(x)}{d\theta}$ 存在。
    (iv) 对任意 $\theta \in \Theta$，$E_\theta[\frac{\partial}{\partial\theta}\log f_\theta(X)] = \int \frac{\partial f_\theta(x)}{\partial\theta}\mu(dx) = 0$。
    (v) $I(\theta) > 0$ 对任意 $\theta \in \Theta$ 成立。
    (vi) 目标函数 $g: \Theta \rightarrow \mathbb{R}$ 存在导数 $\frac{dg(\theta)}{d\theta}$。且存在估计量 $\hat{g}: \mathbb{X} \rightarrow \Theta$，使得 $\hat{g}(X)$ 是 $g(\theta)$ 的**无偏估计 (Unbiased estimator)**。
    (vii) 积分与求导可交换：$\frac{d}{d\theta}\int \hat{g}(x)f_\theta(x)\mu(dx) = \int \hat{g}(x)\frac{df_\theta(x)}{d\theta}\mu(dx)$。

    那么，无偏估计量 $\hat{g}(X)$ 的方差满足 Cramér-Rao 下界：

    \[
    var_\theta(\hat{g}(X)) \ge \frac{[g'(\theta)]^2}{I_n(\theta)}
    \]

    对于多变量情况，有 $var_\theta(\hat{g}(X)) \ge [g'(\theta)]^T I_n^{-1}(\theta) [g'(\theta)]$。

在上述定理中，条件 (iv) 和 (vii) （关于求导和积分号交换）是最具限制性的。下面的引理给出了一组确保其成立的充分条件（控制收敛定理的应用）。

!!! info "引理 8.1 (可交换性的充分条件)"

    在满足前述定理条件 (i)-(iii) 的基础上，如果存在一个包络函数 $G: \mathbb{X} \times \Theta \rightarrow \mathbb{R}$ 使得：
    
    (a) 对任意 $\theta \in \Theta$，$G(x, \theta)$ 是 $\mathcal{X}$-可测的。
    (b) 对任意 $\theta \in \Theta$，$E_\theta G^2(X, \theta) < \infty$。
    (c) 对任意 $\theta \in \Theta$，存在 $\epsilon_\theta > 0$，使得当 $|\theta - \theta'| < \epsilon_\theta$ 且 $x \in A$ 时有：
    
    \[
    \left| \frac{df_{\theta'}(x)}{d\theta'} \right| \le G(x, \theta)f_\theta(x)
    \]

    那么条件 (iv) 必定满足；并且对于任何方差有限（$E_\theta(\hat{g}(X))^2 < \infty$）的无偏估计量 $\hat{g}(X)$，条件 (vii) 也同样成立。

??? proof "引理 8.1 的证明（点击展开）"

    **证明利用均值定理 (MVT) 和控制收敛定理 (DCT)：**

    对于任意 $\theta \in \Theta$ 以及 $\theta' \in \Theta$ 满足 $|\theta - \theta'| < \epsilon_\theta$：
    由于 $\int_{\mathcal{X}} f_\theta(x)\mu(dx) = \int_{\mathcal{X}} f_{\theta'}(x)\mu(dx) = 1$，有：

    \[
    \int_{\mathcal{X}} \frac{f_\theta(x) - f_{\theta'}(x)}{\theta - \theta'} \mu(dx) = 0
    \]

    由微分中值定理 (MVT) 和条件 (c)，存在介于 $\theta$ 与 $\theta'$ 之间的 $\tilde{\theta}$ 使得：

    \[
    \left| \frac{f_\theta(x) - f_{\theta'}(x)}{\theta - \theta'} \right| = \left| \frac{df_{\tilde{\theta}}(x)}{d\tilde{\theta}} \right| \le G(x, \theta)f_\theta(x)
    \]

    注意到控制函数是 $\mu$-可积的（利用 Cauchy-Schwarz）：

    \[
    \int_{\mathcal{X}} G(x, \theta)f_\theta(x)\mu(dx) = E_\theta G(X, \theta) \le E_\theta^{1/2}[G^2(X, \theta)] < \infty
    \]

    因此，利用控制收敛定理 (DCT) 取极限 $\theta' \rightarrow \theta$：

    \[
    \int_{\mathcal{X}} \frac{df_\theta(x)}{d\theta}\mu(dx) = \lim_{\theta' \rightarrow \theta} \int_{\mathcal{X}} \frac{f_\theta(x) - f_{\theta'}(x)}{\theta - \theta'} \mu(dx) = 0
    \]

    这就证明了条件 (iv)。

    另一方面，设 $\hat{g}(X)$ 是 $g(\theta)$ 的无偏估计且 $E_\theta \hat{g}^2(X) < \infty$。考虑增量：

    \[
    \int_{\mathcal{X}} \hat{g}(x)\frac{f_\theta(x) - f_{\theta'}(x)}{\theta - \theta'}\mu(dx) = \frac{g(\theta) - g(\theta')}{\theta - \theta'}
    \]

    利用相同的包络界限：

    \[
    \left| \hat{g}(x)\frac{f_\theta(x) - f_{\theta'}(x)}{\theta - \theta'} \right| \le |\hat{g}(x)|G(x, \theta)f_\theta(x)
    \]

    验证新的控制函数是否可积：

    \[
    \int_{\mathcal{X}} |\hat{g}(x)|G(x, \theta)f_\theta(x)\mu(dx) = E_\theta[|\hat{g}(X)|G(X, \theta)] \le [E_\theta \hat{g}^2(X) \cdot E_\theta G^2(X, \theta)]^{1/2} < \infty
    \]

    再次应用 DCT，取 $\theta' \rightarrow \theta$，我们便得到了条件 (vii)。证明完毕。

### 2.2 Bhattacharya 不等式

C-R 下界有时“太低”了，不足以给出更紧的方差下界。Bhattacharya 不等式通过引入更高阶的导数，是对 C-R 下界的自然推广。

!!! success "定理 8.2 (Bhattacharya 不等式)"

    在满足 C-R 下界定理中条件 (i) 和 (ii) 的基础上，如果我们加强其他条件为：
    
    (iii)* 对于 $i=1, \dots, K$，$f_\theta(x)$ 的 $i$ 阶导数存在且 $\int_{\mathcal{X}} \frac{\partial^i f_\theta(x)}{\partial \theta^i}\mu(dx) = 0$。
    (iv)* 对于 $i=1, \dots, K$，高阶对数导数的方差有限：$\int_{\mathcal{X}} \frac{1}{f_\theta(x)} \left( \frac{\partial^i f_\theta(x)}{\partial \theta^i} \right)^2 \mu(dx) < \infty$。
    (v)* $\hat{g}(X)$ 是 $g(\theta)$ 的方差有限的无偏估计，且对于任意 $i=1, \dots, K$ 满足积分可交换性：
    
    \[
    g^{(i)}(\theta) = \frac{\partial^i}{\partial\theta^i}g(\theta) = \int_{\mathcal{X}} \hat{g}(x)\frac{\partial^i f_\theta(x)}{\partial\theta^i}\mu(dx)
    \]

    那么，估计量的方差满足：

    \[
    var_\theta(\hat{g}(X)) \ge \tilde{g}^T(\theta) V^{-1}(\theta) \tilde{g}(\theta)
    \]

    其中 $V(\theta) = (V_{ij}(\theta))$ 矩阵的元素为 $V_{ij}(\theta) = E_\theta[\frac{1}{f^2_\theta(X)} \frac{\partial^i f_\theta(X)}{\partial\theta^i} \frac{\partial^j f_\theta(X)}{\partial\theta^j}]$，且导数向量 $\tilde{g}(\theta) = (g'(\theta), \dots, g^{(K)}(\theta))^T$。

??? proof "Bhattacharya 不等式的证明（点击展开）"

    令向量 $S = S_\theta(x) = (S_\theta^{(1)}(x), \dots, S_\theta^{(K)}(x))^T$，其中第 $i$ 个分量定义为：

    \[
    S_\theta^{(i)}(x) = \frac{1}{f_\theta(x)} \frac{\partial^i f_\theta(x)}{\partial\theta^i}
    \]

    * 由条件 (iii)* 可知，$E_\theta[S] = 0$。
    * 由条件 (iv)* 和 $V(\theta)$ 的定义可知，$var_\theta(S) = V(\theta)$。
    * 由条件 (v)* 可知，协方差 $cov_\theta(\hat{g}(X), S_\theta^{(i)}(X)) = E_\theta[\hat{g}(X) S_\theta^{(i)}(X)] - 0 = g^{(i)}(\theta)$。

    我们考虑组合向量 $(\hat{g}, S^T)^T$ 的协方差分块矩阵 $A$：

    \[
    A := var_\theta \begin{pmatrix} \hat{g} \\ S \end{pmatrix} = \begin{pmatrix} var_\theta(\hat{g}(X)) & \tilde{g}^T(\theta) \\ \tilde{g}(\theta) & V(\theta) \end{pmatrix}
    \]

    因为协方差矩阵 $A$ 是半正定的，其行列式 $|A| \ge 0$。根据分块矩阵行列式公式：

    \[
    |A| = |V(\theta)| \left[ var_\theta(\hat{g}(X)) - \tilde{g}^T(\theta) V^{-1}(\theta) \tilde{g}(\theta) \right]
    \]

    由于 $|V(\theta)| > 0$，必然有：

    \[
    var_\theta(\hat{g}(X)) - \tilde{g}^T(\theta) V^{-1}(\theta) \tilde{g}(\theta) \ge 0
    \]

    即证毕。显然当 $K=1$ 时，这正好退化为 Cramér-Rao 下界。

---

## 3. Kullback-Leibler 散度与可识别性

为了证明极大似然估计的相合性，我们需要一种衡量两个概率分布之间“距离”的工具，并定义参数的可识别性。

!!! info "定义 8.3 (Kullback-Leibler 散度)"

    从概率测度 $P_\theta$ 到 $P_\eta$ 的 Kullback-Leibler (KL) 散度定义为：

    \[
    D_{KL}(P_\eta || P_\theta) = E_\eta \left[ \log \frac{p_\eta(X)}{p_\theta(X)} \right], \quad X \sim P_\eta
    \]

    其中 $p_\theta$ 和 $p_\eta$ 分别是 $P_\theta$ 和 $P_\eta$ 的密度函数。

* KL 散度**不是一个真正的度量 (metric)**，因为一般来说不对称：$D_{KL}(P || Q) \ne D_{KL}(Q || P)$。
* 根据对数函数的凹性 (Jensen 不等式)，我们始终有 $D_{KL}(P || Q) \ge 0$。当且仅当 $P=Q$（在模型可识别的条件下）时，等号成立。

!!! info "定义 8.4 (可识别性 Identifiability)"

    如果对于任意 $\theta_1 \ne \theta_2$ ($\theta_1, \theta_2 \in \Theta$)，都有：

    \[
    \mu(x : f_{\theta_1}(x) \ne f_{\theta_2}(x)) > 0
    \]

    （即两个不同参数给出的概率分布在底测度 $\mu$ 下不几乎处处相等），则称参数族 $\mathbb{P}_\Theta := \{f_\theta(x) : \theta \in \Theta\}$ 是**可识别的 (Identifiable)**。

可识别性是保证 MLE 相合性的必要前置条件：如果参数不可识别，那么一致估计量根本不可能存在。

!!! success "引理 8.2 (最小化 KL 距离)"

    设 $\mathbb{P}_\Theta$ 是可识别的参数族。如果 $E_{\theta_0} \log f_{\theta_0}(X) < \infty$，那么目标函数 $M(\theta) := E_{\theta_0} [\log f_\theta(X) / f_{\theta_0}(X)]$ 在真实参数 $\theta_0$ 处**唯一地**达到它的最大值，即：

    \[
    E_{\theta_0} \log f_\theta(X) \le E_{\theta_0} \log f_{\theta_0}(X) < \infty
    \]

    **证明概要：** 对于 $\theta \in \Theta$，因为 $-\log(t)$ 是严格凸函数，应用 Jensen 不等式：
    
    \[
    E_{\theta_0} \log \frac{f_\theta}{f_{\theta_0}}(X) \le \log E_{\theta_0} \left[ \frac{f_\theta}{f_{\theta_0}}(X) \right] = \log \int \frac{f_\theta(x)}{f_{\theta_0}(x)} f_{\theta_0}(x) dx = \log \int f_\theta(x) dx = \log 1 = 0
    \]
    
    由可识别性条件，等号成立当且仅当 $\theta = \theta_0$。因此期望对数似然在真实参数处最大。

---

## 4. 极大似然估计的相合性 (Consistency of MLE)

!!! success "定理 8.3 (MLE 相合性根的存在定理)"

    设 $X_1, \dots, X_n$ i.i.d. $\sim P_\theta$ 且满足 Cramér-Rao 条件。假设存在 $\theta$ 的开邻域 $\Theta_\theta$ 使得：
    
    (i) 支撑集 $A := \{x | f_\theta(x) > 0\}$ 不依赖于 $\theta$。
    (ii) 对任意 $x \in A$，$f_{\theta}(x)$ 对任意 $\theta' \in \Theta_\theta$ 可导。
    (iii) 期望 $E_\theta \log f_{\theta'}(X)$ 对所有 $\theta' \in \Theta_\theta$ 存在且有限。
    (iv) 参数族是可识别的。
    
    那么，对于任意 $\epsilon > 0, \delta > 0$，存在 $m_{\epsilon, \delta} > 0$，使得当 $n > m_{\epsilon, \delta}$ 时：

    \[
    P_\theta\left\{ \text{似然方程 } \frac{d}{d\theta'} \sum_{i=1}^n \log f_{\theta'}(X_i) = 0 \text{ 在 } (\theta-\epsilon, \theta+\epsilon) \text{ 内有一个根} \right\} \ge 1 - \delta
    \]

??? proof "定理 8.3 的证明（点击展开）"

    不失一般性，假设 $\epsilon$ 足够小使得 $[\theta-\epsilon, \theta+\epsilon] \subset \Theta_\theta$。由弱大数定律 (WLLN) 及条件 (iii)：

    \[
    \frac{1}{n}\sum_{i=1}^n \log \frac{f_{\theta\pm\epsilon}(X_i)}{f_\theta(X_i)} \xrightarrow{P_\theta} E_\theta \log \frac{f_{\theta\pm\epsilon}(X)}{f_\theta(X)} := -\eta_{\theta\pm\epsilon} < 0
    \]

    （最后的不等号由引理 8.2 保证，因为真实参数处期望对数似然最大）。
    因此，对于任意 $\delta > 0, \xi > 0$，存在 $m = m_{\epsilon, \delta}$，使得对所有 $n > m$：

    \[
    P_\theta \left\{ \left| \frac{1}{n}\sum_{i=1}^n \log \frac{f_{\theta\pm\epsilon}(X_i)}{f_\theta(X_i)} + \eta_{\theta\pm\epsilon} \right| < \xi \right\} \ge 1 - \frac{\delta}{2}
    \]

    通过选取足够小的 $0 < \xi < \min\{\eta_{\theta-\epsilon}, \eta_{\theta+\epsilon}\}$，上述界限意味着当 $n > m$ 时，我们分别有极大概率保证边界点的对数似然值低于中心点：

    \[
    P_\theta(A) := P_\theta \left( \frac{1}{n}\sum_{i=1}^n \log f_\theta(X_i) > \frac{1}{n}\sum_{i=1}^n \log f_{\theta+\epsilon}(X_i) \right) \ge 1 - \frac{\delta}{2}
    \]

    \[
    P_\theta(B) := P_\theta \left( \frac{1}{n}\sum_{i=1}^n \log f_\theta(X_i) > \frac{1}{n}\sum_{i=1}^n \log f_{\theta-\epsilon}(X_i) \right) \ge 1 - \frac{\delta}{2}
    \]

    由于 $P(AB) = P(A) - P(AB^C) \ge P(A) - P(B^C) \ge 1 - \frac{\delta}{2} - \frac{\delta}{2} = 1 - \delta$，我们有：

    \[
    P_\theta \Big( l_n(\theta-\epsilon) < l_n(\theta) \text{ 且 } l_n(\theta+\epsilon) < l_n(\theta) \Big) \ge 1 - \delta
    \]

    既然在端点处函数值都比内部点小，且对数似然函数 $l_n(\theta')$ 是连续可导的，那么在区间 $(\theta-\epsilon, \theta+\epsilon)$ 内部必定存在一个局部极大值点。即导数为 0 的根存在的概率至少为 $1-\delta$。

**注**：上述定理仅保证了在真实参数附近“存在一个一致的根”，但并**没有保证该根就是我们要找的全局极大似然估计 (MLE)**。

!!! success "定理 8.4 (MLE 的相合性)"

    在定理 8.3 的条件下，定义 $\hat{\theta}_n$ 为似然方程的根（当方程恰好有唯一根时）。如果：

    \[
    \lim_{n \rightarrow \infty} P_\theta(\text{似然方程只有单一根}) = 1
    \]

    那么：

    \[
    \hat{\theta}_n \xrightarrow{P_\theta} \theta
    \]

证明非常直接：当 $n$ 足够大时，区间 $(\theta-\epsilon, \theta+\epsilon)$ 内有根的概率趋于 1（定理 8.3），同时全空间内只有单一根的概率也趋于 1，这两个事件的交集概率趋于 1，意味着这个唯一的根必定落在 $(\theta-\epsilon, \theta+\epsilon)$ 这个极小区间内。

**反例与改进：Cauchy 极大似然估计**

Cauchy 分布 $f_\theta(x) = \frac{1}{\pi\{1+(x-\theta)^2\}}$ 的对数似然方程往往存在**多个局部极值点**（Reeds, 1985 证明其根的数目渐近表现为 $2 \times \text{Poisson}(1/\pi) + 1$）。全局最大值点通常能很好地从其他极值点中分离出来 (well-separated)。
实践中，为了避免陷入错误的局部极值，我们通常采用 **One-step update (单步更新) 方法**（Newton-Raphson）：

1. 先计算一个具有 $\sqrt{n}$-相合性的稳健初始估计量 $\hat{\theta}_0$（如样本中位数）。

2. 执行单步 Newton-Raphson 更新：$\hat{\theta}_n := \hat{\theta}_0 - (\mathbb{P}_n \ddot{m}(\hat{\theta}_0))^{-1} \mathbb{P}_n \dot{m}(\hat{\theta}_0)$

这在保证相合性的同时，能提升统计效率。

---

## 5. 极大似然估计的渐近正态性 (Asymptotic Normality of MLE)

!!! success "定理 8.5 (MLE 的渐近正态性)"

    设 $X_1, \dots, X_n$ i.i.d. $\sim P_{\theta_0}$，$\Theta \subset \mathbb{R}$，且在真实参数 $\theta_0$ 的开邻域 $\Theta_0$ 内满足：
    
    (i) $f_{\theta'}(x) > 0, \forall x, \theta' \in \Theta_0$。
    (ii) $\forall x$，$f_{\theta'}(x)$ 在 $\Theta_0$ 内 3 次可导。
    (iii) 存在非负函数 $M(x)$，且 $E_{\theta_0}M(X) < \infty$，使得对任意 $\theta' \in \Theta_0$ 有 $\left| \frac{d^3}{d\theta'^3}\log f_{\theta'}(x) \right| \le M(x)$。
    (iv) 对 $l=1,2$，积分等式 $\int f_{\theta'}(x)\mu(dx)=1$ 可在积分号下求两次导。
    (v) $0 < I_1(\theta') < \infty$。
    (vi) $\lim_{n\rightarrow \infty} P_\theta(\hat{\theta}_n \text{ 是似然方程的根}) = 1$ 且 $E_{\theta} |\log f_{\theta'}(X)| < \infty$。
    (vii) MLE 满足相合性 $\hat{\theta}_n \rightarrow \theta_0$ 且具有可识别性。

    那么，MLE 的渐近分布为：

    \[
    \sqrt{n}(\hat{\theta}_n - \theta_0) \xrightarrow{d} N(0, I_1^{-1}(\theta_0))
    \]

??? proof "定理 8.5 的证明（泰勒展开法，点击展开）"

    令 $l_n(\theta') = \frac{1}{n} \sum_{i=1}^n \log f_{\theta'}(X_i)$。由于 $\hat{\theta}_n$ 是得分方程的根，在 $\theta_0$ 处对 $l'_n(\hat{\theta}_n)$ 进行二阶泰勒展开：

    \[
    0 = l'_n(\hat{\theta}_n) = l'_n(\theta_0) + l''_n(\theta_0)(\hat{\theta}_n - \theta_0) + \frac{1}{2}l'''_n(\theta_1)(\hat{\theta}_n - \theta_0)^2
    \]

    其中 $\theta_1$ 介于 $\hat{\theta}_n$ 与 $\theta_0$ 之间。根据弱大数定律 (WLLN) 及 Fisher 信息量的定义：

    \[
    l''_n(\theta_0) = \frac{1}{n} \sum_{i=1}^n \frac{d^2}{d\theta^2}\log f_\theta(X_i) \xrightarrow{P_\theta} -I_1(\theta_0) \in (-\infty, 0)
    \]

    即可以写为 $l''_n(\theta_0) = -I_1(\theta_0) + o_{P_\theta}(1)$。
    
    对于三阶余项，由于条件 (iii)：

    \[
    |l'''_n(\theta_1)| = \left| \frac{1}{n}\sum_{i=1}^n \frac{d^3}{d\theta^3}\log f_\theta(X_i) \Big|_{\theta=\theta_1} \right| \le \frac{1}{n}\sum_{i=1}^n M(X_i) \xrightarrow{P_\theta} E_{\theta_0} M(X) < \infty
    \]

    因此序列 $\{l'''_n(\theta_1)\}$ 依概率有界，即胎紧 (tight)，写为 $O_{P_\theta}(1)$。
    已知 $\hat{\theta}_n \xrightarrow{P_\theta} \theta_0$，即 $\hat{\theta}_n - \theta_0 = o_{P_\theta}(1)$。所以三阶余项部分整体为：

    \[
    (\hat{\theta}_n - \theta_0)^2 l'''_n(\theta_1) = o_{P_\theta}(\hat{\theta}_n - \theta_0)
    \]

    将上述渐近项代回展开式：

    \[
    0 = l'_n(\theta_0) + (-I_1(\theta_0) + o_{P_\theta}(1))(\hat{\theta}_n - \theta_0) + o_{P_\theta}(\hat{\theta}_n - \theta_0)
    \]

    两边同乘 $\sqrt{n}$，并重排项解出 $\sqrt{n}(\hat{\theta}_n - \theta_0)$：

    \[
    \sqrt{n}(\hat{\theta}_n - \theta_0) = -I_1^{-1}(\theta_0)\sqrt{n}l'_n(\theta_0) + o_{P_\theta}(\sqrt{n}(\hat{\theta}_n - \theta_0))
    \]

    由中心极限定理 (CLT)，得分函数的标准化项趋于正态：

    \[
    \sqrt{n}l'_n(\theta_0) = \frac{1}{\sqrt{n}}\sum_{i=1}^n \frac{\partial \log f_\theta(X_i)}{\partial\theta} \xrightarrow{d} N(0, I_1(\theta_0))
    \]

    最后，应用 Slutsky 定理（常数乘法）：

    \[
    -I_1^{-1}(\theta_0)\sqrt{n}l'_n(\theta_0) \xrightarrow{d} N(0, I_1^{-1}(\theta_0))
    \]

    即证 $\sqrt{n}(\hat{\theta}_n - \theta_0) \xrightarrow{d} N(0, I_1^{-1}(\theta_0))$。

此外，对于**紧致凸参数空间**，我们还有如下定理（Bijma & Jonker & Van der Vaart, 2017）：
如果参数空间 $\Theta$ 是紧致且凸的，模型是可识别的；映射 $\vartheta \mapsto \log p_\vartheta(x)$ 连续可导，且导数 $|l_\vartheta(x)| \le L(x)$ 有平方可积包络 $E_\theta L^2(X_1) < \infty$；并且真实参数 $\theta_0$ 是 $\Theta$ 的内点，Fisher 信息矩阵 $I(\vartheta)$ 连续且正定，那么：

\[
\sqrt{n}(\hat{\theta}_n - \theta_0) \sim N(0, I^{-1}(\theta_0))
\]

---

## 6. 基于似然的统计推断与似然比检验 (Likelihood-based Inference & LRT)

由 $\hat{\theta}_n$ 的渐近正态性，我们可以为 $\theta$ 构造置信区或进行假设检验。利用 Slutsky 定理和连续映射定理，如果 $I_1(\theta)$ 连续：

\[
\sqrt{n} I_1^{1/2}(\hat{\theta}_n) (\hat{\theta}_n - \theta) \xrightarrow{d} N_p(0, I_p)
\]

对于多维假设检验 $H_0: \theta = \theta_0$ VS $H_1: \theta \ne \theta_0$，常用的 **Wald 检验统计量**为：

\[
\mathcal{W}_n \hat{=} n(\hat{\theta}_n - \theta_0)^T I_1(\theta_0) (\hat{\theta}_n - \theta_0) \rightarrow \chi^p_2
\]

如果 $W_n > \chi^2_{p, 1-\alpha}$，则拒绝原假设。

### 似然比检验 (Likelihood Ratio Test, LRT)

根据 Neyman-Pearson 引理，在简单原假设对简单备择假设时，似然比检验是一致最优检验 (UMP)。我们考虑更一般的复合假设检验：

\[
H_0: \theta \in \Theta_0 \quad \text{VS} \quad H_1: \theta \in \Theta_1
\]

其中全空间 $\Theta = \Theta_0 \cup \Theta_1$ 且 $\Theta_0 \cap \Theta_1 = \emptyset$。
设 $\hat{\theta}_n$ 是在 $\Theta$ 上的全局 MLE，$\hat{\theta}_{n,0}$ 是在 $\Theta_0$ 限制下的受限 MLE。

**对数似然比统计量 (Log Likelihood Ratio Statistic)** 定义为：

\[
LR_n = -2 \log \frac{L_n(\hat{\theta}_{n,0})}{L_n(\hat{\theta}_n)} = 2 \left\{ l_n(\hat{\theta}_n) - l_n(\hat{\theta}_{n,0}) \right\}
\]

一般的决策规则是：当 $LR_n$ 过大时，拒绝原假设 $H_0$。

**将 LR 视为渐近距离度量**

如果对数似然函数足够平滑（如 3 次连续可导）且 $\tilde{\theta}_n$ 具有 $\sqrt{n}$-相合性，由全局 MLE $\hat{\theta}_n$ 处的一阶最优性条件（导数为 0）和泰勒展开，我们有：

\[
\log L_n(\tilde{\theta}_n) - \log L_n(\hat{\theta}_n) = 0 + \frac{1}{2}(\tilde{\theta}_n - \hat{\theta}_n)^T \frac{\partial^2 \log L_n(\hat{\theta}_n)}{\partial\theta\partial\theta^T} (\tilde{\theta}_n - \hat{\theta}_n) + o_p(1)
\]

因此，似然比 $LR_n$ 实际上可以看作是**在原假设和全空间下两个极大似然估计量之间的渐近距离**（以 Fisher 信息量为权重的二次型）：

\[
LR_n = (\tilde{\theta}_n - \hat{\theta}_n)^T \left( -\frac{\partial^2 \log L_n(\theta_0)}{\partial\theta\partial\theta^T} \right) (\tilde{\theta}_n - \hat{\theta}_n) + o_p(1)
\]