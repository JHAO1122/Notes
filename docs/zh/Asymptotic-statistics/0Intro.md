# 第十一章：M-估计量与 Z-估计量（一）

在本章中，我们将研究 **M-估计量**（由 Peter J. Huber 提出）和 **Z-估计量** 的相合性 (Consistency) 与渐近正态性 (Asymptotic Normality)。

极大似然估计量 (MLEs) 和矩估计量 (ME) 分别可以看作是 M-估计量和 Z-估计量的特例。

* 假设我们感兴趣的参数（或“泛函”）$\theta$ 附加在观测值的分布上。

* 观测值序列服从分布：

\[
X_n := (X_1, \dots, X_n) \sim p_\theta(X)
\]

---

## 1. 基本定义与示例

### 1.1 M-估计量

!!! info "定义 11.1 (M-估计量 M-estimator)"

    M-估计量的目标是寻找一个估计量 $\hat{\theta}_n := \hat{\theta}_n(X_1, \dots, X_n)$，使得某种类型的随机（经验）准则函数达到最大化：

    \[
    \theta \mapsto M_n(\theta)
    \]

    例如，通常的形式为样本平均：

    \[
    M_n(\theta) = \frac{1}{n} \sum_{i=1}^n m_\theta(X_i)
    \]

### 1.2 Z-估计量

通常，最大化问题可以通过令导数（或梯度）等于零来求解。这就自然引出了 Z-估计量的概念。

!!! info "定义 11.2 (Z-估计量 Z-estimator)"

    Z-估计量满足如下的**估计方程 (Estimating Equations)**：

    \[
    \Psi_n(\theta) = \frac{1}{n} \sum_{i=1}^n \psi_\theta(X_i) = 0
    \]

    更一般地，Z-估计量就是寻找一个估计量 $\hat{\theta}_n$，使其成为上述估计方程的根。

**关于 M- 和 Z-估计量的几点注记：**

* M- 和 Z-估计量**不要求**观测值必须独立同分布 (IID) 或具有独立的结构。

* 目标函数 $-M_n(\theta)$ 的最小化问题可能是**非凸的 (non-convex)**。

* 在实际计算中，Z-估计量通常通过数值优化算法求解，例如（拟）牛顿法、梯度下降法、随机梯度下降法等（尤其在面临非凸问题时）。

---

### 1.3 经典示例

#### 示例 1：极大似然估计与伪极大似然估计

* **MLE**：设 $X_1, \dots, X_n \sim p_\theta$。极大似然估计量最大化似然函数 $\prod_{i=1}^n p_\theta(X_i)$，等价于最大化对数似然函数：

\[
M_n(\theta) := \frac{1}{n} \sum_{i=1}^n \log p_\theta(X_i)
\]

* **伪极大似然估计 (Pseudo-MLE)**：此时 $X_i$ 之间可能是相依的（dependent），但我们仍然借用并最大化上述的对数似然准则函数。

#### 示例 2：位置估计量 (Location Estimators)

样本均值和样本中位数都可以视为 Z-估计量，它们分别满足如下估计方程：

* **样本均值**：

\[
\Psi_n(\theta) = \frac{1}{n} \sum_{i=1}^n (X_i - \theta) = 0
\]

* **样本中位数**：

\[
\Psi_n(\theta) = \frac{1}{n} \sum_{i=1}^n \text{sign}(X_i - \theta) = 0
\]

#### 示例 3：分位数 (Quantiles)

对于分布函数 $F$，其分位数函数定义为：

\[
F^{-1}(p) := \inf\{x \in \mathbb{R} : F(x) \ge p\}
\]

* **中位数**：$\text{med}(X) = F^{-1}(0.5)$

* **一般分位数**：参数可以通过最小化期望损失得到，即 $\theta_0 = \arg\min_\theta E \rho_\tau(X - \theta) = F^{-1}(\tau)$。

引入**检验函数 (Check function)** 作为损失函数：

\[
\rho_\tau(y) = y(\tau - I_{(y<0)})
\]

由此，$\tau$-样本分位数 $\hat{\theta}$ 既可以看作 M-估计量也可以看作 Z-估计量：

* 作为 M-估计量的准则函数：

\[
M_n(\theta) := \frac{1}{n} \sum_{i=1}^n \rho_\tau(X_i - \theta)
\]

* 作为 Z-估计量的估计方程：

\[
\Psi_n(\theta) := \frac{1}{n} \sum_{i=1}^n \left( (1-\tau) 1\{X_i < \theta\} - \tau 1\{X_i > \theta\} \right) = 0
\]

*对于小样本量 $n$，van der Vaart 的书给出了一种分位数 $\hat{\theta}$ 的替代定义，即求解不等式：$-1 < n\Psi_n(\theta) < 1$。*

#### 示例 4：Huber 估计量

Huber 估计量起源于稳健统计学 (Robust Statistics)，旨在控制极端异常数据点对估计结果的影响。
它基于 Huber $\psi$ 函数：

\[
\psi(x) = [x]_{-k}^k := \begin{cases} -k & \text{if } x \le -k \\ x & \text{if } |x| \le k \\ k & \text{if } x \ge k \end{cases}
\]

Huber 估计量求解如下的估计方程：

\[
\Psi_n(\theta) = \frac{1}{n} \sum_{i=1}^n \psi(X_i - \theta) = 0
\]

当 $k$ 较大时，Huber 估计量的行为更接近于非稳健的样本均值；当 $k$ 较小时，其行为更接近于高度稳健的样本中位数。因此，它在均值和中位数之间取得了极好的平衡。

*(注：经验上的 Z-估计量方程图像通常呈现阶梯状（如分位数的经验方程），而其总体期望方程则是平滑的曲线。)*

---

## 2. M-估计量的相合性 (Consistency)

对于估计量 $\hat{\theta}_n$，我们希望证明当样本量趋于无穷时，它能够依概率收敛到真实参数：$\hat{\theta}_n \xrightarrow{p} \theta_0$，其中 $\theta, \theta_0 \in \Theta$，$\Theta$ 赋予了距离度量 $d$。

假设 $\hat{\theta}_n$ 最大化了随机（经验）准则函数：

\[
\hat{\theta}_n = \arg\max_{\theta \in \Theta} M_n(\theta) = \arg\min_{\theta \in \Theta} -M_n(\theta)
\]

其中 $-M_n(\theta) =: L(P_n, P)$ 可以被看作经验损失函数。

!!! info "定义 11.3 (真实参数 True parameter)"

    真实参数 $\theta_0$ 通常定义为确定性（真实）准则函数 $M(\theta) =: E m_\theta(X)$ 的最大化点：

    \[
    \theta_0 = \arg\max_{\theta \in \Theta} E m_\theta(X) = \arg\min_{\theta \in \Theta} E [-m_\theta(X)]
    \]

由大数定律 (LLN)，对任意给定的 $\theta$，有 $M_n(\theta) \xrightarrow{p} M(\theta)$。但需要注意的是，**这种逐点收敛通常并非在整个参数空间 $\Theta$ 上一致成立的！**

为了证明相合性，给定一个任意的随机函数 $\theta \mapsto M_n(\theta)$，我们考虑满足如下**近似最大化条件 (Nearly maximization condition)** 的估计量序列 $\{\hat{\theta}_n\}$：

\[
M_n(\hat{\theta}_n) \ge \sup_{\theta \in \Theta} M_n(\theta) - o_P(1) \ge M_n(\theta_0) - o_P(1)
\]

*(例如：当 $-M_n(\theta)$ 是强凸函数时满足此条件，即当且仅当二阶导 $\ddot{M}_n(\theta) \ge O(1)I_p > 0, \forall \theta \in \Theta$。)*

### 2.1 M-估计量相合性基本定理

!!! success "定理 11.1 (M-估计量的相合性)"

    设 $M_n$ 为一列随机函数，$M$ 为一个固定的确定性函数。若对任意 $\epsilon > 0$，以下三个条件均成立：

    * **C1. 一致收敛 (Uniform convergence)**：

    \[
    \sup_{\theta \in \Theta} |M_n(\theta) - M(\theta)| \xrightarrow{p} 0
    \]

    * **C2. 良好分离/可识别性 (Well-separation / Identifiability)**：

    \[
    \sup_{\theta: d(\theta, \theta_0) \ge \epsilon} M(\theta) < M(\theta_0)
    \]

    * **C3. 近似最大化**：序列 $\{\hat{\theta}_n\}$ 满足上述的近似最大化条件。

    那么，必有 $\hat{\theta}_n \xrightarrow{p} \theta_0$。

??? proof "定理 11.1 的证明（点击展开）"

    由良好分离假设 C2，对于任意 $\epsilon > 0$，存在 $\eta > 0$ 使得：
    
    \[
    M(\theta) < M(\theta_0) - \eta \quad \text{对所有满足 } d(\theta, \theta_0) \ge \epsilon \text{ 的 } \theta \text{ 成立。}
    \]

    将 $\theta$ 替换为 $\hat{\theta}_n$，我们得到事件的包含关系：
    $\{d(\hat{\theta}_n, \theta_0) \ge \epsilon\} \subseteq \{M(\hat{\theta}_n) < M(\theta_0) - \eta\}$。
    
    因此在概率上：
    
    \[
    P\{d(\hat{\theta}_n, \theta_0) \ge \epsilon\} \le P\{M(\theta_0) - M(\hat{\theta}_n) > \eta\} \quad \text{--- (11.2)}
    \]

    接下来，我们要证明 (11.2) 式的右侧趋于 0。综合使用条件 C3 和 C1。
    由 C3（近似最大化）：
    
    \[
    M_n(\hat{\theta}_n) \ge M_n(\theta_0) - o_P(1) = M(\theta_0) - o_P(1) \quad \text{--- (11.3)}
    \]
    
    注意 (11.3) 中最后一步等号的成立，是因为 C1 暗示了单点处的依概率收敛：$M_n(\theta_0) \xrightarrow{p} M(\theta_0)$。
    
    进一步利用 (11.3) 和 C1：
    
    \[
    M(\theta_0) - M(\hat{\theta}_n) \le M_n(\hat{\theta}_n) - M(\hat{\theta}_n) + o_P(1) \le \sup_{\theta \in \Theta} |M_n(\theta) - M(\theta)| + o_P(1)
    \]
    
    由于一致收敛性 C1，上式的上确界项是 $o_p(1)$。
    所以，整体上 $M(\theta_0) - M(\hat{\theta}_n) = o_P(1)$。
    
    令 $n \rightarrow \infty$，代入 (11.2) 式，这意味着极限概率为 0，从而推出：
    
    \[
    d(\hat{\theta}_n, \theta_0) \xrightarrow{p} 0
    \]
    
    证明完毕。 $\square$

*(注：如果目标函数 $M(\theta)$ 在极大值点附近非常平坦，或者存在局部极大值任意逼近全局极大值，则良好分离条件 C2 将不成立。)*

### 2.2 修改条件下的相合性推论

!!! success "推论 11.2 (修改条件下的 M-估计量相合性)"

    假设以下条件成立：

    * (i) 定理 11.1 中的一致收敛条件 C1；

    * (ii) **唯一最大化 (Unique maximization)**：$M(\theta) := E m_\theta(X)$ 在 $\theta_0$ 处唯一取得最大值；

    * (iii) **紧致性 (Compactification)**：参数空间 $\Theta$ 是紧致的；

    * (iv) **连续性 (Continuous M-function)**：确定性函数 $M(\theta)$ 是连续的。

    那么，对于任何满足近似最大化条件 $(v) \ M_n(\hat{\theta}_n) \ge M_n(\theta_0) - o_P(1)$ 的估计量 $\hat{\theta}_n$，均有 $\hat{\theta}_n \xrightarrow{p} \theta_0$。

??? proof "推论 11.2 的证明（点击展开）"

    对于任意给定的 $\delta > 0$，定义邻域的外部集合 $B_\delta^c(\theta_0) := \{\theta : d(\theta, \theta_0) \ge \delta\}$。
    根据条件 (ii)-(iv)（连续函数在紧集上必定取得最大值），存在某个点 $\theta^* \in \Theta \cap B_\delta^c(\theta_0)$ 使得：
    
    \[
    \sup_{\theta \in \Theta \cap B_\delta^c(\theta_0)} M(\theta) = M(\theta^*) < M(\theta_0)
    \]

    接下来我们证明，对于任意 $\epsilon > 0$，当 $n$ 足够大时，以下不等式链以趋于 1 的概率 (with probability approaching 1, wpal) 成立：
    
    \[
    M(\hat{\theta}_n) > M_n(\hat{\theta}_n) - \epsilon/3 > M_n(\theta_0) - 2\epsilon/3 > M(\theta_0) - \epsilon \quad \text{--- (11.4)}
    \]

    我们逐项分析 (11.4) 中的不等关系：
    
    **第一个 ">" (基于 C1)：**
    定义事件 $E_1 := \{M(\hat{\theta}_n) \le M_n(\hat{\theta}_n) - \epsilon/3\}$。由于一致收敛性 $\sup_{\theta \in \Theta} |M_n(\theta) - M(\theta)| \xrightarrow{p} 0$，必然有 $P(E_1) = o(1)$。
    
    **第二个 ">" (基于条件 v)：**
    定义事件 $E_2 := \{M_n(\hat{\theta}_n) \le M_n(\theta_0) - \epsilon/3\}$。由于近似最大化要求 $M_n(\hat{\theta}_n) \ge M_n(\theta_0) - o_P(1)$，故 $P(E_2) = o(1)$。
    
    **第三个 ">" (基于 C1)：**
    类似于第一项，定义事件 $E_3 := \{M_n(\theta_0) - 2\epsilon/3 \le M(\theta_0) - \epsilon\} = \{M_n(\theta_0) \le M(\theta_0) - \epsilon/3\}$。同样由点态收敛性，$P(E_3) = o(1)$。

    现在，令 $\epsilon = M(\theta_0) - M(\theta^*) > 0$。将此 $\epsilon$ 代入上述不等式链的头尾，我们得到以趋向于 1 的概率成立：
    
    \[
    M(\hat{\theta}_n) > M(\theta^*) \quad \text{wpal}
    \]
    
    具体而言：
    
    \[
    P(E_1^c \cap E_2^c \cap E_3^c) \ge 1 - \sum_{i=1}^3 P(E_i) \rightarrow 1
    \]

    需要注意到，因为 $\theta^*$ 是在边界区域 $\Theta \cap B_\delta^c(\theta_0)$ 上的最大值点，所以如果 $M(\hat{\theta}_n) > M(\theta^*)$，则必然意味着 $\hat{\theta}_n$ 落在了中心邻域内：
    
    \[
    \{M(\hat{\theta}_n) > M(\theta^*)\} \subseteq \{d(\hat{\theta}_n, \theta_0) < \delta\}
    \]

    因此，令 $n \rightarrow \infty$ 时：
    
    \[
    1 \leftarrow P\{M(\hat{\theta}_n) > M(\theta^*)\} \le P\{d(\hat{\theta}_n, \theta_0) < \delta\} \le 1
    \]
    
    这就证明了 $d(\hat{\theta}_n, \theta_0) \xrightarrow{p} 0$。 $\square$

---

## 3. Z-估计量的相合性

由于 Z-估计量是通过令目标函数等于零（而不是最大化）来求解的，其相合性定理也是 M-估计量的一个自然推广。

!!! success "定理 11.3 (Z-估计量的相合性 Consistency of Z-Estimator)"

    设 $\Psi_n$ 为一列随机的向量值函数，$M$ 为固定的向量值函数。若对任意 $\epsilon > 0$，以下条件成立：

    * **C1*. 一致收敛**：

    \[
    \sup_{\theta \in \Theta} \|\Psi_n(\theta) - \Psi(\theta)\| \xrightarrow{p} 0
    \]

    * **C2*. 良好分离 (Identifiability)**：

    \[
    \inf_{\theta: d(\theta, \theta_0) \ge \epsilon} \|\Psi(\theta)\| > 0 = \|\Psi(\theta_0)\|
    \]

    * **C3*. 近似零点条件**：估计量序列 $\{\hat{\theta}_n\}$ 满足 $\Psi_n(\hat{\theta}_n) = o_P(1)$。

    则必有 $\hat{\theta}_n \xrightarrow{p} \theta_0$。

**证明思路**：此定理直接由 M-估计量的相合性定理得出。只需构造准则函数 $M_n(\theta) = -\|\Psi_n(\theta)\|$ 以及 $M(\theta) = -\|\Psi(\theta)\|$。可以看到，Z-估计量中的**近似零点条件 (nearly zero condition)** 完美对应于 M-估计量中的**近似最大化条件**。

### 3.1 弱化一致收敛条件

在实际应用中，定理 11.3 中要求的一致收敛性 $C1^*$ 通常过于严格且难以验证（通常需要借助经验过程理论 Empirical Process 才能建立严谨的保证）。
引理 11.4 提供了一个无需一致收敛性即可证明相合性的途径。

!!! success "引理 11.4 (基于单调性/连续性的相合性，vdv p47)"

    设 $\Theta$ 为实数轴的子集，$\Psi_n$ 为随机函数，$\Psi$ 为固定函数，且对每一个给定的 $\theta$，有逐点收敛：$\Psi_n(\theta) \xrightarrow{p} \Psi(\theta)$。

    假设以下条件之一成立：

    * **(a1) 连续且唯一根**：每个映射 $\theta \mapsto \Psi_n(\theta)$ 都是连续的，且恰好有唯一的一个零点 $\hat{\theta}_n$；或者
    * **(a2) 单调有界**：映射是非减的（或非增的），并且 $\Psi_n(\hat{\theta}_n) = o_P(1)$；

    同时假设：

    * **(b) 真实根的符号交替**：设 $\theta_0$ 满足对任意 $\epsilon > 0$，都有 $\Psi(\theta_0 - \epsilon) < 0 < \Psi(\theta_0 + \epsilon)$。

    那么，$\hat{\theta}_n \xrightarrow{p} \theta_0$。

??? proof "引理 11.4 的证明（点击展开）"

    **情况 1：连续且唯一根 (a1 成立时)**

    如果映射 $\theta \mapsto \Psi_n(\theta)$ 是连续的，且在 $\hat{\theta}_n$ 处有唯一的零点。那么，事件 $\{\Psi_n(\theta_0 - \epsilon) < 0 \text{ 且 } \Psi_n(\theta_0 + \epsilon) > 0\}$ 必定蕴含了零点介于这两者之间，即 $\{\theta_0 - \epsilon < \hat{\theta}_n < \theta_0 + \epsilon\}$。
    
    在概率上，这表示为：
    
    \[
    P(\Psi_n(\theta_0 - \epsilon) < 0, \Psi_n(\theta_0 + \epsilon) > 0) \le P(\theta_0 - \epsilon < \hat{\theta}_n < \theta_0 + \epsilon)
    \]
    
    由于对于给定的 $\theta_0 \pm \epsilon$，我们有逐点依概率收敛 $\Psi_n(\theta_0 \pm \epsilon) \xrightarrow{p} \Psi(\theta_0 \pm \epsilon)$。结合条件 (b)，不等式左侧的概率当 $n \rightarrow \infty$ 时趋于 1。
    因此，右侧的概率也必定趋于 1，这就证明了 $\hat{\theta}_n$ 的相合性。

    **情况 2：非减且近似零点 (a2 成立时)**

    如果映射 $\theta \mapsto \Psi_n(\theta)$ 是非减的。首先注意到：
    
    \[
    \{\Psi_n(\theta_0 - \epsilon) < -\eta\} = \{\Psi_n(\theta_0 - \epsilon) < -\eta, \hat{\theta}_n \le \theta_0 - \epsilon\} \cup \{\Psi_n(\theta_0 - \epsilon) < -\eta, \hat{\theta}_n \ge \theta_0 - \epsilon\}
    \]
    
    由于 $\Psi_n$ 的非减性，若 $\hat{\theta}_n \ge \theta_0 - \epsilon$，则 $\Psi_n(\hat{\theta}_n) \ge \Psi_n(\theta_0 - \epsilon) < -\eta$。
    令 $E_- := \{\Psi_n(\hat{\theta}_n) < -\eta\}$。由于已知条件 $\Psi_n(\hat{\theta}_n) = o_P(1)$，对任意 $\eta > 0$ 都有 $P(E_-) = o(1)$。
    所以，上述事件并集可以放缩为：
    
    \[
    \subseteq E_- \cup \{\hat{\theta}_n \ge \theta_0 - \epsilon\}
    \]
    
    同样地，我们可以为右尾定义事件 $E_+ := \{\Psi_n(\hat{\theta}_n) > \eta\}$，也有 $P(E_+) = o(1)$。
    
    因此，对于任意给定的 $\epsilon, \eta > 0$：
    
    \[
    P(\Psi_n(\theta_0 - \epsilon) < -\eta, \Psi_n(\theta_0 + \epsilon) > \eta) \le P(\theta_0 - \epsilon < \hat{\theta}_n < \theta_0 + \epsilon) + P(E_-) + P(E_+)
    \]
    
    整理后得到：
    
    \[
    P(\Psi_n(\theta_0 - \epsilon) < -\eta, \Psi_n(\theta_0 + \epsilon) > \eta) = P(\theta_0 - \epsilon < \hat{\theta}_n < \theta_0 + \epsilon) + o(1)
    \]
    
    选择合适的阈值，令 $2\eta = \min\{-\Psi(\theta_0 - \epsilon), \Psi(\theta_0 + \epsilon)\} > 0$。由于逐点收敛性，左侧的概率趋于 1，因此右侧的概率同样趋于 1，相合性得证。 $\square$

### 3.2 示例：应用引理 11.4 证明中位数的相合性

**样本中位数** $\hat{\theta}_n$ 可以视为方程 $\theta \mapsto \Psi_n(\theta) = n^{-1} \sum_{i=1}^n \text{sign}(X_i - \theta)$ 的零点。

根据大数定律 (LLN)，对于每一个固定的 $\theta$，我们有：

\[
\Psi_n(\theta) \xrightarrow{p} \Psi(\theta) = E[\text{sign}(X - \theta)] = P(X > \theta) - P(X < \theta)
\]

现在我们验证引理 11.4 的条件：

* **条件 (a) 验证**：函数 $\theta \mapsto \Psi_n(\theta)$ 显然是**非增的 (non-increasing)**（此处的引理对称适用于非增情况，只需在证明中取反号）。

* **条件 (b) 验证**：我们需要证明 $-\Psi(\theta_0 - \epsilon) < 0 < -\Psi(\theta_0 + \epsilon)$。
观察极限函数在 $\theta_0$ 附近的值：

\[
\Psi(\theta_0 - \epsilon) = P(X > \theta_0 - \epsilon) - P(X < \theta_0 - \epsilon) = 1 - 2P(X < \theta_0 - \epsilon)
\]

在总体中位数 $\theta_0$ 处，有：

\[
\Psi(\theta_0) = P(X > \theta_0) - P(X < \theta_0) = 0 \implies P(X > \theta_0) = P(X < \theta_0) = 0.5
\]

同理：

\[
\Psi(\theta_0 + \epsilon) = P(X > \theta_0 + \epsilon) - P(X < \theta_0 + \epsilon) = 1 - 2P(X < \theta_0 + \epsilon)
\]

如果 $X$ 的分布是连续的，并且总体的中位数是唯一确定的，那么我们必然有：

\[
P(X < \theta_0 - \epsilon) < 0.5 < P(X < \theta_0 + \epsilon) \quad \text{对所有 } \epsilon > 0 \text{ 成立。}
\]

因此，条件 (a2) 与 (b) 均满足。应用引理 11.4，直接得出结论：样本中位数相合地收敛于总体中位数，即 $\hat{\theta}_n \xrightarrow{p} \theta_0$。这避免了复杂的经验过程理论的推导。