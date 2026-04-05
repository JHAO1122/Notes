# 第一章：数据分布与随机收敛 

在统计推断中，寻找统计量 $\hat{\theta}_n$ 在有限样本下的精确分布 $F_{\hat{\theta}_n}(x)$ 往往是极其困难的。然而，通过让样本量 $n \to \infty$，我们可以利用渐近理论极大地简化问题，并获得质量极高的近似分布。这不仅能帮助我们构建近似的置信区间与假设检验，还能在理论上评估不同推断方法的渐近效率。

---

## 1. 随机收敛的基本定义 (Stochastic Convergence)

设 $\{X_n\}$ 为定义在同一个概率空间 $(\Omega, \mathcal{A}, P)$ 上的 $\mathbb{R}^p$ 维随机向量序列，$d(x,y)$ 为欧氏距离。

!!! abstract "定义：四大随机收敛"

    **1. 几乎必然收敛 (Almost-Sure Convergence, $X_n \xrightarrow{a.s.} X$)**：
    
    $$
    P(\lim_{n \to \infty} d(X_n, X) = 0) = 1
    $$
    
    *(直观理解：100% 确定且 100% 准确)*。

    **2. 依概率收敛 (Convergence in Probability, $X_n \xrightarrow{P} X$)**：
    对于任意给定的 $\epsilon > 0$：
    
    $$
    \lim_{n \to \infty} P(d(X_n, X) < \epsilon) = 1
    $$
    
    *(直观理解：100% 确定，但不一定绝对准确)*。

    **3. $r$ 阶平均收敛 (Convergence in $r_{th}$ mean, $X_n \xrightarrow{L^r} X$)**：
    
    $$
    \lim_{n \to \infty} E[d(X_n, X)^r] = 0
    $$

    **4. 依分布收敛 / 弱收敛 (Convergence in Distribution, $X_n \xrightarrow{d} X$)**：
    设 $F_n$ 和 $F$ 分别为 $X_n$ 和 $X$ 的累积分布函数 (CDF)。如果对于 $F$ 的所有连续点 $x \in \mathcal{C}_F$：
    
    $$
    \lim_{n \to \infty} F_n(x) = F(x)
    $$
    
    则称 $X_n$ 依分布收敛到 $X$。
    
    > **注记**：依分布收敛是最弱的收敛，也是统计推断中最核心的收敛。它**不要求** $X_n$ 和 $X$ 定义在同一个概率空间上。且不连续点集合 $\mathcal{C}_F^c$ 至多只有可数个。

---

## 2. Polya 定理与渐近正态性

在微积分中，闭区间上的连续函数是一致连续的。这可以推广到全空间上的累积分布函数：

!!! info "引理 1.2"
    
    如果 $F$ 是 $\mathbb{R}$ 上的一个连续分布函数，那么 $F$ 在 $\mathbb{R}$ 上是一致连续的。

基于此，我们可以得到强化依分布收敛的一个极其优美的结论：

!!! success "定理 1.3：Polya 定理 (Polya's Theorem)"

    假设 $X_n \xrightarrow{d} X$，且极限随机变量 $X$ 的分布函数 $F(x)$ 是**连续**的。
    那么，这种逐点收敛会自动升级为**一致收敛**：
    
    $$
    \sup_{x \in \mathbb{R}} |F_n(x) - F(x)| \rightarrow 0 \quad \text{as } n \rightarrow \infty
    $$
    
    ??? proof "Polya 定理详细证明（点击展开）"
        
        我们需要证明对于任意 $\epsilon > 0$，当 $n$ 足够大时，$\sup_x |F_n(x) - F(x)| < \epsilon$。

        **1. 构造有限划分（Partitioning）**
        
        由于 $F(x)$ 是连续分布函数，其值域为 $[0, 1]$。对于给定的 $\epsilon > 0$，我们可以找到有限个点 $-\infty = x_0 < x_1 < x_2 < \dots < x_K = \infty$，使得：
        
        $$
        F(x_i) - F(x_{i-1}) < \frac{\epsilon}{2}, \quad \forall i = 1, \dots, K
        $$
        
        （注：这里规定 $F(x_0) = 0$ 且 $F(x_K) = 1$）。

        **2. 利用点点收敛（Pointwise Convergence）**
        
        因为 $X_n \xrightarrow{d} X$，对于上述每一个有限的网格点 $x_i$（当 $1 \le i \le K-1$ 时），根据分布收敛的定义，当 $n \to \infty$ 时：
        
        $$
        F_n(x_i) \rightarrow F(x_i)
        $$
        
        由于网格点是有限的，必存在 $N$，使得当 $n > N$ 时，对所有 $i=1, \dots, K-1$ 均有：
        
        $$
        |F_n(x_i) - F(x_i)| < \frac{\epsilon}{2}
        $$

        **3. 单调性夹逼（The Sandwich Argument）**
        
        对于轴上任意一点 $x$，它必然落在某个区间 $[x_{i-1}, x_i]$ 内。利用 $F_n$ 和 $F$ 的单调不减性：
        
        * **上界：**
        
        $$
        F_n(x) - F(x) \le F_n(x_i) - F(x_{i-1}) = [F_n(x_i) - F(x_i)] + [F(x_i) - F(x_{i-1})]
        $$
        
        当 $n > N$ 时，代入前两步的结论：
        
        $$
        F_n(x) - F(x) < \frac{\epsilon}{2} + \frac{\epsilon}{2} = \epsilon
        $$
        
        * **下界：**
        
        $$
        F_n(x) - F(x) \ge F_n(x_{i-1}) - F(x_i) = [F_n(x_{i-1}) - F(x_{i-1})] - [F(x_i) - F(x_{i-1})]
        $$
        
        当 $n > N$ 时，同理可得：
        
        $$
        F_n(x) - F(x) > -\frac{\epsilon}{2} - \frac{\epsilon}{2} = -\epsilon
        $$

        **4. 结论**
        
        综合上下界，对于所有 $x \in \mathbb{R}$，只要 $n > N$，就有：
        
        $$
        |F_n(x) - F(x)| < \epsilon
        $$
        
        由此证明了：
        
        $$
        \sup_{x \in \mathbb{R}} |F_n(x) - F(x)| \rightarrow 0 \quad (n \rightarrow \infty)
        $$
        
        $\square$

统计学中最常见的一种收敛便是向正态分布的收敛：

!!! abstract "定义 1.4 & 1.5：渐近正态性 (Asymptotic Normality, AN)"

    **1. 一维渐近正态**：
    序列 $\{X_n\}$ 称为具有“均值” $\mu_n$ 和“方差” $\sigma_n^2 > 0$ 的渐近正态分布，记作 $X_n \sim AN(\mu_n, \sigma_n^2)$，如果当 $n$ 足够大时满足：
    
    $$
    \frac{X_n - \mu_n}{\sigma_n} \xrightarrow{d} N(0, 1)
    $$
    
    *(注：这里的 $\mu_n$ 和 $\sigma_n^2$ 不一定真的是 $X_n$ 的均值和方差，有时 $X_n$ 的真实矩甚至不存在！)*
    
    **2. 多维渐近正态**：
    对于随机向量序列 $X_n$，称其服从 $AN(\mu_n, \Sigma_n)$，如果对于任意非零向量 $a \in \mathbb{R}^p$，通过一维化投影均满足：
    
    $$
    a'X_n \sim AN(a'\mu_n, a'\Sigma_n a)
    $$

---

## 3. 随机收敛的基石：Portmanteau 引理

依分布收敛除了通过 CDF 定义，还可以通过期望、开集/闭集等多种拓扑方式等价描述。这构成了渐近理论中最核心的工具。

!!! info "定理 1.6：Portmanteau 综合引理 (Portmanteau Lemma)"

    对于任意随机向量 $X_n$ 和 $X$，以下陈述是**完全等价**的：

    **(i)** $X_n \xrightarrow{d} X$;
    
    **(ii)** 对于任意**有界连续函数** $f \in C_B$，$E[f(X_n)] \to E[f(X)]$;
    
    **(iii)** 对于任意**有界 Lipschitz 连续函数** $f \in C_{B, Lip}$，$E[f(X_n)] \to E[f(X)]$;
    
    **(iv)** 对于任意**非负连续函数** $f$，$\liminf E[f(X_n)] \ge E[f(X)]$;
    
    **(v)** 对于任意**开集** $G$，$\liminf P(X_n \in G) \ge P(X \in G)$;
    
    **(vi)** 对于任意**闭集** $F$，$\limsup P(X_n \in F) \le P(X \in F)$;
    
    **(vii)** 对于任意**边界测度为 0** 的 Borel 集 $B$（即 $P(X \in \partial B) = 0$），$P(X_n \in B) \to P(X \in B)$。

    ??? proof "Portmanteau 引理核心步骤的严格推导（点击展开）"

        **证明 (i) $\Rightarrow$ (ii)**：
        不失一般性，假设 $\sup |f(x)| \le 1$。对于任意 $\epsilon > 0$，选取一个足够大的矩形区域 $I$，使得尾部概率 $P(X \in I^c) < \epsilon$。
        将 $I$ 划分为有限个互不重叠的小矩形 $I = \cup_{j=1}^K I_j$，并在每个小矩形内取代表点 $x_j$。构造简单阶梯函数：
        
        $$
        f_\epsilon(x) = \sum_{j=1}^K f(x_j) \mathbb{I}(x \in I_j)
        $$
        
        由于 $f$ 在有界闭集 $I$ 上一致连续，通过足够细的划分，可以保证在 $I$ 内 $|f(x) - f_\epsilon(x)| < \epsilon$。同时，根据构造有 $\sup |f_\epsilon| \le \sup |f| \le 1$。
        
        利用指示函数，我们将全空间的期望根据 $X_n$ 是否落入区域 $I$ 进行**分解**：
        
        $$
        \begin{aligned}
        |E[f(X_n)] - E[f_\epsilon(X_n)]| &\le E\left[ |f(X_n) - f_\epsilon(X_n)| \right] \\
        &= E\left[ |f(X_n) - f_\epsilon(X_n)| \cdot \mathbb{I}(X_n \in I) \right] + E\left[ |f(X_n) - f_\epsilon(X_n)| \cdot \mathbb{I}(X_n \in I^c) \right]
        \end{aligned}
        $$
        
        对于第一项，由于在 $I$ 内误差受 $\epsilon$ 控制：
        
        $$
        E\left[ |f(X_n) - f_\epsilon(X_n)| \cdot \mathbb{I}(X_n \in I) \right] < \epsilon \cdot P(X_n \in I) \le \epsilon
        $$
        
        对于第二项，利用 $\sup |f| \le 1$ 和 $\sup |f_\epsilon| \le 1$，可知 $|f - f_\epsilon| \le 2$：
        
        $$
        E\left[ |f(X_n) - f_\epsilon(X_n)| \cdot \mathbb{I}(X_n \in I^c) \right] \le 2 \cdot P(X_n \in I^c)
        $$
        
        将两部分合并，即得到该期望误差的上界：
        
        $$
        |E[f(X_n)] - E[f_\epsilon(X_n)]| \le \epsilon + 2P(X_n \in I^c)
        $$
        
        同理对于极限变量：
        
        $$
        |Ef(X) - Ef_\epsilon(X)| \le \epsilon + 2P(X \in I^c) < 3\epsilon
        $$
        
        而对于简单函数部分：
        
        $$
        |Ef_\epsilon(X_n) - Ef_\epsilon(X)| \le \sum_{j=1}^K |f(x_j)| |P(X_n \in I_j) - P(X \in I_j)| \to 0
        $$
        
        由于 $K$ 是有限的，且每个 $I_j$ 是一个连续集，联合这三项即可得证 $E[f(X_n)] \to E[f(X)]$。

        **证明 (iii) $\Rightarrow$ (v)**：
        对于任意开集 $G$，我们构造一列非负的 Lipschitz 函数来逼近其指示函数：设 $f_m(x) = (m \cdot d(x, G^c)) \wedge 1$。
        当 $m \to \infty$ 时，$f_m \uparrow \mathbb{I}_G$。对于固定的 $m$：
        
        $$
        \liminf_{n \to \infty} P(X_n \in G) \ge \liminf_{n \to \infty} E[f_m(X_n)] = E[f_m(X)]
        $$
        
        由单调收敛定理 (Monotone Convergence Theorem)，令 $m \to \infty$，右侧单调增加至 $P(X \in G)$。

        **证明 (v) $\Leftrightarrow$ (vi)**：
        利用开集和闭集的补集对应关系（De Morgan 定律），直接取补集反转不等式方向即可。

        **证明 (v) + (vi) $\Rightarrow$ (vii)**：
        设 $B$ 的内部为 $B^\circ$，闭包为 $\overline{B}$。利用前两个性质：
        
        $$
        P(X \in B^\circ) \le \liminf P(X_n \in B) \le \limsup P(X_n \in B) \le P(X \in \overline{B})
        $$
        
        因为已知边界测度为零，即 $P(X \in \partial B) = 0$，所以 $P(X \in B^\circ) = P(X \in \overline{B})$。由此夹逼出中间的极限存在且等于 $P(X \in B)$。$\square$

        **证明 (vii) $\Rightarrow$ (i)**：
        对于任意实数 $x$，构造一个左侧无限的闭区间 $B = (-\infty, x]$。该集合的边界仅为单点集 $\partial B = \{x\}$。
        如果 $x$ 是累积分布函数 $F(x) = P(X \le x)$ 的连续点，那么在这一点上的概率测度为零，即 $P(X \in \partial B) = P(X = x) = 0$。
        既然该集合的边界测度为零，由条件 (vii) 可知：
        
        $$
        P(X_n \le x) = P(X_n \in B) \to P(X \in B) = P(X \le x)
        $$
        
        由于这个等式对于所有 $F(x)$ 的连续点 $x$ 都成立，这恰好就是依分布收敛 $X_n \xrightarrow{d} X$ 的严格定义。$\square$


!!! success "定理 1.6 补充：Lévy 连续性定理 (Lévy's Continuity Theorem)"

    除了 Portmanteau 引理给出的几种等价拓扑条件外，依分布收敛还有一个极其重要且极具计算价值的等价刻画，即**特征函数 (Characteristic Function)** 的点态收敛。
    
    设 $\{X_n\}$ 和 $X$ 为 $\mathbb{R}^d$ 中的随机向量，$\phi_{X_n}(t)$ 和 $\phi_X(t)$ 分别为它们的特征函数（定义为 $\phi_X(t) = E[e^{i t^\top X}]$）。那么：
    
    $$
    X_n \xrightarrow{d} X \iff \phi_{X_n}(t) \to \phi_X(t), \quad \forall t \in \mathbb{R}^d
    $$
    
    > *(注：特征函数的收敛是证明中心极限定理 (CLT) 等渐近分布时最常用的工具！)*

---

## 4. 连续映射定理 (Continuous Mapping Theorem, CMT)

如果一个随机变量序列是收敛的，那么当它们经过一个“足够好”的函数映射后，收敛性质是否依然保持？映射定理给出了肯定的回答。

!!! success "定理 1.7：连续映射定理 (Mapping Theorem)"

    设函数 $g: \mathbb{R}^k \to \mathbb{R}^m$ 在连续点集 $\mathcal{C}_g$ 上连续，且满足 $P(X \in \mathcal{C}_g) = 1$（即 $X$ 几乎必然落在 $g$ 的连续点上）。
    那么，操作算子 $g(\cdot)$ 会完美地继承并传递以下三种收敛性：
    
    1. 若 $X_n \xrightarrow{a.s.} X$，则 $g(X_n) \xrightarrow{a.s.} g(X)$
    2. 若 $X_n \xrightarrow{P} X$，则 $g(X_n) \xrightarrow{P} g(X)$
    3. 若 $X_n \xrightarrow{d} X$，则 $g(X_n) \xrightarrow{d} g(X)$

    ??? proof "映射定理的严格证明（点击展开）"

        我们在此重点证明**依分布收敛**的情况 $X_n \xrightarrow{d} X$。我们将利用极其巧妙的 Portmanteau 引理的闭集性质 (vi) 来证明。
        
        对于任意闭集 $F$，考虑其原像 $g^{-1}(F) = \{x : g(x) \in F\}$。
        由于 $g$ 并非处处连续，我们需要分析原像闭包 $\overline{g^{-1}(F)}$ 的结构：
        
        $$
        g^{-1}(F) \subset \overline{g^{-1}(F)} \subset g^{-1}(F) \cup \mathcal{C}_g^c
        $$
        
        *(解释：如果一个极限点 $x$ 是连续点，即 $x \in \mathcal{C}_g$，那么序列 $x_m \to x$ 必有 $g(x_m) \to g(x)$。由于 $F$ 是闭集，自然有 $g(x) \in F$，故 $x \in g^{-1}(F)$)。*
        
        对该集合运用 Portmanteau 引理 (vi)：
        
        $$
        \limsup P(g(X_n) \in F) \le \limsup P\left(X_n \in \overline{g^{-1}(F)}\right)
        $$
        
        由引理性质，上式小于等于极限在闭包上的概率：
        
        $$
        \le P\left(X \in \overline{g^{-1}(F)}\right)
        $$
        
        将其拆分为连续点和非连续点两部分：
        
        $$
        \le P(X \in g^{-1}(F)) + P(X \notin \mathcal{C}_g)
        $$
        
        根据定理前提，$P(X \notin \mathcal{C}_g) = 0$。因此：
        
        $$
        \limsup P(g(X_n) \in F) \le P(g(X) \in F)
        $$
        
        再次使用 Portmanteau 引理的反向推导 (vi) $\Rightarrow$ (i)，立刻得证 $g(X_n) \xrightarrow{d} g(X)$。$\square$

        **证明 (ii) 依概率收敛的映射性质**：
        我们需要证明：对于任意给定的 $\epsilon > 0$，$P(|g(X_n) - g(X)| > \epsilon) \to 0$。
        
        固定任意 $\epsilon > 0$。对于任意 $\delta > 0$，我们定义一个“坏集” $B_\delta$，它包含了所有可能使得函数值发生剧烈突变的 $x$ 点：
        
        $$
        B_\delta = \{x : \exists y \text{ 满足 } |x - y| < \delta \text{ 但 } |g(x) - g(y)| > \epsilon\}
        $$
        
        现在考察事件 $\{|g(X_n) - g(X)| > \epsilon\}$。如果该事件发生，且极限变量 $X$ 恰好不在“坏集” $B_\delta$ 中（即 $X \notin B_\delta$），那么必然是因为 $|X_n - X| \ge \delta$。
        利用全概率放缩，我们可以得到：
        
        $$
        P(|g(X_n) - g(X)| > \epsilon) \le P(X \in B_\delta) + P(|X_n - X| \ge \delta)
        $$
        
        接下来对右边两项分别取极限：
        **第一项**：由于 $g$ 在 $\mathcal{C}_g$ 上连续，当 $\delta \downarrow 0$ 时，集合 $B_\delta$ 与连续点集 $\mathcal{C}_g$ 的交集必然为空集。又因为已知 $P(X \in \mathcal{C}_g) = 1$，所以当 $\delta \to 0$ 时，$P(X \in B_\delta) \to 0$。
        **第二项**：对于任何固定的 $\delta > 0$，由于已知 $X_n \xrightarrow{P} X$，当 $n \to \infty$ 时，$P(|X_n - X| \ge \delta) \to 0$。
        
        综合两项，令 $n \to \infty$ 再令 $\delta \downarrow 0$，即可得证原概率趋于 0，即 $g(X_n) \xrightarrow{P} g(X)$。$\square$

!!! tip "映射定理的经典应用示例 (Applications)"

    Mapping Theorem 在推导复杂统计量的渐近分布时堪称“神兵利器”：

    1. **卡方分布的引出**：
       若一维序列 $X_n \xrightarrow{d} X \sim N(0,1)$，取连续映射 $g(x) = x^2$，则立刻得到 $X_n^2 \xrightarrow{d} \chi_1^2$。
       
    2. **柯西分布的引出**：
       若二维序列 $(X_n, Y_n)^\top \xrightarrow{d} N_2(0, I_2)$，取映射 $g(x,y) = x/y$（在 $y=0$ 处不连续，但标准正态分布下 $P(Y=0)=0$，满足几乎必然连续条件），则 $X_n/Y_n \xrightarrow{d} Cauchy$。
       
    3. **样本方差的依概率收敛**：
       由大数定律，$(\overline{X}, \frac{1}{n}\sum X_i^2)^\top \xrightarrow{P} (\mu, \mu_2)^\top$。取连续函数 $g(x,y) = y - x^2$，直接可得样本方差 $S_n^2 = g(\overline{X}, \frac{1}{n}\sum X_i^2) \xrightarrow{P} \mu_2 - \mu^2 = \sigma^2$。
       
    4. **多元正态的仿射变换**：
       若 $X_n \xrightarrow{d} N_p(\mu, \Sigma)$，对于任意常数矩阵 $C \in \mathbb{R}^{m \times p}$，有 $C X_n \xrightarrow{d} N_m(C\mu, C\Sigma C^\top)$。

## 5. 随机收敛的相互关系 (Relations of Stochastic Convergence)

四种随机收敛之间存在着严格的强弱蕴含关系。在实际应用中，我们经常需要利用这些关系在依概率收敛和依分布收敛之间进行转换。

!!! success "定理 1.8：随机收敛关系定理"

    设 $X_n$、$X$ 和 $Y_n$ 为随机向量。那么以下推论成立：

    1. **几乎必然 $\Rightarrow$ 依概率**：若 $X_n \xrightarrow{a.s.} X$，则 $X_n \xrightarrow{P} X$；
    2. **依概率 $\Rightarrow$ 依分布**：若 $X_n \xrightarrow{P} X$，则 $X_n \xrightarrow{d} X$；
    3. **向常数收敛的等价性**：$X_n \xrightarrow{P} c$（$c$ 为常数） 当且仅当 $X_n \xrightarrow{d} c$；
    4. **距离收敛传递 (Convergence Lemma)**：若 $X_n \xrightarrow{d} X$ 且 $d(X_n, Y_n) \xrightarrow{P} 0$，则 $Y_n \xrightarrow{d} X$；
    5. **联合分布收敛 (1)**：若 $X_n \xrightarrow{d} X$ 且 $Y_n \xrightarrow{P} c$（常数），则 $(X_n, Y_n) \xrightarrow{d} (X, c)$；
    6. **联合分布收敛 (2)**：若 $X_n \xrightarrow{P} X$ 且 $Y_n \xrightarrow{P} Y$，则 $(X_n, Y_n) \xrightarrow{P} (X, Y)$。

    > **重要注记**：单一边缘的依概率收敛可以推出联合的依概率收敛。但是，**单一边缘的依分布收敛，通常推不出联合的依分布收敛**（除非用到 Copula 或者其中一个是向常数收敛，如性质 5）。

    ??? proof "核心性质的严格证明（点击展开）"

        **(1) 几乎必然 $\Rightarrow$ 依概率**：
        定义事件序列 $A_n = \cup_{m \ge n} \{||X_m - X|| > \epsilon\}$。该集合序列是单调递减的。
        如果对于所有的 $\omega \in \Omega$ 都有 $X_n(\omega) \to X(\omega)$，那么当 $n \to \infty$ 时 $A_n$ 递减趋于空集。
        若 $X_n \xrightarrow{a.s.} X$，由概率的连续性：
        
        $$
        P(||X_n - X|| > \epsilon) \le P(A_n) \to 0
        $$

        **(4) 距离收敛传递 (证明的基础)**：
        我们利用 Portmanteau 引理 (iii) 证明。对于任意有界 Lipschitz 连续函数 $f \in C_{B, Lip}$，设其 Lipschitz 常数为 $L$，且 $\sup |f| \le M$。
        考察期望差值的绝对值，通过引入指示函数截断：
        
        $$
        \begin{aligned}
        |Ef(X_n) - Ef(Y_n)| &\le E\left[ |f(X_n) - f(Y_n)| \cdot \mathbb{I}(||X_n - Y_n|| \le \epsilon) \right] \\
        &\quad + E\left[ |f(X_n) - f(Y_n)| \cdot \mathbb{I}(||X_n - Y_n|| > \epsilon) \right]
        \end{aligned}
        $$
        
        对于第一部分，利用 Lipschitz 性质放缩为 $L \epsilon \cdot P(||X_n - Y_n|| \le \epsilon) \le L \epsilon$；
        对于第二部分，利用有界性放缩为 $2M \cdot P(||X_n - Y_n|| > \epsilon)$。
        因此：
        
        $$
        |Ef(X_n) - Ef(Y_n)| \le L\epsilon + 2M \cdot P(||X_n - Y_n|| > \epsilon)
        $$
        
        由于 $d(X_n, Y_n) \xrightarrow{P} 0$，当 $n \to \infty$ 时第二项趋于 0。因为 $\epsilon$ 是任意小的，所以 $Ef(X_n) - Ef(Y_n) \to 0$。
        结合前提 $X_n \xrightarrow{d} X$（即 $Ef(X_n) \to Ef(X)$），必然有 $Ef(Y_n) \to Ef(X)$。再次由 Portmanteau 引理，得证 $Y_n \xrightarrow{d} X$。

        **(2) 依概率 $\Rightarrow$ 依分布**：
        将 $X_n$ 拆分为 $X_n = X + (X_n - X)$。
        由于显然有 $X \xrightarrow{d} X$，并且已知 $X_n - X \xrightarrow{P} 0$，我们直接套用刚刚证明的性质 (4)（令性质 (4) 中的 $X_n$ 角色为 $X$，$Y_n$ 角色为 $X_n$），立刻得到 $X_n \xrightarrow{d} X$。

        **(3) 向常数收敛等价性**：
        充分性由性质 (2) 保证。证明必要性 ($X_n \xrightarrow{d} c \Rightarrow X_n \xrightarrow{P} c$)：
        事件 $\{||X_n - c|| \ge \epsilon\}$ 等价于 $X_n \in B(c, \epsilon)^c$。这是一个闭集。
        由 Portmanteau 引理 (vi)：
        
        $$
        \limsup P(||X_n - c|| \ge \epsilon) \le P(c \in B(c, \epsilon)^c) = 0
        $$
        
        因此极限为 0，依概率收敛得证。

        **(5) $(X_n, Y_n) \xrightarrow{d} (X, c)$**：
        将联合变量拆解为 $(X_n, Y_n) = (X_n, c) + (0, Y_n - c)$。
        由于 $Y_n \xrightarrow{P} c$，误差项 $(0, Y_n - c) \xrightarrow{P} (0,0)$。
        由性质 (4)，我们只需证明主项 $(X_n, c) \xrightarrow{d} (X, c)$。对于任意二维连续有界函数 $f(x,y)$，固定 $y=c$ 后边缘函数 $f_m(x) = f(x,c)$ 也是连续有界的。
        因为 $X_n \xrightarrow{d} X$，所以 $E[f(X_n, c)] = E[f_m(X_n)] \to E[f_m(X)] = E[f(X, c)]$。得证！
        
        **(6) $(X_n, Y_n) \xrightarrow{P} (X, Y)$**：
        由三角不等式 $||(X_n, Y_n) - (X, Y)|| \le ||X_n - X|| + ||Y_n - Y||$。
        
        $$
        P(||(X_n, Y_n) - (X, Y)|| > \epsilon) \le P(||X_n - X|| > \epsilon/2) + P(||Y_n - Y|| > \epsilon/2) \to 0
        $$
        
        得证。$\square$

---

## 6. Slutsky 定理 (Slutsky's Theorem)

作为定理 1.8 和连续映射定理 (CMT) 的直接推论，Slutsky 定理为统计量之间的代数运算提供了极为便利的准则。它是构建 $t$ 统计量、Wald 统计量等渐近分布的基石。

!!! success "定理 1.9：Slutsky 定理 (Slutsky, 1925)"

    设 $X_n, X, Y_n$ 为随机向量或标量随机变量。
    如果 $X_n \xrightarrow{d} X$ 且 $Y_n \xrightarrow{P} c$（常数），那么：

    1. **加法法则**：$X_n + Y_n \xrightarrow{d} X + c$
    2. **乘法法则**：$Y_n X_n \xrightarrow{d} cX$
    3. **除法法则**：$Y_n^{-1} X_n \xrightarrow{d} c^{-1}X$ （前提是 $c \neq 0$ 且对于矩阵而言可逆）

    ??? proof "推导极简证明（点击展开）"

        根据定理 1.8 的性质 (5)，已知 $X_n \xrightarrow{d} X$ 且 $Y_n \xrightarrow{P} c$，可以直接推导出它们的联合分布收敛：
        
        $$
        (X_n, Y_n) \xrightarrow{d} (X, c)
        $$
        
        接下来，分别构造二元连续映射函数 $g(x,y) = x+y$、$g(x,y) = yx$ 以及 $g(x,y) = y^{-1}x$。
        直接套用连续映射定理 (CMT) $g(X_n, Y_n) \xrightarrow{d} g(X, c)$，Slutsky 定理的三条法则立刻得证。$\square$

!!! tip "经典应用：$t$-统计量的渐近正态性"

    设 $Y_1, \dots, Y_n$ 是 i.i.d. 的样本，满足 $E[Y_1] = 0, E[Y_1^2] = \sigma^2$。
    我们要推导 $t$-统计量 $t_n := \frac{\sqrt{n}\overline{Y}}{S_n}$ 的极限分布，其中 $S_n^2$ 是样本方差。

    **步骤 1：分析分子**
    由中心极限定理 (CLT)，均值的渐近分布为：
    
    $$
    \sqrt{n}\overline{Y} \xrightarrow{d} N(0, \sigma^2)
    $$
    
    **步骤 2：分析分母**
    由大数定律，样本方差依概率收敛于总体方差：$S_n^2 \xrightarrow{P} \sigma^2$。再由映射定理 ($g(x)=\sqrt{x}$)，得到 $S_n \xrightarrow{P} \sigma$。

    **步骤 3：应用 Slutsky 定理**
    将分子视为 $X_n$，分母视为 $Y_n^{-1}$：
    
    $$
    \frac{\sqrt{n}\overline{Y}}{S_n} \xrightarrow{d} \frac{1}{\sigma} N(0, \sigma^2) \stackrel{d}{=} N(0, 1)
    $$
    
    这从理论上证明了在大样本下，$t$-检验可以近似使用标准正态分布的临界值。

---

## 7. 胎紧性与随机有界 (Tightness & Stochastic Boundedness)

在研究一列随机变量是否收敛时，“它们是否会跑到无穷远处逃逸掉？”是一个核心问题。这就引出了胎紧性（Tightness）的概念。

!!! abstract "定义 1.10：随机有界 (Stochastically Bounded) / 胎紧 (Tight)"

    序列 $\{X_n\}$ 被称为是**随机有界**的（或**胎紧的**），如果对于任意给定的 $\epsilon > 0$，都存在一个有限的常数 $M_\epsilon > 0$，使得对于所有的 $n$：
    
    $$
    \sup_n P(||X_n|| > M_\epsilon) < \epsilon
    $$
    
    在渐近统计中，我们通常将这种性质记为大 O 符号：**$X_n = O_p(1)$**。

一个单独的随机变量 $X$ 本身必定是胎紧的（因为分布函数 $F(\infty)=1, F(-\infty)=0$）。进而，任何**有限个**随机变量的集合也是胎紧的。真正需要警惕的是无限序列 $\{X_n\}$ 的逃逸。

!!! info "定理 1.11：Prohorov 定理 (Prohorov's Theorem)"

    依分布收敛与胎紧性之间有着极其深刻的拓扑联系：

    1. 若 $X_n \xrightarrow{d} X$，那么 $\{X_n\}$ **必然是胎紧的**（即收敛必有界）。
    2. 若 $\{X_n\}$ 是胎紧的，那么它**必然存在一个依分布收敛的子列** $\{X_{n_i}\}$，使得当 $n_i \to \infty$ 时，$X_{n_i} \xrightarrow{d} X$。（这可以看作是实分析中“有界数列必有收敛子列”的 Bolzano-Weierstrass 定理在概率空间上的推广）。

    ??? proof "定理 1.11 (1) 的严格证明（点击展开）"

        因为极限随机变量 $X$ 单独是一个随机变量，所以它是胎紧的。对于任意 $\epsilon > 0$，我们可以找到一个足够大的 $M_\epsilon$，且满足 $P(||X|| = M_\epsilon) = 0$（避开间断点），使得：
        
        $$
        P(||X|| \ge M_\epsilon) < \epsilon
        $$
        
        考虑闭集 $F = \{x : ||x|| \ge M_\epsilon\}$。由 Portmanteau 引理 (vi) 和 $X_n \xrightarrow{d} X$：
        
        $$
        \limsup_{n \to \infty} P(||X_n|| \ge M_\epsilon) \le P(||X|| \ge M_\epsilon) < \epsilon
        $$
        
        既然上极限严格小于 $\epsilon$，那么必定存在一个正整数 $N$，使得对于所有的 $n \ge N$：
        
        $$
        P(||X_n|| \ge M_\epsilon) \le 2\epsilon
        $$
        
        对于前面有限的 $N-1$ 个随机变量 $\{X_1, \dots, X_{N-1}\}$，由于有限集必然胎紧，我们可以适当放大 $M_\epsilon$，使得对于所有的 $n \in \mathbb{N}_+$ 都有 $P(||X_n|| > M_\epsilon) < \epsilon$ 成立。故序列全局胎紧。$\square$

---

## 8. 随机阶符号 $o_p$ 与 $O_p$ (Stochastic Order)

为了在概率论中进行类似微积分的泰勒展开，我们需要一套描述“概率意义下无穷小”和“概率意义下同阶”的符号系统。

!!! abstract "定义 1.12：随机小 $o_p$ 与大 $O_p$"

    设 $\{a_n\}$ 为一个常数序列。对于随机向量序列 $\{X_n\}$：

    * **大 $O_p$ (Stochastically of order $a_n$)**：若 $\frac{X_n}{a_n} = O_p(1)$，即 $\frac{X_n}{a_n}$ 是胎紧的（随机有界的），则记作 $X_n = O_p(a_n)$。
    * **小 $o_p$ (Stochastically of smaller order $a_n$)**：若 $\frac{X_n}{a_n} \xrightarrow{P} 0$，即比 $a_n$ 衰减得更快，则记作 $X_n = o_p(a_n)$。
    
    > *特例*：$X_n = o_p(1)$ 意味着 $X_n \xrightarrow{P} 0$。

在确定统计量的阶时，切比雪夫不等式 (Chebyshev's Inequality) 是最常用的手段。例如，如果 $E[T_n] = \mu_n$ 且 $Var(T_n) = \sigma_n^2$，则必定有 $T_n - \mu_n = O_p(\sigma_n)$。

**$o_p$ 与 $O_p$ 的运算准则：**
在进行渐近展开时，我们可以像处理确定性极限一样处理这些随机符号（从左向右阅读）：

* $o_p(1) + o_p(1) = o_p(1)$
* $O_p(1) + o_p(1) = O_p(1)$
* $O_p(1) \cdot o_p(1) = o_p(1)$
* $(1 + o_p(1))^{-1} = O_p(1)$
* $O_p(1) + O_p(1) = O_p(1)$
* $o_p(O_p(1)) = O_p(o_p(1)) = o_p(1)$

---

## 9. 随机代入引理 (Lemma of Stochastic Plug-in)

这是我们在统计推断中推导 Delta 方法 (Delta Method) 和极大似然估计 (MLE) 渐近正态性的核心预备引理。

!!! success "引理 1.13：随机代入引理"

    设 $R: \mathbb{R}^k \to \mathbb{R}$ 是一个满足 $R(0)=0$ 的实函数。设 $\{X_n\}$ 是一列以 0 为概率极限的随机向量（即 $X_n \xrightarrow{P} 0$）。那么对于任意 $p > 0$：

    1. 若当 $h \to 0$ 时有确定性极限 $R(h) = o(||h||^p)$，则 $R(X_n) = o_p(||X_n||^p)$。
    2. 若当 $h \to 0$ 时有确定性极限 $R(h) = O(||h||^p)$，则 $R(X_n) = O_p(||X_n||^p)$。

    > *注记：这个引理非常强大，它允许我们将确定性微积分中的泰勒展开余项 $o(|x|)$，直接无缝转化为随机变量下的 $o_p(|X_n|)$，而无需担心函数 $R$ 在其他地方是否连续。*

    ??? proof "小 $o_p$ 和 大 $O_p$ 截断证明（点击展开）"

        **证明 (1) 小 $o_p$ 情况**：
        定义辅助函数：
        
        $$
        g(h) = \begin{cases} \frac{R(h)}{||h||^p}, & h \neq 0 \\ 0, & h = 0 \end{cases}
        $$
        
        则可以写成 $R(X_n) = ||X_n||^p g(X_n)$。要证 $R(X_n) = o_p(||X_n||^p)$，只需证明 $g(X_n) \xrightarrow{P} 0$。
        由于前提 $R(h) = o(||h||^p)$，所以当 $h \to 0$ 时 $g(h) \to 0 = g(0)$。即 $g(h)$ 在 0 点是连续的。
        既然已知 $X_n \xrightarrow{P} 0$，通过连续映射定理 (CMT)，立刻得到 $g(X_n) \xrightarrow{P} g(0) = 0$。得证！

        **证明 (2) 大 $O_p$ 情况**：
        同样利用上面的 $g(h)$。由于 $R(h) = O(||h||^p)$，在 $h=0$ 附近 $g(h)$ 是有界的。
        即存在 $M > 0$ 和 $\delta > 0$，使得当 $||h|| < \delta$ 时，$|g(h)| \le M$。
        现在考察事件 $\{|g(X_n)| > M\}$。如果这个事件发生，说明 $X_n$ 必定逃出了半径为 $\delta$ 的安全区域。因此存在集合包含关系：
        
        $$
        \{\omega : |g(X_n(\omega))| > M\} \subset \{\omega : ||X_n(\omega)|| > \delta\}
        $$
        
        取概率测度：
        
        $$
        P(|g(X_n)| > M) \le P(||X_n|| > \delta)
        $$
        
        因为已知 $X_n \xrightarrow{P} 0$，对于任意给定的 $\epsilon > 0$，当 $n$ 足够大时，右侧的 $P(||X_n|| > \delta) < \epsilon$。
        所以：
        
        $$
        P(|g(X_n)| > M) < \epsilon, \quad \forall n \ge N
        $$
        
        这就完全满足了随机有界 $O_p(1)$ 的定义，即 $g(X_n) = O_p(1)$，从而 $R(X_n) = O_p(||X_n||^p)$。$\square$

!!! tip "综合应用示例：样本方差的渐近分布"

    假设 $X_1, \dots, X_n$ i.i.d. $\sim F(\mu, \sigma^2)$，且四阶矩存在 $E[X^4] < \infty$。推导不偏方差 $S_n^2$ 的极限分布。
    
    首先将方差展开：
    
    $$
    S_n^2 = \frac{1}{n}\sum_{i=1}^n (X_i - \overline{X})^2 = \frac{1}{n}\sum_{i=1}^n (X_i - \mu)^2 - (\overline{X} - \mu)^2
    $$
    
    对于第一部分，由于四阶矩存在，利用标准中心极限定理：
    
    $$
    \sqrt{n}\left( \frac{1}{n}\sum_{i=1}^n (X_i - \mu)^2 - \sigma^2 \right) \xrightarrow{d} N(0, Var((X_i - \mu)^2)) := N(0, \nu^2)
    $$
    
    对于第二部分，由 CLT 知 $\overline{X} - \mu = O_p(n^{-1/2})$。利用随机代入引理（取 $p=2$），平方后即为：
    
    $$
    (\overline{X} - \mu)^2 = O_p(n^{-1}) = o_p(n^{-1/2})
    $$
    
    因此，当乘以 $\sqrt{n}$ 后：
    
    $$
    \sqrt{n}(\overline{X} - \mu)^2 = o_p(1)
    $$
    
    利用 Slutsky 定理，高阶误差项在极限下消失：
    
    $$
    \sqrt{n}(S_n^2 - \sigma^2) = \sqrt{n}\left( \frac{1}{n}\sum_{i=1}^n (X_i - \mu)^2 - \sigma^2 \right) - \sqrt{n}(\overline{X} - \mu)^2 \xrightarrow{d} N(0, \nu^2) - 0 = N(0, \nu^2)
    $$