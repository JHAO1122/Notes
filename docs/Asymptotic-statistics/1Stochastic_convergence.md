# 第一章：数据分布与随机收敛 (Data Distributions and Stochastic Convergence)

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
    \sup_x |F_n(x) - F(x)| \rightarrow 0 \quad \text{as } n \rightarrow \infty
    $$
    
    ??? proof "Polya 定理证明思路（点击展开）"
        
        该证明采用了经典的**有限覆盖法 (Covering Method)**。
        
        由于 $F$ 是一致连续的分布函数，我们在无穷区间上截取一个足够大的 $[-M, M]$，使得尾部概率极其微小。
        将 $[-M, M]$ 进行等宽划分，构造网格点 $\{x_i\}_{i=0}^{K+1}$。
        
        定义在这些网格点上的最大误差为：
        
        $$
        \Delta_n = \max_{i=0}^{K+1} |F_n(x_i) - F(x_i)|
        $$
        
        因为 $X_n \xrightarrow{d} X$ 且 $F$ 无间断点，对于这**有限个**网格点 $K$，当 $n \to \infty$ 时必然有 $\Delta_n \to 0$。利用单调性与三角不等式，即可将网格点上的收敛推广到全空间的 $\sup$ 一致收敛。$\square$

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
        
        通过足够细的划分，保证在 $I$ 内 $|f(x) - f_\epsilon(x)| < \epsilon$。
        利用三角不等式拆解期望误差：
        
        $$
        |Ef(X_n) - Ef_\epsilon(X_n)| \le \epsilon + 2P(X_n \in I^c)
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