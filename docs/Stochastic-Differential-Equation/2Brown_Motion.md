# 第二章：布朗运动 (Brownian Motion)

在有了条件期望与鞅论的测度论基础后，我们正式引入一个经典的连续时间随机过程——**布朗运动 (Brownian Motion, BM)**，也被称为维纳过程 (Wiener Process)。它是构建随机微分方程 (SDE) 积分理论（如 Itô 积分）的基础。

## 1. 基础定义与性质

!!! abstract "定义：标准布朗运动 (Standard Brownian Motion)"

    设 $(\Omega, \mathcal{F}, P)$ 为一个概率空间，其上定义了一个实值随机过程 $W = \{W(t), t \ge 0\}$。如果 $W$ 满足以下四个条件，则称其为**标准布朗运动**：

    **1. 初始零点**：$P(W(0) = 0) = 1$ a.s.
    
    **2. 独立增量 (Independent Increments)**：对于任意的时间划分 $0 \le t_1 < t_2 < \dots < t_n$，增量序列：

    $$
    W(t_1), W(t_2) - W(t_1), \dots, W(t_n) - W(t_{n-1})
    $$

    相互独立。
    
    **3. 平稳的高斯增量 (Stationary Gaussian Increments)**：对于任意 $0 \le s < t$，增量服从均值为 0，方差为时间差的正态分布：

    $$
    W(t) - W(s) \sim \mathcal{N}(0, t-s)
    $$
    
    **4. 轨道连续性 (Continuous Paths)**：几乎所有的样本轨道 $t \mapsto W(t, \omega)$ 都是连续的（即 a.s. 连续）。

由以上定义，我们可以立刻得出布朗运动的低阶矩与协方差结构，这是后续推导白噪声性质的关键。

!!! info "基本性质：矩与协方差结构"

    **1. 均值与方差**：
    由于 $W(t) = W(t) - W(0) \sim \mathcal{N}(0, t)$，显然有：

    $$
    E[W(t)] = 0, \quad Var(W(t)) = t
    $$

    **2. 自协方差函数 (Autocovariance)**：
    对于任意的 $s, t \ge 0$，有：

    $$
    R(s, t) = Cov(W(s), W(t)) = E[W(s)W(t)] = s \wedge t = \min(s, t)
    $$

    ??? proof "自协方差的推导（点击展开）"

        不失一般性，假设 $0 \le s \le t$。我们将 $W(t)$ 拆分为增量形式：

        $$
        E[W(s)W(t)] = E\big[W(s) \big( W(t) - W(s) + W(s) \big)\big]
        $$

        展开得到：

        $$
        = E[W(s)(W(t) - W(s))] + E[W(s)^2]
        $$

        由于布朗运动具有独立增量性质，$W(t) - W(s)$ 与 $W(s) = W(s) - W(0)$ 相互独立。又因为增量均值为 0，所以第一项为 0：

        $$
        = E[W(s)] E[W(t) - W(s)] + Var(W(s)) = 0 + s = s
        $$

        由于我们假设了 $s \le t$，因此一般情况可以写为 $s \wedge t$。$\square$

!!! tip "布朗运动的联合概率密度与转移密度"

    布朗运动在时间点 $0 < t_1 < t_2 < \dots < t_n$ 上的取值 $(W(t_1), \dots, W(t_n))$ 服从多元正态分布。
    由于马尔可夫性，其联合概率密度 (Joint PDF) 可以用**转移密度函数 (Transition Density)** 乘积的形式优雅地表达：

    定义高斯转移密度（从空间点 $y$ 经过时间 $t$ 转移到 $x$）：

    $$
    g(x, t | y) = \frac{1}{\sqrt{2\pi t}} e^{-\frac{|x-y|^2}{2t}}
    $$

    则联合概率密度为：

    $$
    p(x_1, t_1; x_2, t_2; \dots; x_n, t_n) = g(x_1, t_1 | 0) g(x_2, t_2 - t_1 | x_1) \dots g(x_n, t_n - t_{n-1} | x_{n-1})
    $$

---

## 2. 引入白噪声与 SDE 雏形

从常微分方程 (ODE) $\frac{dX(t)}{dt} = b(X(t), t)$ 过渡到随机微分方程 (SDE) 时，我们需要加入一个噪声项 $\xi(t)$。
理想的白噪声 $\xi(t)$ 应该在不同时刻完全不相关，即协方差呈现 Dirac $\delta$ 函数的性质：$E[\xi(s)\xi(t)] = \delta(s-t)$。
在数学上，这个所谓的白噪声正是**布朗运动的“形式导数”** $\dot{W}(t)$。

!!! info "定理：布朗运动增量商的协方差极限为 $\delta$ 函数"
    
    考虑布朗运动的差商过程 $\xi_h(t) = \frac{W(t+h) - W(t)}{h}$ ($h > 0$)。
    当 $h \to 0$ 时，其自协方差函数在广义函数（分布）的意义下收敛于 Dirac $\delta$ 函数。

    ??? proof "极限推导过程（点击展开）"

        我们计算差商在不同时刻 $s, t$ 的协方差 $E[\xi_h(s) \xi_h(t)]$。
        利用布朗运动的协方差性质 $E[W(u)W(v)] = u \wedge v$：
        
        $$
        E\left[ \frac{W(s+h)-W(s)}{h} \frac{W(t+h)-W(t)}{h} \right] = \frac{1}{h^2} \Big( (s+h \wedge t+h) - (s \wedge t+h) - (s+h \wedge t) + (s \wedge t) \Big)
        $$

        假设 $s \le t$，我们分析上式的非零区域：
        1. 当时间差 $|t-s| \ge h$ 时，上述四个项相互抵消，结果为 $0$。这表明只要时间间隔大于 $h$，差商过程就是不相关的。
        2. 当时间差 $|t-s| < h$ 时，存在重叠区间。计算可得协方差为：

        $$
        \varphi_h(t-s) = \frac{1}{h^2} (h - |t-s|)
        $$

        这是一个底边宽为 $2h$、高为 $\frac{1}{h}$ 的等腰三角形函数。
        显然，积分 $\int_{-\infty}^\infty \varphi_h(x) dx = 1$。
        当 $h \to 0$ 时，该函数在非零点处趋于 0，在 0 处趋于无穷大，且积分恒为 1。这正是 Dirac $\delta$ 函数的定义：

        $$
        \lim_{h \to 0} E[\xi_h(s) \xi_h(t)] = \delta(t-s)
        $$

        这就解释了为何 SDE 通常写成微分形式 $dX(t) = b(X,t)dt + \sigma(X,t)dW(t)$，因为 $W(t)$ 的真正导数并不存在，只能作为广义函数处理。$\square$

---

## 3. 多维布朗运动与核心性质

在量化金融与多粒子系统中，我们常面临高维随机现象。

!!! abstract "定义：多维布朗运动 (n-dimensional BM)"

    $n$ 维布朗运动定义为一个向量过程 $W(t) = (W^1(t), W^2(t), \dots, W^n(t))^T$，其中：
    1. 每一个分量 $W^k(t)$ 都是一个标准的一维布朗运动。
    2. 分量之间相互独立：即对于任意 $k \neq l$，$\sigma$-代数 $\sigma(W^k(t), t \ge 0)$ 与 $\sigma(W^l(t), t \ge 0)$ 独立。
    
    其分量间的协方差结构为：

    $$
    E[W^k(t) W^l(s)] = (t \wedge s) \delta_{kl}
    $$

    *(此处 $\delta_{kl}$ 为 Kronecker 记号，不是 Dirac $\delta$ 函数)*

布朗运动的经典之处在于是它兼具了**鞅 (Martingale)** 和 **马尔可夫过程 (Markov Process)** 的特征。

!!! success "定理：布朗运动是连续鞅"

    设 $\{\mathcal{F}_t\}_{t \ge 0}$ 为布朗运动自然生成的流 $\mathcal{F}_t = \sigma(W(u), 0 \le u \le t)$。则 $W(t)$ 是一个鞅。

    ??? proof "鞅性质证明（点击展开）"

        对于任意 $s \le t$，我们需要证明 $E[W(t) | \mathcal{F}_s] = W(s)$。
        通过增量拆分构造独立性：

        $$
        E[W(t) | \mathcal{F}_s] = E[W(t) - W(s) + W(s) | \mathcal{F}_s]
        $$

        由条件期望的线性性质：

        $$
        = E[W(t) - W(s) | \mathcal{F}_s] + E[W(s) | \mathcal{F}_s]
        $$

        因为布朗运动具有独立增量，未来的增量 $W(t) - W(s)$ 与历史信息流 $\mathcal{F}_s$ 完全独立，故独立则无关（等于无条件期望）；
        而 $W(s)$ 本身是 $\mathcal{F}_s$-可测的，故已知即常数（直接提出）：

        $$
        = E[W(t) - W(s)] + W(s) = 0 + W(s) = W(s)
        $$

        证毕。$\square$

!!! success "定理：布朗运动是马尔可夫过程"

    布朗运动满足马尔可夫性：即“未来仅依赖于现在，而与过去无关”。
    对于任意 Borel 集 $B \in \mathcal{B}(\mathbb{R}^n)$ 以及 $s \le t$：

    $$
    P(W(t) \in B | \mathcal{F}_s) = P(W(t) \in B | W(s)) \quad a.s.
    $$

    ??? proof "马尔可夫性严格测度论证明（点击展开）"

        利用指示函数 $\chi_B$（或记作 $I_B$），我们可以将概率写为条件期望：

        $$
        P(W(t) \in B | \mathcal{F}_s) = E[\chi_B(W(t)) | \mathcal{F}_s]
        $$

        引入增量，设函数 $f(x, y) = \chi_B(x + y)$。将 $W(t)$ 拆分为 $W(s)$ 和 $W(t) - W(s)$：

        $$
        E[\chi_B(W(t)) | \mathcal{F}_s] = E[f(W(s), W(t) - W(s)) | \mathcal{F}_s]
        $$

        此处应用测度论中独立性与条件期望的核心引理（冻结引理/替换定理）：
        因为 $W(s)$ 是 $\mathcal{F}_s$-可测的，而 $W(t) - W(s)$ 与 $\mathcal{F}_s$ 独立，我们可以将 $W(s)$ 视作“冻结”的常数 $x$，对另一部分求无条件期望，然后再把 $W(s)$ 代回来：

        $$
        = E[f(x, W(t) - W(s))]\Big|_{x = W(s)}
        $$

        设增量 $Z = W(t) - W(s) \sim \mathcal{N}(0, t-s)$，其密度函数为 $g(z)$，上式等于：

        $$
        = \int_{\mathbb{R}^n} \chi_B(x + z) g(z) dz \Bigg|_{x = W(s)}
        $$

        换元令 $y = x + z$，则 $dz = dy$：

        $$
        = \int_B g(y - x) dy \Bigg|_{x = W(s)} = \int_B \frac{1}{\sqrt{2\pi(t-s)}} e^{-\frac{|y - W(s)|^2}{2(t-s)}} dy
        $$

        这个积分结果显然只依赖于随机变量 $W(s)$ 的取值，而不依赖于 $\mathcal{F}_s$ 中 $s$ 时刻之前的任何信息。
        根据条件期望的定义，这恰好等于 $E[\chi_B(W(t)) | W(s)]$，即 $P(W(t) \in B | W(s))$。$\square$

---

## 4. Kolmogorov 连续性定理与轨道性质

布朗运动定义中直接假设了“轨道连续”。但数学上我们需要问：给定了一致的有限维分布后，是否**一定存在**一个具有连续轨道的修正版本？这需要借助极其强大的 Kolmogorov 连续性定理。

为了量化连续的“粗糙程度”，我们引入 Hölder 空间。

!!! abstract "定义：Hölder 连续空间 $C^\gamma$"

    如果一个函数 $f(t)$ 满足存在常数 $K > 0$ 使得：

    $$
    |f(t) - f(s)| \le K |t - s|^\gamma
    $$

    则称其具有指数为 $\gamma$ 的 Hölder 连续性。记为 $f \in C^\gamma$。
    （注：当 $\gamma = 1$ 时即为 Lipschitz 连续。布朗运动的轨道极其粗糙，它处处不可微，因此达不到 Lipschitz 连续）。

!!! info "定理：Kolmogorov 连续性定理 (Kolmogorov Continuity Theorem)"

    设 $X(t)$ 为一定义在区间 $[0, T]$ 上的随机过程。如果存在常数 $\alpha > 0, \beta > 0, C > 0$，使得对于所有的 $s, t \in [0, T]$：

    $$
    E[|X(t) - X(s)|^\beta] \le C |t - s|^{1+\alpha}
    $$

    那么 $X(t)$ 存在一个具有连续样本轨道的修正版本。更进一步，对于任意的 $\gamma \in \left(0, \frac{\alpha}{\beta}\right)$，该修正版本的样本轨道几乎必然 (a.s.) 是局部 $\gamma$-Hölder 连续的。

    ??? proof "核心证明思路 (基于 Borel-Cantelli 引理)（点击展开）"

        定理的严格证明十分繁琐，但其核心机理非常优美。
        
        **1. 考察二进有理数点：** 我们在区间上取二进网络点 $t_i = \frac{i}{2^n}$。
        考虑相邻两点增量过大的概率事件 $A_n$：

        $$
        A_n = \left\{ \omega : \max_{0 \le i < 2^n} \left| X\left(\frac{i+1}{2^n}\right) - X\left(\frac{i}{2^n}\right) \right| > \left(\frac{1}{2^n}\right)^\gamma \right\}
        $$
        
        **2. 利用 Chebyshev/Markov 不等式放缩：**

        $$
        P(A_n) \le \sum_{i=0}^{2^n-1} P\left( |X_{i+1} - X_i| > 2^{-n\gamma} \right) \le \sum_{i=0}^{2^n-1} \frac{E[|X_{i+1} - X_i|^\beta]}{2^{-n\gamma\beta}}
        $$

        代入定理给定的条件 $E[|X_{i+1} - X_i|^\beta] \le C (2^{-n})^{1+\alpha}$，可得：

        $$
        P(A_n) \le 2^n \cdot C 2^{-n(1+\alpha)} \cdot 2^{n\gamma\beta} = C 2^{-n(\alpha - \gamma\beta)}
        $$
        
        **3. 应用 Borel-Cantelli 引理：**
        由于我们选取了 $\gamma < \frac{\alpha}{\beta}$，故 $\alpha - \gamma\beta > 0$。这是一个公比小于 1 的几何级数，因此：

        $$
        \sum_{n=1}^\infty P(A_n) < \infty
        $$

        由 Borel-Cantelli 引理（第一引理）可知，$P(\limsup A_n) = 0$。即几乎必然地，存在 $N(\omega)$，使得对所有 $n \ge N(\omega)$，二进网络上的增量都被紧紧控制在 $2^{-n\gamma}$ 之内。通过一致收敛性，这保证了极限过程的 Hölder 连续性。$\square$

之后将此定理应用到布朗运动的连续性分析上，得出其重要的 $C^{\frac{1}{2}-}$ 性质。

!!! success "结论：布朗运动轨道的 Hölder 粗糙度为 $C^{\frac{1}{2}-}$"

    布朗运动 $W(t)$ 的样本轨道几乎处处是 $\gamma$-Hölder 连续的，对于任意 $\gamma \in \left(0, \frac{1}{2}\right)$。
    （即它无限逼近 $1/2$ 次 Hölder 连续，但不能达到 $1/2$）。

    ??? proof "推导过程（利用高斯矩）（点击展开）"

        由于布朗运动的增量 $W(t) - W(s) \sim \mathcal{N}(0, |t-s|)$，我们知道正态分布的高阶偶数矩有明确的公式。
        对于偶数次幂 $2m$（$m \in \mathbb{N}$）：

        $$
        E[|W(t) - W(s)|^{2m}] = (2m-1)!! \big( Var(W(t)-W(s)) \big)^m = (2m-1)!! |t - s|^m
        $$

        这完美匹配了 Kolmogorov 连续性定理的条件。我们令：
        - $\beta = 2m$
        - $1 + \alpha = m \implies \alpha = m - 1$
        
        根据定理，轨道的 Hölder 指数 $\gamma$ 必须满足：

        $$
        \gamma < \frac{\alpha}{\beta} = \frac{m - 1}{2m} = \frac{1}{2} - \frac{1}{2m}
        $$

        由于这个性质对所有的正整数 $m$ 都成立。我们可以让 $m \to \infty$，此时 $\frac{1}{2m} \to 0$。
        因此，对于任意严格小于 $\frac{1}{2}$ 的 $\gamma$，布朗运动都是 $\gamma$-Hölder 连续的。进而就得出了 **布朗运动轨道是 $C^{\frac{1}{2}-}$ 连续的** 这一结论。$\square$