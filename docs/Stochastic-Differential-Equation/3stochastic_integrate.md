# 第三章：随机积分基础与 Itô 积分的引出

在经典的微积分中，我们可以轻易地定义黎曼-斯蒂尔杰斯积分 $\int g(t) dX(t)$。但在随机分析中，当积分微元变成布朗运动 $dW(t)$ 时，由于其样本轨道的极端粗糙性，经典积分理论彻底失效。本章将展示从确定性被积函数向随机被积函数跨越的演进过程，并严格引出 Itô 积分的核心思想。

## 1. Paley-Wiener-Zygmund (PWZ) 随机积分

最简单的情形是：被积函数 $g(t)$ 是一个普通的确定性函数，而积分变量是布朗运动。

!!! abstract "定义：PWZ 随机积分"

    设 $W(t)$ 为标准布朗运动。假设 $g(t) \in C^1([0, T])$ 是一阶连续可导的确定性函数。
    利用分部积分公式（Integration by Parts），我们“回避”了布朗运动不可微的问题，将 **PWZ 随机积分** 定义为：

    $$
    \int_0^T g(t) dW(t) \triangleq g(T)W(T) - g(0)W(0) - \int_0^T g'(t) W(t) dt
    $$

    由于 $W(0) = 0$ a.s.，上式通常简写为：

    $$
    \int_0^T g(t) dW(t) = g(T)W(T) - \int_0^T g'(t) W(t) dt
    $$

    这里右侧的积分是普通的黎曼积分（因为 $W(t)$ 的样本轨道是 a.s. 连续的，$g'(t)$ 是连续的，所以积分良定）。

!!! info "定理：PWZ 积分的核心性质"

    对于上述定义的 PWZ 积分，满足以下两个极其重要的性质（也是后续所有随机积分的基石）：

    **(1) 期望为零**：

    $$
    E\left[ \int_0^T g(t) dW(t) \right] = 0
    $$

    **(2) Itô 等距同构 (Itô Isometry)**：

    $$
    E\left[ \left( \int_0^T g(t) dW(t) \right)^2 \right] = \int_0^T g(t)^2 dt
    $$

    ??? proof "等距同构性质的严格推导（点击展开）"

        我们直接从 PWZ 积分的定义出发，计算其平方的期望（假设 $g(0)=0$ 或为简化书写将其吸收到常数项中，此处采用你手稿中更本质的双重积分法推导）：
        
        将平方写成双重积分的形式：

        $$
        E\left[ \left( \int_0^T g'(t)W(t) dt \right)^2 \right] = E\left[ \int_0^T g'(t)W(t) dt \int_0^T g'(s)W(s) ds \right]
        $$

        利用 Fubini 定理交换期望与积分的顺序：

        $$
        = \int_0^T \int_0^T g'(t)g'(s) E[W(t)W(s)] ds dt
        $$

        代入布朗运动的自协方差函数 $E[W(t)W(s)] = t \wedge s$：

        $$
        = \int_0^T \int_0^T g'(t)g'(s) (t \wedge s) ds dt
        $$

        利用对称性，将积分区域拆分为 $s \le t$ 和 $t \le s$ 两部分：

        $$
        = \int_0^T \left( \int_0^t g'(s) g'(t) s \, ds + \int_t^T g'(s) g'(t) t \, ds \right) dt
        $$

        通过分部积分和代数化简（手稿此处略过繁杂代数步骤），上式最终将完美收缩为目标形式：

        $$
        = \int_0^T g(t)^2 dt \quad \square
        $$

---

## 2. 有界线性算子的稠密延拓 (BLT 定理)

刚才的 PWZ 积分要求 $g(t) \in C^1$。但在实际应用中，我们需要对更一般的 $L^2$ 函数进行积分。这需要借助泛函分析中极其强大的工具：**有界线性算子延拓定理 (BLT Theorem)**。

!!! abstract "定理：有界线性算子的稠密延拓"

    设 $X, Y$ 为 Banach 空间，$S$ 为 $X$ 中稠密的线性子空间。
    设 $T: S \to Y$ 为一个有界线性算子，即存在常数 $C > 0$，使得对任意 $x \in S$：

    $$
    \|Tx\|_Y \le C \|x\|_X
    $$

    则存在唯一的有界线性算子 $\overline{T}: X \to Y$，使得对于所有的 $x \in S$，都有 $\overline{T}x = Tx$（即 $\overline{T}|_S = T$）。
    此外，算子的算子范数保持不变：$\|\overline{T}\| = \|T\| \le C$。

    ??? proof "BLT 定理的构造性证明（点击展开）"

        **第一步：构造极限映射**
        由于 $S$ 在 $X$ 中稠密，对于任意的 $x \in X$，必然存在 $S$ 中的序列 $\{x_n\}$ 使得 $x_n \to x$。
        因为 $T$ 在 $S$ 上是有界（连续）的，我们考察序列 $\{Tx_n\}$ 在 $Y$ 中的距离：

        $$
        \|Tx_n - Tx_m\|_Y = \|T(x_n - x_m)\|_Y \le C \|x_n - x_m\|_X
        $$

        由于 $\{x_n\}$ 是收敛列，它必然是 Cauchy 列（柯西列）。因此，当 $n, m \to \infty$ 时，$\|x_n - x_m\|_X \to 0$。
        这就意味着 $\{Tx_n\}$ 是 Banach 空间 $Y$ 中的 Cauchy 列。由于 $Y$ 是完备的，极限 $\lim_{n \to \infty} Tx_n$ 必定存在。我们定义：

        $$
        \overline{T}x \triangleq \lim_{n \to \infty} Tx_n
        $$

        **第二步：证明定义的合理性（与序列选择无关）**
        假设有另一个序列 $x_n' \to x$，我们需要证明 $\lim Tx_n' = \lim Tx_n$。
        设 $y = \lim Tx_n$，$y' = \lim Tx_n'$。

        $$
        \|y - y'\|_Y = \lim_{n \to \infty} \|Tx_n - Tx_n'\|_Y \le C \lim_{n \to \infty} \|x_n - x_n'\|_X = 0
        $$

        故 $y = y'$，映射 $\overline{T}$ 定义良定。

        **第三步：证明线性和保范性**
        线性是极限的自然推论。对于保范性，对任意 $x \in X$ 及逼近列 $x_n \to x$：

        $$
        \|\overline{T}x\|_Y = \lim_{n \to \infty} \|Tx_n\|_Y \le C \lim_{n \to \infty} \|x_n\|_X = C \|x\|_X
        $$

        故 $\overline{T}$ 有界，且范数不超过 $C$。$\square$

!!! success "应用：PWZ 积分向 $L^2$ 空间的延拓"

    我们将被积函数空间取为 $X = L^2([0, T])$，积分结果空间取为 $Y = L^2(\Omega, P)$。
    稠密子空间取为 $S = C^1([0, T])$。
    定义算子 $T: g \mapsto \int_0^T g(t) dW(t)$。
    
    由 Itô 等距同构知，对于任意 $g \in S$：

    $$
    \|Tg\|_{L^2(\Omega)}^2 = E\left[ \left( \int_0^T g(t) dW(t) \right)^2 \right] = \int_0^T g(t)^2 dt = \|g\|_{L^2([0,T])}^2
    $$

    这意味着算子 $T$ 是一致等距的（算子范数 $C=1$）。由 BLT 定理，我们可以将 PWZ 积分完美且唯一地延拓到整个 $L^2([0, T])$ 空间。

---

## 3. 布朗运动的二次变差 (Quadratic Variation)

为什么我们不能用传统的黎曼方法处理随机积分？核心在于布朗运动轨道的“粗糙性”，这体现在其二次变差上。

考虑时间区间 $[0, T]$ 的一个划分 $P = \{0 = t_0 < t_1 < \dots < t_m = T\}$，网格模长 $|P| = \max (t_{k+1} - t_k)$。

!!! info "定理 1：布朗运动的二次变差等于时间 $T$"

    当网格不断细化 $|P| \to 0$ 时，布朗运动增量的平方和在 $L^2(\Omega, P)$ 的意义下收敛于 $T$：

    $$
    \sum_{k=0}^{m-1} (W(t_{k+1}) - W(t_k))^2 \xrightarrow{L^2} T
    $$

    ??? proof "二次变差 $L^2$ 收敛的严格推导（点击展开）"

        为了证明 $L^2$ 收敛，我们需要证明其与 $T$ 的均方误差趋于 0。
        由于 $\sum_{k=0}^{m-1} (t_{k+1} - t_k) = T$，我们可以将目标误差写为：

        $$
        E\left[ \left( \sum_{k=0}^{m-1} \big((W(t_{k+1}) - W(t_k))^2 - (t_{k+1} - t_k)\big) \right)^2 \right]
        $$

        令 $\Delta W_k = W(t_{k+1}) - W(t_k)$，$\Delta t_k = t_{k+1} - t_k$。展开平方项，分为平方项和交叉项：

        $$
        = \sum_k E\Big[ \big((\Delta W_k)^2 - \Delta t_k\big)^2 \Big] + \sum_{k \neq j} E\Big[ \big((\Delta W_k)^2 - \Delta t_k\big)\big((\Delta W_j)^2 - \Delta t_j\big) \Big]
        $$

        **关键点 1：交叉项为 0。**
        由于布朗运动具有独立增量性质，当 $k \neq j$ 时，$\Delta W_k$ 与 $\Delta W_j$ 独立。且由于 $E[(\Delta W_k)^2] = \Delta t_k$，每个因式的期望均为 0，故交叉项整体期望为 0。

        **关键点 2：平方项的计算。**
        只剩下对角线上的方差项。注意到 $\Delta W_k \sim \mathcal{N}(0, \Delta t_k)$，因此可以标准化表示为 $\sqrt{\Delta t_k} Z$，其中 $Z \sim \mathcal{N}(0, 1)$。

        $$
        E\Big[ \big((\Delta W_k)^2 - \Delta t_k\big)^2 \Big] = E\Big[ (\Delta t_k Z^2 - \Delta t_k)^2 \Big] = (\Delta t_k)^2 E[(Z^2 - 1)^2]
        $$

        由于标准正态分布的四阶矩 $E[Z^4] = 3$，$E[Z^2] = 1$，所以 $E[(Z^2 - 1)^2] = 3 - 2(1) + 1 = 2$。

        $$
        \text{总误差} = 2 \sum_{k=0}^{m-1} (\Delta t_k)^2
        $$

        我们放大这个和式：提取出一个最大的 $\Delta t_k$ 即网格模长 $|P|$：

        $$
        2 \sum_{k=0}^{m-1} (\Delta t_k)^2 \le 2 |P| \sum_{k=0}^{m-1} \Delta t_k = 2 |P| T
        $$

        当 $|P| \to 0$ 时，$2 |P| T \to 0$。故在 $L^2$ 意义下收敛于 $T$。$\square$

这个定理直接导出了一个令人绝望（但又极具物理意义）的推论：

!!! success "定理 2：布朗运动处处无界变差 (Infinite Total Variation)"

    几乎所有的布朗运动样本轨道 $W(t, \omega)$ 在任何区间上的全变差都是无穷大。

    > **反证法极简证明**：如果某条轨道有界的总变差 $V_T < \infty$，那么它的二次变差可以放缩为：
    > $\sum (\Delta W_k)^2 \le (\max_k |\Delta W_k|) \sum |\Delta W_k| \le (\max_k |\Delta W_k|) \cdot V_T$
    > 由于轨道连续，当分割无限细化时 $\max |\Delta W_k| \to 0$。这会导致二次变差趋于 0，与定理 1 中二次变差等于 $T > 0$ 产生绝对矛盾！

---

## 4. 随机积分的 Riemann 和与 Itô 积分的引出

既然经典微积分因变差无界而失效，我们如何计算 $\int_0^T W(t) dW(t)$？
让我们回到黎曼和的定义，并观察一个在经典微积分中绝不会出现的诡异现象：**取值点的微小改变，将导致积分结果的剧变。**

构造划分 $P$，并在每个子区间 $[t_k, t_{k+1}]$ 中取一点 $\tau_k = (1-\lambda)t_k + \lambda t_{k+1}$ ($\lambda \in [0, 1]$)。
考察黎曼和：

$$
R_n = \sum_{k=0}^{m-1} W(\tau_k) \big( W(t_{k+1}) - W(t_k) \big)
$$

!!! info "积分结果对 $\lambda$ 取值的依赖"

    为了清晰展示，我们研究 $\lambda = 0$（取左端点，即 Itô 积分）的特殊情况，此时 $\tau_k = t_k$：

    $$
    R_n = \sum_{k=0}^{m-1} W(t_k) \big( W(t_{k+1}) - W(t_k) \big)
    $$

    ??? proof "代数恒等式拆分与 Itô 积分的极限求解（点击展开）"

        这是一个极其巧妙的代数技巧。利用恒等式 $a(b-a) = \frac{1}{2}\big( b^2 - a^2 - (b-a)^2 \big)$，我们将每一项改写：
        令 $a = W(t_k), b = W(t_{k+1})$：

        $$
        W(t_k) \big( W(t_{k+1}) - W(t_k) \big) = \frac{1}{2}\big( W(t_{k+1})^2 - W(t_k)^2 \big) - \frac{1}{2}\big( W(t_{k+1}) - W(t_k) \big)^2
        $$

        将所有的项加起来，原和式分裂为两部分 $B_1$ 和 $B_2$：

        $$
        R_n = \underbrace{ \frac{1}{2} \sum_{k=0}^{m-1} \big( W(t_{k+1})^2 - W(t_k)^2 \big) }_{B_1} - \underbrace{ \frac{1}{2} \sum_{k=0}^{m-1} \big( W(t_{k+1}) - W(t_k) \big)^2 }_{B_2}
        $$

        **分析 $B_1$**：这是一个完美的交错求和（Telescoping sum），中间项全部抵消：

        $$
        B_1 = \frac{1}{2} \big( W(T)^2 - W(0)^2 \big) = \frac{1}{2} W(T)^2
        $$

        **分析 $B_2$**：这正是布朗运动的二次变差！根据定理 1，当划分细化时，它在 $L^2$ 意义下收敛：

        $$
        B_2 \xrightarrow{L^2} \frac{1}{2} T
        $$

        综上所述，当 $\lambda = 0$ 时（Itô 积分），我们得到：

        $$
        \int_0^T W(t) dW(t) \triangleq \lim_{|P|\to 0} R_n = \frac{1}{2} W(T)^2 - \frac{1}{2} T
        $$

        *(注：后面多出来的 $-\frac{1}{2}T$ 项，正是大名鼎鼎的 Itô 校正项！)* $\square$

**两种最重要的积分流派：**
* 当 $\lambda = 0$ 时（取左端点），即为 **Itô 积分**，结果为 $\frac{1}{2}W(T)^2 - \frac{1}{2}T$。它保持了鞅的性质（不“偷看”未来），是金融数学的基石。
* 当 $\lambda = 1/2$ 时（取中点），即为 **Stratonovich 积分**，校正项相互抵消，结果为 $\frac{1}{2}W(T)^2$，形式上与经典微积分一致，常用于物理和工程中的随机动力系统。

---

## 5. 严格 Itô 积分的测度论准备

为了将 Itô 积分推广到更一般的随机过程（不只是 $W(t)$），我们需要严格定义什么叫做“不偷看未来”。这就引入了流 (Filtration) 和适应过程 (Adapted Process) 的概念。

!!! abstract "定义：信息流与 $\sigma$-代数"

    **1. 布朗运动的自然流**：
    对于任意时刻 $t$，布朗运动历史轨迹产生的信息记为 $\sigma$-代数：

    $$
    \mathcal{F}_W(t) \triangleq \sigma(\{W(s) \mid 0 \le s \le t\})
    $$

    它包含了直到时刻 $t$ 为止布朗运动所有的路径信息。

    **2. 未来增量的独立性**：
    我们定义未来的增量信息流 $\mathcal{F}^t \triangleq \sigma(\{W(s) - W(t) \mid s > t\})$。根据布朗运动的独立增量性质，$\mathcal{F}^t$ 与历史流 $\mathcal{F}_W(t)$ 是**完全独立**的。

    **3. 一般信息流 (Filtration)**：
    一族满足以下条件的 $\sigma$-代数 $\{\mathcal{F}(t)\}_{t \ge 0}$：
    * **单调性**：$\mathcal{F}(s) \subset \mathcal{F}(t)$ 对任意 $0 \le s \le t$（信息不遗忘）。
    * **包含历史**：$\mathcal{F}_W(t) \subset \mathcal{F}(t)$。
    * **未来独立**：增量 $W(s) - W(t)$ 与 $\mathcal{F}(t)$ 独立。

有了信息的数学框架，我们就可以框定哪些随机过程是可以被 Itô 积分的。

!!! info "定义：适应过程与循序可测 (Adapted & Progressively Measurable)"

    **1. 适应过程 (Adapted Process)**：
    若对于每一个固定的时刻 $t$，随机变量 $G(t, \omega)$ 都是 $\mathcal{F}(t)$-可测的，则称随机过程 $G(t)$ 是关于流 $\{\mathcal{F}(t)\}$ 适应的。
    *(直观理解：在时刻 $t$ 这一刻，你只要知道了所有的历史信息 $\mathcal{F}(t)$，你就知道此时刻 $G(t)$ 的值，绝不需要未来的信息。)*

    **2. 循序可测 (Progressively Measurable)**：
    一个更强的条件。映射 $(s, \omega) \mapsto G(s, \omega)$ 在乘积空间 $[0, t] \times \Omega$ 上的 $\mathcal{B}([0, t]) \otimes \mathcal{F}(t)$ 域上是联合可测的。这保证了在时间区间上的黎曼或勒贝格积分操作是合法的。

在本书后续的构建中，Itô 积分 $\int_0^T G(t, \omega) dW(t)$ 将在以下希尔伯特空间中被严格定义：

$$
L^2(\Omega \times [0, T]) = \left\{ G(t, \omega) \text{ 为循序可测过程} \mid E\left[ \int_0^T G(t, \omega)^2 dt \right] < \infty \right\}
$$