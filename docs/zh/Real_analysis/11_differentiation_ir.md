# 4.3 实轴上的微分与微积分基本定理

在经典的数学分析中，微积分基本定理（Newton-Leibniz 公式）揭示了微分与积分的互逆关系。而在 Lebesgue 积分和测度论的框架下，我们不仅可以将积分推广到抽象的测度空间，同样也可以将“导数”的概念推广为测度关于测度的导数（即 Radon-Nikodym 导数）。本节我们将回归到实数轴 $\mathbb{R}$ 和 $\mathbb{R}^n$ 上，探讨局部可积函数的微分性质、有界变差函数、绝对连续函数，并最终给出 Lebesgue 意义下的微积分基本定理。

---

## 1. 测度的微分与 Hardy-Littlewood 极大函数

在 $\mathbb{R}^n$ 中，我们试图通过考察一个点附近小球内的平均测度密度来定义测度的“导数”。设 $\nu$ 为 $\mathbb{R}^n$ 上的复测度或符号测度，$m$ 为 Lebesgue 测度，定义在点 $x$ 处的微分为：

\[
F(x) := \lim_{r \to 0} \frac{\nu(B(r, x))}{m(B(r, x))}
\]

如果 $d\nu = f dm$，即 $\nu$ 关于 $m$ 绝对连续且密度为 $f$，那么上述极限可以写为 $f$ 在球 $B(r, x)$ 上的平均值：

\[
\frac{\nu(B(r, x))}{m(B(r, x))} = \frac{\int_{B(r, x)} f(y) \, dm(y)}{\int_{B(r, x)} 1 \, dm(y)}
\]

我们自然希望当半径 $r \to 0$ 时，这个平均值几乎处处收敛到 $f(x)$。为了证明这一点，我们需要引入分析学中极具威力的工具——Hardy-Littlewood 极大函数。

!!! info "定义 4.3.1 (局部可积与平均算子)"

    * **局部可积空间 $L_{local}^1(\mathbb{R}^n)$**：如果函数 $f: \mathbb{R}^n \to \mathbb{C}$ 在任意紧集（或有界集）上都是 Lebesgue 可积的，即对于任何有限测度的可测集 $K$ 都有 $\int_K |f| dm < \infty$，则称 $f$ 属于 $L_{local}^1(\mathbb{R}^n)$。

    * **平均算子 $A_r f(x)$**：对于 $f \in L_{local}^1(\mathbb{R}^n)$，定义其在以 $x$ 为心、$r$ 为半径的开球 $B(r, x)$ 上的平均值为：

    \[
    A_r f(x) = \frac{1}{m(B(r, x))} \int_{B(r, x)} f(y) \, dy
    \]

!!! info "定义 4.3.2 (Hardy-Littlewood 极大函数)"

    对于 $f \in L_{local}^1(\mathbb{R}^n)$，定义其极大函数为所有半径 $r > 0$ 下绝对值平均的上确界：

    \[
    Hf(x) = \sup_{r > 0} A_r |f|(x) = \sup_{r > 0} \frac{1}{m(B(r, x))} \int_{B(r, x)} |f(y)| \, dy
    \]

极大函数 $Hf(x)$ 刻画了函数 $f$ 在 $x$ 点附近局部最极端的振荡情况。极大函数最重要的性质是它满足**弱型 (1,1) 不等式**，即它的尾概率受控于原函数的 $L^1$ 范数。

!!! success "定理 4.3.1 (极大值定理 Maximal Theorem)"

    如果 $f \in L^1(\mathbb{R}^n)$，则对于任意 $\alpha > 0$，极大函数水平集 $E_\alpha = \{x : Hf(x) > \alpha\}$ 的 Lebesgue 测度满足：

    \[
    m(\{x : Hf(x) > \alpha\}) \le \frac{3^n}{\alpha} \int_{\mathbb{R}^n} |f| \, dm
    \]

??? proof "极大值定理的证明要点 (覆盖引理)"

    任取 $x \in E_\alpha$，由上确界的定义，必定存在一个半径 $r_x > 0$，使得在该球上的平均值大于 $\alpha$：

    \[
    \frac{1}{m(B(r_x, x))} \int_{B(r_x, x)} |f| \, dm > \alpha \implies m(B(r_x, x)) < \frac{1}{\alpha} \int_{B(r_x, x)} |f| \, dm
    \]

    所有的这些球 $\{B(r_x, x)\}_{x \in E_\alpha}$ 构成对集合 $E_\alpha$ 的一个开覆盖。
    
    利用 **Wiener/Vitali 覆盖引理**，我们可以从中选出一个**互不相交**的有限子集或可数子集 $\{B_j\}_{j=1}^\infty$，使得将这些选出的球的半径放大 3 倍后，能够覆盖住任何一个有界子集乃至整个原覆盖的测度核心部分。根据覆盖引理的几何性质，测度放大倍数为 $3^n$：

    \[
    m(E_\alpha) \le m\left(\bigcup_{j} 3 B_j\right) \le 3^n \sum_j m(B_j)
    \]

    由于球 $B_j$ 是互不相交的，我们可以进行求和放缩：

    \[
    m(E_\alpha) \le 3^n \sum_j \frac{1}{\alpha} \int_{B_j} |f| \, dm \le \frac{3^n}{\alpha} \int_{\bigcup B_j} |f| \, dm \le \frac{3^n}{\alpha} \int_{\mathbb{R}^n} |f| \, dm
    \]

    定理得证。

利用极大值定理，我们可以通过逼近手段证明实分析中最美丽的定理之一——Lebesgue 微分定理。

!!! success "定理 4.3.2 (Lebesgue 微分定理 Lebesgue Differentiation Theorem)"

    如果 $f \in L_{local}^1(\mathbb{R}^n)$，那么对于几乎所有的 $x \in \mathbb{R}^n$（即除了一个 Lebesgue 零测集外），都有：

    \[
    \lim_{r \to 0} A_r f(x) = f(x)
    \]

    更强地，几乎所有的 $x$ 实际上都是 $f$ 的 **Lebesgue 点**，满足：

    \[
    \lim_{r \to 0} \frac{1}{m(B(r, x))} \int_{B(r, x)} |f(y) - f(x)| \, dy = 0
    \]

---

## 2. 有界变差函数 (BV) 与 Borel 测度

在一维实轴 $\mathbb{R}$ 上，符号测度（或复测度）可以通过累积分布函数来描述。但与单调函数对应于非负测度不同，符号测度对应的是**有界变差函数**。

!!! info "定义 4.3.3 (有界变差函数 Bounded Variation, BV)"

    设 $F: \mathbb{R} \to \mathbb{C}$。如果在 $\mathbb{R}$ 上任意取有限个严格递增的分点 $x_0 < x_1 < \dots < x_n$，相应的差值绝对值之和都存在一致上界，即全变差：

    \[
    T_F(x) = \sup \sum_{j=1}^n |F(x_j) - F(x_{j-1})| < \infty
    \]

    则称 $F$ 为**有界变差函数**，记作 $F \in BV$。

任何有界变差函数都可以分解为两个单调递增函数的差（Jordan 分解）。为了与 Borel 测度建立严格的一一对应关系，我们需要对 BV 函数进行规范化处理，类似于测度的右连续性与空集为零性质。

!!! info "定义 4.3.4 (规范有界变差空间 NBV)"

    定义**规范有界变差函数 (Normalized Bounded Variation)** 空间 $NBV$ 为满足以下两个条件的 BV 函数集合：

    * $F$ 是右连续的，即 $F(x^+) = F(x)$。

    * $F$ 在负无穷处消失，即 $F(-\infty) = \lim_{x \to -\infty} F(x) = 0$。

!!! success "定理 4.3.3 (NBV 与 Borel 测度的对应关系)"

    * 若 $\nu$ 是 $\mathbb{R}$ 上的有限复 Borel 测度，定义 $F(x) = \nu((-\infty, x])$，那么必定有 $F \in NBV$。

    * 反之，若 $F \in NBV$，则存在**唯一**的有限复 Borel 测度 $\mu_F$，使得对于任意区间，有 $\mu_F((a, b]) = F(b) - F(a)$。

    * 此时，全变差测度与函数的全变差完美对应：$|\mu_F| = \mu_{T_F}$。

由于任意有界变差函数都可以分解为单调函数，而单调函数根据 Lebesgue 的结论是几乎处处可导的，因此我们有如下推论：

!!! success "命题 4.3.1 (NBV 函数的可导性)"

    * 如果 $F \in NBV$，则 $F$ 几乎处处可导，并且其导数 $F' \in L^1(m)$。

    * 如果 $\mu_F \perp m$（测度 $\mu_F$ 关于 Lebesgue 测度是奇异的），那么 $F'(x) = 0$ a.e.。

---

## 3. 绝对连续函数 (AC) 与微积分基本定理

虽然 NBV 函数几乎处处可导，但这并不意味着牛顿-莱布尼茨公式对它们都成立。为了使微积分基本定理严格成立，我们需要寻找一类比连续函数更强、使得“局部微小波动不能积累成宏大跳跃”的函数类——**绝对连续函数**。

!!! info "定义 4.3.5 (绝对连续函数 Absolutely Continuous, AC)"

    称函数 $F: [a, b] \to \mathbb{C}$ 在区间 $[a, b]$ 上是**绝对连续的**，如果对于任意给定的 $\epsilon > 0$，必定存在 $\delta > 0$，使得对于 $[a, b]$ 中任意**互不相交的有限个开区间** $\{(a_j, b_j)\}_{j=1}^n$，只要它们的总长度满足：

    \[
    \sum_{j=1}^n (b_j - a_j) < \delta
    \]

    就必定有对应的函数值波动总和满足：

    \[
    \sum_{j=1}^n |F(b_j) - F(a_j)| < \epsilon
    \]

显然，绝对连续（AC）蕴含着一致连续，进而蕴含连续。不仅如此，它在测度层面上等价于测度的绝对连续性：

!!! success "命题 4.3.2"

    如果 $F \in NBV$，那么 $F \in AC$ 当且仅当由它诱导的测度绝对连续于 Lebesgue 测度，即 $\mu_F \ll m$。

基于这层关系，我们终于可以给出实分析中的微积分基本定理：

!!! success "定理 4.3.4 (微积分基本定理 The Fundamental Theorem of Calculus)"

    设 $F: [a, b] \to \mathbb{C}$ 是一个函数，则以下命题是等价的：

    * **1. AC 性质**：$F$ 在 $[a, b]$ 上是绝对连续的 ($F \in AC$)。

    * **2. 积分表示**：存在一个局部可积函数 $f \in L^1([a, b])$，使得对于所有的 $x \in [a, b]$ 有：

    \[
    F(x) - F(a) = \int_a^x f(t) \, dt
    \]

    * **3. 几乎处处可导且导数可积恢复**：$F$ 几乎处处可微，$F' \in L^1([a, b])$，并且牛顿-莱布尼茨公式严格成立：

    \[
    F(x) - F(a) = \int_a^x F'(t) \, dt
    \]

---

## 4. 测度的 Lebesgue 分解与奇异测度

结合 Radon-Nikodym 定理与 Lebesgue 分解定理，$\mathbb{R}$ 上任何一个有限的 Borel 测度 $\nu$ 都可以被彻底分解为三个互不相交的纯粹部分：

\[
\nu = \nu_{ac} + \nu_{cs} + \nu_{d}
\]

* **1. 离散测度 (Discrete Measure) $\nu_d$**：
  它集中在有限或可数个点上，形式为 $\nu_d = \sum c_j \delta_{x_j}$，其中 $\delta_{x_j}$ 是狄拉克点测度。其对应的分布函数是一个纯粹的“阶梯函数”，只在可数个点处发生跳跃，除此之外处处导数为 0。

* **2. 绝对连续测度 (Absolutely Continuous Measure) $\nu_{ac}$**：
  它受控于 Lebesgue 测度 $m$，即 $\nu_{ac} \ll m$。必定存在一个密度函数 $f \in L^1$，使得 $d\nu_{ac} = f dm$。其对应的分布函数是绝对连续函数 (AC)。

* **3. 连续奇异测度 (Continuous Singular Measure) $\nu_{cs}$**：
  这是测度论中最违反直觉的存在。它是一个奇异测度（$\nu_{cs} \perp m$），意味着它完全集中在一个 Lebesgue 测度为 0 的集合上；但同时它又是一个连续测度，即它在任何单个点上的测度都为 0（没有跳跃点）。
  
  **经典的例子**：建立在 Cantor 三分集上的 Cantor 测度，其对应的 Cantor 函数（恶魔的阶梯）连续、单调递增，且几乎处处导数为 0（在除了零测集的 Cantor 集之外都是平坦的），但却完成了从 0 到 1 的攀升。它证明了并非所有几乎处处可导且连续的函数都能通过对其导数积分来恢复原函数！