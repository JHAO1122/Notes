# 3.2 三大收敛定理与依测度收敛

在黎曼积分中，要交换极限和积分号通常需要极为苛刻的“一致收敛”条件。而在勒贝格积分的框架下，由于积分是建立在测度基础上的，我们只需要很弱的条件即可实现积分与极限的交换。构筑这一宏大理论的基石便是著名的“三大收敛定理”。

---

## 1. 单调收敛定理 (MCT)

单调收敛定理 (Monotone Convergence Theorem) 是建立积分理论的第一块基石。它表明，对于非负递增的函数列，极限的积分等于积分的极限。

!!! success "定理 3.2.1 (单调收敛定理 MCT)"

    设 $\{f_n\}_{n=1}^\infty \subset L^+$ 为一列非负可测函数。如果它几乎处处单调递增，且逐点收敛于 $f$：

    \[
    0 \le f_1 \le f_2 \le \dots \le f_n \le \dots \text{ a.e.} \quad \text{且} \quad \lim_{n \to \infty} f_n(x) = f(x) \text{ a.e.}
    \]

    那么 $f \in L^+$，并且有：

    \[
    \int_X f \, d\mu = \lim_{n \to \infty} \int_X f_n \, d\mu
    \]

??? proof "单调收敛定理的证明"

    首先，由于 $f_n \le f$，由积分的单调性显然有 $\int f_n \le \int f$，从而 $\lim_{n \to \infty} \int f_n \le \int f$。

    接下来我们需要证明反向的不等式。任取一个满足 $0 \le \phi \le f$ 的非负简单函数 $\phi$，以及任意常数 $\alpha \in (0, 1)$。
    定义可测集序列：

    \[
    E_n = \{x \in X : f_n(x) \ge \alpha \phi(x)\}
    \]

    因为 $f_n$ 单调递增，所以集合序列是递增的，即 $E_1 \subset E_2 \subset \dots$。
    又因为 $f_n \to f \ge \phi > \alpha \phi$ (当 $\phi > 0$ 时)，所以 $\bigcup_{n=1}^\infty E_n = X$。

    根据积分的性质，我们有：

    \[
    \int_X f_n \, d\mu \ge \int_{E_n} f_n \, d\mu \ge \alpha \int_{E_n} \phi \, d\mu
    \]

    令 $n \to \infty$，利用测度（及简单函数积分）的自下连续性：

    \[
    \lim_{n \to \infty} \int_X f_n \, d\mu \ge \alpha \int_X \phi \, d\mu
    \]

    由于 $\alpha \in (0, 1)$ 是任意的，令 $\alpha \to 1$，得到 $\lim \int f_n \ge \int \phi$。再对所有满足 $0 \le \phi \le f$ 的简单函数取上确界，即得：

    \[
    \lim_{n \to \infty} \int_X f_n \, d\mu \ge \int_X f \, d\mu
    \]

    结合两个方向的不等式，定理得证。

利用 MCT，我们可以立刻推导出级数积分的逐项求和性质。

!!! info "推论 (非负级数的逐项积分)"

    设 $\{f_n\}_{n=1}^\infty \subset L^+$ 为非负可测函数列，则：

    \[
    \int_X \left( \sum_{n=1}^\infty f_n \right) d\mu = \sum_{n=1}^\infty \int_X f_n \, d\mu
    \]

---

## 2. Fatou 引理

对于不一定单调的非负函数列，我们无法直接得到等号，但可以得到一个极具价值的不等式。它是连接 MCT 和 DCT 的桥梁。

!!! success "定理 3.2.2 (Fatou 引理)"

    设 $\{f_n\}_{n=1}^\infty \subset L^+$ 是一列非负可测函数，则：

    \[
    \int_X \left( \liminf_{n \to \infty} f_n \right) d\mu \le \liminf_{n \to \infty} \int_X f_n \, d\mu
    \]

??? proof "Fatou 引理的证明"

    定义函数列 $g_n(x) = \inf_{k \ge n} f_k(x)$。

    显然 $\{g_n\}$ 也是非负可测函数列，且由于随着 $n$ 的增大，求下确界的范围在缩小，所以 $\{g_n\}$ 是单调递增的，即 $g_1 \le g_2 \le \dots$。
    由下极限的定义：

    \[
    \lim_{n \to \infty} g_n(x) = \liminf_{n \to \infty} f_n(x)
    \]

    因为对于所有的 $k \ge n$ 都有 $g_n \le f_k$，特别地有 $g_n \le f_n$，所以：

    \[
    \int_X g_n \, d\mu \le \int_X f_n \, d\mu
    \]

    两边取下极限：

    \[
    \liminf_{n \to \infty} \int_X g_n \, d\mu \le \liminf_{n \to \infty} \int_X f_n \, d\mu
    \]

    另一方面，由于 $g_n$ 单调递增，应用单调收敛定理 (MCT) 于 $\{g_n\}$：

    \[
    \liminf_{n \to \infty} \int_X g_n \, d\mu = \lim_{n \to \infty} \int_X g_n \, d\mu = \int_X \left( \lim_{n \to \infty} g_n \right) d\mu = \int_X \left( \liminf_{n \to \infty} f_n \right) d\mu
    \]

    综合以上两式，即得 Fatou 引理。

---

## 3. 控制收敛定理 (DCT)

这是实分析中使用频率最高的收敛定理。它放宽了单调性的要求，代价是需要一个可积的“控制函数”。

!!! success "定理 3.2.3 (Lebesgue 控制收敛定理 DCT)"

    设 $\{f_n\}$ 是一列可测函数，几乎处处收敛于 $f$：

    \[
    \lim_{n \to \infty} f_n(x) = f(x) \quad \text{a.e.}
    \]

    如果存在一个非负的**可积函数** $g \in L^1(\mu)$，使得对于所有的 $n$，几乎处处有：

    \[
    |f_n(x)| \le g(x) \quad \text{a.e.}
    \]

    那么 $f \in L^1(\mu)$，并且：

    \[
    \lim_{n \to \infty} \int_X f_n \, d\mu = \int_X f \, d\mu \quad \text{且} \quad \lim_{n \to \infty} \int_X |f_n - f| \, d\mu = 0
    \]

??? proof "控制收敛定理的证明"

    由于 $|f_n| \le g$ 且 $f_n \to f$，取极限知 $|f| \le g$ a.e.。因为 $g \in L^1$，所以 $f \in L^1$。
    
    我们构造两个非负函数列：$g + f_n \ge 0$ 以及 $g - f_n \ge 0$。

    对 $g + f_n$ 应用 Fatou 引理：

    \[
    \int_X (g + f) \, d\mu \le \liminf_{n \to \infty} \int_X (g + f_n) \, d\mu = \int_X g \, d\mu + \liminf_{n \to \infty} \int_X f_n \, d\mu
    \]

    因为 $g$ 可积，可以消去两边的 $\int g$，得到：

    \[
    \int_X f \, d\mu \le \liminf_{n \to \infty} \int_X f_n \, d\mu
    \]

    同理，对 $g - f_n$ 应用 Fatou 引理：

    \[
    \int_X (g - f) \, d\mu \le \liminf_{n \to \infty} \int_X (g - f_n) \, d\mu = \int_X g \, d\mu - \limsup_{n \to \infty} \int_X f_n \, d\mu
    \]

    再次消去 $\int g$，化简得到：

    \[
    \limsup_{n \to \infty} \int_X f_n \, d\mu \le \int_X f \, d\mu
    \]

    结合上述两式，必有 $\liminf = \limsup = \lim$，从而证明了 $\lim_{n \to \infty} \int f_n = \int f$。
    将上述证明中的 $f_n$ 替换为 $|f_n - f|$（它被 $2g$ 控制），同理可证 $\lim \int |f_n - f| = 0$。

---

## 4. 依测度收敛与经典关系

除了几乎处处收敛 (a.e.) 和 $L^1$ 收敛，测度论中还有一种在概率论中极为重要的收敛性——依测度收敛（对应于概率论中的依概率收敛）。

!!! info "定义 3.2.1 (依测度收敛)"

    设 $f_n, f$ 为可测函数。如果对于任意 $\epsilon > 0$，满足以下条件的集合的测度趋于 0：

    \[
    \lim_{n \to \infty} \mu(\{x \in X : |f_n(x) - f(x)| \ge \epsilon\}) = 0
    \]

    则称 $f_n$ **依测度收敛**于 $f$，记作 $f_n \xrightarrow{\mu} f$。

同理，若 $\mu(\{x : |f_n(x) - f_m(x)| \ge \epsilon\}) \to 0$ ($n, m \to \infty$)，则称 $\{f_n\}$ 是**依测度的 Cauchy 列**。

依测度收敛与几乎处处收敛之间有着微妙的关系。著名的 Riesz 定理指出了这两种收敛的子列关系。

!!! success "定理 3.2.4 (Riesz 定理)"

    * 如果 $\{f_n\}$ 是依测度的 Cauchy 列，则它必然依测度收敛到某个可测函数 $f$。

    * 如果 $f_n \xrightarrow{\mu} f$，那么必定存在 $\{f_n\}$ 的一个子列 $\{f_{n_k}\}$，使得该子列几乎处处收敛于 $f$：

    \[
    f_{n_k}(x) \to f(x) \quad \text{a.e. as } k \to \infty
    \]

在有限测度空间上，几乎处处收敛甚至可以推导出比依测度收敛更强的性质，即**近乎一致收敛**。

!!! success "定理 3.2.5 (Egoroff 定理)"

    设 $\mu(X) < \infty$。如果 $\{f_n\}$ 在 $X$ 上几乎处处收敛于 $f$，那么对于任意 $\epsilon > 0$，存在一个可测子集 $E \subset X$ 满足：

    * $\mu(E) < \epsilon$ （剔除的坏集测度极小）

    * $f_n$ 在保留的集合 $E^c = X \setminus E$ 上**一致收敛**于 $f$。

??? proof "Egoroff 定理的证明核心"

    构造集合 $E_n(k) = \bigcup_{m=n}^\infty \{x : |f_m(x) - f(x)| \ge \frac{1}{k}\}$。
    因为 $f_n \to f$ a.e.，对于固定的 $k$，随着 $n \to \infty$，集合序列 $E_n(k)$ 单调递减且其交集的测度为 0。
    由于全空间测度有限，利用测度的自上连续性，必然有 $\lim_{n \to \infty} \mu(E_n(k)) = 0$。
    
    对于每一个 $k$，我们可以选取足够大的 $N_k$，使得 $\mu(E_{N_k}(k)) < \epsilon \cdot 2^{-k}$。
    令坏集 $E = \bigcup_{k=1}^\infty E_{N_k}(k)$。由次可加性，$\mu(E) < \sum \epsilon 2^{-k} = \epsilon$。
    在 $E^c$ 上，对于任意 $k$ 只要 $n > N_k$ 就有 $|f_n - f| < 1/k$，这恰好是一致收敛的定义。