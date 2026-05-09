# 4.2 绝对连续与 Radon-Nikodym 定理

在上一节中，我们探讨了单一符号测度的内部结构（Hahn-Jordan 分解）。本节我们将视角转向同一可测空间上**两个不同测度之间的关系**。如果说“相互奇异”刻画了两个测度完全不相干的极端情况，那么“绝对连续”则描述了另一个极端：一个测度完全受制于另一个测度。

探索这一关系的顶点，便是现代分析和概率论的核心支柱——**Lebesgue-Radon-Nikodym 定理**。

---

## 1. 测度的绝对连续性

在微积分中，绝对连续函数具有良好的性质。在测度论中，绝对连续性反映了“零测集的遗传性”。

!!! info "定义 4.2.1 (绝对连续 Absolute Continuity)"

    设 $\mu$ 是可测空间 $(X, \mathcal{M})$ 上的一个非负测度，$\nu$ 是同一个空间上的符号测度（或复测度）。
    
    如果对于任意满足 $\mu(E) = 0$ 的可测集 $E \in \mathcal{M}$，都有：

    \[
    \nu(E) = 0
    \]

    则称测度 $\nu$ 关于 $\mu$ **绝对连续**，记作 $\nu \ll \mu$。

* **注记**：因为全变差测度 $|\nu|$ 在 $E$ 上的值由 $E$ 的子集的 $\nu$ 值决定，所以 $\nu \ll \mu$ 当且仅当 $|\nu| \ll \mu$。

绝对连续性有一个极其直观的 $\epsilon-\delta$ 等价刻画，它解释了“连续”一词的由来：只要集合的 $\mu$ 测度足够小，它的 $\nu$ 测度就可以任意小。

!!! success "定理 4.2.1 (绝对连续的 $\epsilon-\delta$ 刻画)"

    设 $\nu$ 为一个**有限**符号测度（即 $|\nu|(X) < \infty$）。那么 $\nu \ll \mu$ 当且仅当：对于任意 $\epsilon > 0$，必定存在 $\delta > 0$，使得对任意满足 $\mu(E) < \delta$ 的可测集 $E \in \mathcal{M}$，都有：

    \[
    |\nu(E)| < \epsilon
    \]

??? proof "定理 4.2.1 的证明"

    **充分性 ($\Leftarrow$)**：若 $\mu(E) = 0$，则对于任意 $\delta > 0$ 都满足 $\mu(E) < \delta$。由条件可知对任意 $\epsilon > 0$ 都有 $|\nu(E)| < \epsilon$，因此必定有 $\nu(E) = 0$，即 $\nu \ll \mu$。

    **必要性 ($\Rightarrow$)**：采用反证法。假设结论不成立，那么存在某个 $\epsilon > 0$，使得对任意的 $n \in \mathbb{N}$（取 $\delta = 1/2^n$），都存在对应的集合 $E_n \in \mathcal{M}$，满足：

    \[
    \mu(E_n) < \frac{1}{2^n} \quad \text{但是} \quad |\nu(E_n)| \ge \epsilon
    \]

    定义集合序列的“上极限”：

    \[
    F_k = \bigcup_{n=k}^\infty E_n, \quad F = \bigcap_{k=1}^\infty F_k
    \]

    由于测度 $\mu$ 的次可加性：

    \[
    \mu(F_k) \le \sum_{n=k}^\infty \mu(E_n) < \sum_{n=k}^\infty \frac{1}{2^n} = \frac{1}{2^{k-1}}
    \]

    显然当 $k \to \infty$ 时，$\mu(F_k) \to 0$。由自上连续性得 $\mu(F) = 0$。
    
    因为已知 $\nu \ll \mu$，所以必须有 $|\nu|(F) = 0$。
    
    然而，由于 $\nu$ 是有限测度，利用自上连续性：

    \[
    |\nu|(F) = \lim_{k \to \infty} |\nu|(F_k) \ge \limsup_{k \to \infty} |\nu(E_k)| \ge \epsilon > 0
    \]

    这导致了严重矛盾！因此假设不成立，原命题得证。

---

## 2. Lebesgue-Radon-Nikodym 定理

在很多应用场景中，给定积分测度 $d\nu = f d\mu$，我们可以轻易判定 $\nu \ll \mu$。那么反过来，如果已知 $\nu \ll \mu$，能否找到一个函数 $f$，使得 $d\nu = f d\mu$ 呢？

Lebesgue-Radon-Nikodym 定理完美地回答了这个问题，并给出了更为一般的分解结构。

!!! success "定理 4.2.2 (Lebesgue-Radon-Nikodym 定理)"

    设 $\mu$ 和 $\nu$ 是可测空间 $(X, \mathcal{M})$ 上的两个 **$\sigma$-有限测度**（$\nu$ 可以是符号测度）。
    
    那么，存在**唯一**的一对测度 $\lambda$ 和 $\rho$，使得：

    \[
    \nu = \lambda + \rho
    \]

    并且满足以下两个核心条件：

    * **奇异性**：$\lambda \perp \mu$ （即存在互补集将它们完全隔开）。

    * **绝对连续性**：$\rho \ll \mu$。而且存在一个 $\mu$-几乎处处唯一的函数 $f \in L^1(\mu)$（如果 $\nu$ 不取负无穷则为扩展的局部可积函数），使得对于所有的可测集 $E \in \mathcal{M}$ 有：

    \[
    \rho(E) = \int_E f \, d\mu
    \]

* **注记**：这个分解被称为 $\nu$ 关于 $\mu$ 的 **Lebesgue 分解**。$\lambda$ 称为奇异部分，$\rho$ 称为绝对连续部分。函数 $f$ 正好充当了两者之间的“密度”。

??? proof "Lebesgue-Radon-Nikodym 定理的核心构造证明（针对有限非负测度情形）"

    为简化，假设 $\mu$ 和 $\nu$ 都是有限非负测度。我们通过构造最大的“下方函数”来寻找 $f$。

    定义如下函数族：

    \[
    \mathcal{F} = \left\{ g : X \to [0, +\infty] \text{ 为可测函数 } \mid \int_E g \, d\mu \le \nu(E), \forall E \in \mathcal{M} \right\}
    \]

    显然常数函数 $g = 0$ 在 $\mathcal{F}$ 中，所以 $\mathcal{F}$ 非空。并且如果 $g_1, g_2 \in \mathcal{F}$，可以证明 $\max(g_1, g_2) \in \mathcal{F}$。

    令 $M = \sup_{g \in \mathcal{F}} \int_X g \, d\mu$。由于 $\nu(X) < \infty$，故 $M < \infty$。
    我们可以取到一个序列 $\{g_n\} \subset \mathcal{F}$，使得 $\int g_n d\mu \to M$。
    令 $f_n = \max(g_1, \dots, g_n)$，则 $\{f_n\}$ 是一个单调递增的函数列，且 $f_n \in \mathcal{F}$。
    
    根据单调收敛定理 (MCT)，令 $f = \lim_{n \to \infty} f_n = \sup_n f_n$，则 $f \in \mathcal{F}$ 并且：

    \[
    \int_X f \, d\mu = M
    \]

    现在，我们定义剩余的测度部分 $\lambda$ 为：

    \[
    \lambda(E) = \nu(E) - \int_E f \, d\mu
    \]

    由 $f \in \mathcal{F}$ 可知，$\lambda$ 是一个非负测度。现在我们只需证明 $\lambda \perp \mu$。
    如果 $\lambda$ 与 $\mu$ 不奇异，那么根据测度论技巧，必然能在某个测度为正的集合上找到一个极小的扰动 $\epsilon \mu \le \lambda$，从而使得 $f + \epsilon \chi_A$ 依然属于 $\mathcal{F}$。
    但这将导致 $\int (f + \epsilon \chi_A) d\mu = M + \epsilon \mu(A) > M$，与 $M$ 是上确界相矛盾！
    
    因此，$\lambda$ 必定关于 $\mu$ 奇异，定理的非负情形得证。符号测度或 $\sigma$-有限测度可通过 Hahn 分解和全空间分割转化为该情形。

---

## 3. Radon-Nikodym 定理与导数

作为 Lebesgue 分解定理的直接推论，当我们提前知道了 $\nu \ll \mu$ 时，奇异部分 $\lambda$ 就会退化为 0。

!!! success "定理 4.2.3 (Radon-Nikodym 定理)"

    设 $\mu, \nu$ 为 $\sigma$-有限测度，且 $\nu \ll \mu$。
    
    那么必定存在一个 $\mu$-几乎处处唯一的函数 $f$，使得对于任意可测集 $E \in \mathcal{M}$：

    \[
    \nu(E) = \int_E f \, d\mu
    \]

这个神奇的函数 $f$ 就被称为 $\nu$ 关于 $\mu$ 的 **Radon-Nikodym 导数 (Radon-Nikodym Derivative)**，通常记作符号分式：

\[
f = \frac{d\nu}{d\mu}
\]

Radon-Nikodym 导数在积分运算中完美模拟了微积分中的导数链式法则。

!!! success "命题 4.2.1 (Radon-Nikodym 导数的性质)"

    设所有的测度均为 $\sigma$-有限测度。

    * **积分法则**：如果 $g \in L^1(\nu)$，那么 $g \frac{d\nu}{d\mu} \in L^1(\mu)$，并且有直观的变量代换公式：

    \[
    \int_X g \, d\nu = \int_X g \frac{d\nu}{d\mu} \, d\mu
    \]

    * **链式法则**：如果 $\nu \ll \mu$ 且 $\mu \ll \lambda$，则 $\nu \ll \lambda$，并且几乎处处有：

    \[
    \frac{d\nu}{d\lambda} = \frac{d\nu}{d\mu} \frac{d\mu}{d\lambda}
    \]

    * **倒数法则**：如果 $\mu \ll \nu$ 且 $\nu \ll \mu$（即它们**等价**，具有完全相同的零测集），那么几乎处处有：

    \[
    \frac{d\nu}{d\mu} = \left( \frac{d\mu}{d\nu} \right)^{-1}
    \]

---

## 4. 复测度与极分解 (Complex Measures)

符号测度的取值限于实数。如果我们允许测度取复数值，就得到了泛函分析中非常重要的**复测度**。

!!! info "定义 4.2.2 (复测度 Complex Measure)"

    设 $(X, \mathcal{M})$ 为可测空间。集函数 $\nu: \mathcal{M} \to \mathbb{C}$ 称为一个**复测度**，如果它对任何互不相交的可测集列 $\{E_j\}$ 满足：

    \[
    \nu\left(\bigcup_{j=1}^\infty E_j\right) = \sum_{j=1}^\infty \nu(E_j)
    \]

    注意：为使复数级数的和与排列顺序无关，该级数必须是**绝对收敛**的。这意味着复测度**天然地只取有限值**，不能取 $\infty$。

任何复测度 $\nu$ 都可以被唯一地分解为实部和虚部：$\nu = \nu_r + i \nu_i$，其中 $\nu_r$ 和 $\nu_i$ 都是有限的符号测度。

!!! info "定义 4.2.3 (全变差 Total Variation)"

    对复测度 $\nu$，定义其全变差测度 $|\nu|$ 为：

    \[
    |\nu|(E) = \sup \sum_{j=1}^n |\nu(E_j)|
    \]

    其中上确界是对集合 $E$ 的所有有限可测划分 $\{E_j\}_{j=1}^n$ 取得的。$|\nu|$ 也是一个有限的非负测度。

有了 Radon-Nikodym 定理的加持，我们可以像表示复数的极坐标形式 $z = r e^{i\theta}$ 一样，对复测度进行绝妙的“极分解”。

!!! success "定理 4.2.4 (复测度的极分解 Polar Decomposition)"

    设 $\nu$ 为一个复测度。由于根据全变差的定义显然有 $\nu \ll |\nu|$，由 Radon-Nikodym 定理，必定存在一个复可测函数 $h = \frac{d\nu}{d|\nu|}$，使得：

    \[
    d\nu = h \, d|\nu|
    \]

    更重要的是，这个函数 $h$ 满足 $|\nu|$-几乎处处有：

    \[
    |h(x)| = 1
    \]

利用极分解，关于复测度的积分可以严格地定义转化为关于其全变差测度（非负测度）的积分。如果 $f \in L^1(|\nu|)$，则定义：

\[
\int_X f \, d\nu = \int_X f \cdot h \, d|\nu|
\]

从这个定义中，我们可以非常直接地推导出复积分绝对值的不等式界限：

\[
\left| \int_X f \, d\nu \right| \le \int_X |f| \, d|\nu|
\]