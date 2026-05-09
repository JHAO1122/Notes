# 1.3 Lebesgue-Stieltjes 测度与 Lebesgue 测度

在前面两节中，我们建立了从外测度到测度空间的抽象扩张理论。本节我们将该理论应用到实数集 $\mathbb{R}$ 上，通过单调增函数构造一类极其重要的 Borel 测度，即 Lebesgue-Stieltjes 测度，并由此引出分析学的基础——Lebesgue 测度。

---

## 1. 分布函数与 Borel 测度

在 $\mathbb{R}$ 的 Borel $\sigma$-代数 $\mathcal{B}_{\mathbb{R}}$ 上定义的测度称为 **Borel 测度**。如果该测度在任何有界区间上取值有限，则它与一类特殊的实函数——分布函数一一对应。

!!! info "定义 1.3.1 (分布函数 Distribution Function)"

    设 $F: \mathbb{R} \to \mathbb{R}$ 为一个函数。如果 $F$ 满足以下性质，则称其为分布函数：

    * **单调性**：$F$ 是单调递增的，即若 $x \le y$，则 $F(x) \le F(y)$。

    * **右连续性**：对于任意 $x_0 \in \mathbb{R}$，$F(x_0^+) = \lim_{x \to x_0^+} F(x) = F(x_0)$。

    通过分布函数，我们可以定义半开区间上的测度值为：

    \[
    \mu((a, b]) = F(b) - F(a)
    \]

---

## 2. Lebesgue-Stieltjes 测度的构造

构造的第一步是在一个代数上定义预测度。令 $\mathcal{A}$ 是由所有形如 $(a, b]$ 的半开区间（包括 $(a, \infty)$ 和 $(-\infty, b]$）的有限不相交并组成的族。易证 $\mathcal{A}$ 是一个代数。

!!! success "定理 1.3.1 (预测度的建立)"

    给定分布函数 $F$，在代数 $\mathcal{A}$ 上定义集函数 $\mu_0$：

    \[
    \mu_0\left( \bigcup_{j=1}^n (a_j, b_j] \right) = \sum_{j=1}^n [F(b_j) - F(a_j)]
    \]

    则 $\mu_0$ 是 $\mathcal{A}$ 上的一个预测度。

??? proof "定理 1.3.1 的证明要点"

    证明预测度的核心在于验证**可数可加性**。设 $I = (a, b]$ 且 $I = \bigcup_{j=1}^\infty I_j$，其中 $I_j = (a_j, b_j]$ 为互不相交的区间。

    **1. 有限可加性**：

    利用 $F$ 的单调性，容易证明对于有限项并集，$\mu_0$ 是可加的。对于 $n$ 项并集，通过对端点的整理即可得到 $\mu_0(\cup_{j=1}^n I_j) = \sum_{j=1}^n \mu_0(I_j)$。

    **2. 可数可加性的证明 (有限覆盖原理)**：

    * **不等式一**：由于有限并包含在可数并中，由单调性可知 $\mu_0(I) \ge \sum_{j=1}^n \mu_0(I_j)$，令 $n \to \infty$ 得 $\mu_0(I) \ge \sum_{j=1}^\infty \mu_0(I_j)$。

    * **不等式二**：任取 $\epsilon > 0$。由于 $F$ 是右连续的，对于 $I = (a, b]$，存在 $\delta > 0$ 使得 $F(a+\delta) - F(a) < \epsilon$，从而紧集 $K = [a+\delta, b] \subset (a, b]$。

    * 对于每个 $I_j = (a_j, b_j]$，存在 $\delta_j > 0$ 使得 $F(b_j+\delta_j) - F(b_j) < \epsilon / 2^j$，从而开区间 $O_j = (a_j, b_j+\delta_j) \supset (a_j, b_j]$。

    * 由于紧集 $K$ 被开区间族 $\{O_j\}$ 覆盖，由 **Heine-Borel 定理**，存在有限个区间 $O_{j_1}, \dots, O_{j_N}$ 覆盖 $K$。利用有限可加性的性质并令 $\epsilon \to 0$，最终证得 $\mu_0(I) \le \sum_{j=1}^\infty \mu_0(I_j)$。

!!! success "定理 1.3.2 (Lebesgue-Stieltjes 测度扩张)"

    由预测度 $\mu_0$ 诱导出的外测度 $\mu_F^*$ 限制在 Carathéodory 可测集族 $\mathcal{M}_F$ 上得到的测度称为 **Lebesgue-Stieltjes 测度**。由扩张的唯一性定理，若 $F$ 确定，该测度在 Borel $\sigma$-代数上的限制是唯一的。

---

## 3. 正则性与近似定理

Lebesgue-Stieltjes 测度具有良好的拓扑逼近性质，这使得我们可以用结构简单的开集或紧集来研究复杂的测度集。

!!! success "定理 1.3.3 (正则性 Regularity)"

    对于任何 $E \in \mathcal{M}_F$，其测度满足：

    * **外正则性**：$\mu_F(E) = \inf \{ \mu_F(U) : E \subset U, U \text{ 为开集} \}$。

    * **内正则性**：$\mu_F(E) = \sup \{ \mu_F(K) : K \subset E, K \text{ 为紧集} \}$。

由此推论，任何可测集 $E$ 具有以下结构描述：

* **$G_\delta$ 逼近**：存在 $G_\delta$ 集 $V$（可数个开集的交）使得 $E \subset V$ 且 $\mu_F(V \setminus E) = 0$。

* **$F_\sigma$ 逼近**：存在 $F_\sigma$ 集 $H$（可数个闭集的并）使得 $H \subset E$ 且 $\mu_F(E \setminus H) = 0$。

---

## 4. Lebesgue 测度 (Lebesgue Measure)

当分布函数取恒等函数，即 $F(x) = x$ 时，所得到的测度即为经典的 **Lebesgue 测度**。

!!! info "定义 1.3.2 (Lebesgue 测度)"

    由 $F(x) = x$ 诱导的测度空间 $(\mathbb{R}, \mathcal{L}, m)$ 称为 Lebesgue 测度空间。对于任何半开区间 $(a, b]$，其测度恰好等于其几何长度：

    \[
    m((a, b]) = b - a
    \]

### Lebesgue 测度的不变性性质

Lebesgue 测度与实数的加法和数乘结构相容，展现出极好的对称性。

!!! success "定理 1.3.4 (平移与伸缩不变性)"

    设 $E$ 为 Lebesgue 可测集，对于任何 $s, r \in \mathbb{R}$：

    * **平移不变性**：$m(E+s) = m(E)$，其中 $E+s = \{x+s : x \in E\}$。

    * **伸缩性质**：$m(rE) = |r|m(E)$，其中 $rE = \{rx : x \in E\}$。

---

## 5. Cantor 集与 Cantor 函数

Cantor 集是测度论中讨论“大小”与“势”最著名的反例。

!!! info "定义 1.3.3 (Cantor 三分集)"

    从闭区间 $[0, 1]$ 出发，第一步挖去中间的开区间 $(\frac{1}{3}, \frac{2}{3})$，第二步挖去剩下两段中间的各三分之一，如此无限循环。最终剩下的点集 $C$ 称为 **Cantor 集**。

!!! success "Cantor 集的性质"

    * **零测性**：Cantor 集的测度为 0。因为挖去的区间总长度为：

    \[
    \sum_{n=1}^\infty \frac{2^{n-1}}{3^n} = \frac{1/3}{1-2/3} = 1
    \]

    * **不可数性**：虽然测度为 0，但 Cantor 集的势与 $[0, 1]$ 相同，即基数为 $c$。

    * **拓扑性质**：$C$ 是紧集、无处稠密、且没有孤立点（完美集）。

!!! info "Cantor 函数 (Cantor Function / Devil's Staircase)"

    基于 Cantor 集的构造，可以定义一个连续单调函数 $f: [0, 1] \to [0, 1]$：

    * $f(x)$ 在 Cantor 集的补集（挖去的区间）上是常数。

    * $f(0) = 0, f(1) = 1$。

    * $f'(x) = 0$ 几乎处处成立，但函数却是连续增长的。这说明了“几乎处处导数为 0”并不一定意味着函数是常数。