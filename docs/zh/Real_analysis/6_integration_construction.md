# 3.1 抽象积分的构造与简单函数

在黎曼积分中，我们通过划分定义域上的区间来构造积分（即竖切法）。而在 Lebesgue 抽象积分论中，我们通过划分值域来构造积分（即横切法）。这种构造的核心桥梁是**简单函数 (Simple Function)**，它是连接测度与积分的基础。

---

## 1. 特征函数与简单函数

!!! info "定义 3.1.1 (特征函数与简单函数)"

    * **特征函数 (Characteristic Function)**：对于测度空间 $X$ 中的子集 $E \in \mathcal{M}$，定义其特征函数为：

    \[
    \chi_E(x) = \begin{cases} 1, & x \in E \\ 0, & x \notin E \end{cases}
    \]

    * **简单函数 (Simple Function)**：如果一个可测函数 $\phi: X \to \mathbb{C}$ 的值域是有限集，则称其为简单函数。任何简单函数都可以表示为特征函数的有限线性组合：

    \[
    \phi(x) = \sum_{j=1}^n a_j \chi_{E_j}(x)
    \]

    其中 $a_j \in \mathbb{C}$，且 $E_j = \{x \in X : \phi(x) = a_j\} \in \mathcal{M}$ 是一组互不相交的可测集分割。这被称为简单函数的**标准表示**。

简单函数的重要性在于，任何可测函数都可以被简单函数完美逼近。

!!! success "定理 3.1.1 (简单函数逼近定理)"

    设 $(X, \mathcal{M})$ 为可测空间。

    * 对于任意非负可测函数 $f: X \to [0, +\infty]$，必定存在一列非负简单函数 $\{\phi_n\}_{n=1}^\infty$，满足单调递增且逐点收敛于 $f$：

    \[
    0 \le \phi_1 \le \phi_2 \le \dots \le f, \quad \lim_{n \to \infty} \phi_n(x) = f(x)
    \]

    * 若 $f$ 是一般复值可测函数，则存在一列复简单函数序列 $\{\psi_n\}$，使得 $0 \le |\psi_1| \le |\psi_2| \le \dots \le |f|$ 且 $\psi_n \to f$。

---

## 2. 简单函数的积分

我们将从最简单的阶梯块开始定义“面积”或积分。

!!! info "定义 3.1.2 (非负简单函数的积分)"

    设 $(X, \mathcal{M}, \mu)$ 为测度空间，$\phi = \sum_{j=1}^n a_j \chi_{E_j}$ 是一非负简单函数的标准表示（$a_j \ge 0$），则定义 $\phi$ 关于测度 $\mu$ 的积分为：

    \[
    \int_X \phi \, d\mu = \sum_{j=1}^n a_j \mu(E_j)
    \]

    **注意**：在测度论中，我们采用约定 $0 \cdot \infty = 0$。这意味着即使某个集合的测度为无穷大，只要函数在该集合上取值为 0，其积分贡献就是 0。对于可测子集 $A \subset X$，局部积分定义为 $\int_A \phi \, d\mu = \int_X \phi \cdot \chi_A \, d\mu$。

简单函数的积分拥有完美的代数结构。

!!! success "定理 3.1.2 (简单函数积分的性质)"

    设 $\phi, \psi$ 为非负简单函数，$c \ge 0$ 为常数，则：

    * **线性性**：$\int c\phi = c \int \phi$ 且 $\int (\phi + \psi) = \int \phi + \int \psi$。

    * **单调性**：若 $\phi \le \psi$，则 $\int \phi \le \int \psi$。

    * **集函数性质**：映射 $A \mapsto \nu(A) = \int_A \phi \, d\mu$ 构成了 $\mathcal{M}$ 上的一个新测度。

---

## 3. 非负可测函数的积分

利用上确界，我们可以将积分从简单函数推广到任意非负可测函数。

!!! info "定义 3.1.3 (非负可测函数的积分)"

    设 $f: X \to [0, +\infty]$ 是一个非负可测函数（记为 $f \in L^+$）。定义 $f$ 的积分为所有位于 $f$ 下方的非负简单函数积分的上确界：

    \[
    \int_X f \, d\mu = \sup \left\{ \int_X \phi \, d\mu : 0 \le \phi \le f, \phi \text{ 为简单函数} \right\}
    \]

对于非负可测函数的积分，有以下与“几乎处处 (a.e.)”相关的核心性质。

!!! success "命题 3.1.1 (非负积分的几乎处处性质)"

    设 $f \in L^+$，则：

    * $\int f \, d\mu = 0$ 当且仅当 $f = 0$ 几乎处处成立 ($f = 0$ a.e.)。

    * 若 $\int f \, d\mu < \infty$，则 $f < \infty$ 几乎处处成立（即可积函数必然几乎处处有限）。

??? proof "证明：$\int f = 0 \iff f = 0 \text{ a.e.}$"

    **充分性 ($\Leftarrow$)**：若 $f = 0$ a.e.，对于任意满足 $0 \le \phi \le f$ 的简单函数 $\phi = \sum a_j \chi_{E_j}$，只要 $a_j > 0$，必有 $\mu(E_j) = 0$。因此 $\int \phi = 0$，由上确界定义得 $\int f = 0$。

    **必要性 ($\Rightarrow$)**：若 $\int f = 0$。定义集合 $E_n = \{x : f(x) \ge \frac{1}{n}\}$，那么显然有 $\frac{1}{n} \chi_{E_n} \le f$。
    
    由积分的单调性：
    
    \[
    \int f \ge \int \frac{1}{n} \chi_{E_n} = \frac{1}{n} \mu(E_n)
    \]
    
    既然 $\int f = 0$，必然推导出 $\mu(E_n) = 0$。
    而 $\{x : f(x) > 0\} = \bigcup_{n=1}^\infty E_n$，由测度的次可加性可知该集合测度为 0，即 $f = 0$ a.e.。

---

## 4. 一般可积函数与 $L^1$ 空间

有了非负函数的积分，我们就可以通过正负部拆分定义一般的实数或复数积分。

!!! info "定义 3.1.4 (实可积函数与 $L^1$ 空间)"

    设 $f: X \to \mathbb{R}$ 是实可测函数，将其分解为正部和负部 $f = f^+ - f^-$（且 $|f| = f^+ + f^-$）。

    如果正部和负部的积分都有限，即 $\int f^+ < \infty$ 且 $\int f^- < \infty$，则称 $f$ 在 $X$ 上是 **Lebesgue 可积的**。全体这样的函数构成的空间记为 $L^1(\mu)$。

    此时，定义 $f$ 的积分为：

    \[
    \int_X f \, d\mu = \int_X f^+ \, d\mu - \int_X f^- \, d\mu
    \]

    显然，$f \in L^1(\mu)$ 当且仅当 $\int |f| \, d\mu < \infty$。

对于复值函数 $f(x) = u(x) + i v(x)$，如果 $|f| \in L^1(\mu)$，则定义其积分为 $\int f = \int u + i \int v$。

!!! success "定理 3.1.3 (积分的绝对不等式)"

    若 $f \in L^1(\mu)$，则有绝对值不等式：

    \[
    \left| \int_X f \, d\mu \right| \le \int_X |f| \, d\mu
    \]

??? proof "绝对值不等式的证明 (复数情形)"

    设积分的值为复数 $z = \int f d\mu$。若 $z = 0$，不等式显然成立。
    若 $z \ne 0$，设 $z = r e^{i\theta}$ ($r > 0$)，令 $\alpha = e^{-i\theta}$，则 $\alpha z = |\int f d\mu|$。
    由于 $\alpha z$ 是实数，有：

    \[
    \left| \int f \right| = \alpha \int f = \int (\alpha f) = \text{Re} \left( \int \alpha f \right) = \int \text{Re}(\alpha f)
    \]

    因为对于任意复数都有 $\text{Re}(\alpha f) \le |\alpha f| = |f|$，由积分的单调性即得：

    \[
    \int \text{Re}(\alpha f) \le \int |f|
    \]

    综合即证得绝对不等式成立。

!!! success "命题 3.1.2 (积分在 a.e. 等价下的不变性)"

    设 $f, g$ 为可测函数。如果 $f = g$ 几乎处处成立 ($f = g$ a.e.)，且 $f \in L^1(\mu)$，那么 $g \in L^1(\mu)$ 且它们具有相同的积分：

    \[
    \int_X f \, d\mu = \int_X g \, d\mu
    \]

这一性质表明，Lebesgue 积分对于零测集上的改变是完全“免疫”的。在泛函分析中，我们通常将彼此几乎处处相等的函数视为 $L^1$ 空间中的同一个元素。