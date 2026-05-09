# 1.1 $\sigma$-代数与测度空间基础

在本章中，我们将建立测度论最基本的框架。与传统的黎曼积分依赖于区间的划分不同，勒贝格积分（以及抽象测度论）建立在对集合的划分之上。因此，我们需要首先定义哪些集合是“可测的”，这就引入了 $\sigma$-代数的概念。

---

## 1. 集合的代数与 $\sigma$-代数

为了使得测度的加法运算有意义，可测集族必须对集合的基本运算封闭。

!!! info "定义 1.1 (代数 Algebra)"

    设 $X$ 为一个非空集合，$\mathcal{A} \subset 2^X$ 为 $X$ 的一个子集族。如果 $\mathcal{A}$ 满足以下条件，则称其为一个**代数**：

    * 对于任意 $E_1, E_2 \in \mathcal{A}$，有 $E_1 \cup E_2 \in \mathcal{A}$。

    * 对于任意 $E \in \mathcal{A}$，其补集 $E^c \in \mathcal{A}$。

由上述定义容易推导，代数对有限次交集运算也是封闭的，并且必然包含全集 $X$ 和空集 $\emptyset$。

!!! info "定义 1.2 ($\sigma$-代数 $\sigma$-Algebra)"

    设 $\mathcal{M}$ 是一个代数。如果 $\mathcal{M}$ 进一步对**可数并**封闭，即对于任意集合列 $\{E_k\}_{k=1}^\infty \subset \mathcal{M}$，都有：

    \[
    \bigcup_{k=1}^\infty E_k \in \mathcal{M}
    \]

    则称 $\mathcal{M}$ 为一个 **$\sigma$-代数**。

结合德摩根定律和补集的封闭性，$\sigma$-代数对**可数交**必然也是封闭的。

!!! success "定理 1.1 (生成的 $\sigma$-代数与 Borel $\sigma$-代数)"

    * **生成的 $\sigma$-代数**：给定 $X$ 上的任意子集族 $\mathcal{E} \subset 2^X$，包含 $\mathcal{E}$ 的所有 $\sigma$-代数的交集仍然是一个 $\sigma$-代数，记为 $\mathcal{M}(\mathcal{E})$。它被称为由 $\mathcal{E}$ 生成的最小 $\sigma$-代数。

    * **Borel $\sigma$-代数**：如果 $X$ 是一个拓扑空间，设 $\mathcal{T}$ 为其上所有开集构成的族，则由 $\mathcal{T}$ 生成的 $\sigma$-代数称为 **Borel $\sigma$-代数**，记为 $\mathcal{B}_X = \mathcal{M}(\mathcal{T})$。

在实数空间 $\mathbb{R}^n$ 中，Borel $\sigma$-代数不仅可以由开集生成，还可以由开区间 $(a,b)$、闭区间 $[a,b]$、半开半闭区间 $[a,b)$，甚至形如 $(-\infty, a)$ 或 $[a, +\infty)$ 的射线族生成。

### 乘积 $\sigma$-代数与基本族

!!! info "定义 1.3 (乘积 $\sigma$-代数)"

    设 $\{(X_\alpha, \mathcal{M}_\alpha)\}_{\alpha \in A}$ 是一族可测空间。定义乘积空间 $X = \prod_{\alpha \in A} X_\alpha$，并设投影映射为 $\pi_\alpha: X \rightarrow X_\alpha$。

    那么，$X$ 上的乘积 $\sigma$-代数定义为使所有投影映射 $\pi_\alpha$ 都可测的最小 $\sigma$-代数：

    \[
    \bigotimes_{\alpha \in A} \mathcal{M}_\alpha = \mathcal{M} \left( \left\{ \pi_\alpha^{-1}(E_\alpha) \mid E_\alpha \in \mathcal{M}_\alpha, \alpha \in A \right\} \right)
    \]

!!! info "定义 1.4 (基本族 Elementary Family)"

    集族 $\mathcal{E} \subset 2^X$ 称为一个**基本族**，如果满足：

    * $\emptyset \in \mathcal{E}$。

    * 若 $E, F \in \mathcal{E}$，则 $E \cap F \in \mathcal{E}$。

    * 若 $E \in \mathcal{E}$，则 $E^c$ 可以表示为 $\mathcal{E}$ 中有限个互不相交集合的并集。

---

## 2. 测度与测度空间

在有了 $\sigma$-代数赋予“可测集”概念后，我们现在可以在其上定义“长度”、“面积”或“体积”的抽象概念，即测度。

!!! info "定义 1.5 (测度 Measure)"

    设 $(X, \mathcal{M})$ 为一个可测空间。如果集函数 $\mu: \mathcal{M} \rightarrow [0, +\infty]$ 满足以下条件，则称 $\mu$ 为 $\mathcal{M}$ 上的一个**测度**：

    * **非负性**：测度的取值在 $[0, +\infty]$ 之间。

    * **空集为零**：$\mu(\emptyset) = 0$。

    * **可数可加性 (Countable Additivity)**：对于 $\mathcal{M}$ 中任意互不相交的集合列 $\{E_k\}_{k=1}^\infty$ （即 $E_i \cap E_j = \emptyset, i \ne j$），有：

    \[
    \mu\left(\bigcup_{k=1}^\infty E_k\right) = \sum_{k=1}^\infty \mu(E_k)
    \]

如果给定了一个空间 $X$、其上的 $\sigma$-代数 $\mathcal{M}$ 以及测度 $\mu$，三元组 $(X, \mathcal{M}, \mu)$ 就构成了一个**测度空间 (Measure Space)**。

**几个常见的简单测度实例：**

* **计数测度 (Counting Measure)**：对于任意集合 $E$，$\mu(E)$ 等于集合中包含的元素个数。

* **狄拉克测度 (Dirac Delta Measure)**：给定空间中固定一点 $x_0$，若 $x_0 \in E$ 则 $\mu(E) = 1$，若 $x_0 \notin E$ 则 $\mu(E) = 0$。

### 测度的基本性质

测度函数由其公理定义，可以自然地推导出以下重要分析性质。

!!! success "定理 1.2 (测度的基本性质)"

    设 $(X, \mathcal{M}, \mu)$ 为测度空间，则 $\mu$ 满足：

    * **单调性 (Monotonicity)**：若 $E, F \in \mathcal{M}$ 且 $E \subset F$，则 $\mu(E) \le \mu(F)$。进一步，如果 $\mu(E) < \infty$，则有 $\mu(F \setminus E) = \mu(F) - \mu(E)$。

    * **次可加性 (Subadditivity)**：对于任意集合列 $\{E_j\}_{j=1}^\infty \subset \mathcal{M}$ （不要求互不相交），恒有：

    \[
    \mu\left(\bigcup_{j=1}^\infty E_j\right) \le \sum_{j=1}^\infty \mu(E_j)
    \]

    * **自下连续性 (Continuity from Below)**：对于单调递增的集合列 $E_1 \subset E_2 \subset \cdots$，有：

    \[
    \mu\left(\bigcup_{j=1}^\infty E_j\right) = \lim_{j \rightarrow \infty} \mu(E_j)
    \]

    * **自上连续性 (Continuity from Above)**：对于单调递减的集合列 $E_1 \supset E_2 \supset \cdots$，且满足 $\mu(E_1) < \infty$，有：

    \[
    \mu\left(\bigcap_{j=1}^\infty E_j\right) = \lim_{j \rightarrow \infty} \mu(E_j)
    \]

---

## 3. 零测集与测度空间的完备化

在实分析的后续讨论（例如“几乎处处”概念）中，零测集扮演着核心角色。然而，一般的测度空间并不能保证零测集的所有子集都可测，这就引出了完备化的概念。

!!! info "定义 1.6 (零测集与完备测度空间)"

    * **零测集 (Null Set)**：若 $E \in \mathcal{M}$ 满足 $\mu(E) = 0$，则称 $E$ 为 $\mu$-零测集。

    * **完备测度空间 (Complete Measure Space)**：如果零测集的任意子集都属于 $\mathcal{M}$ （因此其测度也必然为零），则称测度空间 $(X, \mathcal{M}, \mu)$ 是**完备的**。

一般测度空间中的零测子集不一定是可测的。但我们可以通过向原 $\sigma$-代数中人为添加这些零测子集，将任何测度空间“完备化”。

!!! success "定理 1.3 (测度空间的完备化)"

    对于任意测度空间 $(X, \mathcal{M}, \mu)$，我们可以构造其完备化空间 $(X, \overline{\mathcal{M}}, \overline{\mu})$，其中新的 $\sigma$-代数定义为：

    \[
    \overline{\mathcal{M}} = \{E \cup F \mid E \in \mathcal{M}, F \subset N \text{ 且 } N \in \mathcal{M}, \mu(N)=0\}
    \]

    在这个扩张的 $\sigma$-代数上定义 $\overline{\mu}(E \cup F) = \mu(E)$，则 $\overline{\mathcal{M}}$ 构成一个 $\sigma$-代数，且 $\overline{\mu}$ 是其上唯一的完备测度扩张。

??? proof "完备化定理的证明"

    要证明 $(X, \overline{\mathcal{M}}, \overline{\mu})$ 是一个完备的测度空间，我们需要验证两点：$\overline{\mathcal{M}}$ 是一个 $\sigma$-代数，且 $\overline{\mu}$ 在其上是良定义的完备测度。

    **第一步：验证 $\overline{\mathcal{M}}$ 是 $\sigma$-代数**

    对于可数并的封闭性：设 $A_j = E_j \cup F_j \in \overline{\mathcal{M}}$，其中 $E_j \in \mathcal{M}$，$F_j \subset N_j$ 且 $\mu(N_j) = 0$。

    \[
    \bigcup_{j=1}^\infty A_j = \left(\bigcup_{j=1}^\infty E_j\right) \cup \left(\bigcup_{j=1}^\infty F_j\right)
    \]

    因为 $\mathcal{M}$ 是 $\sigma$-代数，所以 $\bigcup_{j=1}^\infty E_j \in \mathcal{M}$。同时，$\bigcup_{j=1}^\infty F_j \subset \bigcup_{j=1}^\infty N_j$。由测度的次可加性：

    \[
    \mu\left(\bigcup_{j=1}^\infty N_j\right) \le \sum_{j=1}^\infty \mu(N_j) = 0
    \]

    所以并集依然具有 $E \cup F$ 的形式，证明了对可数并封闭。

    接着验证对补集封闭。设 $A = E \cup F \in \overline{\mathcal{M}}$，为了方便取补集可以假设 $E \cap N = \emptyset$。
    由于 $A^c = (E \cup F)^c = E^c \cap F^c$，并结合 $F \subset N$，可以改写为：

    \[
    A^c = (E \cup N)^c \cup (N \setminus F) \cap E^c
    \]

    这里 $(E \cup N)^c \in \mathcal{M}$，而 $(N \setminus F) \cap E^c \subset N$ 为零测集的子集，说明 $A^c \in \overline{\mathcal{M}}$，故 $\overline{\mathcal{M}}$ 是 $\sigma$-代数。

    **第二步：验证完备性**

    假设存在某个子集 $S \subset A \in \overline{\mathcal{M}}$ 且 $\overline{\mu}(A) = \mu(E) = 0$。
    因为 $A = E \cup F \subset E \cup N$，所以 $S \subset E \cup N$。
    由于 $\mu(E \cup N) \le \mu(E) + \mu(N) = 0$，这就说明 $S$ 包含在 $\mathcal{M}$ 的一个零测集内。按定义 $S = \emptyset \cup S \in \overline{\mathcal{M}}$。这证明了测度空间在完备化下的完备性。