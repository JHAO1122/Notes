# 1.2 外测度与 Carathéodory 扩张

在建立测度论的过程中，一个核心难题是如何在足够大的 $\sigma$-代数上构造测度。Carathéodory 扩张理论提供了一种标准方法：先定义一个涵盖所有子集的“外测度”，再通过一种巧妙的可测性准则筛选出那些表现良好的集合，从而形成 $\sigma$-代数。

---

## 1. 外测度 (Outer Measure)

外测度是定义在全集 $X$ 的所有子集上的集函数。与测度不同，它只需要满足次可加性，而不需要满足可数可加性。

!!! info "定义 1.2.1 (外测度)"

    设 $X$ 为一个非空集合。称集函数 $\mu^*: 2^X \to [0, +\infty]$ 为一个**外测度**，如果它满足以下性质：

    * **空集性质**：$\mu^*(\emptyset) = 0$。

    * **单调性**：若 $A \subset B \subset X$，则 $\mu^*(A) \le \mu^*(B)$。

    * **可数次可加性**：对于任意子集序列 $\{A_j\}_{j=1}^\infty \subset 2^X$，有：

    \[
    \mu^*\left( \bigcup_{j=1}^\infty A_j \right) \le \sum_{j=1}^\infty \mu^*(A_j)
    \]

### 外测度的普遍构造方法

我们可以通过一个简单的集族 $\Sigma$ 上的初等函数 $\rho$ 来生成外测度。这种方法在构造勒贝格测度时至关重要。

!!! success "定理 1.2.1 (外测度的生成)"

    设 $\Sigma \subset 2^X$ 为一族子集，满足 $\emptyset \in \Sigma$。设 $\rho: \Sigma \to [0, +\infty]$ 满足 $\rho(\emptyset) = 0$。对于任何 $A \subset X$，定义：

    \[
    \mu^*(A) = \inf \left\{ \sum_{j=1}^\infty \rho(E_j) : A \subset \bigcup_{j=1}^\infty E_j, E_j \in \Sigma \right\}
    \]

    则 $\mu^*$ 是 $X$ 上的一个外测度。

---

## 2. Carathéodory 可测准则

由于外测度在一般集合上不满足可加性，Carathéodory 提出了一种判别条件，用来定义哪些集合是“可测”的。

!!! info "定义 1.2.2 (Carathéodory 可测性)"

    设 $\mu^*$ 是 $X$ 上的一个外测度。称集合 $A \subset X$ 是 **$\mu^*$-可测的**，如果对任何试验集 $E \subset X$，都满足：

    \[
    \mu^*(E) = \mu^*(E \cap A) + \mu^*(E \cap A^c)
    \]

    习惯上，所有 $\mu^*$-可测集的全体记为 $\mathcal{M}$。



**直观理解**：一个集合 $A$ 是可测的，意味着它能够对空间中的任何集合 $E$ 进行“整齐”的分割，使得 $E$ 在 $A$ 内外的两部分测度之和等于 $E$ 整体的测度。验证可测性时，只需证明不等式 $\mu^*(E) \ge \mu^*(E \cap A) + \mu^*(E \cap A^c)$，因为另一半不等式由次可加性自动成立。

---

## 3. Carathéodory 扩张定理

这是测度论的基石定理，它揭示了从外测度筛选出的可测集族具有极好的代数结构。

!!! success "定理 1.2.2 (Carathéodory Theorem)"

    设 $\mu^*$ 是 $X$ 上的一个外测度，$\mathcal{M}$ 为所有 $\mu^*$-可测集构成的族。则：

    * **结构**：$\mathcal{M}$ 是一个 $\sigma$-代数。

    * **测度性**：$\mu^*$ 限制在 $\mathcal{M}$ 上是一个测度（满足可数可加性）。

    * **完备性**：测度空间 $(X, \mathcal{M}, \mu^*)$ 是完备的，即 $\mu^*$ 的零测集的所有子集都属于 $\mathcal{M}$。

??? proof "Carathéodory 定理证明概要"

    **1. 证明 $\mathcal{M}$ 是一个代数**：
    由于定义中 $A$ 与 $A^c$ 是对称的，显然 $A \in \mathcal{M} \implies A^c \in \mathcal{M}$。对于并集的封闭性，可以选取试验集 $E$，利用 $A$ 的可测性进行第一次拆分，再利用 $B$ 的可测性对拆分后的部分进行第二次拆分，通过次可加性组合各项即可。

    **2. 证明有限可加性**：
    设 $A, B \in \mathcal{M}$ 且 $A \cap B = \emptyset$。利用 $A$ 对 $E \cap (A \cup B)$ 的分割作用，可得 $\mu^*(E \cap (A \cup B)) = \mu^*(E \cap A) + \mu^*(E \cap B)$。

    **3. 扩展到可数可加性**：
    设 $\{A_j\}$ 是 $\mathcal{M}$ 中互不相交的序列。令 $B_n = \bigcup_{j=1}^n A_j$，$B = \bigcup_{j=1}^\infty A_j$。
    通过归纳可知 $\mu^*(E \cap B_n) = \sum_{j=1}^n \mu^*(E \cap A_j)$。
    由于 $\mu^*(E) = \mu^*(E \cap B_n) + \mu^*(E \cap B_n^c) \ge \sum_{j=1}^n \mu^*(E \cap A_j) + \mu^*(E \cap B^c)$。
    令 $n \to \infty$，并利用次可加性 $\mu^*(E \cap B) \le \sum \mu^*(E \cap A_j)$，两端夹逼即可证明 $\mu^*(E) = \mu^*(E \cap B) + \mu^*(E \cap B^c)$ 且测度满足可数可加。

---

## 4. 预测度与唯一性扩张

在实际应用中，测度通常起始于一个较小的集族（如代数 $\mathcal{A}$）。

!!! info "定义 1.2.3 (预测度 Pre-measure)"

    设 $\mathcal{A}$ 是 $X$ 上的一个代数。称集函数 $\mu_0: \mathcal{A} \to [0, +\infty]$ 为一个**预测度**，如果它满足 $\mu_0(\emptyset) = 0$，且对 $\mathcal{A}$ 中任何互不相交的集合列 $\{A_j\}$，只要 $\bigcup A_j \in \mathcal{A}$，就有：

    \[
    \mu_0\left(\bigcup_{j=1}^\infty A_j\right) = \sum_{j=1}^\infty \mu_0(A_j)
    \]

!!! success "定理 1.2.3 (扩张的唯一性)"

    设 $\mu_0$ 是代数 $\mathcal{A}$ 上的预测度，$\mu^*$ 是由 $\mu_0$ 诱导的外测度。

    * **存在性**：$\mu^*$ 在由 $\mathcal{A}$ 生成的 $\sigma$-代数 $\sigma(\mathcal{A})$ 上的限制是一个测度，且该测度在 $\mathcal{A}$ 上与 $\mu_0$ 一致。

    * **唯一性**：如果 $\mu_0$ 是 **$\sigma$-有限**的（即 $X$ 可表示为 $\mathcal{A}$ 中测度有限的集合的可数并），则该测度在 $\sigma(\mathcal{A})$ 上的扩张是唯一的。

如果预测度不满足 $\sigma$-有限性，扩张后的测度可能不唯一。这一结论强调了 $\sigma$-有限性在测度论研究中的重要性。