# 📏 测度论基础 (Measure Theory Basics)

本模块涵盖了集合论与测度论的核心内容，包括 Lebesgue 测度、可测函数、Lebesgue 积分理论以及三大控制收敛定理。

---

## 一、 集合代数与基数 (Sets and Cardinality)

!!! info "对等与基数 (Equipotence)"
    若集合 $A$ 到 $B$ 之间存在一个双射，则称 $A$ 与 $B$ 对等，记作 $A \sim B$。
    
    * **可数集**：与正整数集 $\mathbb{N}$ 对等的集合称为可数集，其基数记为 $\aleph_0$。
    
    * **连续统基数**：与实数集 $\mathbb{R}$ 对等的集合其基数记为 $\mathfrak{c}$。已知 $\mathfrak{c} = 2^{\aleph_0}$。

??? note "Bernstein 定理"
    设 $A$ 和 $B$ 是两个集合。若存在 $A$ 到 $B$ 的单射，且存在 $B$ 到 $A$ 的单射，则 $A \sim B$。

!!! success "开集、闭集与 Cantor 集"
    * **开集的结构**：$\mathbb{R}$ 上的任意非空开集都可以唯一地表示为至多可数个互不相交的开区间的并。
    
    * **Cantor 三分集 ($C$)**：通过对闭区间 $[0, 1]$ 不断挖去中间三分之一开区间构造而成。
        * **性质**：$C$ 是紧集、完全集、无处稠密集，其基数为 $\mathfrak{c}$，但其 Lebesgue 测度为 0。

---

## 二、 Lebesgue 测度理论 (Lebesgue Measure)

!!! info "外测度 (Outer Measure)"
    对于 $\mathbb{R}^n$ 中的任意子集 $E$，其 Lebesgue 外测度 $m^*(E)$ 定义为：

    \[
    m^*(E) = \inf \left\{ \sum_{i=1}^\infty |I_i| \ \middle|\  E \subset \bigcup_{i=1}^\infty I_i \right\}
    \]

    其中 $\{I_i\}$ 为覆盖 $E$ 的至多可数个开矩体序列，$|I_i|$ 表示矩体的体积。

!!! success "Carathéodory 可测性定义"
    集合 $E \subset \mathbb{R}^n$ 称为 **Lebesgue 可测的**，若对于任意的“探测集” $T \subset \mathbb{R}^n$，均满足：

    \[
    m^*(T) = m^*(T \cap E) + m^*(T \cap E^c)
    \]

    此时，其外测度称为 Lebesgue 测度，记为 $m(E) = m^*(E)$。可测集的全体构成一个 $\sigma$-代数。

??? note "测度的连续性"
    1. **递增测度列的连续性**：若 $E_1 \subset E_2 \subset \dots \subset E_n \subset \dots$ 且均可测，则：

    \[
    m\left( \bigcup_{n=1}^\infty E_n \right) = \lim_{n \to \infty} m(E_n)
    \]

    2. **递减测度列的连续性**：若 $E_1 \supset E_2 \supset \dots \supset E_n \supset \dots$ 且均可测，**且 $m(E_1) < \infty$**，则：

    \[
    m\left( \bigcap_{n=1}^\infty E_n \right) = \lim_{n \to \infty} m(E_n)
    \]

---

## 三、 可测函数及其收敛性 (Measurable Functions)

!!! abstract "可测函数定义"
    设 $E \subset \mathbb{R}^n$ 为可测集，$f$ 是定义在 $E$ 上的广义实值函数。若对于任意的实数 $a$，集合

    \[
    \{x \in E \mid f(x) > a\}
    \]

    均是 Lebesgue 可测集，则称 $f$ 是 $E$ 上的**可测函数**。

!!! info "几乎处处 (Almost Everywhere, a.e.)"
    若某一性质在集合 $E$ 中除去一个测度为零的子集之外，在其余所有点上都成立，则称该性质在 $E$ 上**几乎处处成立**。

??? tip "四种函数收敛性的关系"
    设 $\{f_n\}$ 与 $f$ 是 $E$ 上的有限可测函数。
    
    1. **几乎处处收敛 ($f_n \xrightarrow{a.e.} f$)**：$m(\{x \in E \mid \lim_{n\to\infty} f_n(x) \neq f(x)\}) = 0$。
    
    2. **依测度收敛 ($f_n \xrightarrow{m} f$)**：对任意 $\epsilon > 0$，有 $\lim_{n \to \infty} m(\{x \in E \mid |f_n(x) - f(x)| > \epsilon\}) = 0$。

    * **蕴含链与转换**：
        * 若 $m(E) < \infty$，则 $f_n \xrightarrow{a.e.} f \implies f_n \xrightarrow{m} f$（叶戈罗夫定理）。
        
        * 若 $f_n \xrightarrow{m} f$，则必定存在**子序列** $\{f_{n_k}\}$ 使得 $f_{n_k} \xrightarrow{a.e.} f$（里斯定理）。

!!! success "三大核心定理 (Egorov & Lusin)"
    * **Egorov (叶戈罗夫) 定理**：设 $m(E) < \infty$。若 $f_n \xrightarrow{a.e.} f$，则对任意 $\delta > 0$，存在可测子集 $E_0 \subset E$ 满足 $m(E \setminus E_0) < \delta$，使得 $\{f_n\}$ 在 $E_0$ 上**一致收敛**于 $f$。
    
    * **Lusin (鲁欣) 定理**：设 $f$ 是 $E$ 上的几乎处处有限的可测函数。则对任意 $\delta > 0$，存在闭集 $F \subset E$ 满足 $m(E \setminus F) < \delta$，使得 $f$ 限制在 $F$ 上是**连续函数**。

---

## 四、 Lebesgue 积分与控制收敛定理 (Lebesgue Integral)

!!! info "非负简单函数的积分"
    设 $\phi(x) = \sum_{i=1}^k c_i \chi_{A_i}(x)$ 是 $E$ 上的非负简单函数，其中 $A_i$ 两两不交且可测，则其 Lebesgue 积分定义为：

    \[
    \int_E \phi(x) dx = \sum_{i=1}^k c_i m(A_i)
    \]

!!! success "一般函数的 Lebesgue 积分"
    * **非负可测函数**：定义 $\int_E f(x) dx = \sup \{ \int_E \phi(x) dx \mid 0 \le \phi \le f, \phi \text{ 为简单函数} \}$。
    
    * **一般可测函数**：引入正部 $f^+ = \max(f, 0)$ 和负部 $f^- = \max(-f, 0)$。若 $\int_E f^+ dx$ 与 $\int_E f^- dx$ 至少有一个有限，则定义：

    \[
    \int_E f(x) dx = \int_E f^+(x) dx - \int_E f^-(x) dx
    \]

    若两部分均有限，则称 $f$ 在 $E$ 上 **Lebesgue 可积**。

??? success "三大核心积分极限控制定理"
    * **1. 列维 (Levi) 单调收敛定理**：若 $0 \le f_1(x) \le f_2(x) \le \dots \le f_n(x) \le \dots$ 且 $f_n(x) \xrightarrow{a.e.} f(x)$，则：

    \[
    \lim_{n \to \infty} \int_E f_n(x) dx = \int_E f(x) dx
    \]

    * **2. 法图 (Fatou) 引理**：若 $\{f_n\}$ 是 $E$ 上的非负可测函数序列，则：

    \[
    \int_E \liminf_{n \to \infty} f_n(x) dx \le \liminf_{n \to \infty} \int_E f_n(x) dx
    \]

    * **3. 勒贝格 (Lebesgue) 控制收敛定理 (LDCT)**：设 $\{f_n\}$ 为可测函数序列且 $f_n(x) \xrightarrow{a.e.} f(x)$。若存在一个**可积**的非负函数 $F(x)$ 满足 $|f_n(x)| \le F(x)$ a.e. 对所有 $n$ 成立，则 $f$ 亦可积，且：

    \[
    \lim_{n \to \infty} \int_E f_n(x) dx = \int_E f(x) dx
    \]

---

## 五、 积分空间与微分理论 (Spaces and Differentiation)

!!! info "积分号下求导与 Fubini 定理"
    * **Fubini (富比尼) 定理**：若 $f(x, y)$ 在 $\mathbb{R}^p \times \mathbb{R}^q$ 的可测集 $E_X \times E_Y$ 上非负可测（或可积），则重积分等于累次积分：

    \[
    \int_{E_X \times E_Y} f(x, y) dxdy = \int_{E_X} \left( \int_{E_Y} f(x, y) dy \right) dx = \int_{E_Y} \left( \int_{E_X} f(x, y) dx \right) dy
    \]

!!! success "有界变差函数与绝对连续函数 (Absolute Continuity)"
    * **有界变差函数 ($BV$)**：若 $f$ 在 $[a, b]$ 上的总变差 $V_a^b(f) < \infty$，则 $f$ 可表示为两个单调递增函数的差（Jordan 分解），且 $f$ 几乎处处可导。
    
    * **绝对连续函数 ($AC$)**：若对任意 $\epsilon > 0$，存在 $\delta > 0$，只要 $[a, b]$ 上任意有限个互不相交的开区间 $\{(a_i, b_i)\}_{i=1}^m$ 满足 $\sum_{i=1}^m (b_i - a_i) < \delta$，就有：

    \[
    \sum_{i=1}^m |f(b_i) - f(a_i)| < \epsilon
    \]

    * **微积分基本定理的充要条件**：公式 $f(x) - f(a) = \int_a^x f'(t) dt$ 成立的**充要条件**是 $f(x)$ 在 $[a, b]$ 上绝对连续。