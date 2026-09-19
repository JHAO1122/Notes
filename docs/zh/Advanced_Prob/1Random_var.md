# 第一章：概率空间与随机变量

本章回顾概率论的基本语言。我们首先从样本空间、事件域和概率测度出发，定义概率空间；随后介绍事件的独立性与乘积概率空间；最后给出随机变量、诱导分布和分布函数的定义及基本性质。

本课程后续将进一步讨论随机变量的数字特征（期望、方差和特征函数）以及经典极限定理。

---

## 1. 概率空间

### 1.1 样本空间与事件

!!! info "定义 1.1（样本空间）"

    随机试验所有可能结果组成的集合称为**样本空间（sample space）**，记为 $\Omega$。样本空间中的元素 $\omega\in\Omega$ 称为样本点。

样本空间的子集称为事件。为了对事件赋予概率，需要先规定一族允许讨论的事件，这便引出了 $\sigma$-代数。

### 1.2 $\sigma$-代数

!!! info "定义 1.2（$\sigma$-代数）"

    设 $\mathcal{A}$ 是 $\Omega$ 的若干子集组成的集合族。如果满足：

    1. $\varnothing,\Omega\in\mathcal{A}$；
    2. 若 $A\in\mathcal{A}$，则 $A^c\in\mathcal{A}$；
    3. 若 $A_1,A_2,\ldots\in\mathcal{A}$，则

       \[
       \bigcup_{n=1}^{\infty}A_n\in\mathcal{A},
       \]

    则称 $\mathcal{A}$ 是 $\Omega$ 上的一个 **$\sigma$-代数（$\sigma$-field）**，其中的元素称为事件。

由补集和可数并的封闭性，根据 De Morgan 公式还可得到：若 $A_1,A_2,\ldots\in\mathcal{A}$，则

\[
\bigcap_{n=1}^{\infty}A_n\in\mathcal{A}.
\]

!!! info "定义 1.3（生成的 $\sigma$-代数）"

    设 $\mathcal{C}$ 是 $\Omega$ 的一个子集族。包含 $\mathcal{C}$ 的最小 $\sigma$-代数称为由 $\mathcal{C}$ **生成的 $\sigma$-代数**，记为

    \[
    \sigma(\mathcal{C}).
    \]

    特别地，由事件 $A_1,\ldots,A_n$ 生成的 $\sigma$-代数记为

    \[
    \sigma(A_1,\ldots,A_n).
    \]

### 1.3 概率测度

!!! info "定义 1.4（概率测度）"

    设 $\mathcal{A}$ 是 $\Omega$ 上的 $\sigma$-代数。映射

    \[
    P:\mathcal{A}\longrightarrow[0,1]
    \]

    如果满足：

    1. $P(\Omega)=1$；
    2. 对任意两两不交的事件 $A_1,A_2,\ldots\in\mathcal{A}$，有

       \[
       P\left(\bigcup_{n=1}^{\infty}A_n\right)
       =\sum_{n=1}^{\infty}P(A_n),
       \]

    则称 $P$ 是 $(\Omega,\mathcal{A})$ 上的一个**概率测度**。

三元组

\[
(\Omega,\mathcal{A},P)
\]

称为一个**概率空间**。

!!! example "例 1.5（古典概率模型）"

    设样本空间有限：

    \[
    \Omega=\{\omega_1,\ldots,\omega_N\},
    \]

    并且每个样本点等可能，即

    \[
    P(\{\omega_i\})=\frac{1}{N},
    \qquad i=1,\ldots,N.
    \]

    则对任意事件 $A\subseteq\Omega$，

    \[
    P(A)=\frac{\# A}{N},
    \]

    其中 $\#A$ 表示事件 $A$ 中样本点的个数。

---

## 2. 事件的独立性

### 2.1 两个事件的独立性

!!! info "定义 2.1（两个事件独立）"

    事件 $A,B\in\mathcal{A}$ 称为相互独立，如果

    \[
    P(A\cap B)=P(A)P(B).
    \]

当 $P(B)>0$ 时，上述条件等价于

\[
P(A\mid B)=P(A).
\]

也就是说，已知 $B$ 发生不会改变事件 $A$ 发生的概率。

### 2.2 三个事件的相互独立

!!! info "定义 2.2（三个事件相互独立）"

    事件 $A,B,C\in\mathcal{A}$ 称为**相互独立**，如果

    \[
    P(A\cap B)=P(A)P(B),
    \]

    \[
    P(A\cap C)=P(A)P(C),
    \]

    \[
    P(B\cap C)=P(B)P(C),
    \]

    并且

    \[
    P(A\cap B\cap C)=P(A)P(B)P(C).
    \]

!!! warning "相互独立与两两独立"

    前三个等式只说明 $A,B,C$ 两两独立。两两独立一般不能推出

    \[
    P(A\cap B\cap C)=P(A)P(B)P(C),
    \]

    因而不能推出三个事件相互独立。

---

## 3. 乘积概率空间

设

\[
(\Omega_1,\mathcal{A}_1,P_1)
\quad\text{和}\quad
(\Omega_2,\mathcal{A}_2,P_2)
\]

是两个概率空间。

### 3.1 乘积样本空间

两个样本空间的笛卡尔积为

\[
\Omega_1\times\Omega_2
=\{(\omega_1,\omega_2):\omega_1\in\Omega_1,\ \omega_2\in\Omega_2\}.
\]

### 3.2 乘积 $\sigma$-代数

!!! info "定义 3.1（乘积 $\sigma$-代数）"

    由所有可测矩形 $A_1\times A_2$ 生成的 $\sigma$-代数称为乘积 $\sigma$-代数，记为

    \[
    \mathcal{A}_1\otimes\mathcal{A}_2
    =\sigma\bigl(\{A_1\times A_2:A_1\in\mathcal{A}_1,\ A_2\in\mathcal{A}_2\}\bigr).
    \]

### 3.3 乘积概率测度

!!! info "定义 3.2（乘积概率测度）"

    乘积概率测度 $P_1\otimes P_2$ 是定义在 $\mathcal{A}_1\otimes\mathcal{A}_2$ 上，并且对任意可测矩形满足

    \[
    (P_1\otimes P_2)(A_1\times A_2)
    =P_1(A_1)P_2(A_2)
    \]

    的概率测度。

因此，两个概率空间的乘积可写为

\[
(\Omega_1\times\Omega_2,
\mathcal{A}_1\otimes\mathcal{A}_2,
P_1\otimes P_2).
\]

---

## 4. 随机变量与诱导分布

### 4.1 Borel $\sigma$-代数

记 $\mathcal{B}(\mathbb{R})$ 为实数轴上的 Borel $\sigma$-代数，即由所有开集生成的 $\sigma$-代数。它也可以由形如 $(-\infty,x]$ 的半直线生成。

例如，开半直线可以表示为

\[
(-\infty,x)
=\bigcup_{n=1}^{\infty}
\left(-\infty,x-\frac{1}{n}\right].
\]

### 4.2 随机变量

!!! info "定义 4.1（随机变量）"

    设 $(\Omega,\mathcal{A},P)$ 是概率空间。映射

    \[
    X:\Omega\longrightarrow\mathbb{R}
    \]

    如果对任意 $B\in\mathcal{B}(\mathbb{R})$ 都满足

    \[
    X^{-1}(B)
    =\{\omega\in\Omega:X(\omega)\in B\}
    \in\mathcal{A},
    \]

    则称 $X$ 是一个**随机变量**，也称 $X$ 是从 $(\Omega,\mathcal{A})$ 到 $(\mathbb{R},\mathcal{B}(\mathbb{R}))$ 的可测映射。

随机变量并不是“随机取值的变量”这一朴素概念，而是把抽象样本点映射为实数的可测函数。

### 4.3 随机变量的诱导分布

!!! info "定义 4.2（诱导分布）"

    随机变量 $X$ 在 $\mathbb{R}$ 上诱导的概率测度定义为

    \[
    P_X(B)
    =(P\circ X^{-1})(B)
    =P\bigl(X^{-1}(B)\bigr),
    \qquad B\in\mathcal{B}(\mathbb{R}).
    \]

    $P_X$ 称为 $X$ 的**分布**或**概率律**。

---

## 5. 分布函数

### 5.1 定义

在诱导分布的定义中取

\[
B=(-\infty,x],
\]

便得到随机变量的分布函数。

!!! info "定义 5.1（分布函数）"

    随机变量 $X$ 的**分布函数（cumulative distribution function）**定义为

    \[
    F_X(x)
    =P_X((-\infty,x])
    =P\{\omega:X(\omega)\leq x\}
    =P(X\leq x),
    \qquad x\in\mathbb{R}.
    \]

    因而

    \[
    F_X:\mathbb{R}\longrightarrow[0,1].
    \]

### 5.2 分布函数的基本性质

!!! success "命题 5.2（分布函数的性质）"

    任意分布函数 $F_X$ 都满足：

    1. **单调不减：**若 $x_1<x_2$，则

       \[
       F_X(x_1)\leq F_X(x_2).
       \]

    2. **右连续：**对任意 $x\in\mathbb{R}$，

       \[
       \lim_{h\downarrow0}F_X(x+h)=F_X(x).
       \]

    3. **两端极限：**

       \[
       \lim_{x\to-\infty}F_X(x)=0,
       \qquad
       \lim_{x\to+\infty}F_X(x)=1.
       \]

??? proof "命题 5.2 的说明（点击展开）"

    若 $x_1<x_2$，则

    \[
    \{X\leq x_1\}\subseteq\{X\leq x_2\},
    \]

    因此由概率测度的单调性可得

    \[
    F_X(x_1)\leq F_X(x_2).
    \]

    对任意满足 $x_n\downarrow x$ 的数列，有

    \[
    \{X\leq x_n\}\downarrow\{X\leq x\}.
    \]

    由概率测度对递减事件列的连续性，得到

    \[
    F_X(x_n)\longrightarrow F_X(x),
    \]

    即 $F_X$ 右连续。

    最后，由

    \[
    \{X\leq x\}\downarrow\varnothing
    \quad (x\to-\infty)
    \]

    以及

    \[
    \{X\leq x\}\uparrow\Omega
    \quad (x\to+\infty),
    \]

    可分别得到两端极限。$\square$

---

## 6. 本节小结

本节建立了概率论的基本框架：

1. 概率空间由样本空间、$\sigma$-代数和概率测度组成；
2. 独立性可以通过联合事件的概率乘法公式刻画；
3. 乘积概率空间用于描述由多个随机试验组成的联合试验；
4. 随机变量是从样本空间到实数轴的可测映射；
5. 随机变量通过 $P\circ X^{-1}$ 在实数轴上诱导分布，其分布函数为 $F_X(x)=P(X\leq x)$。
