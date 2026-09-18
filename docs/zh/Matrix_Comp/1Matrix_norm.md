# 第一章：线性代数基础与矩阵范数

矩阵计算研究如何以高效、稳定的方式完成线性代数运算。本章首先回顾线性无关、子空间、基、值域、零空间和秩等基本概念，随后介绍低秩修正下的 Sherman–Morrison–Woodbury 公式、正交性与正交补，最后系统讨论向量范数和矩阵范数。

本章的重点不仅是记住各种范数的定义，还要理解它们之间的关系，以及如何用范数量化矩阵运算中的误差和扰动。

---

## 1. 线性代数基础

### 1.1 线性无关

!!! info "定义 1.1（线性无关）"

    设 $a_1,\ldots,a_k\in\mathbb{R}^n$。如果

    \[
    \sum_{i=1}^k\beta_i a_i=0
    \quad\Longrightarrow\quad
    \beta_i=0,
    \qquad i=1,\ldots,k,
    \]

    则称向量组 $a_1,\ldots,a_k$ **线性无关（linearly independent）**；否则称其线性相关。

等价地，将这些向量作为列组成矩阵

\[
A=(a_1,\ldots,a_k)\in\mathbb{R}^{n\times k},
\]

则 $a_1,\ldots,a_k$ 线性无关，当且仅当齐次方程 $A\beta=0$ 只有零解，即

\[
\operatorname{Null}(A)=\{0\}.
\]

### 1.2 张成空间与子空间

!!! info "定义 1.2（张成空间）"

    向量组 $a_1,\ldots,a_k$ 的张成空间定义为

    \[
    \operatorname{span}\{a_1,\ldots,a_k\}
    =\left\{\sum_{i=1}^k\beta_i a_i:\beta_i\in\mathbb{R}\right\}.
    \]

    它是 $\mathbb{R}^n$ 的一个线性子空间。

若

\[
\operatorname{span}\{a_1,\ldots,a_k\}
=\operatorname{span}\{b_1,\ldots,b_r\},
\]

则两组向量张成同一个子空间，但它们的向量数目和线性相关性未必相同。

### 1.3 基与维数

!!! info "定义 1.3（基）"

    设 $S\subseteq\mathbb{R}^n$ 是线性子空间。如果向量组 $b_1,\ldots,b_r$ 线性无关，并且

    \[
    S=\operatorname{span}\{b_1,\ldots,b_r\},
    \]

    则称 $b_1,\ldots,b_r$ 是 $S$ 的一组**基（basis）**。

同一个有限维子空间的任意两组基含有相同数目的向量。这个数目称为子空间的维数，记作 $\dim(S)$。

### 1.4 值域与零空间

!!! info "定义 1.4（值域）"

    对 $A\in\mathbb{R}^{m\times n}$，矩阵 $A$ 的值域或列空间定义为

    \[
    \operatorname{Range}(A)
    =\{Ax:x\in\mathbb{R}^n\}
    \subseteq\mathbb{R}^m.
    \]

    若 $A=(a_1,\ldots,a_n)$，则

    \[
    \operatorname{Range}(A)
    =\operatorname{span}\{a_1,\ldots,a_n\}.
    \]

!!! info "定义 1.5（零空间）"

    矩阵 $A$ 的零空间定义为

    \[
    \operatorname{Null}(A)
    =\{x\in\mathbb{R}^n:Ax=0\}.
    \]

值域和零空间都是线性子空间。

### 1.5 秩与秩-零化度定理

!!! info "定义 1.6（秩）"

    矩阵 $A$ 的秩定义为其值域的维数：

    \[
    \operatorname{rank}(A)
    =\dim\operatorname{Range}(A).
    \]

    它也等于 $A$ 的线性无关列的最大数目。

!!! success "定理 1.7（秩-零化度定理）"

    若 $A\in\mathbb{R}^{m\times n}$，则

    \[
    \dim\operatorname{Range}(A)
    +\dim\operatorname{Null}(A)
    =n.
    \]

??? proof "定理 1.7 的证明（点击展开）"

    设 $v_1,\ldots,v_s$ 是 $\operatorname{Null}(A)$ 的一组基，并将其扩充为 $\mathbb{R}^n$ 的一组基

    \[
    v_1,\ldots,v_s,v_{s+1},\ldots,v_n.
    \]

    我们证明 $Av_{s+1},\ldots,Av_n$ 是 $\operatorname{Range}(A)$ 的一组基。

    任取 $x\in\mathbb{R}^n$，可写成

    \[
    x=\sum_{i=1}^n\alpha_i v_i.
    \]

    因为 $Av_i=0$ 对 $i=1,\ldots,s$ 成立，所以

    \[
    Ax=\sum_{i=s+1}^n\alpha_iAv_i.
    \]

    因而这些向量张成 $\operatorname{Range}(A)$。若

    \[
    \sum_{i=s+1}^n\alpha_iAv_i=0,
    \]

    则 $\sum_{i=s+1}^n\alpha_iv_i\in\operatorname{Null}(A)$。由整组 $v_1,\ldots,v_n$ 的线性无关性可知所有 $\alpha_i=0$，所以 $Av_{s+1},\ldots,Av_n$ 线性无关。

    因此

    \[
    \dim\operatorname{Range}(A)=n-s,
    \qquad
    \dim\operatorname{Null}(A)=s,
    \]

    两者相加等于 $n$。$\square$

---

## 2. Sherman–Morrison–Woodbury 公式

在矩阵计算中，我们经常已经知道 $A^{-1}$，但需要求低秩修正矩阵 $A+UV^T$ 的逆。Sherman–Morrison–Woodbury 公式将一个 $n\times n$ 矩阵的求逆问题转化为一个 $k\times k$ 矩阵的求逆问题；当 $k\ll n$ 时，这一转化非常有用。

!!! info "引理 2.1（两个逆矩阵之差）"

    若 $A$ 和 $B$ 均可逆，则

    \[
    A^{-1}-B^{-1}
    =A^{-1}(B-A)B^{-1}
    =B^{-1}(B-A)A^{-1}.
    \]

??? proof "引理 2.1 的证明（点击展开）"

    直接展开可得

    \[
    A^{-1}(B-A)B^{-1}
    =A^{-1}BB^{-1}-A^{-1}AB^{-1}
    =A^{-1}-B^{-1}.
    \]

    另一个等式同理。$\square$

!!! success "定理 2.2（Sherman–Morrison–Woodbury 公式）"

    设 $A\in\mathbb{R}^{n\times n}$ 可逆，$U,V\in\mathbb{R}^{n\times k}$。则 $A+UV^T$ 可逆，当且仅当 $I_k+V^TA^{-1}U$ 可逆。此时

    \[
    (A+UV^T)^{-1}
    =A^{-1}
    -A^{-1}U(I_k+V^TA^{-1}U)^{-1}V^TA^{-1}.
    \]

??? proof "定理 2.2 的证明（点击展开）"

    首先注意

    \[
    A+UV^T
    =A(I_n+A^{-1}UV^T).
    \]

    Sylvester 行列式恒等式给出

    \[
    \det(I_n+A^{-1}UV^T)
    =\det(I_k+V^TA^{-1}U),
    \]

    因而两个矩阵同时可逆或同时不可逆。

    令

    \[
    M=I_k+V^TA^{-1}U.
    \]

    直接验证右侧候选逆矩阵：

    \[
    \begin{aligned}
    &(A+UV^T)
    \left(A^{-1}-A^{-1}UM^{-1}V^TA^{-1}\right)\\
    &\quad=I_n+UV^TA^{-1}
    -UM^{-1}V^TA^{-1}
    -UV^TA^{-1}UM^{-1}V^TA^{-1}\\
    &\quad=I_n+U\left[I_k-M^{-1}-(M-I_k)M^{-1}\right]V^TA^{-1}\\
    &\quad=I_n.
    \end{aligned}
    \]

    同理可验证反方向的乘积也等于 $I_n$，从而公式成立。$\square$

!!! success "推论 2.3（Sherman–Morrison 公式）"

    当 $k=1$ 时，记 $U=u$、$V=v$。若

    \[
    \alpha=1+v^TA^{-1}u\neq0,
    \]

    则

    \[
    (A+uv^T)^{-1}
    =A^{-1}-\frac{A^{-1}uv^TA^{-1}}{1+v^TA^{-1}u}.
    \]

---

## 3. 正交性

### 3.1 正交与标准正交

!!! info "定义 3.1（正交向量组）"

    向量组 $x_1,\ldots,x_k\in\mathbb{R}^n$ 称为正交向量组，如果

    \[
    x_i^Tx_j=0,
    \qquad i\neq j.
    \]

    若进一步满足

    \[
    x_i^Tx_j=\delta_{ij},
    \]

    则称其为标准正交向量组。

任意不含零向量的正交向量组都是线性无关的。

### 3.2 正交补

!!! info "定义 3.2（正交补）"

    设 $S\subseteq\mathbb{R}^n$ 是线性子空间。$S$ 的正交补定义为

    \[
    S^\perp
    =\{y\in\mathbb{R}^n:y^Tx=0\text{ 对所有 }x\in S\}.
    \]

正交补 $S^\perp$ 也是 $\mathbb{R}^n$ 的线性子空间，并且

\[
\mathbb{R}^n=S\oplus S^\perp,
\qquad
\dim(S)+\dim(S^\perp)=n.
\]

!!! success "命题 3.3（值域与零空间的正交关系）"

    对任意 $A\in\mathbb{R}^{m\times n}$，

    \[
    \operatorname{Range}(A)^\perp
    =\operatorname{Null}(A^T),
    \qquad
    \operatorname{Range}(A^T)^\perp
    =\operatorname{Null}(A).
    \]

??? proof "命题 3.3 的证明（点击展开）"

    对 $y\in\mathbb{R}^m$，

    \[
    \begin{aligned}
    y\in\operatorname{Range}(A)^\perp
    &\Longleftrightarrow y^TAx=0
    \text{ 对所有 }x\in\mathbb{R}^n\\
    &\Longleftrightarrow A^Ty=0\\
    &\Longleftrightarrow y\in\operatorname{Null}(A^T).
    \end{aligned}
    \]

    第二个等式对 $A^T$ 应用第一个等式即可。$\square$

### 3.3 正交矩阵与正交基补全

!!! info "定义 3.4（正交矩阵）"

    若方阵 $Q\in\mathbb{R}^{n\times n}$ 满足

    \[
    Q^TQ=QQ^T=I_n,
    \]

    则称 $Q$ 为正交矩阵。等价地，$Q$ 的列向量构成 $\mathbb{R}^n$ 的一组标准正交基。

正交矩阵满足

\[
Q^{-1}=Q^T,
\]

并保持 Euclidean 内积与 $2$-范数：

\[
(Qx)^T(Qy)=x^Ty,
\qquad
\lVert Qx\rVert_2=\lVert x\rVert_2.
\]

!!! success "定理 3.5（标准正交基的补全）"

    设 $U_1\in\mathbb{R}^{n\times r}$，其中 $r<n$，且 $U_1$ 的列向量标准正交，即

    \[
    U_1^TU_1=I_r.
    \]

    则存在 $U_2\in\mathbb{R}^{n\times(n-r)}$，使得

    \[
    Q=[U_1\ U_2]\in\mathbb{R}^{n\times n}
    \]

    是正交矩阵，并且

    \[
    \operatorname{Range}(U_1)^\perp
    =\operatorname{Range}(U_2).
    \]

??? proof "定理 3.5 的证明（点击展开）"

    $U_1$ 的列向量构成子空间 $\operatorname{Range}(U_1)$ 的一组标准正交基。在其正交补中选取一组标准正交基，并将这些向量作为 $U_2$ 的列。由于

    \[
    \mathbb{R}^n
    =\operatorname{Range}(U_1)
    \oplus\operatorname{Range}(U_1)^\perp,
    \]

    $Q=[U_1\ U_2]$ 的全部列构成 $\mathbb{R}^n$ 的标准正交基，因此 $Q$ 是正交矩阵。$\square$

由该定理，每个 $x\in\mathbb{R}^n$ 都具有正交分解

\[
x=U_1U_1^Tx+U_2U_2^Tx,
\]

其中两个分量分别属于 $\operatorname{Range}(U_1)$ 和 $\operatorname{Range}(U_1)^\perp$。

---

## 4. 向量范数

### 4.1 向量范数的定义

!!! info "定义 4.1（向量范数）"

    函数 $\lVert\cdot\rVert:\mathbb{R}^n\to\mathbb{R}$ 称为一个向量范数，如果对任意 $x,y\in\mathbb{R}^n$ 和任意 $\alpha\in\mathbb{R}$，满足：

    1. **正定性：** $\lVert x\rVert\geq0$，且 $\lVert x\rVert=0$ 当且仅当 $x=0$；
    2. **绝对齐次性：** $\lVert\alpha x\rVert=|\alpha|\lVert x\rVert$；
    3. **三角不等式：** $\lVert x+y\rVert\leq\lVert x\rVert+\lVert y\rVert$。

### 4.2 常用的 $p$-范数

对 $1\leq p<\infty$，定义

\[
\lVert x\rVert_p
=\left(\sum_{i=1}^n|x_i|^p\right)^{1/p}.
\]

特别地，

\[
\lVert x\rVert_1=\sum_{i=1}^n|x_i|,
\qquad
\lVert x\rVert_2=\left(\sum_{i=1}^n|x_i|^2\right)^{1/2}.
\]

当 $p=\infty$ 时，定义

\[
\lVert x\rVert_\infty
=\max_{1\leq i\leq n}|x_i|.
\]

### 4.3 Young 不等式与 Hölder 不等式

!!! info "引理 4.2（Young 不等式）"

    设 $a,b\geq0$，$p,q>1$ 且

    \[
    \frac{1}{p}+\frac{1}{q}=1.
    \]

    则

    \[
    ab\leq\frac{a^p}{p}+\frac{b^q}{q}.
    \]

??? proof "引理 4.2 的证明（点击展开）"

    对凸函数 $f(t)=e^t$ 使用加权 Jensen 不等式，并取

    \[
    \lambda=\frac{1}{p},
    \qquad
    1-\lambda=\frac{1}{q},
    \]

    可得

    \[
    ab
    =\exp(\log a+\log b)
    \leq\frac{a^p}{p}+\frac{b^q}{q}.
    \]

    当 $a=0$ 或 $b=0$ 时结论显然成立。$\square$

!!! success "定理 4.3（Hölder 不等式）"

    设 $x,y\in\mathbb{R}^n$，$1\leq p,q\leq\infty$ 且

    \[
    \frac{1}{p}+\frac{1}{q}=1.
    \]

    则

    \[
    |x^Ty|
    \leq\lVert x\rVert_p\lVert y\rVert_q.
    \]

??? proof "定理 4.3 的证明（点击展开）"

    先考虑 $1<p,q<\infty$，且 $x,y\neq0$。令

    \[
    \bar{x}=\frac{x}{\lVert x\rVert_p},
    \qquad
    \bar{y}=\frac{y}{\lVert y\rVert_q}.
    \]

    对每个 $i$ 应用 Young 不等式，得到

    \[
    |\bar{x}_i\bar{y}_i|
    \leq\frac{|\bar{x}_i|^p}{p}
    +\frac{|\bar{y}_i|^q}{q}.
    \]

    求和可得

    \[
    \sum_{i=1}^n|\bar{x}_i\bar{y}_i|
    \leq\frac{1}{p}\sum_{i=1}^n|\bar{x}_i|^p
    +\frac{1}{q}\sum_{i=1}^n|\bar{y}_i|^q
    =\frac{1}{p}+\frac{1}{q}=1.
    \]

    因此

    \[
    |x^Ty|
    \leq\sum_{i=1}^n|x_iy_i|
    \leq\lVert x\rVert_p\lVert y\rVert_q.
    \]

    当 $(p,q)=(1,\infty)$ 或 $(\infty,1)$ 时，结论由

    \[
    \sum_i|x_iy_i|
    \leq\left(\sum_i|x_i|\right)\max_i|y_i|
    \]

    直接得到。$\square$

当 $p=q=2$ 时，Hölder 不等式就是 Cauchy–Schwarz 不等式：

\[
|x^Ty|\leq\lVert x\rVert_2\lVert y\rVert_2.
\]

### 4.4 Minkowski 不等式

!!! success "定理 4.4（Minkowski 不等式）"

    对任意 $1\leq p\leq\infty$ 和 $x,y\in\mathbb{R}^n$，有

    \[
    \lVert x+y\rVert_p
    \leq\lVert x\rVert_p+\lVert y\rVert_p.
    \]

    因此，$\lVert\cdot\rVert_p$ 确实满足三角不等式，是一个向量范数。

??? proof "定理 4.4 的证明（点击展开）"

    当 $p=1$ 或 $p=\infty$ 时，结论由实数的三角不等式直接得到。

    设 $1<p<\infty$，并令 $q=p/(p-1)$。若 $x+y=0$，结论显然成立。否则，由 Hölder 不等式，

    \[
    \begin{aligned}
    \lVert x+y\rVert_p^p
    &=\sum_{i=1}^n|x_i+y_i|^p\\
    &\leq\sum_{i=1}^n|x_i|\,|x_i+y_i|^{p-1}
    +\sum_{i=1}^n|y_i|\,|x_i+y_i|^{p-1}\\
    &\leq(\lVert x\rVert_p+\lVert y\rVert_p)
    \left(\sum_{i=1}^n|x_i+y_i|^{(p-1)q}\right)^{1/q}\\
    &=(\lVert x\rVert_p+\lVert y\rVert_p)
    \lVert x+y\rVert_p^{p-1}.
    \end{aligned}
    \]

    两边除以 $\lVert x+y\rVert_p^{p-1}$，即得结论。$\square$

### 4.5 有限维空间中的范数等价

!!! info "定义 4.5（等价范数）"

    设 $\lVert\cdot\rVert_\alpha$ 和 $\lVert\cdot\rVert_\beta$ 是 $\mathbb{R}^n$ 上的两个范数。若存在与 $x$ 无关的常数 $c_1,c_2>0$，使得对所有 $x\in\mathbb{R}^n$，

    \[
    c_1\lVert x\rVert_\alpha
    \leq\lVert x\rVert_\beta
    \leq c_2\lVert x\rVert_\alpha,
    \]

    则称这两个范数等价。

!!! success "定理 4.6（有限维空间中所有范数等价）"

    $\mathbb{R}^n$ 上任意两个向量范数都等价。

??? proof "定理 4.6 的证明（点击展开）"

    只需证明任意范数 $\lVert\cdot\rVert_\alpha$ 都与 $2$-范数等价。

    设 $e_1,\ldots,e_n$ 是标准基。对任意 $x=\sum_i x_ie_i$，由三角不等式和 Cauchy–Schwarz 不等式，

    \[
    \lVert x\rVert_\alpha
    \leq\sum_{i=1}^n|x_i|\lVert e_i\rVert_\alpha
    \leq\sqrt{n}\max_i\lVert e_i\rVert_\alpha\lVert x\rVert_2.
    \]

    因此，$\lVert\cdot\rVert_\alpha$ 关于 $2$-范数连续。考虑紧集

    \[
    S^{n-1}=\{x\in\mathbb{R}^n:\lVert x\rVert_2=1\}.
    \]

    连续函数 $x\mapsto\lVert x\rVert_\alpha$ 在 $S^{n-1}$ 上取得最大值 $M$ 和最小值 $m$。由范数的正定性以及 $0\notin S^{n-1}$，有 $m>0$。因此，对任意非零 $x$，

    \[
    m
    \leq\left\lVert\frac{x}{\lVert x\rVert_2}\right\rVert_\alpha
    \leq M.
    \]

    利用齐次性，得到

    \[
    m\lVert x\rVert_2
    \leq\lVert x\rVert_\alpha
    \leq M\lVert x\rVert_2.
    \]

    所以任意范数都与 $2$-范数等价，进而任意两个范数彼此等价。$\square$

对于 $1\leq p\leq q\leq\infty$，常用的显式估计为

\[
\lVert x\rVert_q
\leq\lVert x\rVert_p
\leq n^{1/p-1/q}\lVert x\rVert_q.
\]

特别地，

\[
\lVert x\rVert_2
\leq\lVert x\rVert_1
\leq\sqrt{n}\lVert x\rVert_2,
\]

\[
\lVert x\rVert_\infty
\leq\lVert x\rVert_2
\leq\sqrt{n}\lVert x\rVert_\infty,
\]

以及

\[
\lVert x\rVert_\infty
\leq\lVert x\rVert_1
\leq n\lVert x\rVert_\infty.
\]

!!! note "范数等价与收敛"

    在有限维空间中，如果向量序列在某一个范数下收敛，那么它在任意范数下都收敛。因此，有限维向量序列的收敛性与所选范数无关。

    这一结论一般不能推广到无限维空间。例如，令 $y^{(n)}$ 的前 $n$ 个分量均为 $1/n$，其余分量为零，则

    \[
    \lVert y^{(n)}\rVert_2=\frac{1}{\sqrt{n}}\longrightarrow0,
    \qquad
    \lVert y^{(n)}\rVert_1=1.
    \]

---

## 5. 矩阵范数

### 5.1 矩阵范数与相容性

!!! info "定义 5.1（矩阵范数）"

    函数 $\lVert\cdot\rVert:\mathbb{R}^{m\times n}\to\mathbb{R}$ 若满足正定性、绝对齐次性和三角不等式，则是矩阵空间上的一个向量范数。

    对于尺寸允许相乘的矩阵，如果还满足

    \[
    \lVert AB\rVert\leq\lVert A\rVert\lVert B\rVert,
    \]

    则称该矩阵范数是**次乘的（submultiplicative）**，或称具有相容性。

次乘性并不是前三条范数公理的自动推论。例如，定义

\[
\lVert A\rVert_{\max}=\max_{i,j}|a_{ij}|.
\]

它满足向量范数的三条公理，但不满足次乘性。取

\[
A=\begin{pmatrix}1&1\\1&1\end{pmatrix},
\]

则

\[
\lVert A\rVert_{\max}=1,
\qquad
\lVert A^2\rVert_{\max}=2
>\lVert A\rVert_{\max}^2.
\]

### 5.2 诱导矩阵范数

!!! info "定义 5.2（诱导范数）"

    给定 $\mathbb{R}^n$ 和 $\mathbb{R}^m$ 上的向量范数，矩阵 $A\in\mathbb{R}^{m\times n}$ 的诱导范数定义为

    \[
    \lVert A\rVert
    =\sup_{x\neq0}\frac{\lVert Ax\rVert}{\lVert x\rVert}
    =\sup_{\lVert x\rVert=1}\lVert Ax\rVert.
    \]

    当两端使用 $p$-范数时，记作

    \[
    \lVert A\rVert_p
    =\sup_{x\neq0}\frac{\lVert Ax\rVert_p}{\lVert x\rVert_p}.
    \]

根据定义，诱导范数满足从属性

\[
\lVert Ax\rVert_p
\leq\lVert A\rVert_p\lVert x\rVert_p,
\]

以及次乘性

\[
\lVert AB\rVert_p
\leq\lVert A\rVert_p\lVert B\rVert_p.
\]

??? proof "诱导范数次乘性的证明（点击展开）"

    对任意 $x\neq0$，

    \[
    \lVert ABx\rVert_p
    \leq\lVert A\rVert_p\lVert Bx\rVert_p
    \leq\lVert A\rVert_p\lVert B\rVert_p\lVert x\rVert_p.
    \]

    两边除以 $\lVert x\rVert_p$，再对所有 $x\neq0$ 取上确界，即得结论。$\square$

### 5.3 常用矩阵范数

!!! info "定义 5.3（Frobenius 范数）"

    对 $A=(a_{ij})\in\mathbb{R}^{m\times n}$，Frobenius 范数定义为

    \[
    \lVert A\rVert_F
    =\left(\sum_{i=1}^m\sum_{j=1}^na_{ij}^2\right)^{1/2}
    =\sqrt{\operatorname{tr}(A^TA)}.
    \]

Frobenius 范数不是由同维向量范数诱导出的算子范数，但它满足次乘性。

!!! success "命题 5.4（常用诱导范数的显式公式）"

    对 $A=(a_{ij})\in\mathbb{R}^{m\times n}$，有

    \[
    \lVert A\rVert_1
    =\max_{1\leq j\leq n}\sum_{i=1}^m|a_{ij}|,
    \]

    即最大绝对列和；

    \[
    \lVert A\rVert_\infty
    =\max_{1\leq i\leq m}\sum_{j=1}^n|a_{ij}|,
    \]

    即最大绝对行和；并且

    \[
    \lVert A\rVert_2
    =\sqrt{\lambda_{\max}(A^TA)}.
    \]

??? proof "命题 5.4 的证明（点击展开）"

    对 $1$-范数，

    \[
    \begin{aligned}
    \lVert Ax\rVert_1
    &=\sum_{i=1}^m\left|\sum_{j=1}^na_{ij}x_j\right|\\
    &\leq\sum_{j=1}^n\left(\sum_{i=1}^m|a_{ij}|\right)|x_j|\\
    &\leq\left(\max_j\sum_i|a_{ij}|\right)\lVert x\rVert_1.
    \end{aligned}
    \]

    取 $x=e_{j_*}$，其中 $j_*$ 是达到最大列和的列，即可取到等号。

    对 $\infty$-范数，

    \[
    \lVert Ax\rVert_\infty
    \leq\left(\max_i\sum_j|a_{ij}|\right)\lVert x\rVert_\infty.
    \]

    对达到最大行和的第 $i_*$ 行，取 $x_j=\operatorname{sgn}(a_{i_*j})$，即可取到等号。

    最后，$A^TA$ 是对称半正定矩阵，所以

    \[
    \begin{aligned}
    \lVert A\rVert_2^2
    &=\sup_{x\neq0}\frac{\lVert Ax\rVert_2^2}{\lVert x\rVert_2^2}\\
    &=\sup_{x\neq0}\frac{x^TA^TAx}{x^Tx}\\
    &=\lambda_{\max}(A^TA).
    \end{aligned}
    \]

    证明完毕。$\square$

### 5.4 常用矩阵范数之间的关系

对 $A\in\mathbb{R}^{m\times n}$，有

\[
\lVert A\rVert_2
\leq\lVert A\rVert_F
\leq\sqrt{\min\{m,n\}}\lVert A\rVert_2,
\]

\[
\frac{1}{\sqrt{n}}\lVert A\rVert_\infty
\leq\lVert A\rVert_2
\leq\sqrt{m}\lVert A\rVert_\infty,
\]

\[
\frac{1}{\sqrt{m}}\lVert A\rVert_1
\leq\lVert A\rVert_2
\leq\sqrt{n}\lVert A\rVert_1,
\]

以及

\[
\lVert A\rVert_2^2
\leq\lVert A\rVert_1\lVert A\rVert_\infty.
\]

!!! note "维数因子的意义"

    上述不等式中的 $\sqrt{m}$、$\sqrt{n}$ 等因子来自有限维向量范数之间的等价关系。它们说明不同范数虽然描述同一种有限维拓扑，但在进行定量误差估计时，维数因子可能不可忽略。

### 5.5 用范数量化扰动

!!! success "引理 5.5（Neumann 级数）"

    设 $A\in\mathbb{R}^{n\times n}$，且对某个次乘矩阵范数有

    \[
    \lVert A\rVert<1.
    \]

    则 $I-A$ 可逆，并且

    \[
    (I-A)^{-1}=\sum_{k=0}^{\infty}A^k,
    \qquad
    \lVert(I-A)^{-1}\rVert
    \leq\frac{1}{1-\lVert A\rVert}.
    \]

??? proof "引理 5.5 的证明（点击展开）"

    令

    \[
    S_N=\sum_{k=0}^NA^k.
    \]

    由次乘性，

    \[
    \sum_{k=0}^{\infty}\lVert A^k\rVert
    \leq\sum_{k=0}^{\infty}\lVert A\rVert^k
    =\frac{1}{1-\lVert A\rVert}<\infty.
    \]

    因此 $S_N$ 收敛到某个矩阵 $S$。另一方面，

    \[
    S_N(I-A)=(I-A)S_N=I-A^{N+1}.
    \]

    因为 $\lVert A^{N+1}\rVert\leq\lVert A\rVert^{N+1}\to0$，令 $N\to\infty$ 可得

    \[
    S(I-A)=(I-A)S=I.
    \]

    所以 $S=(I-A)^{-1}$。范数估计由几何级数直接得到。$\square$

由此立即得到

\[
\lVert(I-A)^{-1}-I\rVert
\leq\frac{\lVert A\rVert}{1-\lVert A\rVert}.
\]

!!! success "推论 5.6（逆矩阵对扰动的敏感性）"

    设 $A$ 非奇异，并令

    \[
    r=\lVert A^{-1}E\rVert<1.
    \]

    则 $A+E$ 非奇异，且

    \[
    \lVert(A+E)^{-1}-A^{-1}\rVert
    \leq
    \frac{\lVert E\rVert\lVert A^{-1}\rVert^2}{1-r}.
    \]

??? proof "推论 5.6 的证明（点击展开）"

    因为

    \[
    A+E=A(I+A^{-1}E),
    \]

    且 $\lVert A^{-1}E\rVert<1$，由 Neumann 级数可知 $I+A^{-1}E$ 可逆，因此 $A+E$ 可逆。

    利用引理 2.1，

    \[
    (A+E)^{-1}-A^{-1}
    =-(A+E)^{-1}EA^{-1}.
    \]

    同时

    \[
    \lVert(A+E)^{-1}\rVert
    \leq\frac{\lVert A^{-1}\rVert}{1-r}.
    \]

    结合次乘性即得所需估计。$\square$

### 5.6 正交不变性

!!! success "定理 5.7（Frobenius 范数与谱范数的正交不变性）"

    设 $A\in\mathbb{R}^{m\times n}$，$Q\in\mathbb{R}^{m\times m}$、$Z\in\mathbb{R}^{n\times n}$ 均为正交矩阵，则

    \[
    \lVert QAZ\rVert_F=\lVert A\rVert_F,
    \qquad
    \lVert QAZ\rVert_2=\lVert A\rVert_2.
    \]

??? proof "定理 5.7 的证明（点击展开）"

    对 Frobenius 范数，利用迹的循环不变性，

    \[
    \begin{aligned}
    \lVert QAZ\rVert_F^2
    &=\operatorname{tr}(Z^TA^TQ^TQAZ)\\
    &=\operatorname{tr}(Z^TA^TAZ)\\
    &=\operatorname{tr}(A^TAZZ^T)\\
    &=\operatorname{tr}(A^TA)
    =\lVert A\rVert_F^2.
    \end{aligned}
    \]

    对谱范数，正交矩阵保持 $2$-范数，所以

    \[
    \begin{aligned}
    \lVert QAZ\rVert_2
    &=\sup_{x\neq0}\frac{\lVert QAZx\rVert_2}{\lVert x\rVert_2}\\
    &=\sup_{x\neq0}\frac{\lVert AZx\rVert_2}{\lVert Zx\rVert_2}\\
    &=\lVert A\rVert_2.
    \end{aligned}
    \]

    证明完毕。$\square$

---

## 6. 本章小结

1. 值域、零空间和秩描述线性映射的基本结构，秩-零化度定理给出二者维数之间的关系。
2. Sherman–Morrison–Woodbury 公式能够高效处理矩阵的低秩修正。
3. 标准正交向量组可以补全为正交矩阵，子空间与其正交补给出 $\mathbb{R}^n$ 的正交直和分解。
4. Hölder 不等式和 Minkowski 不等式是 $p$-范数理论的基础。
5. 有限维空间中的所有范数等价，因此收敛性与具体范数无关，但定量估计仍会受到维数因子的影响。
6. 诱导矩阵范数具有从属性和次乘性；$1$-范数、$\infty$-范数和 $2$-范数分别对应最大列和、最大行和及谱范数。
7. Neumann 级数可用于判断扰动后矩阵的可逆性，并给出逆矩阵误差界。
