# Chapter 1: Linear Algebra Foundations and Matrix Norms

Matrix computation studies how to perform linear algebra operations in an efficient and stable manner. This chapter first reviews basic concepts such as linear independence, subspaces, bases, ranges, null spaces, and rank, then introduces the Sherman–Morrison–Woodbury formula under low-rank corrections, orthogonality and orthogonal complements, and finally systematically discusses vector norms and matrix norms.

The focus of this chapter is not only to memorize the definitions of various norms, but also to understand their relationships and how to use norms to quantify errors and perturbations in matrix computations.

---

## 1. Linear Algebra Foundations

### 1.1 Linear Independence

!!! info "Definition 1.1 (Linear Independence)"

    Let $a_1,\ldots,a_k\in\mathbb{R}^n$. If

    \[
    \sum_{i=1}^k\beta_i a_i=0
    \quad\Longrightarrow\quad
    \beta_i=0,
    \qquad i=1,\ldots,k,
    \]

    then the vectors $a_1,\ldots,a_k$ are called **linearly independent**; otherwise they are called linearly dependent.

Equivalently, form a matrix with these vectors as columns:

\[
A=(a_1,\ldots,a_k)\in\mathbb{R}^{n\times k}.
\]

Then $a_1,\ldots,a_k$ are linearly independent if and only if the homogeneous equation $A\beta=0$ has only the zero solution, that is,

\[
\operatorname{Null}(A)=\{0\}.
\]

### 1.2 Span and Subspaces

!!! info "Definition 1.2 (Span)"

    The span of the vectors $a_1,\ldots,a_k$ is defined as

    \[
    \operatorname{span}\{a_1,\ldots,a_k\}
    =\left\{\sum_{i=1}^k\beta_i a_i:\beta_i\in\mathbb{R}\right\}.
    \]

    It is a linear subspace of $\mathbb{R}^n$.

If

\[
\operatorname{span}\{a_1,\ldots,a_k\}
=\operatorname{span}\{b_1,\ldots,b_r\},
\]

then the two sets of vectors span the same subspace, but their numbers of vectors and linear dependence relations need not be the same.

### 1.3 Bases and Dimension

!!! info "Definition 1.3 (Basis)"

    Let $S\subseteq\mathbb{R}^n$ be a linear subspace. If the vectors $b_1,\ldots,b_r$ are linearly independent and

    \[
    S=\operatorname{span}\{b_1,\ldots,b_r\},
    \]

    then $b_1,\ldots,b_r$ are called a **basis** of $S$.

Any two bases of the same finite-dimensional subspace contain the same number of vectors. This number is called the dimension of the subspace and is denoted by $\dim(S)$.

### 1.4 Range and Null Space

!!! info "Definition 1.4 (Range)"

    For $A\in\mathbb{R}^{m\times n}$, the range or column space of the matrix $A$ is defined as

    \[
    \operatorname{Range}(A)
    =\{Ax:x\in\mathbb{R}^n\}
    \subseteq\mathbb{R}^m.
    \]

    If $A=(a_1,\ldots,a_n)$, then

    \[
    \operatorname{Range}(A)
    =\operatorname{span}\{a_1,\ldots,a_n\}.
    \]

!!! info "Definition 1.5 (Null Space)"

    The null space of the matrix $A$ is defined as

    \[
    \operatorname{Null}(A)
    =\{x\in\mathbb{R}^n:Ax=0\}.
    \]

Both the range and the null space are linear subspaces.

### 1.5 Rank and the Rank–Nullity Theorem

!!! info "Definition 1.6 (Rank)"

    The rank of a matrix $A$ is defined as the dimension of its range:

    \[
    \operatorname{rank}(A)
    =\dim\operatorname{Range}(A).
    \]

    It is also equal to the maximum number of linearly independent columns of $A$.

!!! success "Theorem 1.7 (Rank–Nullity Theorem)"

    If $A\in\mathbb{R}^{m\times n}$, then

    \[
    \dim\operatorname{Range}(A)
    +\dim\operatorname{Null}(A)
    =n.
    \]

??? proof "Proof of Theorem 1.7 (click to expand)"

    Let $v_1,\ldots,v_s$ be a basis of $\operatorname{Null}(A)$, and extend it to a basis of $\mathbb{R}^n$

    \[
    v_1,\ldots,v_s,v_{s+1},\ldots,v_n.
    \]

    We prove that $Av_{s+1},\ldots,Av_n$ form a basis of $\operatorname{Range}(A)$.

    Take any $x\in\mathbb{R}^n$ and write

    \[
    x=\sum_{i=1}^n\alpha_i v_i.
    \]

    Since $Av_i=0$ for $i=1,\ldots,s$, we have

    \[
    Ax=\sum_{i=s+1}^n\alpha_iAv_i.
    \]

    Thus these vectors span $\operatorname{Range}(A)$. If

    \[
    \sum_{i=s+1}^n\alpha_iAv_i=0,
    \]

    then $\sum_{i=s+1}^n\alpha_iv_i\in\operatorname{Null}(A)$. By the linear independence of the whole set $v_1,\ldots,v_n$, all $\alpha_i=0$, so $Av_{s+1},\ldots,Av_n$ are linearly independent.

    Therefore

    \[
    \dim\operatorname{Range}(A)=n-s,
    \qquad
    \dim\operatorname{Null}(A)=s,
    \]

    and their sum equals $n$. $\square$

---

## 2. The Sherman–Morrison–Woodbury Formula

In matrix computation, we often already know $A^{-1}$ but need the inverse of a low-rank correction matrix $A+UV^T$. The Sherman–Morrison–Woodbury formula transforms the problem of inverting an $n\times n$ matrix into the problem of inverting a $k\times k$ matrix; when $k\ll n$, this transformation is very useful.

!!! info "Lemma 2.1 (Difference of Two Inverse Matrices)"

    If $A$ and $B$ are both invertible, then

    \[
    A^{-1}-B^{-1}
    =A^{-1}(B-A)B^{-1}
    =B^{-1}(B-A)A^{-1}.
    \]

??? proof "Proof of Lemma 2.1 (click to expand)"

    Direct expansion gives

    \[
    A^{-1}(B-A)B^{-1}
    =A^{-1}BB^{-1}-A^{-1}AB^{-1}
    =A^{-1}-B^{-1}.
    \]

    The other equality follows similarly. $\square$

!!! success "Theorem 2.2 (Sherman–Morrison–Woodbury Formula)"

    Let $A\in\mathbb{R}^{n\times n}$ be invertible and let $U,V\in\mathbb{R}^{n\times k}$. Then $A+UV^T$ is invertible if and only if $I_k+V^TA^{-1}U$ is invertible. In this case,

    \[
    (A+UV^T)^{-1}
    =A^{-1}
    -A^{-1}U(I_k+V^TA^{-1}U)^{-1}V^TA^{-1}.
    \]

??? proof "Proof of Theorem 2.2 (click to expand)"

    First note that

    \[
    A+UV^T
    =A(I_n+A^{-1}UV^T).
    \]

    The Sylvester determinant identity gives

    \[
    \det(I_n+A^{-1}UV^T)
    =\det(I_k+V^TA^{-1}U),
    \]

    so the two matrices are either both invertible or both singular.

    Let

    \[
    M=I_k+V^TA^{-1}U.
    \]

    Directly verify the candidate inverse on the right-hand side:

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

    Similarly, the product in the reverse direction also equals $I_n$, so the formula holds. $\square$

!!! success "Corollary 2.3 (Sherman–Morrison Formula)"

    When $k=1$, write $U=u$ and $V=v$. If

    \[
    \alpha=1+v^TA^{-1}u\neq0,
    \]

    then

    \[
    (A+uv^T)^{-1}
    =A^{-1}-\frac{A^{-1}uv^TA^{-1}}{1+v^TA^{-1}u}.
    \]

---

## 3. Orthogonality

### 3.1 Orthogonality and Orthonormality

!!! info "Definition 3.1 (Orthogonal Vectors)"

    Vectors $x_1,\ldots,x_k\in\mathbb{R}^n$ are called orthogonal if

    \[
    x_i^Tx_j=0,
    \qquad i\neq j.
    \]

    If in addition

    \[
    x_i^Tx_j=\delta_{ij},
    \]

    then they are called orthonormal.

Any orthogonal set of vectors containing no zero vector is linearly independent.

### 3.2 Orthogonal Complement

!!! info "Definition 3.2 (Orthogonal Complement)"

    Let $S\subseteq\mathbb{R}^n$ be a linear subspace. The orthogonal complement of $S$ is defined as

    \[
    S^\perp
    =\{y\in\mathbb{R}^n:y^Tx=0\text{ for all }x\in S\}.
    \]

The orthogonal complement $S^\perp$ is also a linear subspace of $\mathbb{R}^n$, and

\[
\mathbb{R}^n=S\oplus S^\perp,
\qquad
\dim(S)+\dim(S^\perp)=n.
\]

!!! success "Proposition 3.3 (Orthogonal Relation Between Range and Null Space)"

    For any $A\in\mathbb{R}^{m\times n}$,

    \[
    \operatorname{Range}(A)^\perp
    =\operatorname{Null}(A^T),
    \qquad
    \operatorname{Range}(A^T)^\perp
    =\operatorname{Null}(A).
    \]

??? proof "Proof of Proposition 3.3 (click to expand)"

    For $y\in\mathbb{R}^m$,

    \[
    \begin{aligned}
    y\in\operatorname{Range}(A)^\perp
    &\Longleftrightarrow y^TAx=0
    \text{ for all }x\in\mathbb{R}^n\\
    &\Longleftrightarrow A^Ty=0\\
    &\Longleftrightarrow y\in\operatorname{Null}(A^T).
    \end{aligned}
    \]

    The second equality follows by applying the first equality to $A^T$. $\square$

### 3.3 Orthogonal Matrices and Orthonormal Basis Completion

!!! info "Definition 3.4 (Orthogonal Matrix)"

    If a square matrix $Q\in\mathbb{R}^{n\times n}$ satisfies

    \[
    Q^TQ=QQ^T=I_n,
    \]

    then $Q$ is called an orthogonal matrix. Equivalently, the columns of $Q$ form an orthonormal basis of $\mathbb{R}^n$.

An orthogonal matrix satisfies

\[
Q^{-1}=Q^T,
\]

and preserves the Euclidean inner product and the $2$-norm:

\[
(Qx)^T(Qy)=x^Ty,
\qquad
\lVert Qx\rVert_2=\lVert x\rVert_2.
\]

!!! success "Theorem 3.5 (Completion of an Orthonormal Basis)"

    Let $U_1\in\mathbb{R}^{n\times r}$, where $r<n$, and suppose the columns of $U_1$ are orthonormal, that is,

    \[
    U_1^TU_1=I_r.
    \]

    Then there exists $U_2\in\mathbb{R}^{n\times(n-r)}$ such that

    \[
    Q=[U_1\ U_2]\in\mathbb{R}^{n\times n}
    \]

    is an orthogonal matrix, and

    \[
    \operatorname{Range}(U_1)^\perp
    =\operatorname{Range}(U_2).
    \]

??? proof "Proof of Theorem 3.5 (click to expand)"

    The columns of $U_1$ form an orthonormal basis of the subspace $\operatorname{Range}(U_1)$. Choose an orthonormal basis in its orthogonal complement, and take these vectors as the columns of $U_2$. Since

    \[
    \mathbb{R}^n
    =\operatorname{Range}(U_1)
    \oplus\operatorname{Range}(U_1)^\perp,
    \]

    all columns of $Q=[U_1\ U_2]$ form an orthonormal basis of $\mathbb{R}^n$, so $Q$ is an orthogonal matrix. $\square$

By this theorem, every $x\in\mathbb{R}^n$ has the orthogonal decomposition

\[
x=U_1U_1^Tx+U_2U_2^Tx,
\]

where the two components belong respectively to $\operatorname{Range}(U_1)$ and $\operatorname{Range}(U_1)^\perp$.

---

## 4. Vector Norms

### 4.1 Definition of a Vector Norm

!!! info "Definition 4.1 (Vector Norm)"

    A function $\lVert\cdot\rVert:\mathbb{R}^n\to\mathbb{R}$ is called a vector norm if for any $x,y\in\mathbb{R}^n$ and any $\alpha\in\mathbb{R}$, it satisfies:

    1. **Positive definiteness:** $\lVert x\rVert\geq0$, and $\lVert x\rVert=0$ if and only if $x=0$;
    2. **Absolute homogeneity:** $\lVert\alpha x\rVert=|\alpha|\lVert x\rVert$;
    3. **Triangle inequality:** $\lVert x+y\rVert\leq\lVert x\rVert+\lVert y\rVert$.

### 4.2 Common $p$-Norms

For $1\leq p<\infty$, define

\[
\lVert x\rVert_p
=\left(\sum_{i=1}^n|x_i|^p\right)^{1/p}.
\]

In particular,

\[
\lVert x\rVert_1=\sum_{i=1}^n|x_i|,
\qquad
\lVert x\rVert_2=\left(\sum_{i=1}^n|x_i|^2\right)^{1/2}.
\]

When $p=\infty$, define

\[
\lVert x\rVert_\infty
=\max_{1\leq i\leq n}|x_i|.
\]

### 4.3 Young's Inequality and Hölder's Inequality

!!! info "Lemma 4.2 (Young's Inequality)"

    Let $a,b\geq0$, $p,q>1$, and

    \[
    \frac{1}{p}+\frac{1}{q}=1.
    \]

    Then

    \[
    ab\leq\frac{a^p}{p}+\frac{b^q}{q}.
    \]

??? proof "Proof of Lemma 4.2 (click to expand)"

    Apply the weighted Jensen inequality to the convex function $f(t)=e^t$, and take

    \[
    \lambda=\frac{1}{p},
    \qquad
    1-\lambda=\frac{1}{q},
    \]

    to obtain

    \[
    ab
    =\exp(\log a+\log b)
    \leq\frac{a^p}{p}+\frac{b^q}{q}.
    \]

    When $a=0$ or $b=0$, the conclusion is obvious. $\square$

!!! success "Theorem 4.3 (Hölder's Inequality)"

    Let $x,y\in\mathbb{R}^n$, $1\leq p,q\leq\infty$, and

    \[
    \frac{1}{p}+\frac{1}{q}=1.
    \]

    Then

    \[
    |x^Ty|
    \leq\lVert x\rVert_p\lVert y\rVert_q.
    \]

??? proof "Proof of Theorem 4.3 (click to expand)"

    First consider $1<p,q<\infty$, with $x,y\neq0$. Let

    \[
    \bar{x}=\frac{x}{\lVert x\rVert_p},
    \qquad
    \bar{y}=\frac{y}{\lVert y\rVert_q}.
    \]

    Applying Young's inequality to each $i$, we obtain

    \[
    |\bar{x}_i\bar{y}_i|
    \leq\frac{|\bar{x}_i|^p}{p}
    +\frac{|\bar{y}_i|^q}{q}.
    \]

    Summing gives

    \[
    \sum_{i=1}^n|\bar{x}_i\bar{y}_i|
    \leq\frac{1}{p}\sum_{i=1}^n|\bar{x}_i|^p
    +\frac{1}{q}\sum_{i=1}^n|\bar{y}_i|^q
    =\frac{1}{p}+\frac{1}{q}=1.
    \]

    Therefore

    \[
    |x^Ty|
    \leq\sum_{i=1}^n|x_iy_i|
    \leq\lVert x\rVert_p\lVert y\rVert_q.
    \]

    When $(p,q)=(1,\infty)$ or $(\infty,1)$, the conclusion follows directly from

    \[
    \sum_i|x_iy_i|
    \leq\left(\sum_i|x_i|\right)\max_i|y_i|.
    \]

    $\square$

When $p=q=2$, Hölder's inequality becomes the Cauchy–Schwarz inequality:

\[
|x^Ty|\leq\lVert x\rVert_2\lVert y\rVert_2.
\]

### 4.4 Minkowski's Inequality

!!! success "Theorem 4.4 (Minkowski's Inequality)"

    For any $1\leq p\leq\infty$ and $x,y\in\mathbb{R}^n$,

    \[
    \lVert x+y\rVert_p
    \leq\lVert x\rVert_p+\lVert y\rVert_p.
    \]

    Therefore, $\lVert\cdot\rVert_p$ indeed satisfies the triangle inequality and is a vector norm.

??? proof "Proof of Theorem 4.4 (click to expand)"

    When $p=1$ or $p=\infty$, the conclusion follows directly from the triangle inequality for real numbers.

    Let $1<p<\infty$ and set $q=p/(p-1)$. If $x+y=0$, the conclusion is obvious. Otherwise, by Hölder's inequality,

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

    Dividing both sides by $\lVert x+y\rVert_p^{p-1}$ gives the result. $\square$

### 4.5 Equivalence of Norms in Finite-Dimensional Spaces

!!! info "Definition 4.5 (Equivalent Norms)"

    Let $\lVert\cdot\rVert_\alpha$ and $\lVert\cdot\rVert_\beta$ be two norms on $\mathbb{R}^n$. If there exist constants $c_1,c_2>0$ independent of $x$ such that for all $x\in\mathbb{R}^n$,

    \[
    c_1\lVert x\rVert_\alpha
    \leq\lVert x\rVert_\beta
    \leq c_2\lVert x\rVert_\alpha,
    \]

    then the two norms are called equivalent.

!!! success "Theorem 4.6 (All Norms Are Equivalent in Finite-Dimensional Spaces)"

    Any two vector norms on $\mathbb{R}^n$ are equivalent.

??? proof "Proof of Theorem 4.6 (click to expand)"

    It suffices to prove that any norm $\lVert\cdot\rVert_\alpha$ is equivalent to the $2$-norm.

    Let $e_1,\ldots,e_n$ be the standard basis. For any $x=\sum_i x_ie_i$, by the triangle inequality and the Cauchy–Schwarz inequality,

    \[
    \lVert x\rVert_\alpha
    \leq\sum_{i=1}^n|x_i|\lVert e_i\rVert_\alpha
    \leq\sqrt{n}\max_i\lVert e_i\rVert_\alpha\lVert x\rVert_2.
    \]

    Therefore, $\lVert\cdot\rVert_\alpha$ is continuous with respect to the $2$-norm. Consider the compact set

    \[
    S^{n-1}=\{x\in\mathbb{R}^n:\lVert x\rVert_2=1\}.
    \]

    The continuous function $x\mapsto\lVert x\rVert_\alpha$ attains a maximum $M$ and a minimum $m$ on $S^{n-1}$. By positive definiteness of the norm and $0\notin S^{n-1}$, we have $m>0$. Hence, for any nonzero $x$,

    \[
    m
    \leq\left\lVert\frac{x}{\lVert x\rVert_2}\right\rVert_\alpha
    \leq M.
    \]

    Using homogeneity, we obtain

    \[
    m\lVert x\rVert_2
    \leq\lVert x\rVert_\alpha
    \leq M\lVert x\rVert_2.
    \]

    Thus every norm is equivalent to the $2$-norm, and hence any two norms are equivalent to each other. $\square$

For $1\leq p\leq q\leq\infty$, a common explicit estimate is

\[
\lVert x\rVert_q
\leq\lVert x\rVert_p
\leq n^{1/p-1/q}\lVert x\rVert_q.
\]

In particular,

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

and

\[
\lVert x\rVert_\infty
\leq\lVert x\rVert_1
\leq n\lVert x\rVert_\infty.
\]

!!! note "Norm Equivalence and Convergence"

    In finite-dimensional spaces, if a sequence of vectors converges in one norm, then it converges in every norm. Therefore, convergence of finite-dimensional vector sequences is independent of the chosen norm.

    This conclusion generally cannot be extended to infinite-dimensional spaces. For example, let $y^{(n)}$ have its first $n$ components all equal to $1/n$ and its remaining components zero. Then

    \[
    \lVert y^{(n)}\rVert_2=\frac{1}{\sqrt{n}}\longrightarrow0,
    \qquad
    \lVert y^{(n)}\rVert_1=1.
    \]

---

## 5. Matrix Norms

### 5.1 Matrix Norms and Consistency

!!! info "Definition 5.1 (Matrix Norm)"

    A function $\lVert\cdot\rVert:\mathbb{R}^{m\times n}\to\mathbb{R}$ is a vector norm on the matrix space if it satisfies positive definiteness, absolute homogeneity, and the triangle inequality.

    For matrices whose dimensions are compatible for multiplication, if it also satisfies

    \[
    \lVert AB\rVert\leq\lVert A\rVert\lVert B\rVert,
    \]

    then the matrix norm is called **submultiplicative**, or said to have consistency.

Submultiplicativity is not an automatic consequence of the first three norm axioms. For example, define

\[
\lVert A\rVert_{\max}=\max_{i,j}|a_{ij}|.
\]

It satisfies the three axioms of a vector norm, but it is not submultiplicative. Take

\[
A=\begin{pmatrix}1&1\\1&1\end{pmatrix},
\]

then

\[
\lVert A\rVert_{\max}=1,
\qquad
\lVert A^2\rVert_{\max}=2
>\lVert A\rVert_{\max}^2.
\]

### 5.2 Induced Matrix Norms

!!! info "Definition 5.2 (Induced Norm)"

    Given vector norms on $\mathbb{R}^n$ and $\mathbb{R}^m$, the induced norm of a matrix $A\in\mathbb{R}^{m\times n}$ is defined as

    \[
    \lVert A\rVert
    =\sup_{x\neq0}\frac{\lVert Ax\rVert}{\lVert x\rVert}
    =\sup_{\lVert x\rVert=1}\lVert Ax\rVert.
    \]

    When $p$-norms are used on both sides, it is denoted by

    \[
    \lVert A\rVert_p
    =\sup_{x\neq0}\frac{\lVert Ax\rVert_p}{\lVert x\rVert_p}.
    \]

By definition, the induced norm satisfies the subordinate property

\[
\lVert Ax\rVert_p
\leq\lVert A\rVert_p\lVert x\rVert_p,
\]

and submultiplicativity

\[
\lVert AB\rVert_p
\leq\lVert A\rVert_p\lVert B\rVert_p.
\]

??? proof "Proof of Submultiplicativity of the Induced Norm (click to expand)"

    For any $x\neq0$,

    \[
    \lVert ABx\rVert_p
    \leq\lVert A\rVert_p\lVert Bx\rVert_p
    \leq\lVert A\rVert_p\lVert B\rVert_p\lVert x\rVert_p.
    \]

    Dividing both sides by $\lVert x\rVert_p$ and taking the supremum over all $x\neq0$ gives the result. $\square$

### 5.3 Common Matrix Norms

!!! info "Definition 5.3 (Frobenius Norm)"

    For $A=(a_{ij})\in\mathbb{R}^{m\times n}$, the Frobenius norm is defined as

    \[
    \lVert A\rVert_F
    =\left(\sum_{i=1}^m\sum_{j=1}^na_{ij}^2\right)^{1/2}
    =\sqrt{\operatorname{tr}(A^TA)}.
    \]

The Frobenius norm is not an operator norm induced by vector norms of the same dimension, but it satisfies submultiplicativity.

!!! success "Proposition 5.4 (Explicit Formulas for Common Induced Norms)"

    For $A=(a_{ij})\in\mathbb{R}^{m\times n}$,

    \[
    \lVert A\rVert_1
    =\max_{1\leq j\leq n}\sum_{i=1}^m|a_{ij}|,
    \]

    that is, the maximum absolute column sum;

    \[
    \lVert A\rVert_\infty
    =\max_{1\leq i\leq m}\sum_{j=1}^n|a_{ij}|,
    \]

    that is, the maximum absolute row sum; and

    \[
    \lVert A\rVert_2
    =\sqrt{\lambda_{\max}(A^TA)}.
    \]

??? proof "Proof of Proposition 5.4 (click to expand)"

    For the $1$-norm,

    \[
    \begin{aligned}
    \lVert Ax\rVert_1
    &=\sum_{i=1}^m\left|\sum_{j=1}^na_{ij}x_j\right|\\
    &\leq\sum_{j=1}^n\left(\sum_{i=1}^m|a_{ij}|\right)|x_j|\\
    &\leq\left(\max_j\sum_i|a_{ij}|\right)\lVert x\rVert_1.
    \end{aligned}
    \]

    Taking $x=e_{j_*}$, where $j_*$ is the column attaining the maximum column sum, achieves equality.

    For the $\infty$-norm,

    \[
    \lVert Ax\rVert_\infty
    \leq\left(\max_i\sum_j|a_{ij}|\right)\lVert x\rVert_\infty.
    \]

    For the $i_*$-th row attaining the maximum row sum, taking $x_j=\operatorname{sgn}(a_{i_*j})$ achieves equality.

    Finally, $A^TA$ is a symmetric positive semidefinite matrix, so

    \[
    \begin{aligned}
    \lVert A\rVert_2^2
    &=\sup_{x\neq0}\frac{\lVert Ax\rVert_2^2}{\lVert x\rVert_2^2}\\
    &=\sup_{x\neq0}\frac{x^TA^TAx}{x^Tx}\\
    &=\lambda_{\max}(A^TA).
    \end{aligned}
    \]

    This completes the proof. $\square$

### 5.4 Relationships Among Common Matrix Norms

For $A\in\mathbb{R}^{m\times n}$,

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

and

\[
\lVert A\rVert_2^2
\leq\lVert A\rVert_1\lVert A\rVert_\infty.
\]

!!! note "Meaning of the Dimension Factors"

    Factors such as $\sqrt{m}$ and $\sqrt{n}$ in the above inequalities come from the equivalence relations among finite-dimensional vector norms. They show that although different norms describe the same finite-dimensional topology, dimension factors may not be negligible when performing quantitative error estimates.

### 5.5 Quantifying Perturbations with Norms

!!! success "Lemma 5.5 (Neumann Series)"

    Let $A\in\mathbb{R}^{n\times n}$ and suppose that for some submultiplicative matrix norm,

    \[
    \lVert A\rVert<1.
    \]

    Then $I-A$ is invertible, and

    \[
    (I-A)^{-1}=\sum_{k=0}^{\infty}A^k,
    \qquad
    \lVert(I-A)^{-1}\rVert
    \leq\frac{1}{1-\lVert A\rVert}.
    \]

??? proof "Proof of Lemma 5.5 (click to expand)"

    Let

    \[
    S_N=\sum_{k=0}^NA^k.
    \]

    By submultiplicativity,

    \[
    \sum_{k=0}^{\infty}\lVert A^k\rVert
    \leq\sum_{k=0}^{\infty}\lVert A\rVert^k
    =\frac{1}{1-\lVert A\rVert}<\infty.
    \]

    Therefore $S_N$ converges to some matrix $S$. On the other hand,

    \[
    S_N(I-A)=(I-A)S_N=I-A^{N+1}.
    \]

    Since $\lVert A^{N+1}\rVert\leq\lVert A\rVert^{N+1}\to0$, letting $N\to\infty$ gives

    \[
    S(I-A)=(I-A)S=I.
    \]

    Hence $S=(I-A)^{-1}$. The norm estimate follows directly from the geometric series. $\square$

This immediately gives

\[
\lVert(I-A)^{-1}-I\rVert
\leq\frac{\lVert A\rVert}{1-\lVert A\rVert}.
\]

!!! success "Corollary 5.6 (Sensitivity of the Inverse Matrix to Perturbations)"

    Let $A$ be nonsingular, and let

    \[
    r=\lVert A^{-1}E\rVert<1.
    \]

    Then $A+E$ is nonsingular, and

    \[
    \lVert(A+E)^{-1}-A^{-1}\rVert
    \leq
    \frac{\lVert E\rVert\lVert A^{-1}\rVert^2}{1-r}.
    \]

??? proof "Proof of Corollary 5.6 (click to expand)"

    Since

    \[
    A+E=A(I+A^{-1}E),
    \]

    and $\lVert A^{-1}E\rVert<1$, the Neumann series shows that $I+A^{-1}E$ is invertible, so $A+E$ is invertible.

    Using Lemma 2.1,

    \[
    (A+E)^{-1}-A^{-1}
    =-(A+E)^{-1}EA^{-1}.
    \]

    Also,

    \[
    \lVert(A+E)^{-1}\rVert
    \leq\frac{\lVert A^{-1}\rVert}{1-r}.
    \]

    Combining this with submultiplicativity gives the desired estimate. $\square$

### 5.6 Orthogonal Invariance

!!! success "Theorem 5.7 (Orthogonal Invariance of the Frobenius Norm and Spectral Norm)"

    Let $A\in\mathbb{R}^{m\times n}$, and let $Q\in\mathbb{R}^{m\times m}$ and $Z\in\mathbb{R}^{n\times n}$ be orthogonal matrices. Then

    \[
    \lVert QAZ\rVert_F=\lVert A\rVert_F,
    \qquad
    \lVert QAZ\rVert_2=\lVert A\rVert_2.
    \]

??? proof "Proof of Theorem 5.7 (click to expand)"

    For the Frobenius norm, using the cyclic invariance of the trace,

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

    For the spectral norm, orthogonal matrices preserve the $2$-norm, so

    \[
    \begin{aligned}
    \lVert QAZ\rVert_2
    &=\sup_{x\neq0}\frac{\lVert QAZx\rVert_2}{\lVert x\rVert_2}\\
    &=\sup_{x\neq0}\frac{\lVert AZx\rVert_2}{\lVert Zx\rVert_2}\\
    &=\lVert A\rVert_2.
    \end{aligned}
    \]

    This completes the proof. $\square$

---

## 6. Summary of This Chapter

1. Range, null space, and rank describe the basic structure of a linear map, and the rank–nullity theorem gives the relationship between their dimensions.
2. The Sherman–Morrison–Woodbury formula can efficiently handle low-rank corrections of matrices.
3. An orthonormal set of vectors can be completed to an orthogonal matrix, and a subspace together with its orthogonal complement gives the orthogonal direct sum decomposition of $\mathbb{R}^n$.
4. Hölder's inequality and Minkowski's inequality are the foundations of $p$-norm theory.
5. All norms on a finite-dimensional space are equivalent, so convergence is independent of the specific norm, but quantitative estimates may still be affected by dimension factors.
6. Induced matrix norms have the subordinate property and submultiplicativity; the $1$-norm, $\infty$-norm, and $2$-norm correspond respectively to the maximum column sum, maximum row sum, and spectral norm.
7. The Neumann series can be used to determine invertibility after perturbation and to give error bounds for the inverse matrix.