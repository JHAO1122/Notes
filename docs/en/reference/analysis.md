# ♾️ Advanced Mathematical Analysis

This module covers the core topological and analytical theorems of multivariate calculus and function sequences, with a focus on hard analysis tools closely related to real variable functions and functional analysis.

---

## I. Multivariate Differential Calculus and the Implicit Function Theorem

!!! info "Multivariate Total Differential and Continuous Differentiability"
    Let \(f: U \subset \mathbb{R}^n \to \mathbb{R}^m\) be a mapping. If at the point \(x_0 \in U\) there exists a linear mapping \(L: \mathbb{R}^n \to \mathbb{R}^m\) such that

    \[
    \lim_{h \to 0} \frac{\|f(x_0 + h) - f(x_0) - L(h)\|}{\|h\|} = 0
    \]

    then \(f\) is said to be **differentiable** at \(x_0\), and we write \(L = Df(x_0)\) (the Jacobian matrix). If \(Df(x)\) exists and is continuous on the whole set \(U\), then \(f\) is said to be **\(C^1\) continuously differentiable**.

!!! success "Implicit Function Theorem"
    Let \(F: \mathbb{R}^n \times \mathbb{R}^m \to \mathbb{R}^m\) be \(C^1\) on an open set \(U\). Suppose that at a point \((x_0, y_0) \in U\) we have:
    
    1. \(F(x_0, y_0) = 0\)
    
    2. The Jacobian matrix with respect to \(y\) is locally invertible, i.e., \(\det \left( \frac{\partial F}{\partial y} (x_0, y_0) \right) \neq 0\)
    
    Then there exists a neighborhood \(V \subset \mathbb{R}^n\) of \(x_0\) and a unique \(C^1\) mapping \(g: V \to \mathbb{R}^m\) satisfying \(g(x_0) = y_0\), such that for every \(x \in V\),

    \[
    F(x, g(x)) = 0
    \]

    Moreover, the derivative (sensitivity) of the implicit function can be obtained directly via the chain rule:

    \[
    Dg(x) = - \left( \frac{\partial F}{\partial y}(x, g(x)) \right)^{-1} \frac{\partial F}{\partial x}(x, g(x))
    \]

---

## II. Uniform Convergence of Function Sequences and Series

!!! abstract "Definition of Uniform Convergence"
    A sequence of functions \(\{f_n\}\) defined on an interval \(I\) **converges uniformly** to a limit function \(f\) if for every \(\epsilon > 0\) there exists an integer \(N\) (independent of the variable \(x\)) such that for all \(n > N\) and **all** \(x \in I\),

    \[
    |f_n(x) - f(x)| < \epsilon
    \]

    * **Functional analysis perspective**: This is equivalent to convergence of \(\{f_n\}\) in the supremum norm on \(I\): \(\lim_{n \to \infty} \|f_n - f\|_\infty = 0\).

!!! success "Continuity Preservation under Uniform Convergence (Foundation for Completeness of \(C^k[a, b]\))"
    If the sequence \(\{f_n\}\) converges uniformly to \(f\) on \([a, b]\), and each \(f_n\) is **continuous** on \([a, b]\), then the limit function \(f\) is **also continuous** on \([a, b]\).

??? note "Extension: Connection to Completeness of the Function Space \(C^k[a, b]\)"
    In functional analysis, the norm on the Banach space \(C^k[a, b]\) is defined by

    \[
    \|f\|_{C^k} = \sum_{j=0}^k \max_{t \in [a, b]} |f^{(j)}(t)|
    \]

    Proving that this space is complete – i.e., that every Cauchy sequence converges to an element inside the space – relies on **the uniform convergence theorem for derivatives** in mathematical analysis:
    
    If the sequences of derivative functions \(\{f_n^{(j)}\}\) each **converge uniformly** on \([a, b]\) to \(g_j\), then the limit \(f\) of the original sequence satisfies \(f^{(j)} = g_j\). Consequently, by the continuity preservation theorem above, the resulting limit function \(f \in C^k[a, b]\).

---

## III. Interchange of Limits with Integrals and Derivatives

!!! info "Passing the Limit under the Riemann Integral"
    Let \(\{f_n\}\) be a sequence of Riemann integrable functions on \([a, b]\). If \(f_n \xrightarrow{\text{uniformly}} f\), then the limit function \(f\) is also Riemann integrable on \([a, b]\), and the limit and integral can be interchanged:

    \[
    \lim_{n \to \infty} \int_a^b f_n(x) dx = \int_a^b \left( \lim_{n \to \infty} f_n(x) \right) dx = \int_a^b f(x) dx
    \]

    *(Note: In real variable theory, when the Riemann integral is upgraded to the Lebesgue integral, this theorem is completely replaced by the more general Lebesgue Dominated Convergence Theorem (LDCT).)*

!!! success "Differentiation under the Riemann Integral (Leibniz Integral Rule)"
    Let \(f(x, y)\) and its partial derivative \(\frac{\partial f}{\partial y}(x, y)\) be continuous on the rectangle \([a, b] \times [c, d]\). Define the parameter‑dependent integral \(F(y) = \int_a^b f(x, y) dx\). Then \(F(y)\) is differentiable and satisfies

    \[
    \frac{d}{dy} \int_a^b f(x, y) dx = \int_a^b \frac{\partial f}{\partial y}(x, y) dx
    \]

---

## IV. Compactness Theory and the Arzelà–Ascoli Theorem

This is one of the most advanced theorems in mathematical analysis and is closely linked to functional analysis (especially the study of compact operators). It provides a criterion for the compactness of subsets in function spaces.

!!! info "Uniform Boundedness and Equicontinuity"
    Let \(\mathcal{F}\) be a family of functions defined on a compact interval \([a, b]\):
    
    * **Uniformly Bounded**: There exists a constant \(M > 0\) such that for all \(f \in \mathcal{F}\) and all \(t \in [a, b]\), \(|f(t)| \le M\).
    
    * **Equicontinuous**: For every \(\epsilon > 0\), there exists a \(\delta > 0\) **depending only on \(\epsilon\)** (and not on the particular function \(f\) or point \(t\)) such that whenever \(|t_1 - t_2| < \delta\),

    \[
    |f(t_1) - f(t_2)| < \epsilon \quad (\forall f \in \mathcal{F})
    \]

!!! success "Arzelà–Ascoli Theorem"
    Let \(\mathcal{F} \subset C[a, b]\). Then \(\mathcal{F}\) is **relatively compact** (i.e., every sequence in \(\mathcal{F}\) has a subsequence that converges uniformly) **if and only if** \(\mathcal{F}\) is **uniformly bounded** and **equicontinuous** on \([a, b]\).