# Chapter 8: Orthogonal Decomposition and Hilbert Space Theory

This chapter will deeply investigate the core topological properties of Hilbert spaces, focusing on the orthogonal decomposition theorem, orthonormal systems, Fourier expansion, and the isometric isomorphism between separable Hilbert spaces and the $l^2$ space.

---

## 1. Orthogonal Decomposition Theorem

!!! success "Theorem 8.1 (Orthogonal Decomposition Theorem)"

    Let $M$ be a closed subspace of a Hilbert space $H$. Then for any element $x$ in $H$, there is the following unique orthogonal decomposition:

    $$
    x = y + z, \quad y \in M, \quad z \in M^\perp
    $$

    $y$ is called the orthogonal projection of $x$ onto $M$.

??? proof "Proof: Existence and Uniqueness of Orthogonal Decomposition"

    **1. Existence:**
    By assumption, $M$ is a closed subspace of $H$, thus it is a closed convex set. 
    By the theorem from the previous chapter, there exists a unique best approximation element $y$ for $x$ in $M$. Let $\|x - y\| = \alpha$.

    * For a real Hilbert space: Since $y \in M$, for any real number $\lambda$ and any element $u \in M$, we have $y + \lambda u \in M$.
    
    * Therefore, the function $\|-x + y + \lambda u\|^2$ attains its minimum at $\lambda = 0$.
    
    * Taking the derivative with respect to $\lambda$:

    $$
    0 = \frac{d}{d\lambda} \| -x + y + \lambda u \|^2 \Big|_{\lambda=0}
    $$

    $$
    = \frac{d}{d\lambda} \left[ \|y - x\|^2 + \lambda(y - x, u) + \lambda(u, y - x) + \lambda^2\|u\|^2 \right] \Big|_{\lambda=0}
    $$

    $$
    = (y - x, u) + (u, y - x) = 2(y - x, u)
    $$
    
    * Since $u$ is any element in $M$, we have $z = x - y \perp M$. Thus $x = y + z$, with $y \in M$ and $z \in M^\perp$.
    
    * For a complex Hilbert space: Choose $\lambda = -\frac{(y - x, u)}{\|u\|^2}$, and noting that $\|x - y\| = \alpha$, we obtain:

    $$
    \alpha^2 \le \|-x + y + \lambda u\|^2
    $$

    $$
    = \|y - x\|^2 + \overline{\lambda}(y - x, u) + \lambda(u, y - x) + |\lambda|^2\|u\|^2
    $$

    $$
    = \|y - x\|^2 + \lambda(u, y - x) + \overline{\lambda}[(y - x, u) + \lambda\|u\|^2]
    $$

    * Substituting $\lambda$ and simplifying gives:

    $$
    \alpha^2 \le \alpha^2 - \frac{|(y - x, u)|^2}{\|u\|^2}
    $$

    * Therefore, it must be that $(y - x, u) = 0$. Hence $z = x - y \in M^\perp$.

    **2. Uniqueness:**
    
    * Suppose there is another decomposition $x = y' + z'$, where $y' \in M$, $z' \in M^\perp$.
    
    * From $y + z = y' + z'$, we get $y - y' = z' - z$.
    
    * Since $y - y' \in M$, $z' - z \in M^\perp$, and $M \cap M^\perp = \{\theta\}$ (if $x \in M \cap M^\perp$, then $(x, x) = 0$ implying $x = \theta$),
    
    * We have $y - y' = z' - z = \theta$. Therefore $y = y'$, $z = z'$. The orthogonal decomposition is unique. $\square$

---

## 2. Orthonormal Systems and Fourier Coefficients

!!! info "Definition 8.1 (Orthonormal System)"

    Let $\{e_n\}$ be a system of elements in an inner product space satisfying:

    $$
    (e_n, e_m) = \begin{cases} 0, & n \neq m \\ 1, & n = m \end{cases}
    $$

    Then $\{e_n\}$ is called an **orthonormal system** in $H$. 
    
    For any element $x \in H$, $c_n = (x, e_n)$ is called the $n$-th **Fourier coefficient** of $x$ with respect to $\{e_n\}$, simply referred to as the Fourier coefficient, and $\{(x, e_n)\}$ is the Fourier coefficient set of $x$.

### 2.1 Orthogonal Projection onto a Finite-Dimensional Subspace

!!! success "Corollary"

    Let $n$ be a given natural number, $\{e_1, \dots, e_n\}$ be an orthonormal system in an inner product space, and $M$ be the subspace spanned by $\{e_1, \dots, e_n\}$. Then for any given $x \in H$, the orthogonal projection of $x$ onto $M$ is:

    $$
    y = \sum_{k=1}^n (x, e_k)e_k
    $$

??? proof "Proof"

    * By the orthogonal decomposition theorem, any element $x$ in $H$ has a unique orthogonal decomposition $x = y + z$, where $y \in M$, $z \in M^\perp$.
    
    * Then we have:

    $$
    y = \sum_{k=1}^n (y, e_k)e_k = \sum_{k=1}^n (x, e_k)e_k - \sum_{k=1}^n (z, e_k)e_k = \sum_{k=1}^n (x, e_k)e_k
    $$

    * (Because $z \in M^\perp$, so $(z, e_k) = 0$) $\square$

### 2.2 Classical Examples of Orthonormal Systems

* **Space $L^2[0, 2\pi]$**: The function system $\{ \frac{1}{\sqrt{2\pi}} e^{int} \}$ ($n = 0, \pm 1, \pm 2, \dots$) is an orthonormal system.
    * Proof: When $n \neq m$:

    $$
    \left( \frac{1}{\sqrt{2\pi}} e^{int}, \frac{1}{\sqrt{2\pi}} e^{imt} \right) = \frac{1}{2\pi} \int_0^{2\pi} e^{i(n-m)t} dt = \frac{1}{2\pi i(n-m)} e^{i(n-m)t} \Big|_0^{2\pi} = 0
    $$

* **Space $l^2$**: The system of elements $e_n = (0, \dots, 0, 1, 0, \dots)$ (where the $n$-th position is 1) is an orthonormal system of $l^2$.

* **Space $L^2([-1,1]; \frac{1}{\sqrt{1-t^2}})$ and Chebyshev Polynomials**: Consider the function system $T_n(t) = \cos(n \arccos t)$, where $t \in [-1,1]$.
    * Substituting $\theta = \arccos t$ into the identity $\cos n\theta + \cos(n-2)\theta = 2 \cos \theta \cos(n-1)\theta$, we get the recurrence formula $T_n(t) = 2t T_{n-1}(t) - T_{n-2}(t)$.
    
    * From $T_0(t)=1, T_1(t)=t$ and the recurrence formula, we have $T_2(t)=2t^2-1, T_3(t)=4t^3-3t$. Using mathematical induction, it can be proven that $T_n(t)$ is indeed a polynomial of degree $n$ (Chebyshev polynomial of the first kind).
    
    * Note that $t = \cos \theta$. When $n \neq m$:

    $$
    \int_{-1}^1 \frac{T_n(t)T_m(t)}{\sqrt{1-t^2}} dt = \int_0^\pi \cos n\theta \cdot \cos m\theta \, d\theta = 0
    $$

    * When $n = m$:

    $$
    \int_{-1}^1 \frac{T_n(t)^2}{\sqrt{1-t^2}} dt = \int_0^\pi \cos^2 n\theta \, d\theta = \begin{cases} \frac{\pi}{2}, & n \neq 0 \\ \pi, & n = 0 \end{cases}
    $$

    * Thus, we can normalize to construct the orthonormal system $\tilde{T}_n(t)$, where the coefficient is $\sqrt{\frac{2}{\pi}}$ for $n \ge 1$ and $\sqrt{\frac{1}{\pi}}$ for $n=0$.

---

## 3. Bessel's Inequality and Riesz-Fischer Theorem

!!! success "Theorem 8.2 (Bessel's Inequality)"

    Let $\{e_n\}$ be an orthonormal system in an inner product space. Then for any $x \in H$, the following holds:

    $$
    \sum_{k=1}^\infty |(x, e_k)|^2 \le \|x\|^2
    $$

??? proof "Proof"

    * For any $x \in H$, by the previous corollary, for any given natural number $n$, $y = \sum_{k=1}^n (x, e_k)e_k$ is the orthogonal projection of $x$ onto the subspace spanned by $\{e_1, \dots, e_n\}$.
    
    * Thus $x = y + z$, where $y \perp z$. This immediately implies $\|y\| \le \|x\|$, meaning:

    $$
    \sum_{k=1}^n |(x, e_k)|^2 \le \|x\|^2
    $$

    * Letting $n \to \infty$, we obtain Bessel's inequality. $\square$

From Bessel's inequality, we know that the sequence composed of the Fourier coefficients of any element with respect to an orthonormal system must belong to $l^2$. Conversely, we have the following theorem:

!!! success "Theorem 8.3 (Riesz-Fischer Theorem)"

    Let $\{e_n\}$ be an orthonormal system in a Hilbert space, and let the sequence $\{c_n\} \in l^2$. Then there exists an element $x$ in $H$ such that $\{c_n\}$ is the Fourier coefficient set of this element with respect to $\{e_n\}$, and the following **Parseval's formula** holds:

    $$
    \sum_{k=1}^\infty |(x, e_k)|^2 = \|x\|^2
    $$

??? proof "Proof"

    * Let $x_n = \sum_{k=1}^n c_k e_k$. Suppose natural numbers $m, n$ satisfy $m > n$.

    $$
    \|x_m - x_n\|^2 = \left\| \sum_{k=n+1}^m c_k e_k \right\|^2 = \sum_{k=n+1}^m |c_k|^2
    $$

    * Since $\{c_n\} \in l^2$, as $m, n \to \infty$, the tail series $\sum_{k=n+1}^m |c_k|^2 \to 0$.
    
    * Therefore $\|x_m - x_n\| \to 0$. Thus $\{x_n\}$ is a Cauchy sequence in $H$.
    
    * Since $H$ is complete, there exists an element $x \in H$ such that $x_n \to x$. By the continuity of the inner product, for any given natural number $k_0$, we have:

    $$
    (x_n, e_{k_0}) \to (x, e_{k_0})
    $$

    * On the other hand, when $n \ge k_0$, $(x_n, e_{k_0}) = \left(\sum_{k=1}^n c_k e_k, e_{k_0}\right) = c_{k_0}$.
    
    * So $(x, e_{k_0}) = c_{k_0}$. Since $k_0$ is arbitrarily given, $\{c_n\}$ is exactly the Fourier coefficient set of the element $x$.
    
    * Furthermore, let:

    $$
    \sum_{k=1}^n |(x, e_k)|^2 = \|x_n\|^2
    $$

    * As $n \to \infty$, Parseval's formula holds. $\square$

---

## 4. Completeness and Totalness of Orthonormal Systems

Parseval's identity may not always hold for all elements because the elements in $\{e_n\}$ might "not be enough". If there are "enough" elements, Parseval's identity might hold for all elements.

!!! info "Definition 8.2 (Complete and Total)"

    * **Complete**: If for every $x \in H$, Parseval's identity $\sum_{k=1}^\infty |(x, e_k)|^2 = \|x\|^2$ holds identically, then $\{e_n\}$ is said to be complete.
    
    * **Total**: If for every $x \in H$, $(x, e_k) = 0$ $(k=1,2,\dots)$ implies $x = \theta$, then $\{e_n\}$ is said to be total.

!!! success "Theorem 8.4"

    Let $\{e_n\}$ be an orthonormal system in a Hilbert space. Then the following properties are equivalent:

    * (i) $\{e_n\}$ is complete;
    
    * (ii) For any element $x$ in $H$, the series converges in norm to $x$, i.e., $x = \sum_{k=1}^\infty (x, e_k)e_k$;
    
    * (iii) For any two elements $x, y$ in $H$, the generalized Parseval's identity holds: $(x, y) = \sum_{k=1}^\infty (x, e_k)\overline{(y, e_k)}$, and the right side converges absolutely;
    
    * (iv) $\{e_n\}$ is total.

??? proof "Proof: Equivalence of the Four Properties"

    * **(i) $\implies$ (ii)**:

    $$
    \left\|x - \sum_{k=1}^n (x, e_k)e_k\right\|^2 = \|x\|^2 - \sum_{k=1}^n |(x, e_k)|^2
    $$

    * Taking the limit based on the completeness assumption (Parseval's identity), we get:

    $$
    \lim_{n \to \infty} \left\|x - \sum_{k=1}^n (x, e_k)e_k\right\|^2 = 0
    $$

    * **(ii) $\implies$ (iii)**: For any two elements $x, y$, let $x_n = \sum_{k=1}^n (x, e_k)e_k$ and $y_n = \sum_{k=1}^n (y, e_k)e_k$.
    
    * Noting the orthogonality of $\{e_n\}$, we have:

    $$
    (x_n, y_n) = \sum_{k=1}^n (x, e_k)\overline{(y, e_k)}
    $$

    * From (ii), $x_n \to x, y_n \to y$. By the continuity of the inner product:

    $$
    (x, y) = \lim_{n \to \infty} (x_n, y_n) = \sum_{k=1}^\infty (x, e_k)\overline{(y, e_k)}
    $$

    * By Bessel's inequality, the sequence formed by the Fourier coefficients belongs to $l^2$, so the series on the right side converges absolutely.
    
    * **(iii) $\implies$ (iv)**: Suppose $x \in H$ satisfies $(x, e_k) = 0$ $(k=1,2,\dots)$. From (iii), for any $y \in H$, we have $(x, y) = 0$.
    
    * Letting $y = x$ gives $(x, x) = 0$. Therefore $x = \theta$.
    
    * **(iv) $\implies$ (i)**: By Bessel's inequality, $\{(x, e_k)\} \in l^2$. By the Riesz-Fischer Theorem, there exists $y$ such that $\{(x, e_k)\}$ is the Fourier coefficient set of $y$, and $\sum |(x, e_k)|^2 = \|y\|^2$.
    
    * Notice that $\{(x, e_k)\}$ are also the Fourier coefficients of $x$, thus $(x-y, e_k) = 0$. By the totality assumption of (iv), this implies $x - y = \theta$, i.e., $x=y$.
    
    * Therefore, Parseval's identity of (i) holds: $\sum_{k=1}^\infty |(x, e_k)|^2 = \|x\|^2$. $\square$

* **Corollary 1**: The element that makes Parseval's identity hold in the theorem is unique.
    * (Proof sketch: If both $x$ and $x'$ satisfy it, the series expansion in (ii) shows $x = x'$).

* **Corollary 2**: If there exists a complete orthonormal system in an inner product space, then the space is separable.
    * (Proof sketch: Due to completeness, the subspace spanned by the orthogonal system is dense. Linear combinations with rational coefficients form a countable and dense set, thus it is separable).

* **Corollary 3**: An equivalent condition to prove that an orthonormal system is complete is: the closed subspace generated by $\{e_k\}$ is the entire space.
    * (Proof sketch: Utilize the property that orthogonal projection is the best approximation element in a finite-dimensional space, and bound the distance to show it approaches 0).

---

## 5. Gram-Schmidt Orthogonalization and Space Isomorphism

!!! success "Theorem 8.5 (Gram-Schmidt Orthogonalization Theorem)"

    Let $\{x_n\}$ be a countable subset in an inner product space. Then a complete orthonormal system $\{e_n\}$ can be constructed from $\{x_n\}$ such that the subspace spanned by it is identical to the subspace spanned by $\{x_n\}$.

??? proof "Proof: Orthogonalization Construction Process"

    * Assume $x_1$ is the first non-zero element. Let $e_1 = \frac{x_1}{\|x_1\|}$, then $\|e_1\| = 1$.
    
    * Let $x_2$ be the first element linearly independent of $e_1$. Let $h_2 = x_2 - (x_2, e_1)e_1$.
    
    * Then $h_2 \neq 0$ (otherwise linearly dependent), and it is easy to prove $h_2 \perp e_1$. Let $e_2 = \frac{h_2}{\|h_2\|}$.
    
    * Assume that mutually orthogonal elements $e_1, \dots, e_{k-1}$ with norm 1 have been constructed. Let $x_k$ be the first element linearly independent of them. Let:

    $$
    h_k = x_k - \sum_{j=1}^{k-1} (x_k, e_j)e_j
    $$

    * Then $h_k \neq 0$ and $h_k \perp e_j$. Let $e_k = \frac{h_k}{\|h_k\|}$.
    
    * By induction, an orthonormal system with at most countably many elements is obtained.
    
    * From the construction process, each $e_k$ can be linearly expressed by $x_1, \dots, x_k$, and conversely, each $x_k$ can be linearly expressed by $e_1, \dots, e_k$. Hence the subspaces spanned by both are exactly the same. $\square$

From the above orthogonalization theorem, an important conclusion follows directly: **Any separable inner product space must possess a complete orthonormal system.** ### 5.1 $L^2(\mathbb{R}; e^{-t^2})$ and Hermite Polynomials

* The function family $\{1, t, t^2, \dots\}$ belongs to the weighted space $L^2(\mathbb{R}; e^{-t^2})$. Denote the weight function as $\omega(t) = e^{-t^2}$.
    * Taking higher-order derivatives of $\omega(t)$, by mathematical induction it can be proven that $\omega^{(n)}(t) = y_n(t)e^{-t^2}$, where $y_n(t)$ is a polynomial of degree $n$ with the highest order term's coefficient being $(-2)^n$.
    
    * For any polynomial $u(t)$ of degree less than $n$, applying integration by parts multiple times yields:

    $$
    \int_{-\infty}^{+\infty} y_n(t)e^{-t^2} u(t) dt = \int_{-\infty}^{+\infty} \omega^{(n)}(t)u(t) dt = \dots = (-1)^n \int_{-\infty}^{+\infty} \omega(t)u^{(n)}(t) dt = 0
    $$

    * This shows that $y_n(t)$ is orthogonal to all polynomials of lower degree, meaning it is the product of the orthogonalization of $\{1, t, t^2, \dots\}$.
    
    * Normalization: Calculating the norm yields $\|y_n\|^2 = 2^n n! \sqrt{\pi}$.
    
    * Letting $H_n(t) = \frac{y_n(t)}{\|y_n\|}$, we obtain the famous **Hermite Polynomial** orthonormal system. It can be expressed as:

    $$
    H_n(t) = \frac{1}{(2^n n! \sqrt{\pi})^{\frac{1}{2}}} e^{t^2} \frac{d^n}{dt^n}(e^{-t^2})
    $$

    * Because the polynomial system is dense, this constitutes a complete orthonormal system.

### 5.2 Isometric Isomorphism of Spaces

!!! success "Theorem 8.6"

    **Every real (or complex) separable Hilbert space is isometrically isomorphic to the real (or complex) $l^2$ space.** Therefore, all separable Hilbert spaces are necessarily isometrically isomorphic to each other.

??? proof "Proof: Construction of the Isomorphic Mapping"

    * Let $H$ be a separable Hilbert space, and $\{e_n\}$ be a complete orthonormal system. Define the mapping $T: H \to l^2$ as $Tx = \{c_n\}$, where $\{c_n\}$ is the Fourier coefficient set of $x$.
    
    * **Linearity**: $T(\alpha x + \beta y) = \{\alpha c_n + \beta c_n'\} = \alpha Tx + \beta Ty$.
    
    * **Isometry**: By Parseval's identity, $\|x\|^2 = \sum |c_n|^2 = \|Tx\|^2$, so $T$ is an isometric homomorphic mapping. 
    
    * **Injectivity**: If $Tx = \{0,0,\dots\}$, totality deduces $x = \theta$.
    
    * **Surjectivity**: For any $\{c_n\} \in l^2$, by the Riesz-Fischer Theorem there exists a corresponding $x \in H$, so $T$ is surjective. $\square$

Because the properties of infinite-dimensional spaces are drastically different from finite-dimensional ones, for example, in $l^2$, a bounded closed ball cannot be covered by a finite number of small open balls (the unit closed ball in infinite-dimensional spaces is not compact). 

---

## 6. Chapter Summary

* The inner product can induce a norm, thereby preserving the geometric and topological properties of Euclidean spaces and normed spaces. The inner product defines orthogonality, which is an extremely important structural feature in functional analysis.

* The orthogonal decomposition theorem guarantees that an element can find a unique projection in any closed subspace, which is the cornerstone for establishing the existence of orthonormal systems and related important conclusions.

* "Completeness" and "Totalness" are equivalent in separable Hilbert spaces; Bessel's inequality holds unconditionally, whereas Parseval's identity holds for all elements only when the orthonormal system is complete.

* Core conclusion: A complete orthonormal system exists in a separable Hilbert space, and all separable Hilbert spaces are isometrically isomorphic to each other. This is one of the most beautiful and fundamental facts in Hilbert space theory.