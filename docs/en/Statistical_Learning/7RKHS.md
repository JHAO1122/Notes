# Chapter 7: Reproducing Kernel Hilbert Space (RKHS) and Kernel Ridge Regression

In Chapter 6, we informally introduced high-dimensional feature mappings via the kernel trick. This chapter rigorously establishes the mathematical foundation of the **Reproducing Kernel Hilbert Space (RKHS)** from the perspective of functional analysis. RKHS constitutes a sufficiently flexible infinite-dimensional function space, while its favorable topological properties ensure that the solution to the infinite-dimensional optimization problem can be reduced to a finite-dimensional space. This is the core theoretical reason for the success of kernel methods in statistical learning.

---

## 1. Hilbert Spaces and Evaluation Functionals

Before delving into reproducing kernels, we need to clarify the completeness of function spaces and the continuity of linear functionals.

!!! info "Definition 1.1 (Hilbert Space)"

    Let \(\mathcal{H}\) be a linear space equipped with an inner product \(\langle \cdot, \cdot \rangle_{\mathcal{H}}\). The norm induced by this inner product is \(\|f\|_{\mathcal{H}} = \sqrt{\langle f, f \rangle_{\mathcal{H}}}\). If \(\mathcal{H}\) is complete under this norm topology (i.e., every Cauchy sequence converges to an element in the space), then \(\mathcal{H}\) is called a **Hilbert space**.

### 1.1 Dirac Evaluation Functional

Let \(\mathcal{H}\) be a space of real-valued functions defined on a set \(\mathcal{X}\). For each fixed point \(x \in \mathcal{X}\), we can define a mapping \(L_x: \mathcal{H} \to \mathbb{R}\) that extracts the value of the function \(f\) at point \(x\).

!!! info "Definition 1.2 (Evaluation Functional)"

    The linear operator \(L_x: \mathcal{H} \to \mathbb{R}\) is called the **Dirac evaluation functional** at point \(x\), defined as:

    \[
    L_x(f) = f(x), \quad \forall f \in \mathcal{H}
    \]

In many traditional Hilbert spaces (e.g., the square-integrable space \(L_2\)), point evaluation functionals are either meaningless or discontinuous (since changing the value at a single point does not affect the \(L_2\) inner product integral). The core feature of RKHS is precisely that point evaluation functionals must be continuous.

---

## 2. Definition of Reproducing Kernel Hilbert Space and the Reproducing Property

!!! info "Definition 2.1 (Reproducing Kernel Hilbert Space, RKHS)"

    A function Hilbert space \(\mathcal{H}\) is called a **Reproducing Kernel Hilbert Space (RKHS)** if and only if for every \(x \in \mathcal{X}\), the corresponding point evaluation functional \(L_x\) is a **bounded (continuous) linear operator** on \(\mathcal{H}\). That is, there exists a constant \(M_x > 0\) such that:

    \[
    |L_x(f)| = |f(x)| \le M_x \|f\|_{\mathcal{H}}, \quad \forall f \in \mathcal{H}
    \]

### 2.1 Introduction of the Reproducing Kernel

By the **Riesz Representation Theorem** in functional analysis, any continuous linear functional defined on a Hilbert space can be represented by a unique element of that space via the inner product.

Since \(L_x\) is continuous in an RKHS, there must exist a unique function (denoted \(K_x \in \mathcal{H}\)) satisfying:

\[
L_x(f) = \langle f, K_x \rangle_{\mathcal{H}}, \quad \forall f \in \mathcal{H}
    \]

Since \(K_x\) itself is a function defined on \(\mathcal{X}\), we can unify the two arguments. Define the bivariate function \(K: \mathcal{X} \times \mathcal{X} \to \mathbb{R}\) as:

\[
K(x, z) = K_x(z)
\]

!!! note "Theorem 2.1 (Reproducing Property of the Kernel)"

    According to the construction above, the bivariate function \(K(\cdot, \cdot)\) satisfies the following two decisive **reproducing properties**:

    1. For any fixed \(x \in \mathcal{X}\), the function \(K(x, \cdot)\) belongs to the space \(\mathcal{H}\).
    2. For any \(x \in \mathcal{X}\) and \(f \in \mathcal{H}\), the following holds:

    \[
    \langle f, K(x, \cdot) \rangle_{\mathcal{H}} = f(x)
    \]

    In particular, replacing \(f\) with \(K(z, \cdot)\) yields:

    \[
    \langle K(x, \cdot), K(z, \cdot) \rangle_{\mathcal{H}} = K(x, z)
    \]

---

## 3. Representer Theorem

In statistical learning, when we seek the optimal function in the infinite-dimensional function space \(\mathcal{H}\), we face an infinite-dimensional optimization problem. The Representer Theorem establishes the computational feasibility of kernel methods.

!!! note "Theorem 3.1 (Representer Theorem)"

    Let \(\mathcal{H}\) be an RKHS induced by a reproducing kernel \(K\). Given a training set \(\mathcal{D}_n = \{(x_i, y_i)\}_{i=1}^n\). Consider the following generic regularized empirical risk minimization problem involving an arbitrary empirical loss function \(L\) and a monotonically increasing regularization term \(\Omega(\|f\|_{\mathcal{H}})\):

    \[
    \min_{f \in \mathcal{H}} \left\{ \frac{1}{n} \sum_{i=1}^n L(y_i, f(x_i)) + \Omega(\|f\|_{\mathcal{H}}^2) \right\}
    \]

    Then any global optimal solution \(f^*\) of this infinite-dimensional optimization problem can be exactly expressed as a finite linear combination of kernel functions evaluated at the training sample points:

    \[
    f^*(x) = \sum_{i=1}^n \alpha_i K(x_i, x)
    \]

    where \(\alpha_i \in \mathbb{R}\) (\(i=1,\dots,n\)) are a finite set of unknown real scalar coefficients.

??? proof "Proof: Proof of the Representer Theorem via Orthogonal Decomposition"

    We can orthogonally decompose the entire Hilbert space \(\mathcal{H}\) along the subspace spanned by the current observed sample points.

    Define the finite-dimensional linear subspace \(\mathcal{H}_0\) spanned by the kernel mappings of all training sample points as follows:

    \[
    \mathcal{H}_0 = \text{span}\{K(x_1, \cdot), K(x_2, \cdot), \dots, K(x_n, \cdot)\}
    \]

    Then any function \(f \in \mathcal{H}\) can be uniquely decomposed into a component in the subspace \(\mathcal{H}_0\) and a component in the orthogonal complement \(\mathcal{H}_0^{\perp}\):

    \[
    f = f_0 + f_{\perp}
    \]

    where \(f_0 \in \mathcal{H}_0\) and \(f_{\perp} \in \mathcal{H}_0^{\perp}\). This means that for every \(i = 1, \dots, n\),

    \[
    \langle f_{\perp}, K(x_i, \cdot) \rangle_{\mathcal{H}} = 0
    \]

    Next, consider the value of the function \(f\) at any training data point \(x_i\). By the reproducing property of the kernel:

    \[
    f(x_i) = \langle f, K(x_i, \cdot) \rangle_{\mathcal{H}} = \langle f_0 + f_{\perp}, K(x_i, \cdot) \rangle_{\mathcal{H}}
    \]

    \[
    f(x_i) = \langle f_0, K(x_i, \cdot) \rangle_{\mathcal{H}} + \langle f_{\perp}, K(x_i, \cdot) \rangle_{\mathcal{H}} = \langle f_0, K(x_i, \cdot) \rangle_{\mathcal{H}} = f_0(x_i)
    \]

    Thus, the perpendicular component \(f_{\perp}\) has no contribution to the predicted values of the function on the training set; therefore, changing \(f_{\perp}\) does not affect the loss function term \(\frac{1}{n} \sum_{i=1}^n L(y_i, f(x_i))\) at all.

    Finally, consider the regularization term. By the Pythagorean theorem and orthogonality:

    \[
    \|f\|_{\mathcal{H}}^2 = \|f_0 + f_{\perp}\|_{\mathcal{H}}^2 = \|f_0\|_{\mathcal{H}}^2 + \|f_{\perp}\|_{\mathcal{H}}^2 \ge \|f_0\|_{\mathcal{H}}^2
    \]

    Since the regularization penalty function \(\Omega\) is monotonically increasing, we have:

    \[
    \Omega(\|f\|_{\mathcal{H}}^2) \ge \Omega(\|f_0\|_{\mathcal{H}}^2)
    \]

    To minimize the overall objective function, we must force the unnecessary residual component \(f_{\perp} = 0\). Therefore, the optimal solution \(f^*\) must lie in the finite-dimensional subspace \(\mathcal{H}_0\), i.e., it can be expressed as a linear combination of the basis:

    \[
    f^*(\cdot) = f_0(\cdot) = \sum_{i=1}^n \alpha_i K(x_i, \cdot)
    \]

    Substituting the argument \(x\), we obtain:

    \[
    f^*(x) = \sum_{i=1}^n \alpha_i K(x_i, x)
    \]

---

## 4. Kernel Ridge Regression

Kernel Ridge Regression (KRR) is a classic application of the Representer Theorem. It seamlessly extends traditional linear ridge regression to nonlinear spaces using the squared loss.

!!! info "Definition 4.1 (Kernel Ridge Regression Functional Optimization Problem)"

    Given samples, the goal of kernel ridge regression is to find a function \(f\) in the RKHS that minimizes the squared loss with an \(L_2\) regularization penalty:

    \[
    \min_{f \in \mathcal{H}} \left\{ \frac{1}{n} \sum_{i=1}^n (y_i - f(x_i))^2 + \lambda \|f\|_{\mathcal{H}}^2 \right\}
    \]

### 4.1 Finite-Dimensional Reformulation in Matrix Form

According to the Representer Theorem, we substitute \(f(x) = \sum_{j=1}^n \alpha_j K(x_j, x)\). Define the kernel matrix (Gram Matrix) \(K \in \mathbb{R}^{n \times n}\) with entries \(K_{ij} = K(x_i, x_j)\), and the coefficient vector \(\alpha = (\alpha_1, \dots, \alpha_n)^{\top}\).

We rewrite the two parts of the objective function separately:

1. The fitted prediction vector is \(\hat{Y} = K\alpha\). Thus, the empirical squared loss term becomes:

\[
\frac{1}{n} \|Y - K\alpha\|_2^2
\]

2. Using the bilinearity of the inner product, compute the squared function norm:

\[
\|f\|_{\mathcal{H}}^2 = \left\langle \sum_{i=1}^n \alpha_i K(x_i, \cdot), \sum_{j=1}^n \alpha_j K(x_j, \cdot) \right\rangle_{\mathcal{H}} = \sum_{i=1}^n \sum_{j=1}^n \alpha_i \alpha_j \langle K(x_i, \cdot), K(x_j, \cdot) \rangle_{\mathcal{H}} = \alpha^{\top}K\alpha
\]

Combining the two terms, the finite-dimensional parameter objective function for kernel ridge regression is:

\[
S(\alpha) = \frac{1}{n} (Y - K\alpha)^{\top}(Y - K\alpha) + \lambda \alpha^{\top}K\alpha
\]

!!! note "Theorem 4.2 (Explicit Coefficient Solution for Kernel Ridge Regression)"

    If the kernel matrix \(K\) is positive definite, the coefficient vector \(\alpha\) that minimizes the above objective function has the following closed-form solution:

    \[
    \alpha = (K + n\lambda I_n)^{-1}Y
    \]

??? proof "Proof: Derivation of the Kernel Ridge Regression Parameter Solution via Matrix Calculus"

    First, fully expand the matrix objective function \(S(\alpha)\):

    \[
    S(\alpha) = \frac{1}{n} \left( Y^{\top}Y - 2\alpha^{\top}K Y + \alpha^{\top}K K \alpha \right) + \lambda \alpha^{\top}K\alpha
    \]

    Since the Gram matrix \(K\) is symmetric (\(K^{\top} = K\)), we take the gradient with respect to the column vector \(\alpha\) using matrix differentiation formulas:

    \[
    \frac{\partial S(\alpha)}{\partial \alpha} = \frac{1}{n} \left( -2KY + 2KK\alpha \right) + 2\lambda K\alpha
    \]

    Set the gradient vector equal to zero to find the extremum:

    \[
    -\frac{2}{n}KY + \frac{2}{n}KK\alpha + 2\lambda K\alpha = 0
    \]

    Multiply both sides by \(\frac{n}{2}\) to simplify the constants:

    \[
    -KY + KK\alpha + n\lambda K\alpha = 0
    \]

    Rearrange terms to isolate the unknown coefficient vector \(\alpha\) on the right:

    \[
    K(K + n\lambda I_n)\alpha = KY
    \]

    Since for any \(\lambda > 0\), the regularized matrix \((K + n\lambda I_n)\) is always strictly positive definite and invertible, we multiply both sides on the left by its inverse \((K + n\lambda I_n)^{-1}\) (or equivalently cancel the positive semi-definite operator \(K\) in the normal equation) to obtain:

    \[
    \alpha = (K + n\lambda I_n)^{-1}Y
    \]

Finally, for any new test point \(x_0\), the nonlinear prediction output of kernel ridge regression is:

\[
\hat{f}(x_0) = \sum_{i=1}^n \alpha_i K(x_i, x_0) = K(x_0, \cdot)^{\top}(K + n\lambda I_n)^{-1}Y
\]