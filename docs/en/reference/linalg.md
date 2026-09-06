# 🔢 Linear Algebra

This module currently focuses on content related to **Linear Algebra**.

## 1. Matrix Calculus

Matrix calculus has always been the most challenging topic for me. I encountered it repeatedly in courses like *Regression Analysis*, *Multivariate Statistical Analysis*, and *Practical Algorithms for Optimization*, yet I never organized it systematically. After losing 15 points on a major matrix derivation problem during the *Practical Algorithms for Optimization* final exam due to time constraints and cumbersome summation notation, I decided to consolidate these methods, specifically focusing on the more elegant **Matrix Differential Method**.

!!! warning "Layout Convention"
    This reference sheet consistently adopts the **Denominator Layout**, which is most common in machine learning and optimization.
    That is: for a scalar $f$ with respect to an $m \times n$ matrix $X$, the derivative $\frac{\partial f}{\partial X}$ has the same dimensions as $X$ ($m \times n$).

---

### (1) Core Method: Matrix Differential Method

The traditional approach of taking partial derivatives element-wise is prone to errors within complex summation symbols. A more elegant alternative is to use the **Total Differential**.

!!! info "Connection Between Differential and Derivative"
    The total differential $df$ of a scalar function $f(X)$ can always be expressed as the trace of the product of a certain matrix and $dX$:
    
    $$
    df = \text{tr}\left( \left(\frac{\partial f}{\partial X}\right)^T dX \right)
    $$
    
    **The Four-Step Method:**

    1. Compute the total differential $df$.
    2. Simplify using trace properties to reach the form $\text{tr}(M^T dX)$.
    3. Remove the $\text{tr}$ and $dX$; the remaining matrix $M$ is the derivative $\frac{\partial f}{\partial X}$.

??? note "Basic Rules of Matrix Differentials"
    * **Addition/Subtraction**: $d(X \pm Y) = dX \pm dY$
    * **Product Rule**: $d(XY) = (dX)Y + X(dY)$
    * **Transpose**: $d(X^T) = (dX)^T$
    * **Trace**: $d(\text{tr}(X)) = \text{tr}(dX)$
    * **Cyclic Property of Trace**: $\text{tr}(ABC) = \text{tr}(CAB) = \text{tr}(BCA)$
    * **Invariance under Transpose**: $\text{tr}(A) = \text{tr}(A^T)$

---

### (2) Trace Derivatives

Trace derivatives are frequently required when deriving OLS (Ordinary Least Squares) or constructing loss functions.

!!! success "Common Trace Derivative Formulas"
    Let $X$ be the variable matrix, and $A, B$ be constant matrices:
    
    1. $$\frac{\partial \text{tr}(AX)}{\partial X} = A^T$$

    2. $$\frac{\partial \text{tr}(AXB)}{\partial X} = A^T B^T$$

    3. $$\frac{\partial \text{tr}(X^TAX)}{\partial X} = (A + A^T)X$$
       
    *(Special Case: If $A$ is symmetric, the result simplifies to $2AX$)*

??? note "Derivation: $\frac{\partial \text{tr}(AX)}{\partial X}$"
    Let $f = \text{tr}(AX)$. Compute the differential:
    
    $$
    df = d(\text{tr}(AX)) = \text{tr}(d(AX)) = \text{tr}(A dX)
    $$
    
    To match the form $\text{tr}(M^T dX)$, observe that $\text{tr}(A dX) = \text{tr}((A^T)^T dX)$. Thus, $M^T = A \implies M = A^T$.
    
    Therefore: $\frac{\partial f}{\partial X} = A^T$.

??? note "Derivation: $\frac{\partial \text{tr}(X^TAX)}{\partial X}$"
    Let $f = \text{tr}(X^TAX)$. Compute the differential (using the product rule $d(UV) = (dU)V + U(dV)$):
    
    $$
    df = \text{tr}(d(X^T) AX + X^T d(AX)) = \text{tr}((dX)^T AX) + \text{tr}(X^T A dX)
    $$
    
    Using the transpose invariance $\text{tr}(M) = \text{tr}(M^T)$ for the first term:
    
    $$
    \text{tr}((dX)^T AX) = \text{tr}(((dX)^T AX)^T) = \text{tr}(X^T A^T dX)
    $$
    
    Combining the terms:
    
    $$
    df = \text{tr}(X^T A^T dX) + \text{tr}(X^T A dX) = \text{tr}((X^T A^T + X^T A) dX)
    $$
    
    Factorizing and using the cyclic property to match $\text{tr}(M^T dX)$:
    
    $$
    df = \text{tr}( (X^T(A^T + A)) dX ) = \text{tr}( ((A^T + A)^T X)^T dX )
    $$
    
    Thus, $M = (A^T + A)^T X = (A + A^T)X$.

---

### (3) Determinant and Inverse Matrix Derivatives

In Maximum Likelihood Estimation (MLE), especially involving the covariance matrix $\Sigma$ of a multivariate normal distribution, derivatives of determinants and inverses are essential.

!!! success "Inverse and Determinant Derivative Formulas"
    Let $X$ be an invertible square matrix:
    
    1. **Inverse Matrix Differential**: $d(X^{-1}) = -X^{-1} (dX) X^{-1}$

    2. **Determinant Differential (Jacobi's Formula)**: $d|X| = |X| \text{tr}(X^{-1} dX)$
    
    3. **Derivative of the Determinant**: $\frac{\partial |X|}{\partial X} = |X| (X^{-1})^T$
    
    4. **Derivative of the Log-Determinant**: $\frac{\partial \ln |X|}{\partial X} = (X^{-1})^T$

??? note "Derivation: Differential of the Inverse $d(X^{-1})$"
    Starting from the identity $XX^{-1} = I$, take the differential of both sides:
    
    $$
    d(XX^{-1}) = d(I) \implies (dX)X^{-1} + X d(X^{-1}) = 0
    $$
    
    Rearranging and left-multiplying by $X^{-1}$:
    
    $$
    X d(X^{-1}) = -(dX)X^{-1} \implies d(X^{-1}) = -X^{-1}(dX)X^{-1}
    $$

??? note "Determinant Derivative $d|X|$ (Cofactor Summation Method)"
    This derivation bridges the gap between scalar summation, the adjugate matrix, and the matrix trace.

    **1. Element-wise derivative: $\frac{\partial |X|}{\partial x_{ij}} = C_{ij}$**
    According to the **Laplace Expansion**, expanding along the $i$-th row:

    $$|X| = \sum_{k=1}^n x_{ik} C_{ik}$$
    
    Since the calculation of the cofactor $C_{ik}$ does not include any elements from the $i$-th row, $C_{ik}$ is treated as a constant with respect to $x_{ij}$. Taking the derivative:

    $$\frac{\partial |X|}{\partial x_{ij}} = \frac{\partial}{\partial x_{ij}} (x_{i1}C_{i1} + \dots + x_{ij}C_{ij} + \dots + x_{in}C_{in}) = C_{ij}$$
    
    **Conclusion**: The partial derivative of the determinant with respect to an element is its corresponding cofactor.

    **2. Converting Differential to Trace: $d|X| = \text{tr}(X^* dX)$**
    Using the definition of the total differential and the adjugate matrix $X^*$ (where $(X^*)_{ji} = C_{ij}$):
    
    $$
    d|X| = \sum_{i=1}^n \sum_{j=1}^n \frac{\partial |X|}{\partial x_{ij}} dx_{ij} = \sum_{i=1}^n \sum_{j=1}^n C_{ij} dx_{ij}
    $$
    
    Rewriting this as the sum of the diagonal elements of a matrix product:
    
    $$
    d|X| = \sum_{j=1}^n \left( \sum_{i=1}^n (X^*)_{ji} (dX)_{ij} \right) = \text{tr}(X^* dX)
    $$

    **3. Deriving Jacobi's Formula: $d|X| = |X| \text{tr}(X^{-1} dX)$**
    Substituting the relationship $X^* = |X| X^{-1}$:
    
    $$d|X| = \text{tr}(|X| X^{-1} dX) = |X| \text{tr}(X^{-1} dX)$$

    This is **Jacobi's Formula**, relating the rate of change of the determinant to the inverse and trace of the matrix.

    **4. Gradient of the Log-Determinant $\frac{\partial \ln |X|}{\partial X}$**
    Compute the total differential for $f = \ln |X|$ and substitute Jacobi's Formula:
    
    $$df = d(\ln |X|) = \frac{1}{|X|} d|X| = \text{tr}(X^{-1} dX)$$

    Comparing this to the standard form $df = \text{tr}\left( \left(\frac{\partial f}{\partial X}\right)^T dX \right)$: $\left( \frac{\partial \ln |X|}{\partial X} \right)^T = X^{-1} \implies \frac{\partial \ln |X|}{\partial X} = (X^{-1})^T$.

    *(Note: If $X$ is symmetric, the result simplifies to $X^{-1}$)*

??? note "Determinant Derivative $d|X|$ (Matrix Differential Method)"
    This derivation is entirely independent of element-wise Laplace expansion, utilizing only the multiplicative properties of determinants and the linearity of the trace.

    **Core Lemma (Infinitesimal around Identity):**
    For an infinitesimal perturbation $dE$ to the identity matrix $I$:
    
    $$|I + dE| = 1 + \text{tr}(dE) + o(\|dE\|)$$
    
    **Intuition**: If $dE$ has eigenvalues $\lambda_i$, then $|I+dE| = \prod (1+\lambda_i) = 1 + \sum \lambda_i + \dots = 1 + \text{tr}(dE)$.

    **General Case Derivation for $X$:**
    To calculate $d|X| = |X + dX| - |X|$, we factor out $X$ to construct the identity form:
    
    $$
    \begin{aligned}
    |X + dX| &= |X(I + X^{-1}dX)| \\
    &= |X| \cdot |I + X^{-1}dX|
    \end{aligned}
    $$
    
    Substituting the lemma with $dE = X^{-1}dX$:
    
    $$
    \begin{aligned}
    |X + dX| &= |X| \cdot (1 + \text{tr}(X^{-1}dX)) \\
    &= |X| + |X|\text{tr}(X^{-1}dX)
    \end{aligned}
    $$
    
    By the definition of the differential $d|X| = |X + dX| - |X|$, we immediately obtain:
    
    $$d|X| = |X|\text{tr}(X^{-1}dX)$$

??? note "Derivation: Gradient of the Determinant $\frac{\partial |X|}{\partial X}$"
    Let $f = |X|$. According to Jacobi's Formula, the differential is:
    
    $$
    df = |X| \text{tr}(X^{-1} dX)
    $$
    
    Rewriting this into the inner product form of the trace:
    
    $$
    df = \text{tr}(|X| X^{-1} dX) = \text{tr}( (|X| (X^{-1})^T)^T dX )
    $$
    
    Comparing to $df = \text{tr}(M^T dX)$, we obtain:
    
    $$
    \frac{\partial |X|}{\partial X} = |X| (X^{-1})^T
    $$

---

### (4) Eigenvalue Derivatives

Highly important for deriving Principal Component Analysis (PCA) or studying perturbations of covariance matrices.

!!! success "Eigenvalue Derivative Formula"
    Let $X$ be a square matrix with a non-zero eigenvalue $\lambda$, its corresponding **right eigenvector** $u$ ($Xu = \lambda u$), and its **left eigenvector** $v$ ($v^T X = \lambda v^T$).
    
    Assuming the eigenvectors are normalized such that $v^T u = 1$, the derivative of the eigenvalue $\lambda$ with respect to $X$ is:
    
    $$
    \frac{\partial \lambda}{\partial X} = v u^T
    $$
    
    *(Note: If $X$ is a real symmetric matrix, the left and right eigenvectors are identical ($u=v$), and the derivative is $uu^T$)*

??? note "Derivation: $\frac{\partial \lambda}{\partial X}$"
    Start from the definition of the right eigenvector:
    
    $$
    Xu = \lambda u
    $$
    
    Take the total differential of both sides:
    
    $$
    (dX)u + X(du) = (d\lambda)u + \lambda(du)
    $$
    
    Left-multiply both sides by the left eigenvector $v^T$:
    
    $$
    v^T(dX)u + v^TX(du) = v^T(d\lambda)u + v^T\lambda(du)
    $$
    
    Since $v^T$ is the left eigenvector ($v^TX = \lambda v^T$), substitute this into the second term on the left:
    
    $$
    v^T(dX)u + \lambda v^T(du) = d\lambda (v^Tu) + \lambda v^T(du)
    $$
    
    The terms $\lambda v^T(du)$ on both sides cancel out:
    
    $$
    v^T(dX)u = d\lambda (v^Tu) \implies d\lambda = \frac{v^T(dX)u}{v^Tu}
    $$
    
    Since the trace of a scalar is itself, we can apply $\text{tr}$ and the cyclic property:
    
    $$
    d\lambda = \text{tr}\left( \frac{v^T(dX)u}{v^Tu} \right) = \text{tr}\left( \frac{u v^T}{v^Tu} dX \right)
    $$
    
    Extracting the scalar $v^T u$ and comparing to the standard form $df = \text{tr}(M^T dX)$:
    
    $$
    M^T = \frac{u v^T}{v^Tu} \implies M = \frac{v u^T}{v^Tu}
    $$
    
    Given the normalization $v^Tu = 1$, the final result is:
    
    $$
    \frac{\partial \lambda}{\partial X} = v u^T
    $$