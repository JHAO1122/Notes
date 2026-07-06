# Chapter 6: Support Vector Machines (SVM)

The Support Vector Machine (SVM) is a classic supervised learning algorithm primarily used for binary classification. Its core idea is to find a hyperplane in the feature space that maximizes the geometric margin between the two classes. This chapter starts from the hard margin maximization under the completely linearly separable case, gradually transitions to the soft margin maximization under the linearly non-separable case, and finally introduces the kernel trick to handle non-linear classification problems, providing a complete derivation of the Lagrange dual problem.

---

## 1. Linearly Separable Support Vector Machines and Hard Margin Maximization

!!! info "Definition 1.1 (Separating Hyperplane and Classification Rule)"

    Given a training dataset \(\mathcal{D}_{n}=\{(x_{i}, y_{i})\}_{i=1}^{n}\), where \(x_{i} \in \mathbb{R}^p\) and the labels \(y_{i} \in \{-1, +1\}\). If the dataset is completely linearly separable, then there exists a hyperplane:

    \[
    f(x) = w^{\top}x + b = 0
    \]

    that can correctly separate the two classes. The corresponding classification decision rule is \(C(x) = \text{sign}(w^{\top}x + b)\). That is, when \(y_i = +1\), we have \(w^{\top}x_i + b > 0\); when \(y_i = -1\), we have \(w^{\top}x_i + b < 0\). These can be uniformly written as the inequality:

    \[
    y_i (w^{\top}x_i + b) > 0, \quad \forall i = 1, \dots, n
    \]

### 1.1 Functional Margin and Geometric Margin

* **Functional Margin**: Defined as \(\hat{\gamma}_i = y_i (w^{\top}x_i + b)\). By proportionally scaling \(w\) and \(b\), the functional margin also scales proportionally, so it cannot objectively reflect the physical distance from a point to the hyperplane.

* **Geometric Margin**: By imposing an \(L_2\) norm constraint on the normal vector \(w\), the true geometric distance from a point to the hyperplane is defined as:

\[
\gamma_i = \frac{y_i (w^{\top}x_i + b)}{\|w\|_2}
\]

The goal of SVM is to maximize the minimum geometric margin over the entire training set \(\gamma = \min_{i} \gamma_i\). Therefore, the original optimization problem can be written as:

\[
\max_{w, b, \gamma} \gamma
\]

\[
\text{subject to } \frac{y_i (w^{\top}x_i + b)}{\|w\|_2} \ge \gamma, \quad \forall i = 1, \dots, n
\]

To eliminate the effect of scaling, set the minimum functional margin \(\hat{\gamma} = 1\). Then the geometric margin becomes \(\gamma = \frac{1}{\|w\|_2}\). Maximizing \(\frac{1}{\|w\|_2}\) is equivalent to minimizing \(\frac{1}{2}\|w\|_2^2\).

!!! note "Theorem 1.1 (Convex Quadratic Programming Problem for Hard Margin SVM)"

    The standard primal optimization problem for the hard margin support vector machine is expressed as follows:

    \[
    \min_{w, b} \frac{1}{2} \|w\|_2^2
    \]

    \[
    \text{subject to } y_i (w^{\top}x_i + b) \ge 1, \quad \forall i = 1, \dots, n
    \]

---

## 2. Derivation of the Dual Problem for Hard Margin Using Lagrange Multipliers

To solve the above constrained optimization problem and to prepare for the introduction of the kernel function, we use Lagrange duality to transform the primal problem into a dual problem.

!!! note "Theorem 2.1 (Hard Margin Dual Optimization Problem)"

    The dual problem of the primal hard margin SVM is a maximization problem in terms of the Lagrange multiplier vector \(\alpha = (\alpha_1, \dots, \alpha_n)^{\top}\):

    \[
    \max_{\alpha \in \mathbb{R}^n} \left\{ \sum_{i=1}^n \alpha_i - \frac{1}{2} \sum_{i=1}^n \sum_{j=1}^n \alpha_i \alpha_j y_i y_j x_i^{\top}x_j \right\}
    \]

    \[
    \text{subject to } \sum_{i=1}^n \alpha_i y_i = 0, \quad \alpha_i \ge 0, \quad \forall i = 1, \dots, n
    \]

??? proof "Proof: Rigorous Derivation of the Hard Margin Lagrange Dual Problem"

    First, introduce non-negative Lagrange multipliers \(\alpha_i \ge 0\) and construct the Lagrangian function as follows:

    \[
    L(w, b, \alpha) = \frac{1}{2} \|w\|_2^2 - \sum_{i=1}^n \alpha_i \left[ y_i (w^{\top}x_i + b) - 1 \right]
    \]

    According to duality theory, we need to first minimize \(L(w, b, \alpha)\) with respect to the primal variables \(w\) and \(b\). To this end, take the partial derivatives and set them equal to zero:

    1. Derivative with respect to \(w\):

    \[
    \frac{\partial L}{\partial w} = w - \sum_{i=1}^n \alpha_i y_i x_i = 0 \implies w = \sum_{i=1}^n \alpha_i y_i x_i
    \]

    2. Derivative with respect to \(b\):

    \[
    \frac{\partial L}{\partial b} = -\sum_{i=1}^n \alpha_i y_i = 0 \implies \sum_{i=1}^n \alpha_i y_i = 0
    \]

    Substitute the expression for \(w\) and the constraint \(\sum_{i=1}^n \alpha_i y_i = 0\) back into the Lagrangian function for elimination:

    \[
    L(w, b, \alpha) = \frac{1}{2} \left( \sum_{i=1}^n \alpha_i y_i x_i \right)^{\top} \left( \sum_{j=1}^n \alpha_j y_j x_j \right) - \sum_{i=1}^n \alpha_i y_i \left( \sum_{j=1}^n \alpha_j y_j x_j \right)^{\top} x_i - b \sum_{i=1}^n \alpha_i y_i + \sum_{i=1}^n \alpha_i
    \]

    Since \(b \sum_{i=1}^n \alpha_i y_i = 0\), the above expression simplifies to:

    \[
    L(w, b, \alpha) = \frac{1}{2} \sum_{i=1}^n \sum_{j=1}^n \alpha_i \alpha_j y_i y_j x_i^{\top}x_j - \sum_{i=1}^n \sum_{j=1}^n \alpha_i \alpha_j y_i y_j x_i^{\top}x_j + \sum_{i=1}^n \alpha_i
    \]

    Combining like terms, we obtain:

    \[
    \min_{w, b} L(w, b, \alpha) = \sum_{i=1}^n \alpha_i - \frac{1}{2} \sum_{i=1}^n \sum_{j=1}^n \alpha_i \alpha_j y_i y_j x_i^{\top}x_j
    \]

    Finally, maximizing this with respect to the multipliers \(\alpha\), along with the constraints on \(\alpha\), yields the dual objective function as stated.

### 2.2 Complementary Slackness and Support Vectors

According to the KKT (Karush-Kuhn-Tucker) conditions, the optimal solution must satisfy the complementary slackness condition:

\[
\alpha_i \left[ y_i (w^{\top}x_i + b) - 1 \right] = 0, \quad \forall i = 1, \dots, n
  \]

* *Definition of Support Vectors*: From the above, if \(\alpha_i > 0\), then it must hold that \(y_i (w^{\top}x_i + b) - 1 = 0\), meaning these sample points lie exactly on the margin boundaries. Such samples are called **support vectors**. The final hyperplane parameter \(w\) is determined solely by these few support vectors with \(\alpha_i > 0\); the majority of the other sample points with \(\alpha_i = 0\) have no contribution to the final decision boundary.

---

## 3. Linearly Non-Separable Support Vector Machines and Soft Margin Maximization

In real-world data, perfect linear separability is rarely satisfied; there is often noise or overlap between the two classes. To improve the model's tolerance, we allow some samples to violate the margin constraint.

!!! info "Definition 3.1 (Slack Variables and Soft Margin)"

    Introduce a slack variable \(\xi_i \ge 0\) for each sample, relaxing the hard constraint to:

    \[
    y_i (w^{\top}x_i + b) \ge 1 - \xi_i, \quad \forall i = 1, \dots, n
    \]

    * When \(0 < \xi_i \le 1\), the sample is correctly classified but lies inside the margin.
    * When \(\xi_i > 1\), the sample crosses the decision boundary and is completely misclassified.

!!! note "Theorem 3.1 (Primal Optimization Problem for Soft Margin SVM)"

    To balance "maximizing the geometric margin" and "minimizing the number of misclassified samples", we introduce a penalty parameter \(C > 0\). The soft margin optimization objective function is expressed as:

    \[
    \min_{w, b, \xi} \left\{ \frac{1}{2} \|w\|_2^2 + C \sum_{i=1}^n \xi_i \right\}
    \]

    \[
    \text{subject to } y_i (w^{\top}x_i + b) \ge 1 - \xi_i, \quad \xi_i \ge 0, \quad \forall i = 1, \dots, n
    \]

### 3.1 Soft Margin Dual Problem and KKT Conditions

Using the same method to construct the Lagrangian function and taking partial derivatives, we can derive that the dual problem for the soft margin is almost identical to that of the hard margin, with the only difference being that the Lagrange multipliers \(\alpha_i\) are bounded above by \(C\) (called the box constraint):

\[
\max_{\alpha \in \mathbb{R}^n} \left\{ \sum_{i=1}^n \alpha_i - \frac{1}{2} \sum_{i=1}^n \sum_{j=1}^n \alpha_i \alpha_j y_i y_j x_i^{\top}x_j \right\}
\]

\[
\text{subject to } \sum_{i=1}^n \alpha_i y_i = 0, \quad 0 \le \alpha_i \le C, \quad \forall i = 1, \dots, n
\]

The corresponding KKT complementary slackness conditions become:

\[
\alpha_i \left[ y_i (w^{\top}x_i + b) - 1 + \xi_i \right] = 0
\]

\[
(C - \alpha_i) \xi_i = 0
\]

* When \(\alpha_i = 0\), then \(\xi_i = 0\), and the sample is correctly classified and lies outside the margin.
* When \(0 < \alpha_i < C\), from the second equation we get \(\xi_i = 0\), meaning the sample lies exactly on the margin boundary.
* When \(\alpha_i = C\), then \(\xi_i > 0\), and the sample lies inside the margin or is misclassified.

---

## 4. Non-linear Support Vector Machines and the Kernel Trick

When the samples exhibit highly non-linear distributions in the original low-dimensional input space, we can use a non-linear mapping \(\phi(x): \mathbb{R}^p \to \mathcal{H}\) to map them into a higher-dimensional Hilbert space \(\mathcal{H}\), where the data become linearly separable.

### 4.1 Kernel Function Definition and Mercer's Theorem

When solving the dual problem in the high-dimensional space, we need to compute the inner product \(\phi(x_i)^{\top}\phi(x_j)\) of the high-dimensional vectors. If we first map and then compute the inner product, the computational complexity would explode. To solve this problem, the model introduces the **kernel function**.

!!! info "Definition 4.1 (Kernel Function)"

    If there exists a mapping \(\phi(x)\) from the input space to a high-dimensional space such that for all \(x, z\),

    \[
    K(x, z) = \langle \phi(x), \phi(z) \rangle
    \]

    holds, then \(K(x, z)\) is called a kernel function.

Commonly used kernel functions include:

* **Polynomial Kernel**:

\[
K(x, z) = (x^{\top}z + 1)^d
\]

* **Gaussian Radial Basis Function (RBF) Kernel**: It can map data into an infinite-dimensional Hilbert space.

\[
K(x, z) = \exp\left( -\gamma \|x - z\|_2^2 \right)
\]

Substituting the kernel function into the dual problem, all inner product terms \(x_i^{\top}x_j\) can be directly replaced by \(K(x_i, x_j)\). The final decision function for the non-linear support vector machine becomes:

\[
f(x) = \sum_{i=1}^n \alpha_i y_i K(x_i, x) + b
\]

---

## 5. Penalty-Loss Perspective: Hinge Loss

In addition to the geometric margin derivation, the support vector machine can also be perfectly integrated into the general machine learning framework of regularized empirical risk minimization (SRM).

!!! note "Theorem 5.1 (Equivalence of SVM to Hinge Loss Formulation)"

    The primal optimization problem of the soft margin support vector machine is mathematically equivalent to the following unconstrained regularization problem using the hinge loss function:

    \[
    \min_{w, b} \left\{ \sum_{i=1}^n \max\left( 0, 1 - y_i(w^{\top}x_i + b) \right) + \frac{1}{2C} \|w\|_2^2 \right\}
    \]

??? proof "Proof: Derivation of the Equivalence to Hinge Loss"

    Let the slack variables in the primal problem satisfy the boundary condition. According to the inequality constraint of the soft margin:

    \[
    y_i (w^{\top}x_i + b) \ge 1 - \xi_i \implies \xi_i \ge 1 - y_i (w^{\top}x_i + b)
    \]

    At the same time, since the primal problem requires \(\xi_i \ge 0\), to minimize the objective function \(\sum_{i=1}^n \xi_i\), \(\xi_i\) should take the smallest possible value that still satisfies the constraint. Therefore, we have:

    \[
    \xi_i = \max\left( 0, 1 - y_i (w^{\top}x_i + b) \right)
    \]

    Substitute this directly into the soft margin objective function:

    \[
    \min_{w, b} \frac{1}{2} \|w\|_2^2 + C \sum_{i=1}^n \max\left( 0, 1 - y_i (w^{\top}x_i + b) \right)
    \]

    Divide the entire objective function by the constant \(C\) and let the new regularization parameter be \(\lambda = \frac{1}{2C}\), we obtain:

    \[
    \min_{w, b} \sum_{i=1}^n \max\left( 0, 1 - y_i(w^{\top}x_i + b) \right) + \lambda \|w\|_2^2
    \]

    This is exactly the traditional "loss term + \(L_2\) regularization term" framework, thus proving the equivalence.
