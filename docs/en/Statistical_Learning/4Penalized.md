# Chapter 4: Advanced Penalized Linear Regression and Selection Consistency

In the previous chapter, we discussed the basic settings of ridge regression and Lasso. Although Lasso achieves variable selection through the \(L_1\) penalty, from a statistical theory perspective, it introduces significant systematic bias when shrinking large coefficients. This chapter starts from asymptotic theory, delving into the "Oracle Properties" of penalized estimators, bias correction methods (Adaptive Lasso, SCAD, MCP), and extensions of the penalty for special data structures.

---

## 1. Theoretical Motivation: High-Dimensional Sparse Setting and Oracle Properties

!!! info "Definition 1.1 (High-Dimensional Sparse Model Setting)"

    Assume the observed data satisfy the multiple linear model \(Y = X\beta + \epsilon\). In the high-dimensional setting where the feature dimension \(p\) is large, even satisfying \(p > n\), we assume that the true regression coefficient vector \(\beta^* = (\beta_1^*, \dots, \beta_p^*)^{\top}\) satisfies a **sparsity** assumption.

    Without loss of generality, let the first \(p_0\) variables be the truly relevant features, and the coefficients corresponding to the remaining \(p - p_0\) variables be zero. Define the index set of true non-zero coefficients as:

    \[
    \mathcal{S} = \{j : \beta_j^* \neq 0\}, \quad |\mathcal{S}| = p_0
    \]

    Correspondingly, the true coefficient vector can be written in block matrix form: \(\beta^* = (\beta_{\mathcal{S}}^{*\top}, \beta_{\mathcal{S}^c}^{*\top})^{\top}\), where \(\beta_{\mathcal{S}^c}^* = 0\).

### 1.1 Oracle Properties

An ideal parameter estimator \(\hat{\beta}\) for a high-dimensional sparse model should perform as if it had an "oracle" (knowing which variables are truly relevant). In statistics, an estimator satisfying the following two conditions is said to possess the **Oracle Properties**:

1. **Selection Consistency (Sparsity/Selection Consistency)**: It can identify the set of truly non-zero variables with probability tending to 1, i.e.:

\[
\lim_{n \to \infty} P(\hat{\mathcal{S}} = \mathcal{S}) = 1
\]

2. **Asymptotic Normality**: For the non-zero coefficient part, the convergence rate of the estimator is exactly the same as that of the OLS estimator when the true sparse structure is known, i.e.:

\[
\sqrt{n}(\hat{\beta}_{\mathcal{S}} - \beta_{\mathcal{S}}^*) \xrightarrow{d} N(0, \Sigma^*)
\]

### 1.2 Limitations of Lasso and the Irreconcilable Conflict

The Lasso estimator uses a uniform penalty strength \(\lambda |\beta_j|\). To eliminate irrelevant variables (making \(\hat{\beta}_{\mathcal{S}^c} = 0\)), \(\lambda\) must be sufficiently large; however, a large penalty strength imposes excessive shrinkage on the truly large coefficients \(\beta_{\mathcal{S}}^*\), leading to significant systematic bias, making it mathematically difficult to simultaneously satisfy both selection consistency and asymptotic normality perfectly.

---

## 2. Bias Correction Method (I): Adaptive Lasso and Non-negative Garrote

To overcome the systematic bias of Lasso, the primary idea is to impose dynamically adjusted, unequal penalty weights on different coefficients.

### 2.1 Non-negative Garrote Estimator

The non-negative garrote is one of the earliest shrinkage and selection methods. It first obtains an initial unbiased estimate \(\hat{\beta}^{\text{init}}\) using OLS (or ridge regression), and then performs a non-negative scaling of the initial estimate by solving the following constrained optimization problem:

\[
\min_{d_1, \dots, d_p} \left\{ \frac{1}{2n} \left\| Y - \sum_{j=1}^p d_j X_j \hat{\beta}_j^{\text{init}} \right\|^2 + \lambda \sum_{j=1}^p d_j \right\}
\]

\[
\text{subject to } d_j \ge 0, \quad \forall j = 1, \dots, p
\]

The final coefficient estimate is \(\hat{\beta}_{j}^{\text{garrote}} = d_j \hat{\beta}_j^{\text{init}}\).

### 2.2 Adaptive Lasso

The adaptive Lasso further generalizes this idea by introducing specific adaptive weights for each feature directly into the \(L_1\) penalty.

!!! info "Definition 2.1 (Adaptive Lasso Optimization Problem)"

    The objective function of the adaptive Lasso is defined as:

    \[
    \min_{\beta} \left\{ \frac{1}{2n} \|Y - X\beta\|^2 + \lambda \sum_{j=1}^p w_j |\beta_j| \right\}
    \]

    The adaptive weight vector \(w = (w_1, \dots, w_p)^{\top}\) is usually constructed using an existing consistent estimator \(\hat{\beta}^{\text{init}}\) (such as the OLS estimator or the ridge regression estimator):

    \[
    w_j = \frac{1}{|\hat{\beta}_j^{\text{init}}|^\gamma}
    \]

    where \(\gamma > 0\) is a constant (usually set to \(\gamma = 1\) or 2).

* *Statistical intuition*: If the true coefficient of a feature \(\beta_j^* \neq 0\), as the sample size increases, its initial estimate \(|\hat{\beta}_j^{\text{init}}|\) will be far from zero, causing the weight \(w_j\) to become small, so the penalty imposed on it is very weak, **eliminating bias**; if the true coefficient \(\beta_j^* = 0\), its initial estimate \(|\hat{\beta}_j^{\text{init}}| \to 0\), causing the weight \(w_j \to \infty\), so it experiences a very strong penalty, thereby **achieving precise variable elimination**. Theoretical results show that the adaptive Lasso perfectly possesses the oracle properties.

---

## 3. Bias Correction Method (II): Non-convex Penalty Operators (SCAD and MCP)

Besides introducing adaptive weights, another more fundamental solution is to directly modify the geometric shape of the penalty function \(P_{\lambda}(\cdot)\) so that its derivative decays to zero for large coefficients.

!!! note "Theorem 3.1 (Three Basic Property Criteria for Oracle Estimators)"

    For a good penalty function \(P_{\lambda}(|\beta|)\) to produce an estimator with oracle properties, its derivative \(P_{\lambda}'(|\beta|)\) should satisfy the following three mathematical calculus criteria simultaneously:

    1. **Sparsity**: \(\min_{t} \{t + P_{\lambda}'(t)\} > 0\), so that small coefficients can be set exactly to zero.
    2. **Unbiasedness**: When \(|t|\) is large, \(P_{\lambda}'(t) = 0\), ensuring that large coefficients do not incur extra bias.
    3. **Continuity**: \(\arg\min_{t} \{ \frac{1}{2}(z - t)^2 + P_{\lambda}(t) \}\) as a function of \(z\) must be continuous, preventing drastic jumps in the predicted values.

### 3.1 SCAD Penalty (Smoothly Clipped Absolute Deviation)

The derivative of the SCAD penalty function is defined by the following piecewise function (where the hyperparameter \(a > 2\), commonly recommended \(a = 3.7\)):

\[
P_{\lambda}'(\theta) = \lambda \mathbb{I}(\theta \le \lambda) + \frac{(a\lambda - \theta)_+}{a - 1} \mathbb{I}(\theta > \lambda)
\]

??? proof "Proof: Derivation of the SCAD Explicit Solution (Threshold Operator) Under Orthogonal Design"

    Under the orthogonal design matrix \(X^{\top}X = I\), the univariate analytical solution of the SCAD estimator for the initial OLS estimate \(z_j = \hat{\beta}_{j}^{\text{OLS}}\) can be derived by integrating the derivative and solving the piecewise extremum problem.

    The piecewise solution yields the following hard-smooth threshold operator:

    \[
    \hat{\beta}_{j}^{\text{SCAD}} = \begin{cases} \text{sign}(z_j)(|z_j| - \lambda)_+, & \text{when } |z_j| \le 2\lambda \\ \frac{(a-1)z_j - \text{sign}(z_j)a\lambda}{a-2}, & \text{when } 2\lambda < |z_j| \le a\lambda \\ z_j, & \text{when } |z_j| > a\lambda \end{cases}
    \]

    From this analytical solution, when the OLS estimate \(|z_j| > a\lambda\), the estimated value is directly equal to \(z_j\) (i.e., the regression itself), **the bias for large coefficients becomes exactly zero**, satisfying the unbiasedness criterion.

### 3.2 MCP Penalty (Minimax Concave Penalty)

Compared to SCAD, MCP has a more rapid change in concavity. Its derivative is defined as (where \(a > 1\)):

\[
P_{\lambda}'(\theta) = \left( \lambda - \frac{\theta}{a} \right)_+
\]

When \(\theta > a\lambda\), the derivative becomes exactly zero. Under the orthogonal design, its univariate analytical solution is:

\[
\hat{\beta}_{j}^{\text{MCP}} = \begin{cases} \text{sign}(z_j) \frac{(|z_j| - \lambda)_+}{1 - 1/a}, & \text{when } |z_j| \le a\lambda \\ z_j, & \text{when } |z_j| > a\lambda \end{cases}
\]

---

## 4. Extensions of the Penalty for Special Data Structures

In many practical scenarios, variables inherently possess some special structure (e.g., group structure or sequential structure). Ordinary Lasso would break these relationships.

### 4.1 Group Lasso

!!! info "Definition 4.1 (Group Lasso Optimization Problem)"

    Suppose the \(p\) predictor variables are divided into \(K\) mutually exclusive groups. Let \(X^{(k)}\) denote the submatrix of features belonging to group \(k\), and \(\beta^{(k)}\) denote the corresponding subvector of coefficients. The group Lasso objective function is defined as:

    \[
    \min_{\beta} \left\{ \frac{1}{2} \left\| Y - \sum_{k=1}^K X^{(k)}\beta^{(k)} \right\|^2 + \lambda \sum_{k=1}^K \sqrt{p_k} \|\beta^{(k)}\|_2 \right\}
    \]

    where \(p_k\) is the number of variables in group \(k\), and \(\|\beta^{(k)}\|_2\) is the Euclidean \(L_2\) norm of the coefficient vector.

* *Statistical implication*: This penalty imposes an \(L_1\) norm between groups (producing between-group sparsity) and an \(L_2\) norm within groups. Consequently, it forces the coefficients of an entire group to be **either all zero (removed as a whole) or all non-zero (retained as a whole)**.

### 4.2 Fused Lasso

When features have a natural ordering (e.g., time series signals, genomic loci arranged along a chromosome), we desire not only sparsity in the coefficients themselves but also smoothness between adjacent coefficients. Its objective function is defined as:

\[
\min_{\beta} \left\{ \frac{1}{2} \sum_{i=1}^n \left( y_i - \sum_{j=1}^p \beta_j x_{ij} \right)^2 + \lambda_1 \sum_{j=1}^p |\beta_j| + \lambda_2 \sum_{j=2}^p |\beta_j - \beta_{j-1}| \right\}
\]

By imposing an \(L_1\) penalty on the differences \(|\beta_j - \beta_{j-1}|\) between adjacent coefficients, the fused Lasso can force adjacent coefficients to produce **exactly equal stepwise block-constant effects**.

---

## 5. Collaborative Regression

In multi-source data integration analysis, such as in biomedicine, we often have features extracted from different sources (e.g., \(X\) representing a DNA methylation matrix, \(Z\) representing an RNA expression matrix). Collaborative regression jointly models them by adding a penalty term that encourages alignment of the predictions from the two sources.

Its basic three-term squared loss optimization objective is:

\[
\min_{\beta, \theta} \left\{ \|Y - X\beta\|_2^2 + \|Y - Z\theta\|_2^2 + \gamma \|X\beta - Z\theta\|_2^2 + \lambda_1 P(\beta) + \lambda_2 P(\theta) \right\}
\]

!!! note "Theorem 5.1 (Design Matrix Augmentation Solution Mechanism for Collaborative Regression)"

    The above collaborative regression optimization problem involving two feature sources, when ignoring the penalty terms, can be perfectly transformed into a standard least squares problem by constructing the following high-dimensional "giant" augmented design matrix and augmented observation vector:

    \[
    \tilde{X} = \begin{pmatrix} X & 0 \\ 0 & Z \\ \sqrt{\gamma}X & -\sqrt{\gamma}Z \end{pmatrix}, \quad \tilde{Y} = \begin{pmatrix} Y \\ Y \\ 0 \end{pmatrix}
    \]

??? proof "Proof: Derivation of the Equivalence of Design Matrix Augmentation"

    Write the matrix expression of the standard residual sum of squares under the augmented system: \(\|\tilde{Y} - \tilde{X} \cdot (\beta^{\top}, \theta^{\top})^{\top} \|_2^2\).

    Using block matrix multiplication, expand it by row blocks:

    \[
    \tilde{Y} - \tilde{X}\begin{pmatrix} \beta \\ \theta \end{pmatrix} = \begin{pmatrix} Y \\ Y \\ 0 \end{pmatrix} - \begin{pmatrix} X\beta + 0 \\ 0 + Z\theta \\ \sqrt{\gamma}X\beta - \sqrt{\gamma}Z\theta \end{pmatrix} = \begin{pmatrix} Y - X\beta \\ Y - Z\theta \\ -\sqrt{\gamma}(X\beta - Z\theta) \end{pmatrix}
    \]

    Taking the square of the Euclidean \(L_2\) norm (i.e., summing the squared Euclidean norms of each block):

    \[
    \left\| \tilde{Y} - \tilde{X}\begin{pmatrix} \beta \\ \theta \end{pmatrix} \right\|_2^2 = \|Y - X\beta\|_2^2 + \|Y - Z\theta\|_2^2 + \left\| -\sqrt{\gamma}(X\beta - Z\theta) \right\|_2^2
    \]

    Removing the constant factor from the norm (squaring it):

    \[
    \left\| \tilde{Y} - \tilde{X}\begin{pmatrix} \beta \\ \theta \end{pmatrix} \right\|_2^2 = \|Y - X\beta\|_2^2 + \|Y - Z\theta\|_2^2 + \gamma \|X\beta - Z\theta\|_2^2
    \]

    This is exactly the main loss function of collaborative regression. Therefore, we can directly use this transformed giant design matrix \(\tilde{X}\) and response vector \(\tilde{Y}\), and apply standard Lasso or ridge regression algorithms to solve the problem.
