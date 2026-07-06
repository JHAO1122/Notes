# Chapter 3: Penalized Linear Regression and Shrinkage Estimation

When the feature dimension \(p\) is large, or even larger than the sample size \(n\) (i.e., the high-dimensional setting \(p > n\)), the traditional Ordinary Least Squares (OLS) estimator suffers from severe overfitting, and may even have non-unique solutions due to the non-invertibility of the design matrix \(X^{\top}X\). This chapter introduces penalized linear regression, which adds a penalty term (regularization) on the coefficients to the loss function, thereby improving the prediction accuracy and interpretability of the model through the bias-variance tradeoff.

---

## 1. Ridge Regression

!!! info "Definition 1.1 (Ridge Regression Optimization Problem)"

    Ridge regression adds an \(L_2\) norm penalty on the regression coefficients to the Residual Sum of Squares (RSS). Its objective function is defined as:

    \[
    \min_{\beta \in \mathbb{R}^p} \left\{ \sum_{i=1}^n \left( y_i - \beta_0 - \sum_{j=1}^p \beta_j x_{ij} \right)^2 + \lambda \sum_{j=1}^p \beta_j^2 \right\}
    \]

    where \(\lambda \ge 0\) is a tuning parameter that controls the complexity.

    * When \(\lambda = 0\), ridge regression reduces to the traditional OLS estimator.
    * When \(\lambda \to \infty\), the penalty term dominates, forcing all regression coefficients \(\beta_j \to 0\) (excluding the intercept term).

### 1.1 Centering of the Data and Matrix Representation

Before solving, it is common to standardize the features and center both the response variable and the features, so that the centered data satisfy \(\sum_{i=1}^n y_i = 0\) and \(\sum_{i=1}^n x_{ij} = 0\). In this case, the intercept estimate \(\hat{\beta}_0 = 0\), and we can directly write the matrix objective function for the parameter vector \(\beta = (\beta_1, \dots, \beta_p)^{\top}\):

\[
S_{\lambda}(\beta) = (Y - X\beta)^{\top}(Y - X\beta) + \lambda \beta^{\top}\beta
\]

!!! note "Theorem 1.1 (Explicit Solution of the Ridge Regression Estimator)"

    For any given \(\lambda > 0\), the ridge regression parameter estimator \(\hat{\beta}_{\text{ridge}}\) has the following unique closed-form solution:

    \[
    \hat{\beta}_{\text{ridge}} = (X^{\top}X + \lambda I_p)^{-1}X^{\top}Y
    \]

    where \(I_p\) is the \(p \times p\) identity matrix.

??? proof "Proof: Derivation of the Ridge Regression Explicit Solution"

    First, expand the matrix objective function \(S_{\lambda}(\beta)\):

    \[
    S_{\lambda}(\beta) = Y^{\top}Y - 2\beta^{\top}X^{\top}Y + \beta^{\top}X^{\top}X\beta + \lambda \beta^{\top}\beta
    \]

    Take the partial derivative with respect to the parameter vector \(\beta\) using matrix calculus:

    \[
    \frac{\partial S_{\lambda}(\beta)}{\partial \beta} = -2X^{\top}Y + 2X^{\top}X\beta + 2\lambda \beta
    \]

    Set the derivative vector equal to zero to obtain the modified normal equations:

    \[
    (X^{\top}X + \lambda I_p)\beta = X^{\top}Y
    \]

    Since when \(\lambda > 0\), even if the base matrix \(X^{\top}X\) is singular (e.g., when \(p > n\)), the matrix \((X^{\top}X + \lambda I_p)\) is strictly positive definite and therefore always invertible. Multiplying both sides on the left by the inverse matrix yields the result:

    \[
    \hat{\beta}_{\text{ridge}} = (X^{\top}X + \lambda I_p)^{-1}X^{\top}Y
    \]

### 1.2 Derivation of the Expectation and Variance of Ridge Regression

Now assume the true model is \(Y = X\beta + \epsilon\), where \(\mathbb{E}(\epsilon) = 0\) and \(\text{Var}(\epsilon) = \sigma^2 I_n\).

!!! note "Theorem 1.2 (Statistical Properties of Ridge Regression)"

    The ridge regression estimator \(\hat{\beta}_{\text{ridge}}\) is a biased estimator. Its expectation vector and covariance matrix are respectively:

    \[
    \mathbb{E}(\hat{\beta}_{\text{ridge}}) = (X^{\top}X + \lambda I_p)^{-1}X^{\top}X\beta
    \]

    \[
    \text{Var}(\hat{\beta}_{\text{ridge}}) = \sigma^2 (X^{\top}X + \lambda I_p)^{-1}X^{\top}X(X^{\top}X + \lambda I_p)^{-1}
    \]

??? proof "Proof: Mathematical Derivation of the Expectation and Variance of Ridge Regression"

    **1. Derivation of the Expectation:**

    Substitute \(Y = X\beta + \epsilon\) into the explicit solution:

    \[
    \hat{\beta}_{\text{ridge}} = (X^{\top}X + \lambda I_p)^{-1}X^{\top}(X\beta + \epsilon) = (X^{\top}X + \lambda I_p)^{-1}X^{\top}X\beta + (X^{\top}X + \lambda I_p)^{-1}X^{\top}\epsilon
    \]

    Taking conditional expectation, since \(\mathbb{E}(\epsilon) = 0\), we obtain:

    \[
    \mathbb{E}(\hat{\beta}_{\text{ridge}}) = (X^{\top}X + \lambda I_p)^{-1}X^{\top}X\beta
    \]

    Clearly, unless \(\lambda = 0\), we have \(\mathbb{E}(\hat{\beta}_{\text{ridge}}) \neq \beta\), meaning the estimator is biased.

    **2. Derivation of the Variance Matrix:**

    Using the variance property of linear transformations, \(\text{Var}(AZ) = A\text{Var}(Z)A^{\top}\), and noting that \(\beta\) is a non-random constant:

    \[
    \text{Var}(\hat{\beta}_{\text{ridge}}) = \text{Var}\left( (X^{\top}X + \lambda I_p)^{-1}X^{\top}\epsilon \right)
    \]

    Extract the constant matrix and expand using \(\text{Var}(\epsilon) = \sigma^2 I_n\) to obtain:

    \[
    \text{Var}(\hat{\beta}_{\text{ridge}}) = (X^{\top}X + \lambda I_p)^{-1}X^{\top} \left( \sigma^2 I_n \right) \left[ (X^{\top}X + \lambda I_p)^{-1}X^{\top} \right]^{\top}
    \]

    \[
    \text{Var}(\hat{\beta}_{\text{ridge}}) = \sigma^2 (X^{\top}X + \lambda I_p)^{-1}X^{\top}X(X^{\top}X + \lambda I_p)^{-1}
    \]

---

## 2. Lasso Regression (Least Absolute Shrinkage and Selection Operator)

!!! info "Definition 2.1 (Lasso Optimization Problem)"

    Lasso regression achieves coefficient shrinkage by adding an \(L_1\) norm penalty on the regression coefficients. The objective function under a centered design matrix is defined as:

    \[
    \min_{\beta \in \mathbb{R}^p} \left\{ \sum_{i=1}^n \left( y_i - \sum_{j=1}^p \beta_j x_{ij} \right)^2 + \lambda \sum_{j=1}^p |\beta_j| \right\}
    \]

### 2.1 Geometric Interpretation of Sparsity

Because the absolute value function \(|\beta_j|\) in the \(L_1\) penalty is non-differentiable at zero, its corresponding constraint region is a polytope (a diamond in two dimensions), while the \(L_2\) constraint region for ridge regression is a hypersphere. When the elliptical contours of the RSS first intersect the constraint region, the intersection point for Lasso is very likely to lie on a vertex or an edge of the polytope (i.e., on some axes), thereby forcing some of the regression coefficients to be **exactly equal to zero**. This endows Lasso with the dual functionality of parameter shrinkage and feature selection.

---

## 3. Analytical Solutions Under an Orthogonal Design Matrix and Special Case Analysis

To gain deeper insight into the shrinkage effects of different penalty terms on the coefficients, consider a special case where the design matrix is strictly orthogonal, i.e., \(X^{\top}X = I_p\). In this case, the traditional OLS estimator is \(\hat{\beta}_{\text{OLS}} = X^{\top}Y\).

!!! note "Theorem 3.1 (Special Estimator Solutions Under Orthogonal Design)"

    Under the orthogonal design \(X^{\top}X = I_p\), the solutions for ridge regression and Lasso can be expressed as explicit univariate functions of \(\hat{\beta}_{\text{OLS}}\):

    * **Ridge Regression Solution**:

    \[
    \hat{\beta}_{j,\text{ridge}} = \frac{1}{1 + \lambda} \hat{\beta}_{j,\text{OLS}}
    \]

    * **Lasso Solution (Soft Thresholding)**:

    \[
    \hat{\beta}_{j,\text{Lasso}} = \text{sign}(\hat{\beta}_{j,\text{OLS}}) \max \left( |\hat{\beta}_{j,\text{OLS}}| - \frac{\lambda}{2}, 0 \right)
    \]

??? proof "Proof: Derivation of Explicit Shrinkage Solutions Under Orthogonal Design"

    **1. Ridge Regression Orthogonal Solution:**

    Substitute \(X^{\top}X = I_p\) into the ridge regression analytical solution:

    \[
    \hat{\beta}_{\text{ridge}} = (I_p + \lambda I_p)^{-1}X^{\top}Y = \frac{1}{1 + \lambda} X^{\top}Y
    \]

    Since \(\hat{\beta}_{\text{OLS}} = X^{\top}Y\), we have for each component:

    \[
    \hat{\beta}_{j,\text{ridge}} = \frac{1}{1 + \lambda} \hat{\beta}_{j,\text{OLS}}
    \]

    **2. Lasso Orthogonal Solution:**

    Under the orthogonal design, the Lasso objective function decomposes into \(p\) independent univariate optimization problems. For the \(j\)-th component we have:

    \[
    Q(\beta_j) = (\hat{\beta}_{j,\text{OLS}} - \beta_j)^2 + \lambda |\beta_j|
    \]

    We solve \(\min Q(\beta_j)\) by discussing the sign of \(\beta_j\):

    * When \(\beta_j > 0\), \(|\beta_j| = \beta_j\). Taking the derivative:

    \[
    \frac{dQ}{d\beta_j} = -2(\hat{\beta}_{j,\text{OLS}} - \beta_j) + \lambda = 0 \implies \hat{\beta}_j = \hat{\beta}_{j,\text{OLS}} - \frac{\lambda}{2}
    \]

    A necessary condition for this solution to be valid is \(\hat{\beta}_{j,\text{OLS}} > \frac{\lambda}{2}\).

    * When \(\beta_j < 0\), \(|\beta_j| = -\beta_j\). Taking the derivative:

    \[
    \frac{dQ}{d\beta_j} = -2(\hat{\beta}_{j,\text{OLS}} - \beta_j) - \lambda = 0 \implies \hat{\beta}_j = \hat{\beta}_{j,\text{OLS}} + \frac{\lambda}{2}
    \]

    A necessary condition for this solution to be valid is \(\hat{\beta}_{j,\text{OLS}} < -\frac{\lambda}{2}\).

    * When \(-\frac{\lambda}{2} \le \hat{\beta}_{j,\text{OLS}} \le \frac{\lambda}{2}\), the minimum is attained at the boundary point \(\beta_j = 0\).

    Combining the three cases and using the sign function \(\text{sign}(\cdot)\) and the truncation function \(\max(\cdot, 0)\), we obtain:

    \[
    \hat{\beta}_{j,\text{Lasso}} = \text{sign}(\hat{\beta}_{j,\text{OLS}}) \max \left( |\hat{\beta}_{j,\text{OLS}}| - \frac{\lambda}{2}, 0 \right)
    \]

---

## 4. Extensions of the Penalty: \(L_q\) Norms and the Elastic Net

To combine the advantages of ridge regression and Lasso, the form of the penalty can be generalized.

### 4.1 Generalized \(L_q\) Penalty (Bridge Penalty)

We consider a general form of the regression coefficient penalty:

\[
\lambda \sum_{j=1}^p |\beta_j|^q, \quad q \ge 0
\]

* When \(q = 2\), this corresponds to ridge regression (\(L_2\) norm).
* When \(q = 1\), this corresponds to Lasso regression (\(L_1\) norm).
* When \(q = 0\), the penalty becomes the number of non-zero coefficients \(\sum \mathbb{I}(\beta_j \neq 0)\), corresponding to best subset selection.
* *Geometric property*: Sparsity (feature selection) is possible if and only if \(q \le 1\), because the constraint region has a sharp corner at zero. When \(q < 1\), the constraint region is non-convex, making the optimization problem much harder to solve.

### 4.2 Elastic Net

The elastic net linearly combines the \(L_1\) and \(L_2\) penalties to overcome a limitation of Lasso when dealing with groups of highly correlated features (Lasso tends to randomly select one from a group of highly correlated features, whereas the elastic net tends to keep them together as a group). Its objective function is defined as:

\[
\min_{\beta} \left\{ \sum_{i=1}^n \left( y_i - \sum_{j=1}^p \beta_j x_{ij} \right)^2 + \lambda_1 \sum_{j=1}^p |\beta_j| + \lambda_2 \sum_{j=1}^p \beta_j ^2 \right\}
\]