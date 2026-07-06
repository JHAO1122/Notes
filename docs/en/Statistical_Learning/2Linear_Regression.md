# Chapter 2: Linear Regression Model and Its Statistical Inference

Linear regression is one of the most fundamental and widely used parametric methods in statistical learning. It assumes a linear relationship between the input variables and the output variable. This chapter starts with the matrix formulation of ordinary least squares estimation, derives its statistical properties, and explores model selection methods in high-dimensional settings.

---

## 1. Basic Setup and Matrix Representation of the Linear Regression Model

!!! info "Definition 1.1 (Linear Regression Model)"

    We observe a training dataset containing \(n\) independent and identically distributed (i.i.d.) samples:

    \[
    \mathcal{D}_{n} = \{(x_{i}, y_{i})\}_{i=1}^{n}
    \]

    where each input vector \(x_{i}\) is a \(p\)-dimensional feature vector:

    \[
    x_{i} = (x_{i1}, x_{i2}, \dots, x_{ip})^{\top}
    \]

    The response variable \(y_{i} \in \mathbb{R}\) is continuous. We assume that the underlying true data-generating mechanism follows the multiple linear model:

    \[
    Y = \beta_0 + \beta_1 X_1 + \beta_2 X_2 + \dots + \beta_p X_p + \epsilon
    \]

    where \(\beta_0\) is the intercept term, \(\beta_j\) (\(j=1,\dots,p\)) are the regression coefficients, and \(\epsilon\) is a random error term satisfying the classical Gauss-Markov assumptions:

    \[
    \mathbb{E}(\epsilon) = 0, \quad \text{Var}(\epsilon) = \sigma^2
    \]

    and the error terms of different samples are independent of each other.

### 1.1 Matrix Representation of the Data

To facilitate mathematical derivations, we usually absorb the intercept term into the feature vector. We let the first column of features be all ones, i.e., \(x_{i} = (1, x_{i1}, \dots, x_{ip})^{\top}\). Define the design matrix \(X \in \mathbb{R}^{n \times (p+1)}\) and the observation vector \(Y \in \mathbb{R}^n\) as follows:

\[
X = \begin{pmatrix} 1 & x_{11} & \dots & x_{1p} \\ 1 & x_{21} & \dots & x_{2p} \\ \vdots & \vdots & \ddots & \vdots \\ 1 & x_{n1} & \dots & x_{np} \end{pmatrix}, \quad Y = \begin{pmatrix} y_1 \\ y_2 \\ \vdots \\ y_n \end{pmatrix}
\]

Let the regression coefficient vector be \(\beta = (\beta_0, \beta_1, \dots, \beta_p)^{\top}\). Then the multiple linear regression model for the entire sample system can be compactly expressed in matrix form:

\[
Y = X\beta + \epsilon
\]

---

## 2. Ordinary Least Squares (OLS) Estimation

To estimate the unknown parameter \(\beta\), we typically employ the method of Ordinary Least Squares (OLS), whose core idea is to minimize the Residual Sum of Squares (RSS).

!!! note "Theorem 2.1 (Explicit Solution of the OLS Estimator)"

    If the design matrix \(X\) has full column rank (i.e., \(\text{Rank}(X) = p+1 < n\)), then the parameter estimator \(\hat{\beta}\) that minimizes the residual sum of squares has the following unique closed-form solution:

    \[
    \hat{\beta} = (X^{\top}X)^{-1}X^{\top}Y
    \]

??? proof "Proof: Derivation of the OLS Explicit Solution Based on Matrix Calculus"

    First, write the matrix expression of the residual sum of squares \(RSS(\beta)\) according to its definition:

    \[
    RSS(\beta) = \sum_{i=1}^n (y_i - x_i^{\top}\beta)^2 = (Y - X\beta)^{\top}(Y - X\beta)
    \]

    Expand it using matrix multiplication:

    \[
    RSS(\beta) = Y^{\top}Y - Y^{\top}X\beta - \beta^{\top}X^{\top}Y + \beta^{\top}X^{\top}X\beta
    \]

    Since \(Y^{\top}X\beta\) is a scalar, its transpose equals itself, i.e., \((Y^{\top}X\beta)^{\top} = \beta^{\top}X^{\top}Y\). Hence the above expression simplifies to:

    \[
    RSS(\beta) = Y^{\top}Y - 2\beta^{\top}X^{\top}Y + \beta^{\top}X^{\top}X\beta
    \]

    To find the \(\beta\) that minimizes \(RSS(\beta)\), we take the partial derivative with respect to the vector \(\beta\). Using the matrix differentiation formulas \(\frac{\partial (\beta^{\top}A)}{\partial \beta} = A\) and \(\frac{\partial (\beta^{\top}B\beta)}{\partial \beta} = 2B\beta\) (where \(B\) is a symmetric matrix), we obtain:

    \[
    \frac{\partial RSS(\beta)}{\partial \beta} = -2X^{\top}Y + 2X^{\top}X\beta
    \]

    Setting the derivative equal to zero gives the well-known normal equations:

    \[
    X^{\top}X\beta = X^{\top}Y
    \]

    Because we assume that \(X\) has full column rank, the matrix \(X^{\top}X\) is invertible. Multiplying both sides of the normal equations on the left by \((X^{\top}X)^{-1}\) yields:

    \[
    \hat{\beta} = (X^{\top}X)^{-1}X^{\top}Y
    \]

---

## 3. Statistical Properties of the OLS Estimator

After obtaining the parameter estimator \(\hat{\beta}\), we need to quantitatively analyze its statistical properties as a random vector (mean and variance-covariance matrix).

!!! note "Theorem 3.1 (Unbiasedness and Variance of the OLS Estimator)"

    Under the Gauss-Markov assumptions, the ordinary least squares estimator \(\hat{\beta}\) is an unbiased estimator of the true parameter \(\beta\), and its covariance matrix is:

    \[
    \text{Var}(\hat{\beta}) = \sigma^2 (X^{\top}X)^{-1}
    \]

??? proof "Proof: Derivation of Unbiasedness and Variance Matrix"

    **1. Proof of Unbiasedness:**

    Substitute the true model \(Y = X\beta + \epsilon\) into the expression for \(\hat{\beta}\):

    \[
    \hat{\beta} = (X^{\top}X)^{-1}X^{\top}(X\beta + \epsilon) = (X^{\top}X)^{-1}X^{\top}X\beta + (X^{\top}X)^{-1}X^{\top}\epsilon
    \]

    Since \((X^{\top}X)^{-1}X^{\top}X = I\) (the identity matrix), we have:

    \[
    \hat{\beta} = \beta + (X^{\top}X)^{-1}X^{\top}\epsilon
    \]

    Conditional on the design matrix \(X\) (which is usually treated as non-random constants), take expectation on both sides:

    \[
    \mathbb{E}(\hat{\beta}) = \mathbb{E}\left[ \beta + (X^{\top}X)^{-1}X^{\top}\epsilon \right] = \beta + (X^{\top}X)^{-1}X^{\top}\mathbb{E}(\epsilon)
    \]

    Because \(\mathbb{E}(\epsilon) = 0\), we obtain:

    \[
    \mathbb{E}(\hat{\beta}) = \beta
    \]

    Unbiasedness is proved.

    **2. Proof of the Covariance Matrix:**

    By the variance property of linear transformations, for a constant matrix \(A\) and a random vector \(Z\), we have \(\text{Var}(AZ) = A\text{Var}(Z)A^{\top}\). Using this, we compute the variance of \(\hat{\beta}\):

    \[
    \text{Var}(\hat{\beta}) = \text{Var}\left( \beta + (X^{\top}X)^{-1}X^{\top}\epsilon \right) = \text{Var}\left( (X^{\top}X)^{-1}X^{\top}\epsilon \right)
    \]

    Treat \((X^{\top}X)^{-1}X^{\top}\) as the constant matrix \(A\) and \(\epsilon\) as the random vector \(Z\). Since the sample errors are independent and identically distributed, their covariance matrix is \(\text{Var}(\epsilon) = \sigma^2 I_n\). Therefore:

    \[
    \text{Var}(\hat{\beta}) = \left[ (X^{\top}X)^{-1}X^{\top} \right] \cdot \text{Var}(\epsilon) \cdot \left[ (X^{\top}X)^{-1}X^{\top} \right]^{\top}
    \]

    \[
    \text{Var}(\hat{\beta}) = (X^{\top}X)^{-1}X^{\top} \cdot (\sigma^2 I_n) \cdot X \left[ (X^{\top}X)^{-1} \right]^{\top}
    \]

    Since \((X^{\top}X)\) is symmetric, its inverse is also symmetric, i.e., \(\left[ (X^{\top}X)^{-1} \right]^{\top} = (X^{\top}X)^{-1}\). Simplifying:

    \[
    \text{Var}(\hat{\beta}) = \sigma^2 (X^{\top}X)^{-1}X^{\top}X(X^{\top}X)^{-1} = \sigma^2 (X^{\top}X)^{-1}
    \]

    The variance matrix is proved.

### 3.1 Unbiased Estimation of the Overall Variance \(\sigma^2\)

Since the true error variance \(\sigma^2\) is unknown, we frequently estimate it using the residual sum of squares. Define the \(n\)-dimensional fitted vector as \(\hat{Y} = X\hat{\beta}\). Then the unbiased estimator \(s^2\) (or \(\hat{\sigma}^2\)) of \(\sigma^2\) is defined as:

\[
s^2 = \frac{RSS(\hat{\beta})}{n - p - 1} = \frac{(Y - \hat{Y})^{\top}(Y - \hat{Y})}{n - p - 1}
\]

The denominator \(n - p - 1\) represents the degrees of freedom, which corrects for the loss of degrees of freedom due to estimating \(p+1\) regression coefficients.

---

## 4. Training Error, Test Error, and Model Selection Criteria

In machine learning and statistics, we must distinguish between the model's performance on the training set and its performance on unseen test data.

* **Training Error**: The average loss computed by the model on the training dataset \(\mathcal{D}_n\). As model complexity increases (i.e., more features are introduced), the training error always decreases monotonically.

* **Test Error**: The expected loss of the model on independently sampled new test points. It exhibits a U-shaped curve that first decreases and then increases.

When the feature dimension \(p\) is large, the full model is prone to serious overfitting. To select the subset of features that yields the best generalization ability, we need to use specific statistical criteria that balance model complexity and goodness of fit.

### 4.1 Mallows' \(C_p\) Criterion

For a submodel containing \(d\) features, the \(C_p\) statistic is defined as:

\[
C_p = \frac{1}{n} \left( RSS_d + 2d\hat{\sigma}^2 \right)
\]

where \(RSS_d\) is the residual sum of squares of the model with \(d\) features, and \(\hat{\sigma}^2\) is the estimate of the overall variance obtained from the full model containing all features.

### 4.2 AIC (Akaike Information Criterion) and BIC (Bayesian Information Criterion)

Under the assumption that the random errors follow a normal distribution, the model selection criteria based on maximum likelihood estimation are expressed as:

* **AIC**:

\[
AIC = \frac{1}{n\hat{\sigma}^2} \left( RSS_d + 2d\hat{\sigma}^2 \right)
\]

* **BIC**: BIC imposes a heavier penalty on complex models compared to AIC (when \(\ln n > 2\)):

\[
BIC = \frac{1}{n\hat{\sigma}^2} \left( RSS_d + \ln(n)d\hat{\sigma}^2 \right)
\]

We should select the feature subset that minimizes \(C_p\), AIC, or BIC.

---

## 5. Model Selection Algorithms

When the total number of features \(p\) is large, exhaustively enumerating all \(2^p\) possible feature combinations is computationally infeasible. Therefore, we need to adopt computationally efficient search algorithms.

### 5.1 Best Subset Selection

* **Basic Steps**: For each fixed model size \(k = 1, 2, \dots, p\), exhaustively examine all \(\binom{p}{k}\) feature combinations and identify the model that minimizes \(RSS\) among those of size \(k\). Then, using criteria such as \(C_p\), AIC, or BIC, select the final global best model from these optimal models across different sizes.

* **Limitations**: Because it requires checking \(2^p\) models, this algorithm is feasible only when \(p\) is small (typically \(p < 50\)).

### 5.2 Stepwise Regression Algorithms

Stepwise regression is a class of greedy algorithms. Although they do not guarantee finding the global optimum, they are extremely fast and perform well in practice.

* **Forward Stepwise Selection**: Start from the null model containing no features. At each step, test the effect of adding each remaining feature individually, and choose the feature that yields the most significant reduction in \(RSS\) to add to the model. Continue until a stopping criterion is met.

* **Backward Elimination**: Start from the full model containing all \(p\) features. At each step, attempt to remove one feature, and select the feature whose removal causes the smallest increase in \(RSS\) (i.e., the least significant feature) to be removed. Continue until all remaining features meet a significance threshold.