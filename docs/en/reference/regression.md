# 📈 Regression Analysis

This module covers the core theory of linear regression models, including Ordinary Least Squares (OLS), the Gauss-Markov theorem, and corrections when model assumptions are violated.

---

## I. Multiple Linear Regression Model and OLS Estimation

!!! info "Matrix Form of the Model"
    The multiple linear regression model can be expressed as:

    \[
    Y = X\beta + \varepsilon
    \]

    where \(Y\) is an \(n \times 1\) vector of observations, \(X\) is an \(n \times k\) design matrix of independent variables (including the constant term), \(\beta\) is a \(k \times 1\) vector of unknown parameters, and \(\varepsilon\) is an \(n \times 1\) vector of random errors.

!!! success "Ordinary Least Squares (OLS)"
    Estimate the parameters by minimizing the Residual Sum of Squares (RSS):

    \[
    \min_{\beta} \text{RSS} = \hat{\varepsilon}^T \hat{\varepsilon} = (Y - X\beta)^T (Y - X\beta)
    \]

    Taking the derivative and setting it to zero yields the **Normal Equations**: \(X^TX\beta = X^TY\). If \(X\) has full column rank, the optimal estimator is:

    \[
    \hat{\beta} = (X^T X)^{-1} X^T Y
    \]

---

## II. Classical Gauss-Markov Assumptions and Properties

!!! info "Classical Linear Regression Model (CLRM) Assumptions"
    * **1. Linearity**: The population model is linear in parameters.
    
    * **2. Strict Exogeneity**: The conditional expectation of the error term is zero, i.e., \(E[\varepsilon \mid X] = 0\).
    
    * **3. No Perfect Multicollinearity**: \(\text{rank}(X) = k < n\), i.e., there is no exact linear relationship among the independent variables.
    
    * **4. Homoskedasticity and No Autocorrelation**: The conditional variance matrix satisfies:

    \[
    \text{Var}(\varepsilon \mid X) = \sigma^2 I_n
    \]

!!! success "Gauss-Markov Theorem"
    Under the assumptions of the classical linear regression model, the OLS estimator \(\hat{\beta}\) has the smallest variance among all linear unbiased estimators; i.e., \(\hat{\beta}\) is **BLUE (Best Linear Unbiased Estimator)**.

??? note "Statistical Properties of the Estimator"
    * **Unbiasedness**: \(E[\hat{\beta}] = \beta\).
    
    * **Covariance matrix**: \(\text{Var}(\hat{\beta}) = \sigma^2 (X^T X)^{-1}\).
    
    * **Unbiased estimator of the error variance \(\sigma^2\)**:

    \[
    s^2 = \hat{\sigma}^2 = \frac{\text{RSS}}{n - k}
    \]

---

## III. Goodness of Fit and Hypothesis Testing

!!! info "Decomposition of Sum of Squares (TSS = ESS + RSS)"
    * **Total Sum of Squares (TSS)**: \(\sum (y_i - \bar{y})^2\).
    
    * **Explained Sum of Squares (ESS)**: \(\sum (\hat{y}_i - \bar{y})^2\).
    
    * **Residual Sum of Squares (RSS)**: \(\sum (y_i - \hat{y}_i)^2\).
    
    * **Coefficient of Determination (\(R^2\))**: A measure of the model's explanatory power:

    \[
    R^2 = \frac{\text{ESS}}{\text{TSS}} = 1 - \frac{\text{RSS}}{\text{TSS}}
    \]

!!! success "General Linear Constraints and the F-test"
    Consider \(m\) linear constraints on the parameters \(\beta\), i.e., the null hypothesis \(H_0: R\beta = r\), where \(R\) is an \(m \times k\) matrix with full row rank.
    
    * **\(\text{RSS}_H\)**: The residual sum of squares under the restricted model that satisfies the null hypothesis \(H_0\) (imposing constraints reduces the model's fit, so \(\text{RSS}_H \ge \text{RSS}\)).
    
    * **\(\text{RSS}\)**: The residual sum of squares under the unrestricted full model.
    
    Under the null hypothesis \(H_0\) and the assumption of normally distributed errors, since \((\text{RSS}_H - \text{RSS})/\sigma^2 \sim \chi^2(m)\) and \(\text{RSS}/\sigma^2 \sim \chi^2(n - \text{rank}(X))\), and they are independent, the following **F-statistic** can be constructed:

    \[
    F = \frac{(\text{RSS}_H - \text{RSS}) / m}{\text{RSS} / (n - \text{rank}(X))} \sim F(m, n - \text{rank}(X))
    \]

    *(Intuition: If the null hypothesis \(H_0\) is true, imposing constraints should not increase the residual sum of squares significantly, so the F-statistic will be small; conversely, if the F-value exceeds the critical value, we reject \(H_0\).)*

??? note "t-test Commonly Used for Single Coefficient Significance"
    * **t-test**: Tests the significance of a single coefficient \(H_0: \beta_j = 0\):

    \[
    t = \frac{\hat{\beta}_j}{\text{se}(\hat{\beta}_j)} \sim t(n - k)
    \]

---

## IV. Violations of Classical Assumptions and Corrections

| Violation | Definition / Effect | Detection Method | Correction Method |
| :--- | :--- | :--- | :--- |
| **Heteroskedasticity** | \(\text{Var}(\varepsilon_i) = \sigma_i^2 \neq \sigma^2\). Leads to inefficient OLS estimators and biased standard errors. | White test, Breusch-Pagan test | Weighted Least Squares (WLS), robust standard errors (White SE) |
| **Autocorrelation** | \(\text{Cov}(\varepsilon_i, \varepsilon_j) \neq 0 \ (i \neq j)\). Common in time series data; leads to inefficient OLS estimators. | Durbin-Watson (DW) test, Breusch-Godfrey test | Generalized Least Squares (GLS), Newey-West robust standard errors |
| **Multicollinearity** | High linear correlation among independent variables. Leads to inflated variance and insignificant individual coefficients. | Variance Inflation Factor (VIF > 10) | Variable removal, Ridge regression, Principal Component Analysis (PCA) |