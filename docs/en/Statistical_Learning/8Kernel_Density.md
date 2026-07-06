# Chapter 8: Kernel Density Estimation and Nonparametric Kernel Regression

In the previous chapters, models often relied on parametric assumptions about the underlying distribution of the data (e.g., multivariate normal distribution) or the functional relationship (e.g., linear model). When these assumptions fail, nonparametric estimation methods can learn the structure directly from the data. This chapter focuses on nonparametric estimation methods based on the theory of local smoothing, systematically discussing the mathematical principles and statistical properties of **Kernel Density Estimation (KDE)** and **Nadaraya-Watson Kernel Regression**.

---

## 1. Kernel Density Estimation

!!! info "Definition 1.1 (Inspiration from the Histogram Estimator)"

    Assume we have independent and identically distributed continuous random variable samples \(X_1, X_2, \dots, X_n \sim f(\cdot)\), where \(f(x)\) is the probability density function to be estimated. According to calculus, the density at a point \(x\) can be expressed as the limit of the probability within a small interval:

    \[
    f(x) = \lim_{\lambda \to 0} \frac{1}{\lambda} P\left( x - \frac{\lambda}{2} \le X \le x + \frac{\lambda}{2} \right)
    \]

    The traditional histogram estimator divides the sample space into fixed independent bins. Instead of using fixed bins, we construct a sliding interval of length \(\lambda\) centered at the target point \(x\), and count the number of samples falling into this interval, yielding a naive density estimator:

    \[
    \hat{f}_{\text{naive}}(x) = \frac{1}{n\lambda} \sum_{i=1}^n \mathbb{I}\left( x - \frac{\lambda}{2} \le X_i \le x + \frac{\lambda}{2} \right) = \frac{1}{n\lambda} \sum_{i=1}^n \mathbb{I}\left( \left| \frac{x - X_i}{\lambda} \right| \le \frac{1}{2} \right)
    \]

### 1.1 Abstraction and Smoothing of Kernel Functions

The weight function corresponding to the naive estimator above is a step function (Uniform Kernel), which is discontinuous at the boundaries and assigns equal weight to all points inside the interval. To obtain smoother density estimates, we introduce a general **kernel function** \(K(t)\).

A valid kernel function must be a legitimate probability density function, satisfying the following basic mathematical properties:

1. **Non-negativity**: \(K(t) \ge 0, \quad \forall t \in \mathbb{R}\)
2. **Integral normalization**: \(\int_{-\infty}^{\infty} K(t) dt = 1\)
3. **Symmetry (usually satisfied)**: \(K(-t) = K(t)\), which also implies its first moment is 0: \(\int_{-\infty}^{\infty} t K(t) dt = 0\)

Commonly used smooth kernel functions include the Gaussian Kernel:

\[
K(t) = \frac{1}{\sqrt{2\pi}} \exp\left( -\frac{t^2}{2} \right)
\]

!!! info "Definition 1.2 (Kernel Density Estimator KDE)"

    Introducing a bandwidth \(h > 0\) (corresponding to \(\lambda\) above), the kernel density estimator is defined as:

    \[
    \hat{f}_h(x) = \frac{1}{nh} \sum_{i=1}^n K\left( \frac{x - X_i}{h} \right)
    \]

---

## 2. Statistical Properties of KDE: Asymptotic Bias and Variance Derivation

To evaluate the quality of the estimator \(\hat{f}_h(x)\), we perform an asymptotic analysis. Assume the true density function \(f(x)\) has a continuous second derivative. Let \(\sigma_K^2 = \int t^2 K(t) dt\) and abbreviate \(R(K) = \int K(t)^2 dt\).

!!! note "Theorem 2.1 (Asymptotic Expectation, Bias, and Variance of KDE)"

    When \(n \to \infty, h \to 0\) and \(nh \to \infty\), the asymptotic bias and asymptotic variance of the kernel density estimator \(\hat{f}_h(x)\) satisfy:

    \[
    \text{Bias}(\hat{f}_h(x)) = \frac{1}{2} h^2 \sigma_K^2 f''(x) + o(h^2)
    \]

    \[
    \text{Var}(\hat{f}_h(x)) = \frac{1}{nh} f(x) R(K) + o\left( \frac{1}{nh} \right)
    \]

??? proof "Proof: Rigorous Integral Derivation of KDE Bias and Variance"

    **1. Derivation of Expectation and Bias:**

    Using the linearity of expectation and the i.i.d. property of the samples:

    \[
    \mathbb{E}[\hat{f}_h(x)] = \mathbb{E}\left[ \frac{1}{nh} \sum_{i=1}^n K\left( \frac{x - X_i}{h} \right) \right] = \frac{1}{h} \mathbb{E}\left[ K\left( \frac{x - X_1}{h} \right) \right]
    \]

    Write it as an integral form for a continuous random variable:

    \[
    \mathbb{E}[\hat{f}_h(x)] = \frac{1}{h} \int_{-\infty}^{\infty} K\left( \frac{x - u}{h} \right) f(u) du
    \]

    Introduce the change of variable \(t = \frac{x - u}{h}\), then \(u = x - ht\), \(du = -h dt\). After transforming the integration limits:

    \[
    \mathbb{E}[\hat{f}_h(x)] = \int_{-\infty}^{\infty} K(t) f(x - ht) dt
    \]

    Perform a second-order Taylor expansion of the true density \(f(x - ht)\) around point \(x\):

    \[
    f(x - ht) = f(x) - ht f'(x) + \frac{1}{2} h^2 t^2 f''(x) + o(h^2)
    \]

    Substitute the expansion back into the integral and use the properties of the kernel function (integral normalization, symmetric first moment zero):

    \[
    \mathbb{E}[\hat{f}_h(x)] = f(x) \int K(t) dt - h f'(x) \int t K(t) dt + \frac{1}{2} h^2 f''(x) \int t^2 K(t) dt + o(h^2)
    \]

    \[
    \mathbb{E}[\hat{f}_h(x)] = f(x) \cdot 1 - 0 + \frac{1}{2} h^2 \sigma_K^2 f''(x) + o(h^2)
    \]

    Therefore, the bias term is:

    \[
    \text{Bias}(\hat{f}_h(x)) = \mathbb{E}[\hat{f}_h(x)] - f(x) = \frac{1}{2} h^2 \sigma_K^2 f''(x) + o(h^2)
    \]

    **2. Derivation of Variance:**

    Since the samples are independent, the variance is additive:

    \[
    \text{Var}(\hat{f}_h(x)) = \text{Var}\left( \frac{1}{nh} \sum_{i=1}^n K\left( \frac{x - X_i}{h} \right) \right) = \frac{1}{n h^2} \text{Var}\left( K\left( \frac{x - X_1}{h} \right) \right)
    \]

    Using the variance formula \(\text{Var}(Z) = \mathbb{E}[Z^2] - (\mathbb{E}[Z])^2\):

    \[
    \text{Var}\left( K\left( \frac{x - X_1}{h} \right) \right) = \int K\left( \frac{x - u}{h} \right)^2 f(u) du - \left( \int K\left( \frac{x - u}{h} \right) f(u) du \right)^2
    \]

    For the first squared integral term, again use the change of variable \(t = \frac{x - u}{h}\):

    \[
    \int K\left( \frac{x - u}{h} \right)^2 f(u) du = h \int K(t)^2 f(x - ht) dt
    \]

    Perform a first-order Taylor expansion on \(f(x - ht)\): \(f(x - ht) = f(x) + o(1)\). Substituting:

    \[
    h \int K(t)^2 [f(x) + o(1)] dt = h f(x) R(K) + o(h)
    \]

    As for the second term \((\mathbb{E}[\dots])^2\), its order is \((h f(x) + o(h))^2 = \mathcal{O}(h^2)\). When \(h \to 0\), it is negligible compared to the first term of order \(\mathcal{O}(h)\). Substituting back:

    \[
    \text{Var}(\hat{f}_h(x)) = \frac{1}{n h^2} \left[ h f(x) R(K) + o(h) \right] = \frac{1}{nh} f(x) R(K) + o\left( \frac{1}{nh} \right)
    \]

    This completes the proof.

### 1.2 Asymptotic Mean Squared Error and Optimal Bandwidth

Combining bias and variance, we can compute the Asymptotic Mean Squared Error (AMSE) at point \(x\):

\[
\text{AMSE}(\hat{f}_h(x)) = \frac{1}{4} h^4 (\sigma_K^2)^2 [f''(x)]^2 + \frac{1}{nh} f(x) R(K)
\]

By taking the derivative \(\frac{d \text{AMSE}}{dh} = 0\), the theoretical optimal bandwidth can be found as:

\[
h_{\text{opt}} = \left( \frac{R(K)}{n (\sigma_K^2)^2 \int [f''(x)]^2 dx} \right)^{1/5} \propto n^{-1/5}
  \]

This indicates that the convergence rate of nonparametric estimation is \(n^{-2/5}\), slower than the traditional parametric rate of \(n^{-1/2}\).

---

## 3. Nonparametric Kernel Regression

Now consider the regression model \(Y = m(X) + \epsilon\), where the goal is to estimate the conditional expectation function (i.e., the regression function) without presupposing its form:

\[
m(x) = \mathbb{E}[Y \mid X = x] = \int y f(y \mid x) dy = \frac{\int y f(x, y) dy}{f_X(x)}
\]

### 3.1 Nadaraya-Watson Estimator

To estimate the numerator and denominator in the expression above, we use a two-dimensional kernel function to perform kernel density estimation for the joint density \(f(x, y)\). Assume the two-dimensional kernel takes a product form \(K_h(x) K_b(y)\).

!!! note "Theorem 3.1 (Form of the Nadaraya-Watson Estimator)"

    By applying KDE to the joint density and marginal density terms respectively and simplifying through integration, we obtain the famous N-W kernel regression estimator:

    \[
    \hat{m}_{\text{NW}}(x) = \sum_{i=1}^n W_i(x) Y_i
    \]

    where the weight function (called the Nadaraya-Watson weight operator) is defined as:

    \[
    W_i(x) = \frac{K\left( \frac{x - X_i}{h} \right)}{\sum_{j=1}^n K\left( \frac{x - X_j}{h} \right)}
    \]

??? proof "Proof: Mathematical Derivation of the Nadaraya-Watson Estimator"

    According to the definition of conditional probability, the numerator term can be written as \(g(x) = \int y f(x, y) dy\). Using multivariate KDE theory, we substitute the joint density estimate from the observed samples \((X_i, Y_i)\):

    \[
    \hat{f}(x, y) = \frac{1}{n h_x h_y} \sum_{i=1}^n K_x\left( \frac{x - X_i}{h_x} \right) K_y\left( \frac{y - Y_i}{h_y} \right)
    \]

    Now estimate the numerator \(g(x)\) by substituting \(\hat{f}(x, y)\) into the integral:

    \[
    \hat{g}(x) = \int y \left[ \frac{1}{n h_x h_y} \sum_{i=1}^n K_x\left( \frac{x - X_i}{h_x} \right) K_y\left( \frac{y - Y_i}{h_y} \right) \right] dy
    \]

    Use the additivity of the integral to swap the summation and integral signs:

    \[
    \hat{g}(x) = \frac{1}{n h_x} \sum_{i=1}^n K_x\left( \frac{x - X_i}{h_x} \right) \left[ \int \frac{y}{h_y} K_y\left( \frac{y - Y_i}{h_y} \right) dy \right]
    \]

    For the inner integral with respect to \(y\), introduce the change of variable \(t = \frac{y - Y_i}{h_y}\), then \(y = Y_i + h_y t\), \(dy = h_y dt\):

    \[
    \int \frac{Y_i + h_y t}{h_y} K_y(t) h_y dt = \int (Y_i + h_y t) K_y(t) dt = Y_i \int K_y(t) dt + h_y \int t K_y(t) dt
    \]

    Since the kernel function \(K_y\) satisfies integral normalization and symmetry (first moment zero), the above simplifies to:

    \[
    Y_i \cdot 1 + 0 = Y_i
    \]

    Substituting this core result back into the numerator estimator \(\hat{g}(x)\):

    \[
    \hat{g}(x) = \frac{1}{n h_x} \sum_{i=1}^n K_x\left( \frac{x - X_i}{h_x} \right) Y_i
    \]

    Meanwhile, the denominator term is the marginal kernel density estimate of \(X\):

    \[
    \hat{f}_X(x) = \frac{1}{n h_x} \sum_{j=1}^n K_x\left( \frac{x - X_j}{h_x} \right)
    \]

    Dividing the two expressions and canceling the common factor \(\frac{1}{n h_x}\), we obtain:

    \[
    \hat{m}_{\text{NW}}(x) = \frac{\hat{g}(x)}{\hat{f}_X(x)} = \frac{\sum_{i=1}^n K_x\left( \frac{x - X_i}{h_x} \right) Y_i}{\sum_{j=1}^n K_x\left( \frac{x - X_j}{h_x} \right)} = \sum_{i=1}^n W_i(x) Y_i
    \]

---

## 4. Local Linear Regression

Although the N-W kernel regression is intuitive, it suffers from severe systematic bias when dealing with **boundary points** or data with **non-uniform design**. This is because the N-W estimator essentially computes a locally weighted constant (zero-order approximation). To remedy this deficiency, we introduce local linear regression.

Local linear regression does not directly compute the local mean; instead, it assumes that the function can be locally approximated by a straight line in the neighborhood of the target point \(x_0\):

\[
m(x) \approx m(x_0) + m'(x_0)(x - x_0) = \alpha + \beta (x - x_0)
\]

Based on this assumption, we solve a locally weighted least squares problem at point \(x_0\):

\[
\min_{\alpha, \beta} \sum_{i=1}^n K\left( \frac{X_i - x_0}{h} \right) \left[ Y_i - \alpha - \beta(X_i - x_0) \right]^2
\]

The solution \(\hat{\alpha}(x_0)\) is the estimate of the regression function value \(m(x_0)\) at that point. Through matrix calculus, it can be shown that local linear regression perfectly cancels the boundary bias effect, providing more robust convergence performance across the entire region.