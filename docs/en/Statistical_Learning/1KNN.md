# Chapter 1: K-Nearest Neighbors and the Bias-Variance Tradeoff

In statistical learning, the K-Nearest Neighbor (KNN) algorithm is a classic and intuitive nonparametric method. Its core idea is to predict a target point by using the labels of the \(K\) known samples that are closest to it in the feature space. This chapter derives in detail the mathematical formulation of KNN for regression and classification tasks, and delves into the crucial bias-variance tradeoff and concepts such as degrees of freedom in machine learning.

---

## 1. KNN Regression Model and Its Mathematical Definition

!!! info "Definition 1.1 (Basic Regression Model Setup)"

    Let the input variable be \(X \in \mathbb{R}^p\) and the output variable be \(Y \in \mathbb{R}\). We assume that the real data-generating process follows the regression model:

    \[
    Y = f(X) + \epsilon
    \]

    where \(f(X)\) is the true regression function to be estimated, and \(\epsilon\) is a random error term satisfying:

    \[
    \mathbb{E}(\epsilon) = 0, \quad \text{Var}(\epsilon) = \sigma^2
    \]

    and \(\epsilon\) is independent of \(X\). We now have a training dataset of \(n\) independent and identically distributed (i.i.d.) samples:

    \[
    \mathcal{D}_n = \{(x_1, y_1), (x_2, y_2), \dots, (x_n, y_n)\}
    \]

### 1.1 K-Nearest Neighbor Regression Estimator

For a given test target point \(x_0\), the KNN algorithm estimates the regression function value at that point by finding the \(k\) samples in the training set \(\mathcal{D}_n\) that are geometrically closest to it.

Let \(N_k(x_0)\) be the set containing the feature values of these \(k\) nearest neighbors. Then the KNN prediction \(\hat{y}\) (or \(\hat{f}(x_0)\)) at the target point \(x_0\) is defined as:

\[
\hat{f}(x_0) = \frac{1}{k} \sum_{x_i \in N_k(x_0)} y_i
\]

---

## 2. Bias-Variance Decomposition and Prediction Error

To evaluate the generalization ability of the KNN model, we use the squared error loss function. The expected prediction error (EPE) at the target point \(x_0\) can be decomposed into three independent parts.

!!! note "Theorem 2.1 (Bias-Variance Decomposition Theorem)"

    Let \(\hat{f}(x_0)\) be the estimator obtained from the training set \(\mathcal{D}_n\). Then the expected error of a new observation \(Y_0 = f(x_0) + \epsilon_0\) at the test point \(x_0\) has the following decomposition:

    \[
    \text{Err}(x_0) = \mathbb{E}_{\mathcal{D}_n, Y_0} \left[ (Y_0 - \hat{f}(x_0))^2 \right] = \sigma^2 + \left[ \text{Bias}(\hat{f}(x_0)) \right]^2 + \text{Var}(\hat{f}(x_0))
    \]

??? proof "Proof: Detailed Derivation of Bias-Variance Decomposition"

    To simplify notation, denote \(\hat{f}(x_0)\) by \(\hat{f}\) and \(f(x_0)\) by \(f\). First, transform \(Y_0 - \hat{f}\) identically and expand:

    \[
    Y_0 - \hat{f} = (f + \epsilon_0) - \hat{f} = (f - \mathbb{E}[\hat{f}]) + (\mathbb{E}[\hat{f}] - \hat{f}) + \epsilon_0
    \]

    Square both sides and take expectation \(\mathbb{E}_{\mathcal{D}_n, Y_0}[\cdot]\):

    \[
    \mathbb{E}[(Y_0 - \hat{f})^2] = \mathbb{E}\left[ \left( (f - \mathbb{E}[\hat{f}]) + (\mathbb{E}[\hat{f}] - \hat{f}) + \epsilon_0 \right)^2 \right]
    \]

    Expand into sums of squares and cross terms:

    \[
    \mathbb{E}[(Y_0 - \hat{f})^2] = \mathbb{E}[(f - \mathbb{E}[\hat{f}])^2] + \mathbb{E}[(\mathbb{E}[\hat{f}] - \hat{f})^2] + \mathbb{E}[\epsilon_0^2] + 2\mathbb{E}[(f - \mathbb{E}[\hat{f}])(\mathbb{E}[\hat{f}] - \hat{f})] + 2\mathbb{E}[(f - \mathbb{E}[\hat{f}])\epsilon_0] + 2\mathbb{E}[(\mathbb{E}[\hat{f}] - \hat{f})\epsilon_0]
    \]

    Now analyze each term one by one:

    1. *First term*: Since \(f - \mathbb{E}[\hat{f}]\) is constant with respect to sampling of the training set:

    \[
    \mathbb{E}[(f - \mathbb{E}[\hat{f}])^2] = (f - \mathbb{E}[\hat{f}])^2 = \left[ \text{Bias}(\hat{f}) \right]^2
    \]

    2. *Second term*: By definition of variance:

    \[
    \mathbb{E}[(\mathbb{E}[\hat{f}] - \hat{f})^2] = \text{Var}(\hat{f})
    \]

    3. *Third term*: Since the new test error \(\epsilon_0\) satisfies \(\mathbb{E}[\epsilon_0] = 0\) and \(\text{Var}(\epsilon_0) = \sigma^2\),

    \[
    \mathbb{E}[\epsilon_0^2] = \sigma^2
    \]

    4. *Cross term 1*: Because \(\mathbb{E}[\mathbb{E}[\hat{f}] - \hat{f}] = \mathbb{E}[\hat{f}] - \mathbb{E}[\hat{f}] = 0\), and \((f - \mathbb{E}[\hat{f}])\) is constant:

    \[
    \mathbb{E}[(f - \mathbb{E}[\hat{f}])(\mathbb{E}[\hat{f}] - \hat{f})] = (f - \mathbb{E}[\hat{f}]) \cdot \mathbb{E}[\mathbb{E}[\hat{f}] - \hat{f}] = 0
    \]

    5. *Cross terms 2 and 3*: Since the test error \(\epsilon_0\) is independent of the training set \(\mathcal{D}_n\), and \(\mathbb{E}[\epsilon_0] = 0\),

    \[
    \mathbb{E}[(f - \mathbb{E}[\hat{f}])\epsilon_0] = (f - \mathbb{E}[\hat{f}]) \cdot \mathbb{E}[\epsilon_0] = 0
    \]

    \[
    \mathbb{E}[(\mathbb{E}[\hat{f}] - \hat{f})\epsilon_0] = \mathbb{E}[\mathbb{E}[\hat{f}] - \hat{f}] \cdot \mathbb{E}[\epsilon_0] = 0
    \]

    In summary, all cross terms are zero. Summing the nonzero terms yields:

    \[
    \text{Err}(x_0) = \sigma^2 + \left[ \text{Bias}(\hat{f}(x_0)) \right]^2 + \text{Var}(\hat{f}(x_0))
    \]

### 2.2 Statistical Meaning of the Three Errors

* **Irreducible Error**: \(\sigma^2\). This is the inherent random noise of the system and cannot be eliminated by improving the algorithm.

* **Squared Bias (Bias\(^2\))**: Reflects the systematic deviation between the model's expected prediction and the true function \(f(x_0)\).

* **Variance**: Reflects the fluctuation of the model's predictions under repeated sampling of different training sets.

### 2.3 Asymptotic Properties of Bias and Variance in KNN

In the KNN algorithm, the hyperparameter \(k\) determines model complexity:

* When \(k=1\) (1-nearest neighbor), the estimator forces a fit to the closest training sample. As \(n \to \infty\), the nearest neighbor converges to \(x_0\), so that

\[
\text{Bias}^2(\hat{f}(x_0)) \to 0
\]

However, because only one observation is used, the prediction variance in the limit is approximately

\[
\text{Var}(\hat{f}(x_0)) = \text{Var}(y_i) = \sigma^2
\]

* As \(k\) increases, because response values from farther points are included, *the bias of the model increases* (underfitting); but because the average is taken over \(k\) samples, *the prediction variance decreases*, approximately

\[
\text{Var}(\hat{f}(x_0)) \approx \frac{\sigma^2}{k}
\]

---

## 3. Model Complexity and Effective Degrees of Freedom

In generalized additive or nonparametric models, model complexity is not easily measured by the number of parameters. Therefore, we introduce a general definition of degrees of freedom.

!!! info "Definition 3.1 (Effective Degrees of Freedom)"

    Assume the feature values \(X_i = x_i\) are fixed constants (non‑random). Let the true observation vector be \(Y = (Y_1, \dots, Y_n)^T\) and the fitted vector be \(\hat{Y} = (\hat{Y}_1, \dots, \hat{Y}_n)^T\). Then the effective degrees of freedom of the model is defined as

    \[
    \text{df}(\hat{f}) = \frac{1}{\sigma^2} \sum_{i=1}^n \text{Cov}(\hat{Y}_i, Y_i) = \frac{1}{\sigma^2} \text{Trace}\left( \text{Cov}(\hat{Y}, Y) \right)
    \]

### 3.1 Verification of Degrees of Freedom for Typical Models

* **1-nearest neighbor model**: The prediction for each point is its own observed value (i.e., \(\hat{Y}_i = Y_i\)). Hence

\[
\text{Cov}(\hat{Y}_i, Y_i) = \text{Cov}(Y_i, Y_i) = \sigma^2
\]

Substituting into the definition gives

\[
\text{df} = \frac{1}{\sigma^2} \sum_{i=1}^n \sigma^2 = n
\]

* **\(n\)-nearest neighbor model**: The prediction for each point is the overall sample mean (i.e., \(\hat{Y}_i = \frac{1}{n}\sum_{j=1}^n Y_j\)). Then

\[
\text{Cov}(\hat{Y}_i, Y_i) = \text{Cov}\left( \frac{1}{n} \sum_{j=1}^n Y_j, Y_i \right) = \frac{1}{n} \sigma^2
\]

Summing over \(i\) yields

\[
\text{df} = \frac{1}{\sigma^2} \sum_{i=1}^n \left( \frac{1}{n}\sigma^2 \right) = 1
\]

* **\(k\)-nearest neighbor model**: Similarly, each sample contributes \(\frac{1}{k}\sigma^2\) to the degrees of freedom, so the total effective degrees of freedom is

\[
\text{df} = \frac{n}{k}
\]

---

## 4. KNN Classification Model and 1NN Error Upper Bound

For classification tasks, KNN uses the majority vote rule.

!!! info "Definition 4.1 (KNN Classification Decision)"

    Let the class labels be \(c \in \mathcal{C}\). Given a test point \(x_0\), its neighborhood \(N_k(x_0)\) contains \(k\) samples. The hard classification prediction of KNN is

    \[
    \hat{y} = \arg\max_{c \in \mathcal{C}} \sum_{x_i \in N_k(x_0)} \mathbb{I}(y_i = c)
    \]

    where \(\mathbb{I}(\cdot)\) is the indicator function.

### 4.1 Proof of Asymptotic Error Rate Upper Bound for 1NN

In a binary classification problem (assuming labels \(Y \in \{0, 1\}\)), as \(n \to \infty\), the asymptotic error rate of 1NN has a classic upper bound governed by the Bayes error rate.

!!! note "Theorem 4.1 (Cover–Hart 1NN Error Upper Bound)"

    Let \(P(Y=1 \mid x)\) be the true conditional probability over the feature space. Denote the optimal Bayes error rate (the smallest possible error) at \(x_0\) as

    \[
    P^* = \min\{P(Y=1 \mid x_0), 1 - P(Y=1 \mid x_0)\}
    \]

    If the probability density and conditional probability are smooth, then as \(n \to \infty\), the expected error rate of the 1-nearest neighbor algorithm \(P_{1NN}\) satisfies

    \[
    P_{1NN} \le 2P^*(1 - P^*) \le 2P^*
    \]

??? proof "Proof: Derivation of the 1NN Asymptotic Error Upper Bound"

    Let \(x_{1NN}\) be the sample point in the training set closest to the target point \(x_0\). As \(n \to \infty\), due to the density of the feature space, the nearest neighbor approaches the target point geometrically:

    \[
    d(x_0, x_{1NN}) \to 0
    \]

    By the continuity assumption of the conditional probability function,

    \[
    P(Y=1 \mid x_{1NN}) \to P(Y=1 \mid x_0)
    \]

    For notational simplicity, let \(p = P(Y=1 \mid x_0)\). At the test point \(x_0\), the 1NN model makes an error in two cases: either the true label is 1 but the neighbor predicts 0, or the true label is 0 but the neighbor predicts 1. Therefore, as \(n \to \infty\), the asymptotic conditional error rate of 1NN is

    \[
    P_{1NN} = P(Y_0=1 \mid x_0)P(Y_{1NN}=0 \mid x_{1NN}) + P(Y_0=0 \mid x_0)P(Y_{1NN}=1 \mid x_{1NN})
    \]

    Substituting the limiting values yields

    \[
    P_{1NN} = p(1 - p) + (1 - p)p = 2p(1 - p)
    \]

    Now introduce the optimal Bayes error rate \(P^* = \min\{p, 1-p\}\). By this definition, clearly

    \[
    p(1 - p) = P^*(1 - P^*)
    \]

    Hence the 1NN error rate can be expressed exactly as

    \[
    P_{1NN} = 2P^*(1 - P^*)
    \]

    Since \(P^* \in [0, 0.5]\), we have \((1 - P^*) \le 1\). Immediately we obtain the upper bound

    \[
    P_{1NN} \le 2P^*
    \]

    That is, in the limit, the error rate of 1NN is at most twice the optimal Bayes error rate.

---

## 5. Distance Metrics and the Curse of Dimensionality

### 5.1 Mathematical Expressions of Common Distance Metrics

To define “neighbors”, distances between points must be computed. Let \(u, v \in \mathbb{R}^p\):

* **Euclidean Distance (Minkowski \(L_2\) norm)**:

\[
d_2(u, v) = \sqrt{\sum_{j=1}^p (u_j - v_j)^2}
\]

* **Mahalanobis Distance**: This eliminates scaling effects and accounts for correlations among features. Let \(\Sigma\) be the covariance matrix of the features:

\[
d_{\Sigma}(u, v) = \sqrt{(u - v)^T \Sigma^{-1} (u - v)}
\]

### 5.2 Geometric Derivation of the Curse of Dimensionality

The curse of dimensionality refers to the exponential increase in data sparsity as the dimension grows.

!!! note "Derivation of Neighborhood Edge Length in a High‑Dimensional Hypercube"

    Suppose \(n\) sample points are uniformly distributed in a \(p\)-dimensional unit hypercube \([0, 1]^p\). To capture a fixed proportion \(r = \frac{k}{n}\) of the local samples, we need to construct a sub‑hypercube whose volume is also a proportion \(r\) of the total volume.

    Let the target edge length of this sub‑hypercube be \(l\). Its volume is \(l^p\). Since the data are uniformly distributed, we have

    \[
    l^p = r \quad \Longrightarrow \quad l = r^{\frac{1}{p}}
    \]

    Plug in specific numeric values to see how the edge length \(l\) changes with dimension \(p\) (set \(r = 0.1\), i.e., attempt to capture \(10\%\) of the data):

    1. When \(p = 1\) (one‑dimensional line):

    \[
    l = 0.1^1 = 0.10
    \]

    2. When \(p = 2\) (two‑dimensional plane):

    \[
    l = 0.1^{0.5} \approx 0.32
    \]

    3. When \(p = 10\) (ten‑dimensional space):

    \[
    l = 0.1^{0.1} \approx 0.63
    \]

    4. When \(p = 100\) (one‑hundred‑dimensional space):

    \[
    l = 0.1^{0.01} \approx 0.955
    \]

    This shows that in 100‑dimensional space, to capture the nearest \(10\%\) of samples, we must cover as much as \(95.5\%\) of the range along each coordinate axis. At this point, the so‑called “local neighbors” are no longer geometrically “local”, causing distance‑based neighborhood estimation to completely break down.