# Chapter 5: Classification Problems and Discriminant Analysis

In the previous chapters, we mainly discussed regression problems where the response variable \(Y\) is continuous. This chapter focuses on classification problems, where the response variable \(Y\) takes a finite set of discrete labels (e.g., \(\{0, 1\}\) or \(\{1, 2, \dots, K\}\)). We will systematically derive the two cornerstones of soft classifiers: Logistic Regression and Linear/Quadratic Discriminant Analysis (LDA/QDA), and provide the relevant mathematical derivations and property proofs rigorously.

---

## 1. Basic Setup of Classification Problems and Decision Theory

!!! info "Definition 1.1 (Classifier and 0-1 Loss)"

    Let the training dataset be \(\mathcal{D}_{n}=\{(x_{i}, y_{i})\}_{i=1}^{n}\), where the feature vector \(x_{i} \in \mathbb{R}^p\) and the response variable \(y_{i} \in \{0, 1\}\) (or \(\{-1, +1\}\) in binary classification). The goal of a classification task is to learn a mapping function (i.e., a classifier):

    \[
    f: \mathbb{R}^p \longrightarrow \{0, 1\}
    \]

    Its performance is typically measured by the 0-1 loss function:

    \[
    L(f(X), Y) = \mathbb{I}(f(X) \neq Y) = \begin{cases} 0 & \text{if } f(X) = Y \\ 1 & \text{if } f(X) \neq Y \end{cases}
    \]

### 1.1 Bayes Optimal Classifier

To minimize the overall expected risk, we need to find the optimal classification rule.

!!! note "Theorem 1.1 (Bayes Decision Criterion)"

    Under the 0-1 loss, the decision function that minimizes the expected risk \(R(f) = \mathbb{E}_{X,Y}[L(f(X), Y)]\) is the Bayes optimal classifier \(f^*(x)\), given by:

    \[
    f^*(x) = \arg\max_{c \in \{0, 1\}} P(Y = c \mid X = x)
    \]

??? proof "Proof: Derivation of the Optimality of the Bayes Decision"

    Using the law of total expectation, the expected risk can be expanded by conditioning on the feature \(X\):

    \[
    R(f) = \mathbb{E}_X \left[ \mathbb{E}_{Y \mid X} [ \mathbb{I}(f(X) \neq Y) \mid X ] \right]
    \]

    To minimize the overall expectation, it suffices to minimize the conditional expectation (conditional risk) for each fixed \(X = x\):

    \[
    \mathbb{E}_{Y \mid X} [ \mathbb{I}(f(x) \neq Y) \mid X = x ] = P(f(x) \neq Y \mid X = x)
    \]

    For binary classification, the above expression can be written as:

    \[
    P(f(x) \neq Y \mid X = x) = 1 - P(Y = f(x) \mid X = x)
    \]

    Clearly, minimizing \(1 - P(Y = f(x) \mid X = x)\) is equivalent to maximizing \(P(Y = f(x) \mid X = x)\). Therefore, the predicted value \(f(x)\) of the classifier at point \(x\) should be the class that maximizes the conditional probability:

    \[
    f^*(x) = \arg\max_{c \in \{0, 1\}} P(Y = c \mid X = x)
    \]

    This completes the proof that the posterior probability maximization criterion (Bayes criterion) is optimal under the 0-1 loss.

---

## 2. Logistic Regression

Logistic regression is a parametric method that directly models the posterior conditional probability \(P(Y=1 \mid X=x)\).

!!! info "Definition 2.1 (Logistic Model and Logit Transformation)"

    Let \(\eta(x) = P(Y = 1 \mid X = x)\). To map a linear predictor defined on \(\mathbb{R}^p\) to the probability interval \([0, 1]\), the model introduces the Logit transformation (log-odds):

    \[
    \ln \left( \frac{\eta(x)}{1 - \eta(x)} \right) = \beta_0 + \beta^{\top}x = \beta_0 + \sum_{j=1}^p \beta_j X_j
    \]

    By algebraic inversion, the conditional probability can be directly expressed in the Sigmoid function form:

    \[
    \eta(x) = \frac{\exp(\beta_0 + \beta^{\top}x)}{1 + \exp(\beta_0 + \beta^{\top}x)}
    \]

### 2.1 Maximum Likelihood Estimation (MLE)

Since the response variable \(Y_i \mid x_i \sim \text{Bernoulli}(\eta(x_i))\), we can write the likelihood function for the parameter \(\theta = (\beta_0, \beta^{\top})^{\top}\). For brevity, let the augmented feature vector be \(x_i = (1, x_{i1}, \dots, x_{ip})^{\top}\).

!!! note "Theorem 2.2 (Log-likelihood Function and Gradient for Logistic Regression)"

    The maximum likelihood estimation of the logistic regression parameters is achieved by minimizing the negative log-likelihood function (i.e., the cross-entropy loss):

    \[
    \ell(\theta) = -\sum_{i=1}^n \left[ y_i \ln \eta(x_i) + (1 - y_i) \ln (1 - \eta(x_i)) \right]
    \]

    Its gradient vector (first derivative) with respect to the parameter \(\theta\) is:

    \[
    \nabla \ell(\theta) = \sum_{i=1}^n (\eta(x_i) - y_i) x_i
    \]

??? proof "Proof: Detailed Calculus Derivation of the Gradient Formula"

    First, compute the derivative of the Sigmoid function \(\eta_i = \eta(x_i) = \frac{1}{1 + e^{-\theta^{\top}x_i}}\) with respect to \(\theta\). By the chain rule:

    \[
    \frac{\partial \eta_i}{\partial \theta} = \frac{e^{-\theta^{\top}x_i}}{(1 + e^{-\theta^{\top}x_i})^2} \cdot x_i = \eta_i (1 - \eta_i) x_i
    \]

    Next, take the partial derivative of the log-likelihood term \(L_i = y_i \ln \eta_i + (1 - y_i) \ln (1 - \eta_i)\) with respect to \(\theta\):

    \[
    \frac{\partial L_i}{\partial \theta} = \frac{y_i}{\eta_i} \frac{\partial \eta_i}{\partial \theta} - \frac{1 - y_i}{1 - \eta_i} \frac{\partial \eta_i}{\partial \theta}
    \]

    Substitute \(\frac{\partial \eta_i}{\partial \theta} = \eta_i (1 - \eta_i) x_i\) into the above expression:

    \[
    \frac{\partial L_i}{\partial \theta} = \left[ \frac{y_i}{\eta_i} \eta_i (1 - \eta_i) - \frac{1 - y_i}{1 - \eta_i} \eta_i (1 - \eta_i) \right] x_i
    \]

    \[
    \frac{\partial L_i}{\partial \theta} = \left[ y_i (1 - \eta_i) - (1 - y_i) \eta_i \right] x_i = (y_i - \eta_i) x_i
    \]

    Therefore, the gradient of the negative log-likelihood function \(\ell(\theta) = -\sum_{i=1}^n L_i\) is:

    \[
    \nabla \ell(\theta) = -\sum_{i=1}^n (y_i - \eta_i) x_i = \sum_{i=1}^n (\eta(x_i) - y_i) x_i
    \]

    This completes the proof. Since this objective function is strictly log-concave with respect to \(\theta\), it is typically solved iteratively using the Newton-Raphson method or gradient descent in practice.

---

## 3. Generative Models: Discriminant Analysis

Unlike logistic regression, discriminant analysis belongs to the family of **generative models**. It first models the conditional probability density of the features within each class \(k\), \(f_k(x) = f(X=x \mid Y=k)\), and the class prior probabilities \(\pi_k = P(Y=k)\), and then uses Bayes' theorem to derive the posterior probabilities.

By Bayes' theorem, given \(X=x\), the posterior probability of belonging to class \(k\) is:

\[
P(Y = k \mid X = x) = \frac{\pi_k f_k(x)}{\sum_{l=1}^K \pi_l f_l(x)}
\]

We assume that the data within each class follow a multivariate normal distribution:

\[
f_k(x) = \frac{1}{(2\pi)^{p/2} |\Sigma_k|^{1/2}} \exp \left( -\frac{1}{2} (x - \mu_k)^{\top} \Sigma_k^{-1} (x - \mu_k) \right)
\]

### 3.1 Linear Discriminant Analysis (LDA)

!!! info "Definition 3.1 (Basic Assumption of LDA)"

    Linear Discriminant Analysis (LDA) assumes that the covariance matrices of all classes are equal, i.e.:

    \[
    \Sigma_k = \Sigma, \quad \forall k = 1, 2, \dots, K
    \]

!!! note "Theorem 3.1 (Linear Discriminant Function of LDA)"

    Under the assumption of equal covariance matrices across classes, maximizing the Bayes posterior probability is equivalent to maximizing the following linear discriminant function \(\delta_k(x)\) in terms of \(x\):

    \[
    \delta_k(x) = x^{\top} \Sigma^{-1} \mu_k - \frac{1}{2} \mu_k^{\top} \Sigma^{-1} \mu_k + \ln \pi_k
    \]

??? proof "Proof: Mathematical Derivation of the LDA Linear Discriminant Function"

    To maximize \(P(Y=k \mid X=x)\), since the denominator \(\sum_{l=1}^K \pi_l f_l(x)\) of Bayes' formula is the same for all classes \(k\), we only need to maximize the numerator term:

    \[
    \arg\max_k P(Y=k \mid X=x) = \arg\max_k \left[ \pi_k f_k(x) \right]
    \]

    Since the logarithm function is monotonically increasing, we take the natural logarithm of the numerator term:

    \[
    \ln \left( \pi_k f_k(x) \right) = \ln \pi_k + \ln f_k(x)
    \]

    Substituting and expanding the multivariate normal density \(f_k(x)\) (with \(\Sigma_k = \Sigma\)):

    \[
    \ln \left( \pi_k f_k(x) \right) = \ln \pi_k - \frac{p}{2}\ln(2\pi) - \frac{1}{2}\ln|\Sigma| - \frac{1}{2}(x - \mu_k)^{\top} \Sigma^{-1} (x - \mu_k)
    \]

    Note that the terms \(-\frac{p}{2}\ln(2\pi) - \frac{1}{2}\ln|\Sigma|\) are constants independent of the class \(k\) and can be removed when making the maximization decision over \(k\). Expand the remaining part:

    \[
    - \frac{1}{2}(x - \mu_k)^{\top} \Sigma^{-1} (x - \mu_k) + \ln \pi_k = -\frac{1}{2} \left[ x^{\top}\Sigma^{-1}x - 2x^{\top}\Sigma^{-1}\mu_k + \mu_k^{\top}\Sigma^{-1}\mu_k \right] + \ln \pi_k
    \]

    Again, note that the term \(x^{\top}\Sigma^{-1}x\) is also independent of the class \(k\) and can be ignored. After simplification and consolidation, the core part relevant to \(k\) remains:

    \[
    \delta_k(x) = x^{\top} \Sigma^{-1} \mu_k - \frac{1}{2} \mu_k^{\top} \Sigma^{-1} \mu_k + \ln \pi_k
    \]

    Since the highest-order term in the unknown variable \(x\) in the above expression is linear, the decision boundary between different classes (i.e., the hyperplane where \(\delta_k(x) = \delta_l(x)\)) is geometrically a **linear hyperplane**. This completes the proof.

### 3.2 Quadratic Discriminant Analysis (QDA)

!!! info "Definition 3.2 (Basic Assumption of QDA)"

    Quadratic Discriminant Analysis (QDA) relaxes the assumption of equal covariance matrices, allowing different classes to have their own distinct feature covariance matrices \(\Sigma_k\).

In this case, the quadratic terms and determinant terms cannot be eliminated from the log-likelihood expression. Preserving the full feature parameters, the quadratic discriminant function of QDA is expressed as:

\[
\delta_k(x) = -\frac{1}{2} \ln |\Sigma_k| - \frac{1}{2} (x - \mu_k)^{\top} \Sigma_k^{-1} (x - \mu_k) + \ln \pi_k
\]

Since the term \(x^{\top}\Sigma_k^{-1}x\) cannot be eliminated, the decision boundary between different classes will be a **quadratic surface** (e.g., hyperbola, ellipse, or parabola).

---

## 4. Deep Statistical Connection between Logistic Regression and LDA

If we carefully examine the log-odds of the two classes under LDA in a binary classification problem, we arrive at a striking conclusion.

!!! note "Theorem 4.1 (Implicit Logit Form of LDA)"

    For the LDA model, the log-odds of class 1 versus class 0 given features \(x\) have an exact linear form:

    \[
    \ln \left( \frac{P(Y=1 \mid X=x)}{P(Y=0 \mid X=x)} \right) = \alpha_0 + \alpha^{\top}x
    \]

    where the parameters correspond as follows:

    \[
    \alpha_0 = \ln \frac{\pi_1}{\pi_0} - \frac{1}{2}\mu_1^{\top}\Sigma^{-1}\mu_1 + \frac{1}{2}\mu_0^{\top}\Sigma^{-1}\mu_0
    \]

    \[
    \alpha = \Sigma^{-1}(\mu_1 - \mu_0)
    \]

### 4.1 Core Difference: Conditional Likelihood vs. Full Likelihood

Although from the above it is evident that LDA and logistic regression both lead to a linear decision boundary for the posterior odds in form, their fundamental difference lies in the **basis of parameter estimation and the strength of assumptions**:

* **Logistic Regression (Conditional Likelihood)**: It only models the conditional probability \(P(Y \mid X)\) in a discriminative manner. It does not require any normality or independence assumptions about the marginal distribution of the features \(X\). Therefore, logistic regression is more robust when the feature data do not follow a normal distribution (e.g., when they include discrete dummy variables).

* **LDA (Full Likelihood)**: It relies on the strong generative assumption that the features jointly follow a multivariate normal distribution. When this assumption is truly satisfied, LDA can fully exploit the information in the feature space, and its parameter estimation efficiency (variance) in the limiting case is significantly better than that of logistic regression.

---

## 5. Naive Bayes Classifier

When the feature dimension \(p\) is very large, accurately estimating the full covariance matrix \(\Sigma\) in LDA or multiple \(\Sigma_k\) in QDA with limited samples suffers from the severe "curse of dimensionality" (too many parameters leading to high variance in estimates). The Naive Bayes algorithm introduces an extreme conditional independence assumption to simplify the computation.

!!! info "Definition 5.1 (Naive Bayes Conditional Independence Assumption)"

    The Naive Bayes classifier assumes that, given the class label \(Y=k\), the individual feature variables \(X_1, X_2, \dots, X_p\) are conditionally independent of each other.

Under this assumption, the joint conditional density function can be factorized as the product of the marginal densities of each variable:

\[
f_k(x) = \prod_{j=1}^p f_{kj}(x_j)
\]

Substituting into the Bayes decision rule, the final discriminant criterion for Naive Bayes becomes:

\[
f_{\text{NB}}(x) = \arg\max_{k} \left\{ \pi_k \prod_{j=1}^p f_{kj}(x_j) \right\}
\]

* *Algorithm Advantage*: By decomposing the multi-dimensional density estimation into \(p\) independent one-dimensional density estimates \(f_{kj}(x_j)\), Naive Bayes reduces the order of magnitude of parameter estimation from \(\mathcal{O}(p^2)\) directly to \(\mathcal{O}(p)\). For continuous features, the one-dimensional density \(f_{kj}\) can be fitted by a univariate normal distribution or univariate kernel density estimation (KDE); for discrete features, simple frequency histograms are used for counting estimates.