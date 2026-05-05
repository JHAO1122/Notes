---
tags:
  - Group Meeting Report
  - Knockoff
  - FDR Control
  - Statistical Machine Learning
---

# 🎙️ 2025 Fall Group Meeting: Derandomized Knockoff

> **"Randomness in Model-X Knockoffs provides exact FDR control, but Derandomized Knockoffs restore selection stability without sacrificing statistical guarantees."**

!!! abstract "Seminar Abstract"
    * **Topic**: Derandomized Knockoff: Method and Theoretical Analysis
    * **Core Problem**: Traditional Model-X Knockoff involves randomness in the generation process, leading to instability in feature selection results. How can we obtain a stable set of selected features through a derandomization procedure?
    * **Theoretical Tool**: By introducing $e$-values and $e$-processes, we reveal the intrinsic equivalence between the Knockoff method and the e-BH procedure.
    * **Main Result**: By repeatedly generating Knockoff variables and aggregating $e$-values, we can substantially improve the stability of feature selection while maintaining rigorous False Discovery Rate (FDR) control in theory.

---

## 1. Review of the Model-X Method

### 1.1 Background

In classical statistics and high-dimensional feature selection models, the nature of the data often leads to different approaches:

* **Fixed-X method**: The design matrix $X$ is assumed to be known and non-random, i.e., treated as a fixed constant matrix.
* **Model-X method**: In practice, most observed data are randomly generated. Therefore, we assume $X$, as a data matrix, follows some known or unknown joint probability distribution.
* Both methods share a similar goal: select features truly associated with the response variable $Y$ while controlling the **False Discovery Rate (FDR)**.

### 1.2 Core of the Method

!!! info "Definition 1.1 (Model-X Knockoff Variables)"

    Let $\tilde{X} = (\tilde{X}_1, \tilde{X}_2, \dots, \tilde{X}_p)$ be a new set of variables generated from the original variables $X = (X_1, X_2, \dots, X_p)$ via the Model-X Knockoff transformation. Then they must satisfy the following properties:

    * **Swap property (exchangeability)**: For any subset $S \subseteq \{1, \dots, p\}$, swapping the features in $S$ with their corresponding Knockoff variables leaves the joint distribution unchanged:

    $$
    (X, \tilde{X})_{\operatorname{swap}(S)} \stackrel{d}{=} (X, \tilde{X})
    $$

    * **Conditional independence**: For the response variable $Y$, given the original features $X$, the Knockoff variables $\tilde{X}$ must be independent of $Y$:

    $$
    \tilde{X} \perp\!\!\!\perp Y \mid X
    $$

For feature selection statistics, the **Lasso Coefficient-Difference (LCD)** statistic is commonly used:

$$
W_j = Z_j - \tilde{Z}_j
$$

where feature importance is obtained via the regularization path:

$$
Z_j = |\hat{b}_j(\lambda)|, \quad \tilde{Z}_j = |\hat{b}_{j+p}(\lambda)|
$$

Under the **Model-X Gauss** setting, a typical explicit construction takes the form:

$$
\tilde{X} = X(I - \Sigma^{-1} S) + EC
$$

where we require $X \sim N(0, \Sigma)$, and $S = \text{diag}(s_1, s_2, \dots, s_p)$ satisfies $0 \le S \le 2\Sigma$. The coefficient matrix $C$ for the error term must satisfy:

$$
CC^T = 2S - S\Sigma^{-1}S
$$

The statistic $W_j$ satisfies the important **coin-flip property**, which provides the symmetry basis for FDR control.

* (Supplementary) For multivariate Gaussian random variables, mature R packages for Knockoff are available:

```
Xk <- create.gaussian(X, mu, Sigma, diag_s = diags)
```


### 2. Introduction to the Derandomized Knockoff Method

### 2.1 Background and Motivation

In the Model-X method, the generation of Knockoff variables is inherently random. This randomness causes the final set of selected features to vary each time the algorithm is run, leading to instability in the feature selection results. This instability can be intuitively seen from the generation formula for Knockoff variables (taking the Gaussian case as an example): $$\widetilde{X} = X\bigl(I - \Sigma^{-1} S\bigr) + E C$$

Compared to the $p$-value used in traditional hypothesis testing, the $e$-value has excellent composability properties (e.g., it can be directly averaged or multiplied). Therefore, the paper designs a derandomization procedure based on $e$-values.

### 2.2 Core Definitions

!!! info "Definition 2.1 (p-value)"
    * **Definition**: Under the null hypothesis $H_0$, the probability of observing a test statistic as extreme as or more extreme than the current data.
    * Under $H_0$ and when the test statistic is continuous, the $p$-value follows a uniform distribution $U[0,1]$.
    * **Decision rule**: If the $p$-value is less than a pre-specified significance level $\alpha$, reject the null hypothesis.

!!! info "Definition 2.2 (e-value and e-process)"
    
    * **Definition (e-value)**: On a probability space equipped with a filtration, if a non-negative random variable $E \ge 0$ (a.s.) satisfies $\mathbb{E}_{H_0}[E] \le 1$ under $H_0$, i.e.,
    
    $$
    \mathbb{E}_{H_0}[E] \le 1
    $$
    
    then $E$ is called an **e-value**.

    * **Common constructions of e-values**:
        1. Based on p-value: $E = \frac{1(p \le \alpha)}{\alpha}$.
        2. Likelihood ratio: $E = \frac{f_1(x)}{f_0(x)}$, where $f_0, f_1$ are the densities under $H_0$ and $H_1$, respectively.

    * **e-process**: In sequential testing, a sequence of random variables $E_t \ge 0$ is constructed, with $E_t$ adapted to a filtration $\mathcal{F}_t$. If it forms a **supermartingale** under $H_0$, i.e.,
    
    $$
    \mathbb{E}_{H_0}(E_t \mid \mathcal{F}_{t-1}) \le E_{t-1}, \quad E_0 \le 1
    $$
    
    then it is called an e-process. In this case, the Optional Stopping theorem can be used to control error rates.

Based on the properties of e-BH above, we can now safely aggregate multiple Knockoff e-values without breaking FDR control.

### 2.3 Main Theory and Properties

!!! success "Theorem 2.1 (Equivalence Theorem, THM1)"
    
    The Knockoff method is equivalent to the e-BH method based on the above $e$-values. That is, the sets of selected features are exactly the same:

    $$
    S_{kn} = S_{ebh}
    $$

??? proof "Proof sketch for THM 1 (click to expand)"

    Based on the construction of $e_j$, for any $j \in S_{kn}$, since $W_j \ge T$, the indicator function equals 1, yielding:
    $e_j = \frac{p}{1 + \sum_{k \in [p]} \mathbb{1}(W_k \le -T)} \ge \frac{p}{\alpha \hat{k}}$.
    Therefore, if a feature belongs to $S_{kn}$, it must also satisfy the e-BH cutoff condition, i.e., $j \in S_{ebh}$. Conversely, if $j \notin S_{kn}$, then $e_j = 0$, so it will never be selected.


!!! success "Theorem 2.2 (FDR control of e-BH procedure, THM2)"

    Assume that any generated $e_1, e_2, \dots, e_p$ satisfy the following expectation and boundary constraint under the null:

    $$
    \sum_{j \in \mathcal{H}_0} \mathbb{E}[e_j] \le p
    $$

    Then the set of features $S_{ebh}$ selected by the e-BH procedure satisfies strict FDR control:

    $$
    \text{FDR} < \alpha
    $$

    *(Note: Here $\text{FDR} = \mathbb{E}\left[ \frac{\sum_{j \in \mathcal{H}_0} \mathbb{1}\left(j \in S_{ebh}\right)}{|\mathcal{S}_{\text{ebh}}| \vee 1} \right]$)*

??? proof "Proof sketch for THM 2 (click to expand)"

    For any feature $j \in H_0$ under the null hypothesis, if $j \in S_{ebh}$, then necessarily $e_j \ge \frac{p}{\alpha |S_{ebh}|}$.
    Hence we can use the following bound on the indicator:
    $\mathbb{1}\{j \in S_{ebh}\} = \mathbb{1}\left\{ e_j \ge \frac{p}{\alpha |S_{ebh}|} \right\} \le \frac{\alpha |S_{ebh}| e_j}{p}$.
    Substituting this into the definition of FDR and taking expectation, the term $|S_{ebh}|$ cancels, leading to $\text{FDR} \le \alpha$.


Based on the properties of e-BH above, we can now safely aggregate multiple Knockoff e-values without breaking FDR control.

!!! success "Theorem 2.3 (FDR control of Derandomized Knockoff, THM3)"
    
    For any $\alpha_{\text{kn}}, \alpha_{\text{ebh}} \in (0, 1)$, and any given number of Knockoff runs $M \ge 1$, the final selection set $\mathcal{S}_{\text{kn-derand}}$ computed by the derandomized algorithm always satisfies:

    $$
    \text{FDR} \le \alpha_{\text{ebh}}
    $$

    * **Asymptotic property under the strong law of large numbers**: As the number of aggregations $M \to \infty$, the average $e$-value converges almost surely (a.s.) to its conditional expectation:

    $$
    e_j^{\text{avg}} \xrightarrow[M \to \infty]{\text{a.s.}} e_j^{\infty} := \mathbb{E}[e_j^{(1)} \mid X, Y]
    $$

### 2.4 Derandomized Knockoff Algorithm

Using the above theoretical support, the core of the algorithm is to generate $M$ different Knockoff copies, convert them to $e$-values, average them, and finally apply the e-BH procedure for selection.

!!! abstract "Algorithm: Derandomized Knockoff Procedure"

    **Input**: Design matrix $X$, response variable $Y$; target e-BH FDR level $\alpha_{ebh}$, base Knockoff threshold level $\alpha_{kn}$; number of iterations $M$.

    **Steps 1 to 4 (Generation and Transformation)**:
    For $m = 1, 2, \dots, M$:

    * Sample the $m$-th Knockoff matrix $\tilde{X}^{(m)}$.
    
    * Compute feature statistics:

    $$
    W^{(m)} = W([X, \tilde{X}^{(m)}], Y)
    $$

    * Compute the $m$-th threshold $T^{(m)}$ using $\alpha_{kn}$:

    $$
    T^{(m)} = \inf\left\{ t > 0 : \frac{1 + \sum_{j} \mathbb{1}(W_j^{(m)} \le -t)}{\sum_{j} \mathbb{1}(W_j^{(m)} \ge t)} \le \alpha_{kn} \right\}
    $$

    * For all features $j \in [p]$, compute the single-run $e$-value:

    $$
    e_j^{(m)} = p \cdot \frac{\mathbb{1}(W_j^{(m)} \ge T^{(m)})}{1 + \sum_{k=1}^p \mathbb{1}(W_k^{(m)} \le -T^{(m)})}
    $$

    **Step 5 (Aggregation)**:
    After the loop, compute the average $e$-value for each feature:

    $$
    e_j^{avg} = \frac{1}{M}\sum_{m=1}^{M}e_j^{(m)}, \quad \forall j \in [p]
    $$

    **Step 6 (e-BH Cutoff)**:
    Find the cutoff rank $\hat{k}$:

    $$
    \hat{k} = \max\left\{ k \in [p] : e_{(k)}^{avg} \ge \frac{p}{\alpha_{ebh} k} \right\}
    $$

    **Output**: The final derandomized feature selection set:

    $$
    S_{kn-derand} := \left\{ j \in [p] : e_j^{avg} \ge \frac{p}{\alpha_{ebh} \hat{k}} \right\}
    $$

    * **Recommended parameter choice**: In practice, the original paper derives and recommends using the following ratio of parameters, which (under the guarantee of Theorem 3) achieves an optimal balance between Power and FDR:

    $$
    \alpha_{ebh} = 2\alpha_{kn}
    $$


## 3. Extension: Copula-Model-X Knockoff

In practice, if the joint distribution $F(x_1, \dots, x_p)$ of the data is not known to be Gaussian, but only the marginal distributions $F_j(x_j)$ are known, we can use a **Copula model** to transform the problem into Gaussian Knockoff.

**Implementation steps**:

1. **Marginal distribution estimation**: Estimate the empirical CDF of each feature, denoted $\hat{F}_j$.

2. **Normal quantile transform**: Apply the probability integral transform to project the original data into the $[0, 1]$ space, and then map to the latent Gaussian space:
   
   $\hat{U}_{ij} = \hat{F}_j(X_{ij}), \quad \tilde{Z}_{ij} = \Phi^{-1}(\hat{U}_{ij})$

3. **Generate latent Knockoff variables**: In the transformed Gaussian space, generate Gaussian Knockoff variables corresponding to the latent variables $\tilde{Z}$.

4. **Inverse transform back to the original space**: Finally, apply the inverse transform to map the Knockoff variables back to the original feature scale:

   $U_{ij}^{\tilde{Z}} = \Phi(\tilde{Z}_{ij}), \quad \tilde{X}_{ij} = \hat{F}_j^{-1}(U_{ij}^{\tilde{Z}})$
