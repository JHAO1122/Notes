# Chapter 9: Tree Models and Random Forests

This chapter discusses a class of extremely important nonparametric learning methods—decision trees and their ensemble extensions. Tree models, represented by Classification and Regression Trees (CART), recursively partition the feature space into several disjoint hyper-rectangular regions and make simple constant predictions within each region. We will systematically derive their splitting criteria, pruning algorithms, and deeply analyze the asymptotic normality and U-statistic theory of Random Forests.

---

## 1. Basic Setup of Classification and Regression Trees (CART)

!!! info "Definition 1.1 (Recursive Partitioning of Feature Space)"

    Given a training dataset \(\mathcal{D}_{n}=\{(x_{i}, y_{i})\}_{i=1}^{n}\), where the feature vectors \(x_{i} \in \mathbb{R}^p\). A decision tree of depth \(M\) partitions the entire input feature space \(X\) into \(M\) disjoint hyper-rectangular regions \(R_1, R_2, \dots, R_M\).

    The corresponding prediction function has the following general piecewise constant form:

    \[
    f(x) = \sum_{m=1}^M c_m \mathbb{I}(x \in R_m)
    \]

### 1.1 Optimal Parameter Selection for Regression Trees

For regression problems, if the commonly used squared error loss is adopted, the optimal prediction constant \(c_m\) within each partitioned region \(R_m\) is the sample mean of the response variable for all training samples in that region:

\[
\hat{c}_m = \text{avg}(y_i \mid x_i \in R_m)
\]

---

## 2. Splitting Criteria and Local Optimization (Greedy Algorithm)

Finding a partition of the feature space that minimizes the global residual sum of squares is an NP-Hard problem. In practice, CART employs a greedy algorithm for top-down, layer-by-layer optimal binary splitting.

!!! info "Definition 2.1 (Optimal Binary Split Optimization Objective)"

    Consider a current node region \(R\) to be split. We select the \(j\)-th feature as the splitting variable and choose \(s\) as the split point, thereby dividing the region into two sub-regions:

    \[
    R_1(j, s) = \{x \mid x_j \le s\}, \quad R_2(j, s) = \{x \mid x_j > s\}
    \]

    For regression trees, the local optimization objective for the splitting variable \(j\) and split point \(s\) can be expressed as:

    \[
    \min_{j, s} \left[ \min_{c_1} \sum_{x_i \in R_1(j, s)} (y_i - c_1)^2 + \min_{c_2} \sum_{x_i \in R_2(j, s)} (y_i - c_2)^2 \right]
    \]

### 2.1 Impurity Measures for Classification Trees

For classification tree problems, since the response variable \(y_i \in \{1, 2, \dots, K\}\) is a discrete label, we first define the empirical probability of class \(k\) in node \(m\) as:

\[
\hat{p}_{mk} = \frac{1}{|R_m|} \sum_{x_i \in R_m} \mathbb{I}(y_i = k)
\]

To measure the impurity of a node, the following two nonlinear metrics are commonly used instead of the misclassification rate:

* **Gini Index**:

\[
H_m = \sum_{k=1}^K \hat{p}_{mk} (1 - \hat{p}_{mk}) = 1 - \sum_{k=1}^K \hat{p}_{mk}^2
\]

* **Cross-Entropy / Information Entropy**:

\[
H_m = -\sum_{k=1}^K \hat{p}_{mk} \ln \hat{p}_{mk}
\]

---

## 3. Cost-Complexity Pruning Algorithm

If the tree model is allowed to grow without restriction until each leaf node contains only one sample, it will lead to overfitting and extremely large estimation variance. To control model complexity, it is necessary to prune the fully grown original tree \(T_0\).

!!! info "Definition 3.1 (Cost-Complexity Objective Function)"

    Let \(T \subset T_0\) be any subtree obtained by pruning \(T_0\). Let \(|T|\) denote the total number of leaf nodes in subtree \(T\), \(R_m\) the region covered by the \(m\)-th leaf node, and \(Q_m(T)\) the total impurity loss corresponding to that node. Introducing a regularization tuning parameter \(\alpha \ge 0\), the cost-complexity is defined as:

    \[
    C_\alpha(T) = \sum_{m=1}^{|T|} |R_m| Q_m(T) + \alpha |T|
    \]

!!! note "Theorem 3.1 (Weakest Link Pruning Property)"

    For any given complexity penalty factor \(\alpha\), there exists a unique minimal optimal subtree \(T_\alpha = \arg\min_{T \subset T_0} C_\alpha(T)\). As \(\alpha\) increases continuously from 0 to \(\infty\), the sequence of optimal subtrees is **finite and nested**.

??? proof "Proof: Derivation of the Critical Pruning Threshold \(\alpha\) for a Single Node"

    Consider a particular internal node \(t\) in tree \(T_0\).
  
    1. Let \(T_t\) be the subtree rooted at \(t\) (containing multiple leaf nodes). If not pruned, the contribution to the cost-complexity of this local subtree is:

    \[
    C_\alpha(T_t) = R(T_t) + \alpha |T_t|
    \]

    where \(R(T_t) = \sum_{m \in \text{leaves}(T_t)} |R_m| Q_m(T)\).
  
    2. If we prune the subtree rooted at \(t\) and collapse it into a single leaf node \(t\), its cost-complexity contribution becomes:

    \[
    C_\alpha(\{t\}) = R(t) + \alpha
    \]

    When \(\alpha\) is small, the accuracy gain from not pruning is large, so \(C_\alpha(T_t) < C_\alpha(\{t\})\). However, as \(\alpha\) increases to a certain critical point, the two become equal. By setting \(C_\alpha(T_t) = C_\alpha(\{t\})\), we can solve:

    \[
    R(T_t) + \alpha |T_t| = R(t) + \alpha \implies \alpha (|T_t| - 1) = R(t) - R(T_t)
    \]

    This yields the critical pruning threshold function:

    \[
    \alpha = \frac{R(t) - R(T_t)}{|T_t| - 1}
    \]

    In the entire large tree structure, by sequentially traversing and calculating the \(\alpha\) values for all internal nodes, and preferentially pruning the node with the smallest \(\alpha\), we can construct the discrete optimal nested subtree chain from bottom to top.

---

## 4. Random Forest and Out-of-Bag Estimation

Random Forest further reduces the variance of a single decision tree through Bootstrap Aggregation (Bagging) and random feature selection.

### 4.1 Core Mechanism of the Algorithm

1. **Bootstrap Sampling**: From the original dataset containing \(n\) samples, perform \(B\) rounds of resampling with replacement, each with sample size \(n\).
2. **Random Node Construction**: When splitting at each node of each tree, instead of searching over all \(p\) features, randomly select without replacement a subset of \(m \ll p\) features (typically \(m = \sqrt{p}\)), and then choose the optimal split point from this subset.

### 4.2 Out-of-Bag Error Estimation

Since bootstrap sampling is with replacement, for any sample \(i\), the probability that it is never selected in a particular resampling is:

\[
P(\text{sample } i \text{ is not selected}) = \left( 1 - \frac{1}{n} \right)^n \xrightarrow{n \to \infty} e^{-1} \approx 0.368
\]

This means that on average, about 36.8% of the training data do not participate in the construction of each tree; this portion of data is called the **Out-of-Bag (OOB)** data for that tree. We can directly use the OOB samples as a validation set to compute the generalization error, thus avoiding the need for tedious cross-validation.

---

## 5. Asymptotic Normality of Random Forests and U-Statistics Theory

In recent years, the asymptotic statistical properties of random forests have been rigorously established through **U-statistics** theory. This theory provides a mathematical foundation for estimating the prediction variance of random forests and constructing confidence intervals.

!!! info "Definition 5.1 (U-Statistic Predictor Based on Incomplete Sub-samples)"

    Let \(Z_i = (X_i, Y_i)\). Let \(h(Z_{i_1}, \dots, Z_{i_r})\) denote a single tree estimator constructed from an incomplete sub-sample of size \(r\) (\(r < n\)) drawn without replacement. If we exhaustively traverse all possible \(\binom{n}{r}\) combinations, the minimum-variance unbiased estimator (complete U-statistic) takes the form:

    \[
    U_n = \frac{1}{\binom{n}{r}} \sum_{(i)} h(Z_{i_1}, \dots, Z_{i_r})
    \]

### 5.1 Asymptotic Variance Theorem with Shared Single Observation

The prediction of a random forest can be approximated as a variant of such U-statistics.

!!! note "Theorem 5.1 (Mentch & Hooker Asymptotic Normality Theorem)"

    Under a high-dimensional large-sample setting, as the total sample size \(n \to \infty\) and the sub-sample size \(r\) for a single tree satisfies \(r = o(\sqrt{n})\), the random forest prediction under the incomplete U-statistic framework satisfies asymptotic normality:

    \[
    \frac{\sqrt{n}(U_n - \theta)}{\sqrt{r^2 \zeta_{1,r}}} \xrightarrow{d} N(0, 1)
    \]

    Here, the key covariance metric \(\zeta_{1,r}\) is defined as the covariance between two independent trees that **share only one identical observation**:

    \[
    \zeta_{1,r} = \text{Cov}\left( h(Z_1, Z_2, \dots, Z_r), h(Z_1, Z_2', \dots, Z_r') \right)
    \]

Through this theory, we only need to extract the prediction results from those tree pairs in the random forest that contain exactly one common sample and compute their sample covariance, thereby enabling the quantification of prediction variance at a given point.