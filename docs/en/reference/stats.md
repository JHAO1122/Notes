# 📖 Mathematical Statistics

This module continuously compiles core foundational concepts related to mathematical statistics.

---

## 1. Hypothesis Testing and Confusion Matrix

Under the frequentist Neyman-Pearson framework, hypothesis testing is essentially a process of making a binary decision about the true parameter space based on a sample.

!!! abstract "Hypothesis Testing Confusion Matrix"
    
    | True State \ Model Prediction | Accept $H_0$ (Do not reject) | Reject $H_0$ (Reject, Discover) |
    | :--- | :--- | :--- |
    | **$H_0$ is True** (Null is True) | ✅ Correct Inference (True Negative, $1-\alpha$) | ❌ **Type I Error** (False Positive, $\alpha$) |
    | **$H_1$ is True** (Alt is True) | ❌ **Type II Error** (False Negative, $\beta$) | ✅ **Statistical Power** (True Positive, Power: $1-\beta$) |

### Core Concepts Quick Reference

!!! info "Type I Error (False Positive)"
    * **Definition**: The null hypothesis $H_0$ is true, but it is incorrectly rejected. 
    * **Mathematical Expression**: $\alpha = P(\text{Reject } H_0 \mid H_0 \text{ is True}) = \frac{FP}{FP+TN}$
    * **Alias**: The **Significance Level** (or size) of the test. In classical statistics, we usually strictly control $\alpha$ first (e.g., set to $0.05$).

!!! info "Type II Error (False Negative)"
    * **Definition**: The alternative hypothesis $H_1$ is true, but we fail to reject $H_0$. 
    * **Mathematical Expression**: $\beta = P(\text{Accept } H_0 \mid H_1 \text{ is True}) = \frac{FN}{TP+FN}$

!!! success "Statistical Power"
    * **Definition**: The probability of correctly rejecting $H_0$ when the alternative hypothesis $H_1$ is true. It represents the test's ability to detect a true effect.
    * **Mathematical Expression**: $\text{Power} = 1 - \beta = P(\text{Reject } H_0 \mid H_1 \text{ is True}) = \frac{TP}{TP+FN}$
    * **Intuition**: Given a fixed $\alpha$, we aim to find a testing method that maximizes the Power (i.e., the Uniformly Most Powerful test, UMP).

!!! info "p-value"
    * **Definition**: The probability of observing the current sample statistic (or a more extreme one) given that the null hypothesis $H_0$ is true.
    * **Pitfall**: A p-value is **absolutely not** the "probability that the null hypothesis is true" (i.e., $p \ne P(H_0 \mid \text{Data})$). It reflects the **degree of inconsistency between the data and the null hypothesis**. The smaller the p-value, the stronger the evidence to reject $H_0$.

---

### Modern Frontiers: Multiple Testing and FDR Control

In modern high-dimensional statistics, we often need to perform thousands or even millions of tests simultaneously. At this point, traditional $\alpha$ control completely fails.

??? note "The Multiple Testing Problem"
    Suppose we independently test $m$ completely ineffective factors (i.e., all $m$ null hypotheses $H_0$ hold), with the significance level of a single test set to $\alpha = 0.05$.
    
    Then the probability of making at least one Type I error (Family-Wise Error Rate, FWER) is:
    
    \[
    \text{FWER} = P(V \ge 1) = 1 - (1 - \alpha)^m
    \]
    
    When $m = 100$, $\text{FWER} \approx 0.994$. This means that as long as enough factors are tested, seemingly significant "false signals" will inevitably be unearthed. This introduces the concept of FDR.

??? success "False Discovery Rate (FDR)"
    Introduced by Benjamini and Hochberg in 1995, it is the core concept of high-dimensional inference.
    
    * **Definition**: The **expected value** of the proportion of incorrect rejections (i.e., Type I errors $V$) among all rejected null hypotheses (i.e., all discoveries claimed to have a significant effect $R$).
    
    \[
    \text{FDR} = E\left[ \frac{V}{\max(R, 1)} \right] = \frac{FP}{TP+FP}
    \]
    
    * **Benjamini-Hochberg (BH) Procedure**:
        1. Sort the p-values of the $m$ hypotheses in ascending order: $p_{(1)} \le p_{(2)} \le \dots \le p_{(m)}$.
        2. Find the largest integer $k$ such that $p_{(k)} \le \frac{k}{m} \alpha$.
        3. Reject the first $k$ null hypotheses (i.e., $H_{(1)}, \dots, H_{(k)}$).
        4. **Theorem**: Under the assumption of independence (or positive dependence), this procedure strictly controls the FDR at $\le \alpha$.