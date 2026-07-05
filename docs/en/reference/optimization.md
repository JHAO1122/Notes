# 🎯 Convex Optimization

This module covers the core concepts in convex optimization theory, with a focus on Duality Theory and KKT Optimality Conditions.

---

## I. Standard Form of a Convex Optimization Problem

!!! info "Definition of Standard Form"
    A classic convex optimization problem can be written in the following standard form:

    \[
    \begin{aligned}
    \min_{x} \quad & f_0(x) \\
    \text{s.t.} \quad & f_i(x) \le 0, \quad i = 1, \dots, m \\
    & h_i(x) = 0, \quad i = 1, \dots, p
    \end{aligned}
    \]

    where the objective function \(f_0(x)\) and inequality constraint functions \(f_i(x)\) must be **convex functions**; the equality constraint functions \(h_i(x) = a_i^T x - b_i\) must be **affine functions** (i.e., linear).

---

## II. Lagrangian Duality

!!! info "Lagrangian Function"
    Construct the function \(L: \mathbb{R}^n \times \mathbb{R}^m \times \mathbb{R}^p \to \mathbb{R}\) that incorporates the constraints into the objective:

    \[
    L(x, \lambda, \nu) = f_0(x) + \sum_{i=1}^m \lambda_i f_i(x) + \sum_{i=1}^p \nu_i h_i(x)
    \]

    where \(\lambda_i \ge 0\) are the Lagrange multipliers for the inequality constraints, and \(\nu_i\) are the Lagrange multipliers for the equality constraints.

!!! abstract "Lagrange Dual Function"
    The dual function \(g(\lambda, \nu)\) is defined as the infimum of the Lagrangian over the primal variable \(x\):

    \[
    g(\lambda, \nu) = \inf_{x \in \mathcal{D}} L(x, \lambda, \nu)
    \]

    * **Key property**: Regardless of whether the primal problem is convex, the dual function \(g(\lambda, \nu)\) is always a **concave function** (since it is the pointwise infimum of a family of affine functions in \((\lambda, \nu)\)).

??? success "Weak Duality and Strong Duality (Duality Gap)"
    Let \(p^*\) be the optimal value of the primal problem, and \(d^*\) the optimal value of the dual problem. The dual problem can be written as \(\max_{\lambda \ge 0, \nu} g(\lambda, \nu)\).

    * **Weak Duality**: For any optimization problem, we always have:

    \[
    d^* \le p^*
    \]

    * **Strong Duality**: We have \(d^* = p^*\) (i.e., the duality gap is zero).
  
    * **Slater's condition (sufficient condition for strong duality)**: If the primal problem is convex and there exists a **strictly feasible point** \(x\) such that \(f_i(x) < 0\) for all inequality constraints, then strong duality holds.

---

## III. KKT Optimality Conditions (Karush–Kuhn–Tucker Conditions)

!!! success "Formal Description of KKT Conditions"
    Assume strong duality holds. Let \(x^*\) be a primal optimal point, and \((\lambda^*, \nu^*)\) a dual optimal point. If the relevant functions are differentiable, then they must satisfy the following **KKT conditions**:

    * **1. Stationarity**: The gradient of the Lagrangian with respect to \(x\) vanishes at the optimum:

    \[
    \nabla f_0(x^*) + \sum_{i=1}^m \lambda_i^* \nabla f_i(x^*) + \sum_{i=1}^p \nu_i^* \nabla h_i(x^*) = 0
    \]

    * **2. Complementary Slackness**: The multipliers and inequality constraints are tightly coupled:

    \[
    \lambda_i^* f_i(x^*) = 0, \quad i = 1, \dots, m
    \]

    *(Intuition: if a constraint is inactive, i.e., \(f_i(x^*) < 0\), then its multiplier must be zero \(\lambda_i^* = 0\); if a multiplier is positive, then we must have \(f_i(x^*) = 0\).)*

    * **3. Primal Feasibility**: The optimal point must satisfy all primal constraints:

    \[
    f_i(x^*) \le 0, \quad h_i(x^*) = 0
    \]

    * **4. Dual Feasibility**: The Lagrange multipliers for inequality constraints must be non‑negative:

    \[
    \lambda_i^* \ge 0
    \]

??? note "Necessity and Sufficiency of KKT Conditions"
    * **Necessity**: For any differentiable optimization problem, if strong duality holds, then the optimal primal‑dual pair \((x^*, \lambda^*, \nu^*)\) **must** satisfy the KKT conditions.
  
    * **Sufficiency**: If the primal problem is **convex**, then any triple \((x^*, \lambda^*, \nu^*)\) that satisfies the KKT conditions is **automatically** a globally optimal primal‑dual pair.