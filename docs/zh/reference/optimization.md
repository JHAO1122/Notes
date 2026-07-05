# 🎯 凸优化 (Convex Optimization) 

本模块涵盖了凸优化理论中的核心概念，重点关注对偶理论（Dual Theory）与 KKT 最优化条件。

---

## 一、 凸优化问题标准型 (Standard Form)

!!! info "标准形式定义"
    一个经典的凸优化问题可以写为如下标准形式：

    \[
    \begin{aligned}
    \min_{x} \quad & f_0(x) \\
    \text{s.t.} \quad & f_i(x) \le 0, \quad i = 1, \dots, m \\
    & h_i(x) = 0, \quad i = 1, \dots, p
    \end{aligned}
    \]

    其中，目标函数 $f_0(x)$ 和不等式约束函数 $f_i(x)$ 必须是**凸函数**；等式约束函数 $h_i(x) = a_i^T x - b_i$ 必须是**仿射函数**（即线性的）。

---

## 二、 拉格朗日对偶理论 (Lagrangian Duality)

!!! info "拉格朗日函数 (Lagrangian)"
    引入构造函数 $L: \mathbb{R}^n \times \mathbb{R}^m \times \mathbb{R}^p \to \mathbb{R}$，将约束条件融合进目标函数：

    \[
    L(x, \lambda, \nu) = f_0(x) + \sum_{i=1}^m \lambda_i f_i(x) + \sum_{i=1}^p \nu_i h_i(x)
    \]

    其中 $\lambda_i \ge 0$ 称为对应不等式约束的拉格朗日乘子，$\nu_i$ 为对应等式约束的拉格朗日乘子。

!!! abstract "拉格朗日对偶函数 (Lagrange Dual Function)"
    对偶函数 $g(\lambda, \nu)$ 定义为拉格朗日函数关于原变量 $x$ 的下确界：

    \[
    g(\lambda, \nu) = \inf_{x \in \mathcal{D}} L(x, \lambda, \nu)
    \]

    * **核心性质**：无论原问题是否为凸问题，对偶函数 $g(\lambda, \nu)$ 永远是一个**凹函数**（因为它是关于 $(\lambda, \nu)$ 的一系列仿射函数的逐点下确界）。

??? success "弱对偶性与强对偶性 (Duality Gap)"
    设原问题的最优解为 $p^*$，对偶问题的最优解为 $m^*$。对偶问题可写为 $\max_{\lambda \ge 0, \nu} g(\lambda, \nu)$。

    * **弱对偶性 (Weak Duality)**：对于任意优化问题，恒有：

    \[
    d^* \le p^*
    \]

    * **强对偶性 (Strong Duality)**：满足 $d^* = p^*$（此时对偶间隙 Duality Gap 为 0）。
    
    * **Slater 条件（强对偶充分条件）**：若原问题是凸问题，且存在一个**严格可行点** $x$ 使得对所有不等式约束均有 $f_i(x) < 0$，则强对偶性成立。

---

## 三、 KKT 最优化条件 (Karush-Kuhn-Tucker Conditions)

!!! success "KKT 条件公式化描述"
    设强对偶性成立，$x^*$ 是原问题的最优解，$(\lambda^*, \nu^*)$ 是对偶问题的最优解。若相关函数均可微，则它们必须满足以下 **KKT 条件**：

    * **1. 定常条件 (Stationarity)**：拉格朗日函数在最优点处的梯度为零：

    \[
    \nabla f_0(x^*) + \sum_{i=1}^m \lambda_i^* \nabla f_i(x^*) + \sum_{i=1}^p \nu_i^* \nabla h_i(x^*) = 0
    \]

    * **2. 互补松弛性 (Complementary Slackness)**：不等式乘子与约束强绑定：

    \[
    \lambda_i^* f_i(x^*) = 0, \quad i = 1, \dots, m
    \]

    *(直觉：若约束不激活 $f_i(x^*) < 0$，则对应的乘子必为零 $\lambda_i^* = 0$；若乘子大于零，则必有 $f_i(x^*) = 0$)*

    * **3. 原问题可行性 (Primal Feasibility)**：最优点必须满足原问题的全部约束：

    \[
    f_i(x^*) \le 0, \quad h_i(x^*) = 0
    \]

    * **4. 对偶问题可行性 (Dual Feasibility)**：不等式约束对应的拉格朗日乘子非负：

    \[
    \lambda_i^* \ge 0
    \]

??? note "KKT 条件的充要性结论"

    * **必要性**：对于任意可微的优化问题，只要强对偶性成立，最优解对 $(x^*, \lambda^*, \nu^*)$ 就**必然**满足 KKT 条件。
    
    * **充分性**：若原问题是**凸问题**，则满足 KKT 条件的任意一组解 $(x^*, \lambda^*, \nu^*)$ **必定**分别是原问题和对偶问题的全局最优解。
