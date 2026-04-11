---
tags:
  - Probability Theory
  - Stochastic Processes
  - Stochastic Differential Equations
  - Exercise Solutions
---

# 📝 Selected Exercise Solutions: Multivariate Itô Formula & SDE Solving

!!! abstract "About This Page"
    This page contains solutions to core exercises from Chapter 5 (Multivariate Itô Formula) and Chapter 6 (Exact Solutions of SDEs) of the *Stochastic Differential Equations* course. The content covers the calculus rules for multi-dimensional Brownian motion, martingale construction, Bessel processes, and advanced techniques such as the integrating factor method for solving complex SDEs.

---

## Part I: Chapter 5 Multivariate Itô Formula & Martingale Verification

### Exercise 1
**Problem**
Let $W = (W^1, \dots, W^n)$ be an $n$-dimensional Brownian motion. Prove that $Y(t) = |W(t)|^2 - nt$ ($t \ge 0$) is a martingale.

??? success "Solution (Click to expand)"

    This problem can be solved directly by applying the multivariate Itô formula to prove that the drift term is zero.
    
    First, expand the squared norm of the $n$-dimensional Brownian motion as the sum of squares of its components:
    
    $$
    |W(t)|^2 = \sum_{i=1}^n (W^i(t))^2
    $$
    
    Apply the one-dimensional Itô formula to a single component $(W^i(t))^2$:
    
    $$
    d((W^i)^2) = 2W^i dW^i + \frac{1}{2} \cdot 2 (dW^i)^2 = 2W^i dW^i + dt
    $$
    
    Summing over all components, according to the linearity of the differential:
    
    $$
    d(|W(t)|^2) = \sum_{i=1}^n d((W^i)^2) = \sum_{i=1}^n (2W^i dW^i + dt) = 2\sum_{i=1}^n W^i dW^i + n dt
    $$
    
    Now return to the original process $Y(t) = |W(t)|^2 - nt$, and differentiate it:
    
    $$
    dY(t) = d(|W(t)|^2) - d(nt) = \left( 2\sum_{i=1}^n W^i dW^i + n dt \right) - n dt = 2\sum_{i=1}^n W^i dW^i
    $$
    
    Writing this in integral form:
    
    $$
    Y(t) - Y(0) = 2 \int_0^t \sum_{i=1}^n W^i(s) dW^i(s)
    $$
    
    Since the right-hand side contains only Itô integrals with respect to Brownian motion, and the integrand $2W^i(s)$ is square-integrable over compact intervals, this integral process is a martingale. Thus, $Y(t)$ is also a martingale. The proof is complete.

---

### Exercise 2
**Problem**
Let $W = (W^1, \dots, W^n)^T$ be an $n$-dimensional Brownian motion. Denote $R = |W|$. Prove that $R$ satisfies the following Stochastic Bessel Equation:
$$dR = \frac{n-1}{2R} dt + \sum_{i=1}^n \frac{W^i}{R} dW^i$$

??? success "Solution (Click to expand)"

    Here we let the multivariate function be $f(x) = |x| = \left(\sum_{i=1}^n (x^i)^2\right)^{1/2}$. Calculate its partial derivatives:
    
    First-order partial derivative:
    
    $$
    f_{x^i} = \frac{1}{2} \left(\sum_{j=1}^n (x^j)^2\right)^{-1/2} \cdot 2x^i = \frac{x^i}{|x|} = \frac{x^i}{R}
    $$
    
    Second-order partial derivative (using the quotient rule):
    
    $$
    f_{x^i x^i} = \frac{1 \cdot R - x^i \cdot f_{x^i}}{R^2} = \frac{R - x^i (x^i / R)}{R^2} = \frac{1}{R} - \frac{(x^i)^2}{R^3}
    $$
    
    The multi-dimensional Itô formula requires the Laplacian (the sum of all second-order pure partial derivatives). Summing them up:
    
    $$
    \sum_{i=1}^n f_{x^i x^i} = \sum_{i=1}^n \left( \frac{1}{R} - \frac{(x^i)^2}{R^3} \right) = \frac{n}{R} - \frac{\sum_{i=1}^n (x^i)^2}{R^3} = \frac{n}{R} - \frac{R^2}{R^3} = \frac{n-1}{R}
    $$
    
    Substitute $X_t = W_t$ into the $n$-dimensional Itô formula $df(W_t) = \sum_{i=1}^n f_{x^i} dW^i + \frac{1}{2} \sum_{i=1}^n f_{x^i x^i} dt$ (Note: since Brownian motions in different dimensions are independent, the cross terms $dW^i dW^j = 0$):
    
    $$
    dR = d(|W|) = \sum_{i=1}^n \frac{W^i}{R} dW^i + \frac{1}{2} \left( \frac{n-1}{R} \right) dt
    $$
    
    Rearranging this gives the classical form of the high-dimensional Bessel process:
    
    $$
    dR = \frac{n-1}{2R} dt + \sum_{i=1}^n \frac{W^i}{R} dW^i
    $$
    
    The proof is complete.

---

### Exercise 3
**Problem**
(1) Verify that $X = (\cos W, \sin W)$ is a solution to $dX^1 = -\frac{1}{2}X^1 dt - X^2 dW, \quad dX^2 = -\frac{1}{2}X^2 dt + X^1 dW$.
(2) Prove that if $X = (X^1, X^2)$ is a solution to the above system, then $|X|$ is a constant independent of time.

??? success "Solution (Click to expand)"

    **(1) Verification of the Solution**
    Given $X^1 = \cos W(t)$, applying the one-dimensional Itô formula to it:
    
    $$
    dX^1 = d(\cos W) = -\sin W dW - \frac{1}{2} \cos W dt = -X^2 dW - \frac{1}{2} X^1 dt
    $$
    
    Given $X^2 = \sin W(t)$, similarly applying Itô's formula:
    
    $$
    dX^2 = d(\sin W) = \cos W dW - \frac{1}{2} \sin W dt = X^1 dW - \frac{1}{2} X^2 dt
    $$
    
    The calculated results match the SDEs given in the problem exactly, thus $X = (\cos W, \sin W)$ is a solution to the system.

    <br>

    **(2) Proof of Conservation of Norm**
    We examine the differential of $|X|^2 = (X^1)^2 + (X^2)^2$. By the product rule (or Itô's formula):
    
    $$
    d((X^1)^2) = 2X^1 dX^1 + \langle dX^1, dX^1 \rangle
    $$
    
    Substitute the system equations and note the quadratic variation $\langle dX^1, dX^1 \rangle = (-X^2)^2 dt = (X^2)^2 dt$:
    
    $$
    d((X^1)^2) = 2X^1 \left( -\frac{1}{2}X^1 dt - X^2 dW \right) + (X^2)^2 dt = -(X^1)^2 dt - 2X^1 X^2 dW + (X^2)^2 dt
    $$
    
    Similarly for $X^2$:
    
    $$
    d((X^2)^2) = 2X^2 dX^2 + \langle dX^2, dX^2 \rangle = 2X^2 \left( -\frac{1}{2}X^2 dt + X^1 dW \right) + (X^1)^2 dt
    $$
    
    $$
    d((X^2)^2) = -(X^2)^2 dt + 2X^1 X^2 dW + (X^1)^2 dt
    $$
    
    Summing the two equations:
    
    $$
    d(|X|^2) = d((X^1)^2) + d((X^2)^2) = \left[ -(X^1)^2 + (X^2)^2 - (X^2)^2 + (X^1)^2 \right] dt + [-2X^1 X^2 + 2X^1 X^2] dW = 0
    $$
    
    Since $d(|X|^2) = 0$, this means the squared norm does not change with time. Therefore, $|X|$ is a constant independent of time. The proof is complete.

---

### Exercise 4
**Problem**
Prove that $X(t) = (W(t)+t)\exp\left(-W(t)-\frac{1}{2}t\right)$ is a martingale.

??? success "Solution (Click to expand)"

    The standard approach for verifying a martingale is: let $X(t) = u(W(t), t)$, expand it using Itô's formula, and prove that the drift term (the $dt$ term) is identically 0.
    
    Let the function $u(x, t) = (x+t)\exp\left(-x-\frac{1}{2}t\right)$. We calculate its partial derivatives respectively:
    
    Derivative with respect to $t$:
    
    $$
    u_t = \exp\left(-x-\frac{1}{2}t\right) - \frac{1}{2}(x+t)\exp\left(-x-\frac{1}{2}t\right)
    $$
    
    First-order derivative with respect to $x$:
    
    $$
    u_x = \exp\left(-x-\frac{1}{2}t\right) - (x+t)\exp\left(-x-\frac{1}{2}t\right)
    $$
    
    Second-order derivative with respect to $x$:
    
    $$
    u_{xx} = -\exp\left(-x-\frac{1}{2}t\right) - \left[ \exp\left(-x-\frac{1}{2}t\right) - (x+t)\exp\left(-x-\frac{1}{2}t\right) \right]
    $$
    
    $$
    u_{xx} = -2\exp\left(-x-\frac{1}{2}t\right) + (x+t)\exp\left(-x-\frac{1}{2}t\right)
    $$
    
    Now we substitute into the Itô drift term formula: $\text{drift} = u_t + \frac{1}{2}u_{xx}$
    
    $$
    u_t + \frac{1}{2}u_{xx} = \left[ 1 - \frac{1}{2}(x+t) \right] \exp(\dots) + \frac{1}{2} \left[ -2 + (x+t) \right] \exp(\dots)
    $$
    
    Extracting the common factor $\exp\left(-x-\frac{1}{2}t\right)$:
    
    $$
    = \left( 1 - \frac{1}{2}(x+t) - 1 + \frac{1}{2}(x+t) \right) \exp\left(-x-\frac{1}{2}t\right) = 0 \cdot \exp(\dots) = 0
    $$
    
    Since the drift term $u_t + \frac{1}{2}u_{xx} \equiv 0$, this implies that $dX(t) = u_x dW(t)$. Because there is no $dt$ term, this integral process constitutes a martingale. The proof is complete.

---

### Exercise 5
**Problem**
Prove that the solution to the stochastic differential equation $dX(t) = \frac{1}{3}X(t)^{1/3} dt + X(t)^{2/3} dW(t)$ with initial value $X(0) = x_0 > 0$ is $X(t) = \left( x_0^{1/3} + \frac{1}{3}W(t) \right)^3$.

??? success "Solution (Click to expand)"

    This is a classic problem to test the "guess and verify" rule.
    We directly treat the given solution $X(t)$ as a composite function and find its stochastic differential via Itô's formula to see if it matches the SDE in the problem.
    
    Let the auxiliary process be $Y(t) = x_0^{1/3} + \frac{1}{3}W(t)$, so the original solution can be written as $X(t) = Y(t)^3$.
    
    First, the differential of the auxiliary process $Y(t)$ is:
    
    $$
    dY(t) = \frac{1}{3} dW(t)
    $$
    
    Its quadratic variation is:
    
    $$
    (dY(t))^2 = \left( \frac{1}{3} dW(t) \right)^2 = \frac{1}{9} dt
    $$
    
    Next, use Itô's formula to expand $f(Y) = Y^3$. Calculate the derivatives: $f'(Y) = 3Y^2$, $f''(Y) = 6Y$.
    
    $$
    dX(t) = d(Y^3) = 3Y^2 dY(t) + \frac{1}{2}(6Y) (dY(t))^2
    $$
    
    Substitute $dY(t)$ and $(dY(t))^2$:
    
    $$
    dX(t) = 3Y^2 \left( \frac{1}{3} dW(t) \right) + 3Y \left( \frac{1}{9} dt \right)
    $$
    
    $$
    dX(t) = Y^2 dW(t) + \frac{1}{3}Y dt
    $$
    
    Finally, since $Y(t) = X(t)^{1/3}$, substitute this back to eliminate the auxiliary variable $Y$:
    
    $$
    dX(t) = (X(t)^{1/3})^2 dW(t) + \frac{1}{3} X(t)^{1/3} dt = X(t)^{2/3} dW(t) + \frac{1}{3} X(t)^{1/3} dt
    $$
    
    This is completely identical to the SDE given in the problem! Furthermore, when $t=0$, $X(0) = (x_0^{1/3} + 0)^3 = x_0$, so the initial condition is also satisfied. Therefore, the proof is complete.

---
<br>

## Part II: Chapter 6 Solving Stochastic Differential Equations

### Exercise 1
**Problem**
Solve the stochastic differential equations:
(1) $dX = X dt + e^{-t} dW$
(2) $dX_1 = dt + dW_1, \quad dX_2 = X_1 dW_2$

??? success "Solution (Click to expand)"

    **(1) Solve $dX_t = X_t dt + e^{-t} dW_t$**
    This is a linear SDE, so we can use the integrating factor method. Rearranging gives $dX_t - X_t dt = e^{-t} dW_t$.
    Consider multiplying by the integrating factor $F_t = e^{-t}$. We examine the differential of the process $Y_t = e^{-t} X_t$:
    
    $$
    dY_t = d(e^{-t}X_t) = -e^{-t}X_t dt + e^{-t}dX_t
    $$
    
    Substitute $dX_t$ from the original equation:
    
    $$
    dY_t = -e^{-t}X_t dt + e^{-t}(X_t dt + e^{-t} dW_t) = -e^{-t}X_t dt + e^{-t}X_t dt + e^{-2t} dW_t = e^{-2t} dW_t
    $$
    
    Integrate both sides:
    
    $$
    Y_t - Y_0 = \int_0^t e^{-2s} dW_s
    $$
    
    Substitute back $Y_t = e^{-t}X_t$ to get the final solution:
    
    $$
    e^{-t}X_t = X_0 + \int_0^t e^{-2s} dW_s \implies X(t) = e^t X(0) + \int_0^t e^{t-2s} dW(s)
    $$

    <br>

    **(2) Solve $dX_1 = dt + dW_1, \quad dX_2 = X_1 dW_2$**
    This is a hierarchical SDE; solve the first layer first. Integrate the first equation directly:
    
    $$
    X_1(t) = X_1(0) + \int_0^t ds + \int_0^t dW_1(s) = X_1(0) + t + W_1(t)
    $$
    
    Substitute the explicit expression found for $X_1(t)$ into the second equation:
    
    $$
    dX_2(t) = (X_1(0) + t + W_1(t)) dW_2(t)
    $$
    
    Integrate directly to get the final solution:
    
    $$
    X_2(t) = X_2(0) + X_1(0) W_2(t) + \int_0^t s dW_2(s) + \int_0^t W_1(s) dW_2(s)
    $$

---

### Exercise 2
**Problem**
Show that $X(t) = (a\cos W(t), b\sin W(t))^T$ (where $a, b$ are positive constants) is a solution to the following stochastic differential equation:
$$dX = -\frac{1}{2}X dt + \begin{pmatrix} 0 & -\frac{a}{b} \\ \frac{b}{a} & 0 \end{pmatrix} X dW$$

??? success "Solution (Click to expand)"

    This problem introduces a matrix form to SDEs. We expand the given matrix into a system of equations.
    Let the matrix $M = \begin{pmatrix} 0 & -a/b \\ b/a & 0 \end{pmatrix}$, then $MX = \begin{pmatrix} 0 & -a/b \\ b/a & 0 \end{pmatrix} \begin{pmatrix} X^1 \\ X^2 \end{pmatrix} = \begin{pmatrix} -\frac{a}{b}X^2 \\ \frac{b}{a}X^1 \end{pmatrix}$.
    
    Therefore, the original SDE is equivalent to the following two scalar equations:
    
    $$
    dX^1 = -\frac{1}{2}X^1 dt - \frac{a}{b}X^2 dW
    $$
    
    $$
    dX^2 = -\frac{1}{2}X^2 dt + \frac{b}{a}X^1 dW
    $$
    
    Now we verify the given solution $X^1 = a\cos W, X^2 = b\sin W$.
    Apply Itô's formula to $X^1$:
    
    $$
    dX^1 = d(a\cos W) = -a\sin W dW - \frac{1}{2}a\cos W dt
    $$
    
    Substitute the definitions of $X^1$ and $X^2$ into the right side: since $X^2 = b\sin W$, we have $\sin W = X^2/b$, thus $-a\sin W = -a(X^2/b) = -\frac{a}{b}X^2$.
    
    $$
    dX^1 = -\frac{a}{b}X^2 dW - \frac{1}{2}X^1 dt
    $$
    
    This matches the first scalar equation perfectly!
    
    Verify $X^2$ similarly:
    
    $$
    dX^2 = d(b\sin W) = b\cos W dW - \frac{1}{2}b\sin W dt
    $$
    
    Since $X^1 = a\cos W$, we have $\cos W = X^1/a$, thus $b\cos W = \frac{b}{a}X^1$. Substituting this gives:
    
    $$
    dX^2 = \frac{b}{a}X^1 dW - \frac{1}{2}X^2 dt
    $$
    
    This matches the second scalar equation perfectly. The verification is complete.

---

### Exercise 3
**Problem**
Show that $X_t = e^{W_t}$ is a solution to the following stochastic differential equation:
$$dX_t = \frac{1}{2}X_t dt + X_t dW_t$$

??? success "Solution (Click to expand)"

    This problem is the simplest degenerate version of the Black-Scholes geometric Brownian motion model.
    
    Let the function $f(x) = e^x$. Substitute its derivatives $f'(x) = e^x, f''(x) = e^x$ and $W_t$ into Itô's formula:
    
    $$
    df(W_t) = f'(W_t)dW_t + \frac{1}{2}f''(W_t)(dW_t)^2
    $$
    
    Since the quadratic variation $(dW_t)^2 = dt$, substitute it into the exponential function:
    
    $$
    d(e^{W_t}) = e^{W_t} dW_t + \frac{1}{2} e^{W_t} dt
    $$
    
    Then substitute the definition $X_t = e^{W_t}$ back into the coefficients on the right side:
    
    $$
    dX_t = X_t dW_t + \frac{1}{2} X_t dt
    $$
    
    Swapping the order gives $dX_t = \frac{1}{2}X_t dt + X_t dW_t$, which perfectly matches the equation in the problem. The proof is complete.

---

### Exercise 4
**Problem**
Solve the stochastic differential equation:
$$dX_t = e^t(1+W_t^2)dt + (1+2e^t W_t)dW_t, \quad X_0 = 0$$

??? success "Solution (Click to expand)"

    Since the coefficients of the equation explicitly contain time $t$ and Brownian motion $W_t$, we use the **Guess and Verify** (method of undetermined functions).
    Assume the solution has the form $X_t = f(t, W_t)$. Expand the multi-variable function $f$ using Itô's formula:
    
    $$
    dX_t = \left( f_t + \frac{1}{2}f_{ww} \right) dt + f_w dW_t
    $$
    
    Compare this term by term with the given SDE:
    
    **1. Match the coefficient of $dW_t$:**
    
    $$
    f_w(t, w) = 1 + 2e^t w
    $$
    
    Integrate partially with respect to $w$ to get the structure of the function $f$:
    
    $$
    f(t, w) = w + e^t w^2 + g(t)
    $$
    
    where $g(t)$ is an undetermined integration constant function depending only on time.
    
    **2. Match the coefficient of $dt$:**
    First calculate the remaining partial derivatives of the assumed function $f$:
    $f_t(t, w) = e^t w^2 + g'(t)$
    $f_{ww}(t, w) = 2e^t$
    
    Substitute them into the drift term, requiring it to equal the $dt$ coefficient of the original equation:
    
    $$
    f_t + \frac{1}{2}f_{ww} = \left( e^t w^2 + g'(t) \right) + \frac{1}{2}(2e^t) = e^t(1+w^2) + g'(t)
    $$
    
    And the $dt$ coefficient of the original equation is $e^t(1+w^2)$. Setting them equal:
    
    $$
    e^t(1+w^2) + g'(t) = e^t(1+w^2) \implies g'(t) = 0
    $$
    
    This means $g(t) = C$ (a constant).
    
    **3. Substitute the initial condition:**
    The form of the solution is determined as $X_t = W_t + e^t W_t^2 + C$.
    Given $X_0 = 0$, and $W_0 = 0$:
    
    $$
    0 = 0 + e^0(0)^2 + C \implies C = 0
    $$
    
    Therefore, the exact solution to this stochastic differential equation is:
    
    $$
    X(t) = W(t) + e^t W^2(t)
    $$

---

### Exercise 5
**Problem**
Solve the stochastic differential equation:
$$dX_t = \frac{1}{X_t}dt + \alpha X_t dW_t, \quad X_0 = x_0 > 0$$

??? success "Solution (Click to expand)"

    This is a highly challenging nonlinear SDE. Since the equation contains $1/X_t$, it suggests that we first use a square transformation to eliminate the denominator.
    
    **Step 1: Variable substitution (Linearization)**
    Let $Y_t = X_t^2$. Apply Itô's formula to it:
    
    $$
    dY_t = 2X_t dX_t + (dX_t)^2
    $$
    
    Substitute the original equation $dX_t = \frac{1}{X_t}dt + \alpha X_t dW_t$, and note the quadratic variation $(dX_t)^2 = \alpha^2 X_t^2 dt = \alpha^2 Y_t dt$:
    
    $$
    dY_t = 2X_t \left( \frac{1}{X_t}dt + \alpha X_t dW_t \right) + \alpha^2 Y_t dt = 2 dt + 2\alpha Y_t dW_t + \alpha^2 Y_t dt
    $$
    
    Rearranging gives a linear SDE with respect to $Y_t$:
    
    $$
    dY_t = (2 + \alpha^2 Y_t) dt + 2\alpha Y_t dW_t
    $$
    
    **Step 2: Construct an integrating factor to solve the linear SDE**
    To eliminate the proportional term of $Y_t$, we construct an integrating factor $F_t$ satisfying $dF_t = F_t ( -\alpha^2 dt - 2\alpha dW_t )$.
    From the knowledge of geometric Brownian motion, this integrating factor is:
    
    $$
    F_t = \exp\left( -2\alpha W_t - \frac{1}{2}(-2\alpha)^2 dt - \alpha^2 dt \right) = \exp\left( -2\alpha W_t + \alpha^2 t \right)
    $$
    
    (Note: According to Itô expansion, it can be verified that $dF_t = F_t (3\alpha^2 dt - 2\alpha dW_t)$).
    We now calculate the differential of the product $d(Y_t F_t)$, using the product rule $d(YF) = Y dF + F dY + \langle dY, dF \rangle$:
    The cross-variation term is $\langle dY_t, dF_t \rangle = (2\alpha Y_t)(-2\alpha F_t) dt = -4\alpha^2 Y_t F_t dt$.
    
    Expanding the calculation:
    
    $$
    d(Y_t F_t) = Y_t F_t(3\alpha^2 dt - 2\alpha dW_t) + F_t\left[ (2 + \alpha^2 Y_t) dt + 2\alpha Y_t dW_t \right] - 4\alpha^2 Y_t F_t dt
    $$
    
    Observe the coefficient of $dW_t$: $Y_t F_t(-2\alpha) + F_t(2\alpha Y_t) = 0$ (perfect cancellation).
    Observe the coefficient of $dt$: $Y_t F_t(3\alpha^2) + F_t(2 + \alpha^2 Y_t) - 4\alpha^2 Y_t F_t = Y_t F_t (3\alpha^2 + \alpha^2 - 4\alpha^2) + 2F_t = 2F_t$ (the $Y_t$ term perfectly cancels out).
    
    Thus, it elegantly simplifies to:
    
    $$
    d(Y_t F_t) = 2F_t dt
    $$
    
    **Step 3: Integration and recovery**
    Integrate the above equation from $0$ to $t$:
    
    $$
    Y_t F_t - Y_0 F_0 = 2 \int_0^t F_s ds
    $$
    
    Substitute the initial conditions $Y_0 = x_0^2$ and $F_0 = \exp(0) = 1$, and rearrange to get $Y_t$:
    
    $$
    Y_t = F_t^{-1} \left( x_0^2 + 2 \int_0^t F_s ds \right)
    $$
    
    Substitute the explicit expressions for $F_t$ and $F_s$:
    
    $$
    Y_t = \exp(2\alpha W_t - \alpha^2 t) \left( x_0^2 + 2 \int_0^t \exp(-2\alpha W_s + \alpha^2 s) ds \right)
    $$
    
    Finally, take the square root to recover $X_t$ (since the problem specifies $X_0 > 0$ and the interior of the integral is always positive, we take the positive root):
    
    $$
    X(t) = \left[ e^{2\alpha W(t) - \alpha^2 t} \left( x_0^2 + 2 \int_0^t e^{-2\alpha W(s) + \alpha^2 s} ds \right) \right]^{1/2}
    $$