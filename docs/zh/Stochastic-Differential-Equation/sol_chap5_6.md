---
tags:
  - 随机微分方程
  - 课后习题
---

# 📝 课后习题精解：多元 Itô 公式与 SDE 求解实战

!!! abstract "关于本页"
    本页面收录了《随机微分方程》课程第五章（多元 Itô 公式）与第六章（随机微分方程的精确求解）的核心课后习题解答。内容涵盖了多维布朗运动的微积分法则、鞅的构造与验证、Bessel 过程，以及利用积分因子法和待定函数法求解复杂 SDE 的高级实战技巧。

---

## Part I: 第五章 多元 Itô 公式与鞅的验证

### 习题 1
**题目**
设 $W = (W^1, \cdots, W^n)$ 是 $n$ 维 Brown 运动。证明 $Y(t) = |W(t)|^2 - nt$ ($t \ge 0$) 是鞅。

??? success "解答 (点击展开)"

    这道题可以通过多维 Itô 公式直接证明漂移项为 0。
    
    首先，将 $n$ 维布朗运动的模长平方展开为各分量的平方和：
    
    $$
    |W(t)|^2 = \sum_{i=1}^n (W^i(t))^2
    $$
    
    对单个分量 $(W^i(t))^2$ 应用一维 Itô 公式：
    
    $$
    d((W^i)^2) = 2W^i dW^i + \frac{1}{2} \cdot 2 (dW^i)^2 = 2W^i dW^i + dt
    $$
    
    将所有分量求和，根据微分的线性性质：
    
    $$
    d(|W(t)|^2) = \sum_{i=1}^n d((W^i)^2) = \sum_{i=1}^n (2W^i dW^i + dt) = 2\sum_{i=1}^n W^i dW^i + n dt
    $$
    
    现在回到原过程 $Y(t) = |W(t)|^2 - nt$，对其求微分：
    
    $$
    dY(t) = d(|W(t)|^2) - d(nt) = \left( 2\sum_{i=1}^n W^i dW^i + n dt \right) - n dt = 2\sum_{i=1}^n W^i dW^i
    $$
    
    将其写成积分形式：
    
    $$
    Y(t) - Y(0) = 2 \int_0^t \sum_{i=1}^n W^i(s) dW^i(s)
    $$
    
    由于等式右侧只包含关于布朗运动的 Itô 积分，且被积函数 $2W^i(s)$ 在紧区间上是平方可积的，因此该积分过程是一个鞅。故 $Y(t)$ 也是一个鞅。得证。

---

### 习题 2
**题目**
设 $W = (W^1, \cdots, W^n)^T$ 是 $n$ 维 Brown 运动。记 $R = |W|$。证明 $R$ 满足如下随机贝塞尔 (Bessel) 方程

$$dR = \frac{n-1}{2R} dt + \sum_{i=1}^n \frac{W^i}{R} dW^i$$

??? success "解答 (点击展开)"

    这里我们令多元函数 $f(x) = |x| = \left(\sum_{i=1}^n (x^i)^2\right)^{1/2}$。计算其偏导数：
    
    一阶偏导数：
    
    $$
    f_{x^i} = \frac{1}{2} \left(\sum_{j=1}^n (x^j)^2\right)^{-1/2} \cdot 2x^i = \frac{x^i}{|x|} = \frac{x^i}{R}
    $$
    
    二阶偏导数（利用商的求导法则）：
    
    $$
    f_{x^i x^i} = \frac{1 \cdot R - x^i \cdot f_{x^i}}{R^2} = \frac{R - x^i (x^i / R)}{R^2} = \frac{1}{R} - \frac{(x^i)^2}{R^3}
    $$
    
    多维 Itô 公式中需要用到拉普拉斯算子（所有二阶纯偏导之和），对其求和：
    
    $$
    \sum_{i=1}^n f_{x^i x^i} = \sum_{i=1}^n \left( \frac{1}{R} - \frac{(x^i)^2}{R^3} \right) = \frac{n}{R} - \frac{\sum_{i=1}^n (x^i)^2}{R^3} = \frac{n}{R} - \frac{R^2}{R^3} = \frac{n-1}{R}
    $$
    
    将 $X_t = W_t$ 代入 $n$ 维 Itô 公式 $df(W_t) = \sum_{i=1}^n f_{x^i} dW^i + \frac{1}{2} \sum_{i=1}^n f_{x^i x^i} dt$ （注：由于不同维度的布朗运动独立，交叉项 $dW^i dW^j = 0$）：
    
    $$
    dR = d(|W|) = \sum_{i=1}^n \frac{W^i}{R} dW^i + \frac{1}{2} \left( \frac{n-1}{R} \right) dt
    $$
    
    整理后即得高维 Bessel 过程的经典形式：
    
    $$
    dR = \frac{n-1}{2R} dt + \sum_{i=1}^n \frac{W^i}{R} dW^i
    $$
    
    得证。

---

### 习题 3
**题目**
(1) 验证 $X = (\cos W, \sin W)$ 是 $dX^1 = -\frac{1}{2}X^1 dt - X^2 dW, \quad dX^2 = -\frac{1}{2}X^2 dt + X^1 dW$ 的解
(2) 证明 若 $X = (X^1, X^2)$ 是上述方程组的解，则 $|X|$ 是与时间无关的常数

??? success "解答 (点击展开)"

    **(1) 验证解的形式**
    已知 $X^1 = \cos W(t)$，对其应用一维 Itô 公式：
    
    $$
    dX^1 = d(\cos W) = -\sin W dW - \frac{1}{2} \cos W dt = -X^2 dW - \frac{1}{2} X^1 dt
    $$
    
    已知 $X^2 = \sin W(t)$，同样应用 Itô 公式：
    
    $$
    dX^2 = d(\sin W) = \cos W dW - \frac{1}{2} \sin W dt = X^1 dW - \frac{1}{2} X^2 dt
    $$
    
    计算结果与题干给出的 SDE 完全一致，故 $X = (\cos W, \sin W)$ 是方程组的解。

    <br>

    **(2) 证明模长守恒**
    我们考察 $|X|^2 = (X^1)^2 + (X^2)^2$ 的微分。由乘积法则（或 Itô 公式）：
    
    $$
    d((X^1)^2) = 2X^1 dX^1 + \langle dX^1, dX^1 \rangle
    $$
    
    将方程组代入，并注意二次变差 $\langle dX^1, dX^1 \rangle = (-X^2)^2 dt = (X^2)^2 dt$：
    
    $$
    d((X^1)^2) = 2X^1 \left( -\frac{1}{2}X^1 dt - X^2 dW \right) + (X^2)^2 dt = -(X^1)^2 dt - 2X^1 X^2 dW + (X^2)^2 dt
    $$
    
    同理，对于 $X^2$：
    
    $$
    d((X^2)^2) = 2X^2 dX^2 + \langle dX^2, dX^2 \rangle = 2X^2 \left( -\frac{1}{2}X^2 dt + X^1 dW \right) + (X^1)^2 dt
    $$
    
    $$
    d((X^2)^2) = -(X^2)^2 dt + 2X^1 X^2 dW + (X^1)^2 dt
    $$
    
    将两式相加：
    
    $$
    d(|X|^2) = d((X^1)^2) + d((X^2)^2) = \left[ -(X^1)^2 + (X^2)^2 - (X^2)^2 + (X^1)^2 \right] dt + [-2X^1 X^2 + 2X^1 X^2] dW = 0
    $$
    
    由于 $d(|X|^2) = 0$，这意味着模长平方不随时间变化，因此 $|X|$ 是与时间无关的常数。得证。

---

### 习题 4
**题目**
证明 $X(t) = (W(t)+t)\exp\left(-W(t)-\frac{1}{2}t\right)$ 是鞅

??? success "解答 (点击展开)"

    验证鞅的标准化做法是：设 $X(t) = u(W(t), t)$，利用 Itô 公式展开并证明漂移项（$dt$ 项）恒为 0。
    
    令函数 $u(x, t) = (x+t)\exp\left(-x-\frac{1}{2}t\right)$。我们分别计算其偏导数：
    
    对 $t$ 求导：
    
    $$
    u_t = \exp\left(-x-\frac{1}{2}t\right) - \frac{1}{2}(x+t)\exp\left(-x-\frac{1}{2}t\right)
    $$
    
    对 $x$ 求一阶导：
    
    $$
    u_x = \exp\left(-x-\frac{1}{2}t\right) - (x+t)\exp\left(-x-\frac{1}{2}t\right)
    $$
    
    对 $x$ 求二阶导：
    
    $$
    u_{xx} = -\exp\left(-x-\frac{1}{2}t\right) - \left[ \exp\left(-x-\frac{1}{2}t\right) - (x+t)\exp\left(-x-\frac{1}{2}t\right) \right]
    $$
    
    $$
    u_{xx} = -2\exp\left(-x-\frac{1}{2}t\right) + (x+t)\exp\left(-x-\frac{1}{2}t\right)
    $$
    
    现在我们代入 Itô 漂移项公式：$drift = u_t + \frac{1}{2}u_{xx}$
    
    $$
    u_t + \frac{1}{2}u_{xx} = \left[ 1 - \frac{1}{2}(x+t) \right] \exp(\dots) + \frac{1}{2} \left[ -2 + (x+t) \right] \exp(\dots)
    $$
    
    提取公因式 $\exp\left(-x-\frac{1}{2}t\right)$：
    
    $$
    = \left( 1 - \frac{1}{2}(x+t) - 1 + \frac{1}{2}(x+t) \right) \exp\left(-x-\frac{1}{2}t\right) = 0 \cdot \exp(\dots) = 0
    $$
    
    因为漂移项 $u_t + \frac{1}{2}u_{xx} \equiv 0$，这就意味着 $dX(t) = u_x dW(t)$，由于没有 $dt$ 项，该积分过程构成一个鞅。得证。

---

### 习题 5
**题目**
证明 随机微分方程 $dX(t) = \frac{1}{3}X(t)^{1/3} dt + X(t)^{2/3} dW(t)$ 满足初值 $X(0) = x_0 > 0$ 的解是 $X(t) = \left( x_0^{1/3} + \frac{1}{3}W(t) \right)^3$

??? success "解答 (点击展开)"

    这道题是检验“猜解验证”法则的经典题目。
    我们直接将给定的解 $X(t)$ 视为复合函数形式，并通过 Itô 公式求出它的随机微分，看是否与题目 SDE 吻合。
    
    令辅助过程 $Y(t) = x_0^{1/3} + \frac{1}{3}W(t)$，那么原解可写为 $X(t) = Y(t)^3$。
    
    首先，辅助过程 $Y(t)$ 的微分为：
    
    $$
    dY(t) = \frac{1}{3} dW(t)
    $$
    
    其二次变差为：
    
    $$
    (dY(t))^2 = \left( \frac{1}{3} dW(t) \right)^2 = \frac{1}{9} dt
    $$
    
    接下来，利用 Itô 公式对 $f(Y) = Y^3$ 进行展开。计算导数：$f'(Y) = 3Y^2$，$f''(Y) = 6Y$。
    
    $$
    dX(t) = d(Y^3) = 3Y^2 dY(t) + \frac{1}{2}(6Y) (dY(t))^2
    $$
    
    代入 $dY(t)$ 和 $(dY(t))^2$：
    
    $$
    dX(t) = 3Y^2 \left( \frac{1}{3} dW(t) \right) + 3Y \left( \frac{1}{9} dt \right)
    $$
    
    $$
    dX(t) = Y^2 dW(t) + \frac{1}{3}Y dt
    $$
    
    最后，因为 $Y(t) = X(t)^{1/3}$，将其代回上式以消除辅助变量 $Y$：
    
    $$
    dX(t) = (X(t)^{1/3})^2 dW(t) + \frac{1}{3} X(t)^{1/3} dt = X(t)^{2/3} dW(t) + \frac{1}{3} X(t)^{1/3} dt
    $$
    
    这与题目中给定的随机微分方程完全一致！且当 $t=0$ 时，$X(0) = (x_0^{1/3} + 0)^3 = x_0$，初值也满足。因此证明完毕。

---
<br>

## Part II: 第六章 随机微分方程的求解实战

### 习题 1
**题目**
求解随机微分方程
(1) $dX = X dt + e^{-t} dW$
(2) $dX_1 = dt + dW_1, \quad dX_2 = X_1 dW_2$

??? success "解答 (点击展开)"

    **(1) 求解 $dX_t = X_t dt + e^{-t} dW_t$**
    这是一个线性 SDE，我们可以利用积分因子法。移项得到 $dX_t - X_t dt = e^{-t} dW_t$。
    考虑乘上积分因子 $F_t = e^{-t}$，我们考察过程 $Y_t = e^{-t} X_t$ 的微分：
    
    $$
    dY_t = d(e^{-t}X_t) = -e^{-t}X_t dt + e^{-t}dX_t
    $$
    
    代入原方程的 $dX_t$：
    
    $$
    dY_t = -e^{-t}X_t dt + e^{-t}(X_t dt + e^{-t} dW_t) = -e^{-t}X_t dt + e^{-t}X_t dt + e^{-2t} dW_t = e^{-2t} dW_t
    $$
    
    对两边积分：
    
    $$
    Y_t - Y_0 = \int_0^t e^{-2s} dW_s
    $$
    
    代回 $Y_t = e^{-t}X_t$，得到最终解：
    
    $$
    e^{-t}X_t = X_0 + \int_0^t e^{-2s} dW_s \implies X(t) = e^t X(0) + \int_0^t e^{t-2s} dW(s)
    $$

    <br>

    **(2) 求解 $dX_1 = dt + dW_1, \quad dX_2 = X_1 dW_2$**
    这是一个层级 SDE，先解第一层。对第一式直接积分：
    
    $$
    X_1(t) = X_1(0) + \int_0^t ds + \int_0^t dW_1(s) = X_1(0) + t + W_1(t)
    $$
    
    将求出的 $X_1(t)$ 显式表达代入第二式：
    
    $$
    dX_2(t) = (X_1(0) + t + W_1(t)) dW_2(t)
    $$
    
    直接积分得到最终解：
    
    $$
    X_2(t) = X_2(0) + X_1(0) W_2(t) + \int_0^t s dW_2(s) + \int_0^t W_1(s) dW_2(s)
    $$

---

### 习题 2
**题目**
证明 $X(t) = (a\cos W(t), b\sin W(t))^T$ (其中 $a, b$ 是正的常数) 是如下随机微分方程的解

$$dX = -\frac{1}{2}X dt + \begin{pmatrix} 0 & -\frac{a}{b} \\ \frac{b}{a} & 0 \end{pmatrix} X dW$$

??? success "解答 (点击展开)"

    这道题将矩阵形式引入了 SDE。我们将给定的矩阵展开为方程组。
    令矩阵 $M = \begin{pmatrix} 0 & -a/b \\ b/a & 0 \end{pmatrix}$，那么 $MX = \begin{pmatrix} 0 & -a/b \\ b/a & 0 \end{pmatrix} \begin{pmatrix} X^1 \\ X^2 \end{pmatrix} = \begin{pmatrix} -\frac{a}{b}X^2 \\ \frac{b}{a}X^1 \end{pmatrix}$。
    
    因此，原 SDE 等价于如下两个标量方程：
    
    $$
    dX^1 = -\frac{1}{2}X^1 dt - \frac{a}{b}X^2 dW
    $$
    
    $$
    dX^2 = -\frac{1}{2}X^2 dt + \frac{b}{a}X^1 dW
    $$
    
    现在我们验证给定的解 $X^1 = a\cos W, X^2 = b\sin W$。
    对 $X^1$ 应用 Itô 公式：
    
    $$
    dX^1 = d(a\cos W) = -a\sin W dW - \frac{1}{2}a\cos W dt
    $$
    
    将 $X^1$ 和 $X^2$ 的定义代入右侧观察：由于 $X^2 = b\sin W$，即 $\sin W = X^2/b$，因此 $-a\sin W = -a(X^2/b) = -\frac{a}{b}X^2$。
    
    $$
    dX^1 = -\frac{a}{b}X^2 dW - \frac{1}{2}X^1 dt
    $$
    
    这与第一条标量方程完全吻合！
    
    同理验证 $X^2$：
    
    $$
    dX^2 = d(b\sin W) = b\cos W dW - \frac{1}{2}b\sin W dt
    $$
    
    由于 $X^1 = a\cos W$，即 $\cos W = X^1/a$，因此 $b\cos W = \frac{b}{a}X^1$。代入得：
    
    $$
    dX^2 = \frac{b}{a}X^1 dW - \frac{1}{2}X^2 dt
    $$
    
    这与第二条标量方程完全吻合。验证通过。

---

### 习题 3
**题目**
证明 $X_t = e^{W_t}$ 是如下随机微分方程的一个解：

$$dX_t = \frac{1}{2}X_t dt + X_t dW_t$$

??? success "解答 (点击展开)"

    这道题是黑尔斯-斯科尔斯 (Black-Scholes) 几何布朗运动模型的最简退化版本。
    
    令函数 $f(x) = e^x$。将其导数 $f'(x) = e^x, f''(x) = e^x$ 以及 $W_t$ 代入 Itô 公式：
    
    $$
    df(W_t) = f'(W_t)dW_t + \frac{1}{2}f''(W_t)(dW_t)^2
    $$
    
    由于二次变差 $(dW_t)^2 = dt$，代入指数函数：
    
    $$
    d(e^{W_t}) = e^{W_t} dW_t + \frac{1}{2} e^{W_t} dt
    $$
    
    再将 $X_t = e^{W_t}$ 的定义代换回等式右侧的系数：
    
    $$
    dX_t = X_t dW_t + \frac{1}{2} X_t dt
    $$
    
    交换顺序即为 $dX_t = \frac{1}{2}X_t dt + X_t dW_t$，与题目中的方程完美一致。得证。

---

### 习题 4
**题目**
求解随机微分方程

$$dX_t = e^t(1+W_t^2)dt + (1+2e^t W_t)dW_t, \quad X_0 = 0$$

??? success "解答 (点击展开)"

    由于方程的系数显式地包含了时间 $t$ 和布朗运动 $W_t$，我们采用**待定函数法** (Guess and Verify)。
    假设解具有形式 $X_t = f(t, W_t)$。对多元函数 $f$ 展开 Itô 公式：
    
    $$
    dX_t = \left( f_t + \frac{1}{2}f_{ww} \right) dt + f_w dW_t
    $$
    
    将它与题目给定的 SDE 进行逐项系数对比：
    
    **1. 匹配 $dW_t$ 的系数：**
    
    $$
    f_w(t, w) = 1 + 2e^t w
    $$
    
    对 $w$ 偏积分，得到函数 $f$ 的结构：
    
    $$
    f(t, w) = w + e^t w^2 + g(t)
    $$
    
    其中 $g(t)$ 是仅依赖于时间的待定积分常数函数。
    
    **2. 匹配 $dt$ 的系数：**
    先计算假设函数 $f$ 的其余偏导：
    $f_t(t, w) = e^t w^2 + g'(t)$
    $f_{ww}(t, w) = 2e^t$
    
    将它们代入漂移项，要求其必须等于原方程的 $dt$ 系数：
    
    $$
    f_t + \frac{1}{2}f_{ww} = \left( e^t w^2 + g'(t) \right) + \frac{1}{2}(2e^t) = e^t(1+w^2) + g'(t)
    $$
    
    而原方程的 $dt$ 系数为 $e^t(1+w^2)$。令两者相等：
    
    $$
    e^t(1+w^2) + g'(t) = e^t(1+w^2) \implies g'(t) = 0
    $$
    
    说明 $g(t) = C$（常数）。
    
    **3. 代入初始条件：**
    解的形式定为 $X_t = W_t + e^t W_t^2 + C$。
    已知 $X_0 = 0$，且 $W_0 = 0$：
    
    $$
    0 = 0 + e^0(0)^2 + C \implies C = 0
    $$
    
    因此，该随机微分方程的精确解为：
    
    $$
    X(t) = W(t) + e^t W^2(t)
    $$

---

### 习题 5
**题目**
求解随机微分方程

$$dX_t = \frac{1}{X_t}dt + \alpha X_t dW_t, \quad X_0 = x_0 > 0$$

??? success "解答 (点击展开)"

    这是一个非常具有挑战性的非线性 SDE。由于方程包含 $1/X_t$，提示我们先利用平方变换来消去分母。
    
    **第一步：变量代换（线性化）**
    令 $Y_t = X_t^2$。对其应用 Itô 公式：
    
    $$
    dY_t = 2X_t dX_t + (dX_t)^2
    $$
    
    代入原方程 $dX_t = \frac{1}{X_t}dt + \alpha X_t dW_t$，并注意二次变差 $(dX_t)^2 = \alpha^2 X_t^2 dt = \alpha^2 Y_t dt$：
    
    $$
    dY_t = 2X_t \left( \frac{1}{X_t}dt + \alpha X_t dW_t \right) + \alpha^2 Y_t dt = 2 dt + 2\alpha Y_t dW_t + \alpha^2 Y_t dt
    $$
    
    整理后得到关于 $Y_t$ 的线性 SDE：
    
    $$
    dY_t = (2 + \alpha^2 Y_t) dt + 2\alpha Y_t dW_t
    $$
    
    **第二步：构造积分因子解线性 SDE**
    为了消去 $Y_t$ 的比例项，我们构造积分因子 $F_t$ 满足 $dF_t = F_t ( -\alpha^2 dt - 2\alpha dW_t )$。
    根据几何布朗运动的知识，该积分子为：
    
    $$
    F_t = \exp\left( -2\alpha W_t - \frac{1}{2}(-2\alpha)^2 dt - \alpha^2 dt \right) = \exp\left( -2\alpha W_t + \alpha^2 t \right)
    $$
    
    （注：根据 Itô 展开，可以验证 $dF_t = F_t (3\alpha^2 dt - 2\alpha dW_t)$）。
    我们现在计算乘积 $d(Y_t F_t)$ 的微分，利用乘积法则 $d(YF) = Y dF + F dY + \langle dY, dF \rangle$：
    交叉变差项 $\langle dY_t, dF_t \rangle = (2\alpha Y_t)(-2\alpha F_t) dt = -4\alpha^2 Y_t F_t dt$。
    
    展开计算：
    
    $$
    d(Y_t F_t) = Y_t F_t(3\alpha^2 dt - 2\alpha dW_t) + F_t\left[ (2 + \alpha^2 Y_t) dt + 2\alpha Y_t dW_t \right] - 4\alpha^2 Y_t F_t dt
    $$
    
    观察 $dW_t$ 的系数：$Y_t F_t(-2\alpha) + F_t(2\alpha Y_t) = 0$（完美抵消）。
    观察 $dt$ 的系数：$Y_t F_t(3\alpha^2) + F_t(2 + \alpha^2 Y_t) - 4\alpha^2 Y_t F_t = Y_t F_t (3\alpha^2 + \alpha^2 - 4\alpha^2) + 2F_t = 2F_t$（$Y_t$ 项完美抵消）。
    
    因此，极其优雅地化简为了：
    
    $$
    d(Y_t F_t) = 2F_t dt
    $$
    
    **第三步：积分还原**
    对上式从 $0$ 到 $t$ 积分：
    
    $$
    Y_t F_t - Y_0 F_0 = 2 \int_0^t F_s ds
    $$
    
    代入初始条件 $Y_0 = x_0^2$ 且 $F_0 = \exp(0) = 1$，并移项得到 $Y_t$：
    
    $$
    Y_t = F_t^{-1} \left( x_0^2 + 2 \int_0^t F_s ds \right)
    $$
    
    代入 $F_t$ 和 $F_s$ 的显式表达：
    
    $$
    Y_t = \exp(2\alpha W_t - \alpha^2 t) \left( x_0^2 + 2 \int_0^t \exp(-2\alpha W_s + \alpha^2 s) ds \right)
    $$
    
    最后开平方根还原回 $X_t$（由于题目指定 $X_0 > 0$ 且积分内部恒正，故取正根）：
    
    $$
    X(t) = \left[ e^{2\alpha W(t) - \alpha^2 t} \left( x_0^2 + 2 \int_0^t e^{-2\alpha W(s) + \alpha^2 s} ds \right) \right]^{1/2}
    $$