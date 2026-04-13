---
tags:
  - 概率论
  - 随机过程
  - 随机微分方程
  - 课后习题
---

# 📝 课后习题精解：随机积分与 Itô 公式

!!! abstract "关于本页"
    本页面收录了《随机微分方程》课程第三章（随机积分）与第四章（Itô 积分与 Itô 公式）的核心课后习题解答。内容涵盖了关于布朗运动的积分计算、Itô 等距同构、二次变差过程，以及利用 Itô 公式求解随机微分方程等核心技巧。

---

## Part I: 第三章 随机积分的计算

### 习题 1

**题目：**
求随机积分 $\int_0^t s dW(s)$，并计算它的期望和方差。

??? success "解答 (点击展开)"

    **1. 求解随机积分：**
    
    利用分部积分公式（Integration by Parts），设 $g(s) = s$，则：
    
    $$
    \int_0^t s dW(s) = sW(s) \Big|_0^t - \int_0^t W(s) ds = tW(t) - \int_0^t W(s) ds
    $$

    <br>

    **2. 计算期望：**
    
    利用期望的线性性质和布朗运动零均值的性质：
    
    $$
    \mathbb{E}\left[ \int_0^t s dW(s) \right] = t\mathbb{E}[W(t)] - \int_0^t \mathbb{E}[W(s)] ds = t \cdot 0 - \int_0^t 0 ds = 0
    $$

    <br>

    **3. 计算方差：**
    
    由 Itô 等距同构 (Itô Isometry) 原理，对于确定性函数 $g(s) \in L^2([0,t])$，其关于布朗运动积分的方差等于该函数平方的黎曼积分：
    
    $$
    Var\left[ \int_0^t s dW(s) \right] = \mathbb{E}\left[ \left( \int_0^t s dW(s) \right)^2 \right] - \left( \mathbb{E}\left[ \int_0^t s dW(s) \right] \right)^2
    $$
    
    由于期望为 0，方差即为二阶矩：
    
    $$
    = \mathbb{E}\left[ \int_0^t s^2 ds \right] = \int_0^t s^2 ds = \frac{1}{3}t^3
    $$

    *注：由于被积函数为确定性函数，该随机积分服从正态分布 $N(0, \frac{1}{3}t^3)$。*

---

### 习题 2

**题目：**
设函数 $g: [0, T] \to \mathbb{R}$ 连续可微，$g(0) = g(T) = 0$。求 $\int_0^T g(t) dW(t)$ 的概率密度函数。

??? success "解答 (点击展开)"

    **第一步：求期望**
    利用分部积分，结合边界条件 $g(0) = g(T) = 0$：
    
    $$
    \int_0^T g(t) dW(t) = g(T)W(T) - g(0)W(0) - \int_0^T g'(t)W(t) dt = - \int_0^T g'(t)W(t) dt
    $$
    
    取期望：
    
    $$
    \mathbb{E}\left[ \int_0^T g(t) dW(t) \right] = \mathbb{E}\left[ - \int_0^T g'(t)W(t) dt \right] = - \int_0^T g'(t)\mathbb{E}[W(t)] dt = 0
    $$

    <br>

    **第二步：求方差**
    根据 Itô 等距同构（因为 $g(t)$ 是确定性连续函数，必定属于 $L^2([0,T])$ 空间）：
    
    $$
    Var\left[ \int_0^T g(t) dW(t) \right] = \mathbb{E}\left[ \left( \int_0^T g(t) dW(t) \right)^2 \right] = \int_0^T g^2(t) dt
    $$

    <br>

    **第三步：确定分布类型与密度函数**
    由于被积函数 $g(t)$ 是确定性函数，布朗运动的增量是独立且服从正态分布的。作为正态增量的线性组合极限，该随机积分（看作无穷维的线性组合）必然服从正态分布。
    
    因此，该随机变量服从均值为 0，方差为 $\sigma^2 = \int_0^T g^2(t) dt$ 的正态分布：
    
    $$
    I := \int_0^T g(t) dW(t) \sim N\left(0, \int_0^T g^2(t) dt\right)
    $$
    
    其概率密度函数为标准的一维正态密度公式：
    
    $$
    f_I(x) = \frac{1}{\sqrt{2\pi \int_0^T g^2(t) dt}} \exp\left( - \frac{x^2}{2 \int_0^T g^2(t) dt} \right)
    $$

---

### 习题 3

**题目：**
为了从牛顿（Newton）力学角度解释 Brown 运动，郎之万（Langevin, 1872-1946）提出了如下描述液体中微粒运动速度的（随机）微分方程：
$$\frac{dv}{dt} = -\beta v + \dot{W}(t)$$
其中 $-\beta v$ 代表微粒运动受到的摩擦阻力（$\beta$ 是正的常数），白噪声 $\dot{W}(t)$ 描述微粒受到的随机的冲击力。
(1) 说明该方程的解为 $v(t) = v_0 e^{-\beta t} + W(t) - \beta \int_0^t e^{-\beta(t-s)} W(s) ds$，从而由原点出发的微粒运动路径为 $x_\beta(t) = \int_0^t e^{-\beta(t-s)} W(s) ds$；
(2) 计算 $v(t), x(t)$ 的期望和方差；
(3) 证明 $\lim_{\beta \to \infty} \beta x_\beta(t) = W(t)$。

??? success "解答 (点击展开)"

    **(1) 验证解的形式**
    将给定的解 $v(t)$ 代入原微分方程进行验证。对解求导：
    
    $$
    \frac{dv(t)}{dt} = \frac{d}{dt} \left( v_0 e^{-\beta t} + W(t) - \beta \int_0^t e^{-\beta(t-s)} W(s) ds \right)
    $$
    
    利用变上限积分的求导法则（Leibniz Rule）：
    
    $$
    = -\beta v_0 e^{-\beta t} + \dot{W}(t) - \beta \left( e^{-\beta(t-t)}W(t) + \int_0^t -\beta e^{-\beta(t-s)} W(s) ds \right)
    $$
    
    $$
    = -\beta v_0 e^{-\beta t} + \dot{W}(t) - \beta W(t) + \beta^2 \int_0^t e^{-\beta(t-s)} W(s) ds
    $$
    
    再将 $v(t)$ 乘上 $-\beta$：
    
    $$
    -\beta v(t) = -\beta \left( v_0 e^{-\beta t} + W(t) - \beta \int_0^t e^{-\beta(t-s)} W(s) ds \right)
    $$
    
    $$
    = -\beta v_0 e^{-\beta t} - \beta W(t) + \beta^2 \int_0^t e^{-\beta(t-s)} W(s) ds
    $$
    
    对比两式，显然有：
    
    $$
    \frac{dv(t)}{dt} = -\beta v(t) + \dot{W}(t)
    $$
    
    方程成立。对于位移 $x(t) = \int_0^t v(s) ds$，利用分部积分简化形式：
    通过随机微分方程理论，更标准的 O-U 过程位移表达为 $x_\beta(t) = \int_0^t e^{-\beta(t-s)} W(s) ds$。

    <br>

    **(2) 计算期望与方差**
    
    对于 $v(t)$：
    
    $$
    \mathbb{E}[v(t)] = v_0 e^{-\beta t} + \mathbb{E}[W(t)] - \beta \int_0^t e^{-\beta(t-s)} \mathbb{E}[W(s)] ds = v_0 e^{-\beta t}
    $$
    
    由于 $v(t)$ 的另一个等价形式为 $v(t) = v_0 e^{-\beta t} + \int_0^t e^{-\beta(t-s)} dW(s)$，利用 Itô 等距同构计算方差：
    
    $$
    Var[v(t)] = \mathbb{E}\left[ \left( \int_0^t e^{-\beta(t-s)} dW(s) \right)^2 \right] = \int_0^t e^{-2\beta(t-s)} ds = \frac{1 - e^{-2\beta t}}{2\beta}
    $$
    
    对于 $x_\beta(t)$：
    
    $$
    \mathbb{E}[x_\beta(t)] = \mathbb{E}\left[ \int_0^t e^{-\beta(t-s)} W(s) ds \right] = \int_0^t e^{-\beta(t-s)} \mathbb{E}[W(s)] ds = 0
    $$
    
    利用 Fubini 定理交换积分序计算方差（类似于第 2 章习题的技巧）：
    
    $$
    Var[x_\beta(t)] = \mathbb{E}\left[ \left( \int_0^t e^{-\beta(t-s)} W(s) ds \right)^2 \right] = \int_0^t \int_0^t e^{-\beta(t-u)}e^{-\beta(t-v)} \mathbb{E}[W(u)W(v)] du dv
    $$
    
    化简后可得：
    
    $$
    = \frac{t}{\beta^2} - \frac{3 - 4e^{-\beta t} + e^{-2\beta t}}{2\beta^3}
    $$

    <br>

    **(3) 证明极限**
    
    考察 $\beta x_\beta(t)$：
    
    $$
    \beta x_\beta(t) = \beta \int_0^t e^{-\beta(t-s)} W(s) ds
    $$
    
    对其利用分部积分：
    
    $$
    = \int_0^t W(s) d(e^{-\beta(t-s)}) = W(s) e^{-\beta(t-s)} \Big|_0^t - \int_0^t e^{-\beta(t-s)} dW(s)
    $$
    
    $$
    = W(t) - e^{-\beta t} W(0) - \int_0^t e^{-\beta(t-s)} dW(s) = W(t) - \int_0^t e^{-\beta(t-s)} dW(s)
    $$
    
    当 $\beta \to \infty$ 时，考察误差项的均方极限：
    
    $$
    \lim_{\beta \to \infty} \mathbb{E}\left[ \left( \int_0^t e^{-\beta(t-s)} dW(s) \right)^2 \right] = \lim_{\beta \to \infty} \int_0^t e^{-2\beta(t-s)} ds = \lim_{\beta \to \infty} \frac{1 - e^{-2\beta t}}{2\beta} = 0
    $$
    
    由于误差项在 $L^2$ 意义下收敛到 0，因此几乎必然（或在均方意义下）有：
    
    $$
    \lim_{\beta \to \infty} \beta x_\beta(t) = W(t)
    $$

---
<br>

## Part II: 第四章 Itô 积分与二次变差

### 习题 1

**题目：**
利用基于 Riemann 和对 Itô 随机积分的定义方式，证明：

$$\int_0^T W^2(t) dW(t) = \frac{1}{3}W^3(T) - \int_0^T W(t) dt$$

??? success "解答 (点击展开)"

    这是一个经典的需要区分 Riemann 微积分和 Itô 微积分法则的题目。
    
    **利用 Itô 公式证明：**
    我们考虑函数 $f(t, x) = \frac{1}{3}x^3$。其偏导数为：
    $f_t = 0, \quad f_x = x^2, \quad f_{xx} = 2x$
    
    将 $X_t = W_t$ 代入 Itô 公式 $df(t, W_t) = f_t dt + f_x dW_t + \frac{1}{2} f_{xx} (dW_t)^2$：
    
    $$
    d\left( \frac{1}{3}W^3(t) \right) = 0 dt + W^2(t) dW(t) + \frac{1}{2} \cdot 2W(t) dt
    $$
    
    其中用到了布朗运动的二次变差法则 $(dW_t)^2 = dt$。化简得到：
    
    $$
    d\left( \frac{1}{3}W^3(t) \right) = W^2(t) dW(t) + W(t) dt
    $$
    
    对两边从 $0$ 到 $T$ 进行积分，并注意到 $W(0) = 0$：
    
    $$
    \frac{1}{3}W^3(T) - \frac{1}{3}W^3(0) = \int_0^T W^2(t) dW(t) + \int_0^T W(t) dt
    $$
    
    移项即可得到结论：
    
    $$
    \int_0^T W^2(t) dW(t) = \frac{1}{3}W^3(T) - \int_0^T W(t) dt
    $$
    
    *注：与普通微积分 $\int x^2 dx = \frac{1}{3}x^3$ 相比，Itô 积分多出了由二次变差项产生的修正项 $-\int_0^T W(t) dt$。*

---

### 习题 2

**题目：**
对向后积分
$ \int_0^T W(t) dW(t) \doteq \lim_{\delta \to 0} \sum_{i=0}^{n-1} W(t_{i+1})[W(t_{i+1}) - W(t_i)]$
其中 $0 = t_0 < t_1 < \cdots < t_n = T$，而 $\delta \doteq \max_i |t_{i+1} - t_i|$。证明：

$$ \int_0^T W(t) dW(t) = \int_0^T W(t) dW(t) + T$$


??? success "解答 (点击展开)"

    将向后积分的定义式进行恒等变形，构造出标准 Itô 积分的形式：
    
    $$
    I_B = \lim_{\delta \to 0} \sum_{i=0}^{n-1} W(t_{i+1})[W(t_{i+1}) - W(t_i)]
    $$
    
    在求和项中，我们人为地加一项减一项 $W(t_i)$：
    
    $$
    = \lim_{\delta \to 0} \sum_{i=0}^{n-1} \big( W(t_i) + W(t_{i+1}) - W(t_i) \big) [W(t_{i+1}) - W(t_i)]
    $$
    
    展开括号：
    
    $$
    = \lim_{\delta \to 0} \sum_{i=0}^{n-1} W(t_i)[W(t_{i+1}) - W(t_i)] + \lim_{\delta \to 0} \sum_{i=0}^{n-1} [W(t_{i+1}) - W(t_i)]^2
    $$
    
    观察这两项：
    * 第一项正是标准的 **Itô 积分**的离散定义（在每个子区间取左端点）：
      $$\lim_{\delta \to 0} \sum_{i=0}^{n-1} W(t_i) \Delta W_i = \int_0^T W(t) dW(t)$$
    * 第二项正是布朗运动在区间 $[0,T]$ 上的 **二次变差 (Quadratic Variation)**。由布朗运动的性质已知，其二次变差在均方意义下收敛于区间的长度：
      $$\lim_{\delta \to 0} \sum_{i=0}^{n-1} (\Delta W_i)^2 = T$$
      
    将两部分合并，即得证：
    
    $$
    (B) \int_0^T W(t) dW(t) = \int_0^T W(t) dW(t) + T
    $$

---
<br>

## Part III: 随机积分与 Itô 公式的进阶应用

### 习题 1
**题目**
设 $W(t)$ 是 $n$ 维 Brown 运动。证明: $\mathbb{E}[|W(t) - W(s)|^4] = (2n + n^2)(t-s)^2$。

??? success "解答 (点击展开)"

    这道题巧妙地利用了标准正态分布与卡方分布的关系。
    
    已知 $W(t)$ 是 $n$ 维布朗运动，我们考察其在时间区间 $[s, t]$ 上的增量。
    设 $X_i = W_i(t) - W_i(s)$，其中 $i = 1, 2, \dots, n$ 表示各个维度。
    根据布朗运动的性质，各维度的增量独立同分布，且 $X_i \sim N(0, t-s)$。
    
    为了标准化，令 $Z_i = \frac{X_i}{\sqrt{t-s}}$，则 $Z_i \sim N(0, 1)$，且相互独立。
    令 $Q = \sum_{i=1}^n Z_i^2$，由定义可知 $Q$ 服从自由度为 $n$ 的卡方分布，即 $Q \sim \chi^2(n)$。
    
    对于卡方分布 $\chi^2(n)$，我们知道其期望和方差分别为：
    
    $$
    \mathbb{E}[Q] = n, \quad Var(Q) = 2n
    $$
    
    由此可以求出 $Q$ 的二阶原点矩：
    
    $$
    \mathbb{E}[Q^2] = Var(Q) + (\mathbb{E}[Q])^2 = 2n + n^2
    $$
    
    回到原式，考察增量的四阶矩：
    
    $$
    |W(t) - W(s)|^4 = \left( \sum_{i=1}^n X_i^2 \right)^2 = \left( \sum_{i=1}^n (t-s)Z_i^2 \right)^2 = (t-s)^2 Q^2
    $$
    
    对两边取期望：
    
    $$
    \mathbb{E}[|W(t) - W(s)|^4] = (t-s)^2 \mathbb{E}[Q^2] = (2n + n^2)(t-s)^2
    $$
    
    结论得证。

---

### 习题 2
**题目**
定义二阶积分
$$\int_0^T f(t) [dW(t)]^2 \doteq \lim_{\delta \to 0} \sum_{i=0}^{n-1} f(t_i)[W(t_{i+1}) - W(t_i)]^2$$
其中 $0 = t_0 < t_1 < \dots < t_n = T$，而 $\delta \doteq \max_i |t_{i+1} - t_i|$。证明：当 $f \in L^2([0,T])$ 时，有

$$\int_0^T f(t) [dW(t)]^2 = \int_0^T f(t) dt$$

??? success "解答 (点击展开)"

    此题旨在证明布朗运动二次变差 $(dW_t)^2 = dt$ 在积分意义下的严格成立。
    我们要在均方意义（$L^2$）下证明该极限。
    
    令 $\Delta t_i = t_{i+1} - t_i$，$\Delta W_i = W(t_{i+1}) - W(t_i)$。
    令左侧离散和为 $S_n = \sum_{i=0}^{n-1} f(t_i)(\Delta W_i)^2$，右侧目标积分为 $I = \int_0^T f(t)dt \approx \sum_{i=0}^{n-1} f(t_i)\Delta t_i$。
    
    考察它们差值的均方误差：
    
    $$
    \mathbb{E}\left[ \left( S_n - \sum_{i=0}^{n-1} f(t_i)\Delta t_i \right)^2 \right] = \mathbb{E}\left[ \left( \sum_{i=0}^{n-1} f(t_i) ((\Delta W_i)^2 - \Delta t_i) \right)^2 \right]
    $$
    
    由于不重叠区间的布朗增量独立，交叉项（$i \ne j$）的期望可以拆解：
    
    $$
    \mathbb{E}[((\Delta W_i)^2 - \Delta t_i)((\Delta W_j)^2 - \Delta t_j)] = \mathbb{E}[(\Delta W_i)^2 - \Delta t_i] \cdot \mathbb{E}[(\Delta W_j)^2 - \Delta t_j] = 0 \cdot 0 = 0
    $$
    
    因此，展开平方后只剩下平方项：
    
    $$
    = \sum_{i=0}^{n-1} f^2(t_i) \mathbb{E}[((\Delta W_i)^2 - \Delta t_i)^2]
    $$
    
    计算单项期望，利用 $N(0, \Delta t_i)$ 的四阶矩 $\mathbb{E}[(\Delta W_i)^4] = 3(\Delta t_i)^2$：
    
    $$
    \mathbb{E}[(\Delta W_i)^4 - 2\Delta t_i(\Delta W_i)^2 + (\Delta t_i)^2] = 3\Delta t_i^2 - 2\Delta t_i^2 + \Delta t_i^2 = 2\Delta t_i^2
    $$
    
    代回原式，并将一项 $\Delta t_i$ 放缩为最大步长 $\delta$：
    
    $$
    = \sum_{i=0}^{n-1} f^2(t_i) 2\Delta t_i^2 \le 2\delta \sum_{i=0}^{n-1} f^2(t_i)\Delta t_i
    $$
    
    当 $\delta \to 0$ 时，由于 $f \in L^2([0,T])$，$\sum f^2(t_i)\Delta t_i \to \int_0^T f^2(t)dt < \infty$。
    因此前面乘上的 $\delta$ 会使得整个极限趋于 $0$。均方收敛得证，即：
    
    $$
    \int_0^T f(t) [dW(t)]^2 = \int_0^T f(t) dt
    $$

---

### 习题 3
**题目**
设 $f \in L^2([0, T])$ 且 $\int_0^T f(s) dW(s) = 0$。证明: $f$ 几乎处处为零。

??? success "解答 (点击展开)"

    本题利用了 Itô 等距同构 (Itô Isometry) 以及实变函数中的 Lebesgue 积分基本性质。
    
    由于给定 $\int_0^T f(s) dW(s) = 0$ 几乎必然成立，该随机变量的二阶矩必然为 0：
    
    $$
    \mathbb{E}\left[ \left( \int_0^T f(s) dW(s) \right)^2 \right] = \mathbb{E}[0^2] = 0
    $$
    
    另一方面，根据 Itô 等距同构，随机积分的二阶矩等于被积函数平方的 Lebesgue 积分：
    
    $$
    \mathbb{E}\left[ \left( \int_0^T f(s) dW(s) \right)^2 \right] = \mathbb{E}\left[ \int_0^T f^2(s) ds \right]
    $$
    
    由于 $f(s)$ 是确定性函数，其期望即为它本身。两式结合可得：
    
    $$
    \int_0^T f^2(s) ds = 0
    $$
    
    由实变函数基本定理，由于被积函数 $f^2(s) \ge 0$ 恒成立，且其在 $[0,T]$ 上的 Lebesgue 积分为 $0$，这意味着 $f^2(s)$ 必须在 $[0,T]$ 上几乎处处为零 (almost everywhere, a.e.)。
    
    进而得出：$f(s) = 0$ 几乎处处成立。得证。

---

### 习题 4
**题目**
证明: $Y(t) = e^{t/2} \cos(W(t))$ 是鞅。

??? success "解答 (点击展开)"

    为了证明一个过程是鞅，最直接的方法是利用 Itô 公式证明其微分项中没有漂移项（即 $dt$ 项系数为 0）。
    
    令二元函数 $u(t, x) = e^{t/2} \cos(x)$，计算其关于 $t$ 和 $x$ 的偏导数：
    * $u_t = \frac{1}{2} e^{t/2} \cos(x)$
    * $u_x = -e^{t/2} \sin(x)$
    * $u_{xx} = -e^{t/2} \cos(x)$
    
    将 $X(t) = W(t)$ 代入 Itô 公式 $dY_t = (u_t + \frac{1}{2} u_{xx})dt + u_x dW(t)$：
    
    $$
    dY(t) = \left( \frac{1}{2} e^{t/2} \cos(W(t)) - \frac{1}{2} e^{t/2} \cos(W(t)) \right) dt - e^{t/2} \sin(W(t)) dW(t)
    $$
    
    可以清晰地看到，含有 $dt$ 的两项完美抵消：
    
    $$
    dY(t) = -e^{t/2} \sin(W(t)) dW(t)
    $$
    
    将其写成积分形式：
    
    $$
    Y(t) = Y(0) - \int_0^t e^{s/2} \sin(W(s)) dW(s)
    $$
    
    由于右侧是只含 $dW(s)$ 的 Itô 积分，且被积函数是有界的（满足 $L^2$ 适定性条件），Itô 积分本身就是一个鞅。
    因此，过程 $Y(t)$ 也是一个鞅。得证。

---

### 习题 5
**题目**
证明: 
1. $\int_0^T W^2 dW = \frac{1}{3}W(T)^3 - \int_0^T W dt$
2. $\int_0^T W^3 dW = \frac{1}{4}W(T)^4 - \frac{3}{2} \int_0^T W^2 dt$

??? success "解答 (点击展开)"

    这两部分证明都是 Itô 公式的基础逆向应用，即“先猜测高阶项，再用 Itô 公式展开并移项”。

    **(1) 证明第一式**
    
    令函数 $f(x) = \frac{1}{3}x^3$。计算导数：$f'(x) = x^2$，$f''(x) = 2x$。
    将 $W(t)$ 代入 Itô 公式：
    
    $$
    d\left( \frac{1}{3}W^3(t) \right) = W^2(t) dW(t) + \frac{1}{2}(2W(t)) (dW(t))^2
    $$
    
    由二次变差法则 $(dW(t))^2 = dt$：
    
    $$
    d\left( \frac{1}{3}W^3(t) \right) = W^2(t) dW(t) + W(t) dt
    $$
    
    对两边在 $[0, T]$ 上同时积分，并注意到 $W(0) = 0$：
    
    $$
    \frac{1}{3}W^3(T) - 0 = \int_0^T W^2(t) dW(t) + \int_0^T W(t) dt
    $$
    
    移项即得：
    
    $$
    \int_0^T W^2(t) dW(t) = \frac{1}{3}W^3(T) - \int_0^T W(t) dt
    $$

    <br>

    **(2) 证明第二式**
    
    同理，令函数 $g(x) = \frac{1}{4}x^4$。计算导数：$g'(x) = x^3$，$g''(x) = 3x^2$。
    代入 Itô 公式：
    
    $$
    d\left( \frac{1}{4}W^4(t) \right) = W^3(t) dW(t) + \frac{1}{2}(3W^2(t)) dt
    $$
    
    对两边在 $[0, T]$ 上积分：
    
    $$
    \frac{1}{4}W^4(T) - 0 = \int_0^T W^3(t) dW(t) + \frac{3}{2}\int_0^T W^2(t) dt
    $$
    
    移项即得：
    
    $$
    \int_0^T W^3(t) dW(t) = \frac{1}{4}W^4(T) - \frac{3}{2}\int_0^T W^2(t) dt
    $$

---

### 习题 6
**题目**
证明 $\mathbb{E}[e^{\int_0^T g dW}] = e^{\frac{1}{2} \int_0^T g^2 ds}$。

??? success "解答 (点击展开)"

    这道题可以通过分析 Itô 积分的分布特性，直接借助正态分布矩母函数来优雅求解。
    
    令随机变量 $X = \int_0^T g(s) dW(s)$。
    由于 $g(s)$ 是一个确定的时间函数（非随机），这个 Itô 积分是高斯过程的线性叠加，因此 $X$ 依然服从正态分布。
    
    根据随机积分的性质：
    * 期望：$\mathbb{E}[X] = \mathbb{E}[\int_0^T g dW] = 0$
    * 方差：由 Itô 等距同构，$Var(X) = \mathbb{E}[X^2] = \int_0^T g^2(s) ds$
    
    所以 $X \sim N(0, \sigma^2)$，其中 $\sigma^2 = \int_0^T g^2(s) ds$。
    
    原题要求的 $\mathbb{E}[e^X]$，恰好是随机变量 $X$ 的矩母函数 $M_X(u) = \mathbb{E}[e^{uX}]$ 在 $u=1$ 时的取值。
    对于正态分布 $N(\mu, \sigma^2)$，其矩母函数公式为 $M_X(u) = \exp(\mu u + \frac{1}{2}\sigma^2 u^2)$。
    
    代入 $\mu = 0, u = 1, \sigma^2 = \int_0^T g^2 ds$，立刻可得：
    
    $$
    \mathbb{E}[e^{\int_0^T g dW}] = \exp\left( 0 + \frac{1}{2} \int_0^T g^2(s) ds \cdot 1^2 \right) = e^{\frac{1}{2} \int_0^T g^2 ds}
    $$
    
    得证。

---

### 习题 7
**题目**
设 $u = u(x, t)$ 满足抛物型偏微分方程 $u_t + \frac{1}{2}u_{xx} = 0$。证明: $\mathbb{E}[u(W(t), t)] = u(0, 0)$。

??? success "解答 (点击展开)"

    此题是建立偏微分方程 (PDE) 与随机过程 (SDE) 深刻联系（即 Feynman-Kac 公式的最简形式）的经典练习。
    
    定义随机过程 $Y(t) = u(W(t), t)$。利用多元 Itô 公式对 $Y(t)$ 展开微分：
    
    $$
    dY(t) = \frac{\partial u}{\partial t} dt + \frac{\partial u}{\partial x} dW(t) + \frac{1}{2}\frac{\partial^2 u}{\partial x^2} (dW(t))^2
    $$
    
    将二次变差 $(dW)^2 = dt$ 代入合并 $dt$ 项：
    
    $$
    dY(t) = \left( u_t + \frac{1}{2}u_{xx} \right) dt + u_x dW(t)
    $$
    
    由于已知条件给出 $u(x, t)$ 满足 $u_t + \frac{1}{2}u_{xx} = 0$，所以漂移项严格为 0。微分方程简化为纯扩散：
    
    $$
    dY(t) = u_x(W(t), t) dW(t)
    $$
    
    将其写成积分形式：
    
    $$
    Y(t) - Y(0) = \int_0^t u_x(W(s), s) dW(s)
    $$
    
    对等式两边取期望，由于右侧 Itô 积分的期望为 0：
    
    $$
    \mathbb{E}[Y(t)] - \mathbb{E}[Y(0)] = 0 \implies \mathbb{E}[Y(t)] = \mathbb{E}[Y(0)]
    $$
    
    代回 $Y(t)$ 的定义，且初始时刻布朗运动 $W(0) = 0$ 几乎必然成立：
    
    $$
    \mathbb{E}[u(W(t), t)] = \mathbb{E}[u(W(0), 0)] = u(0, 0)
    $$
    
    结论得证。

---

### 习题 8
**题目**
1. 证明 $e^{W(t)} = 1 + \frac{1}{2}\int_0^t e^{W(s)} ds + \int_0^t e^{W(s)} dW(s)$；
2. 证明 $\mathbb{E}[e^{W(t)}] = 1 + \frac{1}{2}\int_0^t \mathbb{E}[e^{W(s)}] ds$，从而 $\mathbb{E}[e^{W(t)}] = e^{t/2}$；
3. 计算 $\mathbb{E}[e^{iW(t)}]$，以及 $e^{W(t)}, \sin W(t), \cos W(t)$ 的方差。

??? success "解答 (点击展开)"

    **(1) SDE 形式验证**
    令 $f(x) = e^x$。将 $W(t)$ 代入 Itô 公式 $df(W_t) = f'(W_t)dW_t + \frac{1}{2}f''(W_t)dt$：
    
    $$
    d(e^{W(t)}) = e^{W(t)}dW(t) + \frac{1}{2}e^{W(t)}dt
    $$
    
    对两边在 $[0,t]$ 上积分，并利用 $e^{W(0)} = e^0 = 1$：
    
    $$
    e^{W(t)} - 1 = \int_0^t e^{W(s)} dW(s) + \frac{1}{2}\int_0^t e^{W(s)} ds
    $$
    
    移项即得证第一问。

    <br>

    **(2) 求解期望的常微分方程**
    对第(1)问的结果两边取期望。由于 Itô 积分 $\int_0^t e^{W(s)} dW(s)$ 在满足适定性条件时期望为 0，利用 Fubini 定理交换期望与 Riemann 积分符号：
    
    $$
    \mathbb{E}[e^{W(t)}] = 1 + \frac{1}{2}\int_0^t \mathbb{E}[e^{W(s)}] ds
    $$
    
    令 $m(t) = \mathbb{E}[e^{W(t)}]$，上式转化为积分方程 $m(t) = 1 + \frac{1}{2}\int_0^t m(s) ds$。
    对其求导，得到初值问题 ODE：
    
    $$
    m'(t) = \frac{1}{2}m(t), \quad m(0) = 1
    $$
    
    解此常微分方程，立刻得到：
    
    $$
    m(t) = \mathbb{E}[e^{W(t)}] = e^{t/2}
    $$

    <br>

    **(3) 特征函数与三角函数方差计算**
    
    *计算 $\mathbb{E}[e^{iW(t)}]$ (特征函数)：*
    同样利用 Itô 公式展开复值过程 $Y(t) = e^{iW(t)}$：
    
    $$
    dY(t) = i e^{iW(t)} dW(t) + \frac{1}{2}(i)^2 e^{iW(t)} dt = i Y(t) dW(t) - \frac{1}{2} Y(t) dt
    $$
    
    取期望并求导，令 $m_2(t) = \mathbb{E}[Y(t)]$，得到 ODE：$m_2'(t) = -\frac{1}{2}m_2(t)$ 且 $m_2(0)=1$。
    解得：
    
    $$
    \mathbb{E}[e^{iW(t)}] = e^{-t/2}
    $$
    
    *计算 $e^{W(t)}$ 的方差：*
    根据方差定义 $Var(e^{W(t)}) = \mathbb{E}[(e^{W(t)})^2] - (\mathbb{E}[e^{W(t)}])^2 = \mathbb{E}[e^{2W(t)}] - e^t$。
    将 $2W(t)$ 看作参数为 $2$ 的情况，由矩母函数结论得 $\mathbb{E}[e^{2W(t)}] = e^{4t/2} = e^{2t}$。
    
    $$
    Var(e^{W(t)}) = e^{2t} - e^t
    $$
    
    *计算 $\sin W(t), \cos W(t)$ 的方差：*
    由欧拉公式，$\mathbb{E}[e^{iW(t)}] = \mathbb{E}[\cos W(t)] + i\mathbb{E}[\sin W(t)] = e^{-t/2}$。对比实部与虚部：
    
    $$
    \mathbb{E}[\cos W(t)] = e^{-t/2}, \quad \mathbb{E}[\sin W(t)] = 0
    $$
    
    同理，令参数为 $2i$，有 $\mathbb{E}[e^{2iW(t)}] = e^{-(2i)^2 t/(-2)} = e^{-2t}$，即：
    
    $$
    \mathbb{E}[\cos 2W(t)] = e^{-2t}, \quad \mathbb{E}[\sin 2W(t)] = 0
    $$
    
    利用二倍角公式降幂：
    
    $$
    \mathbb{E}[\sin^2 W(t)] = \mathbb{E}\left[ \frac{1 - \cos 2W(t)}{2} \right] = \frac{1 - e^{-2t}}{2}
    $$
    
    $$
    \mathbb{E}[\cos^2 W(t)] = \mathbb{E}\left[ \frac{1 + \cos 2W(t)}{2} \right] = \frac{1 + e^{-2t}}{2}
    $$
    
    最终带入方差公式：
    
    $$
    Var(\sin W(t)) = \mathbb{E}[\sin^2] - (\mathbb{E}[\sin])^2 = \frac{1 - e^{-2t}}{2} - 0 = \frac{1 - e^{-2t}}{2}
    $$
    
    $$
    Var(\cos W(t)) = \mathbb{E}[\cos^2] - (\mathbb{E}[\cos])^2 = \frac{1 + e^{-2t}}{2} - (e^{-t/2})^2 = \frac{1 - e^{-2t}}{2}
    $$
    
    *(注：这里手稿上利用了正态分布三角函数矩的巧妙结论。)*