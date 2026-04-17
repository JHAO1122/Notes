---
tags:
  - 随机微分方程
  - 课后习题
---

# 📝 课后习题精解：期权定价与停时

!!! abstract "关于本页"
    本页面收录了《随机微分方程》课程中关于停时理论、高维布朗运动击中概率以及欧式期权定价（Black-Scholes 模型）的核心课后习题解答。内容涵盖了停时的严格证明、利用鞅的性质导出 Dirichlet 边值问题，以及利用热传导方程基本解（高斯核）推导看跌期权公式的完整过程。

---

## Part I: 停时与击中概率

### 习题 1
**题目**
设 $W$ 为一维 Brown 运动，定义 $\tau = \inf\{t: W(t) \in (a, b]\}$ 为首次触碰区间 $(a, b]$ 的时间。证明 $\tau$ 是一个停时。

??? success "解答 (点击展开)"

    为了严谨证明 $\tau$ 是一个停时，我们根据区间 $(a, b]$ 的性质，将其拆分为两个首次击中时间来分别讨论。
    
    定义两个辅助变量：
    
    \[
    \tau_a = \inf\{t: W(t) > a\}
    \]
    
    \[
    \tau_b = \inf\{t: W(t) \le b\}
    \]

    **第一步：证明 $\tau_a$ 和 $\tau_b$ 的可测性**

    对于 $\tau_b$（首次进入闭集 $(-\infty, b]$）：由于布朗运动路径连续，闭集的首次击中时间满足：

    \[
    \{\tau_b \le t\} = \inf_{0 \le s \le t} W(s) \le b
    \]

    利用有理数的稠密性，这等价于在有理数集 $\mathbb{Q}$ 上取下确界：

    \[
    \{\tau_b \le t\} = \bigcap_{n=1}^\infty \bigcup_{s \in \mathbb{Q} \cap [0, t]} \left\{ W_s \le b + \frac{1}{n} \right\}
    \]

    因为右侧是 $\mathcal{F}_t$ 可测事件的可数交并，故 $\tau_b$ 是停时。

    对于 $\tau_a$（首次进入开集 $(a, +\infty)$）：同理，利用开集的性质和有理数逼近：

    \[
    \{\tau_a < t\} = \bigcup_{s \in \mathbb{Q} \cap [0, t)} \{ W_s > a \} \in  \mathcal{F}_t
    \]

    结合信息流的右连续性，$\tau_a$ 也是一个停时。

    **第二步：组合逻辑证明 $\tau$ 是停时**
    
    现在我们考察目标停时 $\tau = \inf\{t: W(t) \in (a, b]\}$。根据初始状态 $W(0)$ 的不同，我们可以将 $\{\tau \le t\}$ 拆解：
    
    1. **若 $W(0) \le a$**：必须先向上穿过 $a$ 才能进入 $(a, b]$，此时事件等价于：$\{W(0) \le a\} \cap \{\tau_a \le t, \inf_{\tau_a \le s \le t} W(s) \le b\}$

    2. **若 $W(0) > b$**：必须先向下穿过 $b$ 才能进入 $(a, b]$，此时事件等价于：$\{W(0) > b\} \cap \{\tau_b \le t\}$
       
    3. **若 $W(0) \in (a, b]$**：此时 $\tau = 0 \le t$ 必然成立，事件等价于 $\{W(0) \in (a, b]\}$。

    将上述三种互斥且完备的情况取并集，得到：
    
    \[
    \{\tau \le t\} = \left( \{W(0) \le a\} \cap \{\tau_a \le t, \inf_{\tau_a \le s \le t} W_s \le b\} \right) \cup \left( \{W(0) > b\} \cap \{\tau_b \le t\} \right) \cup \{W(0) \in (a, b]\}
    \]
    
    由于上述表达式中涉及的所有的基本事件（包括初始状态集合、$\tau_a$ 和 $\tau_b$ 的水平穿越事件）均是 $\mathcal{F}_t$ 可测的，且 $\sigma$ 域对有限并集、交集封闭，因此整体事件 $\{\tau \le t\} \in \mathcal{F}_t$。
    
    故 $\tau$ 满足停时的严格定义。得证。

---

### 习题 2
**题目**
设 $W$ 为 $n$ 维 Brown 运动 ($n \ge 3$)。记 $X = W + x_0$，其中 $x_0$ 落在球壳 $R_1 < |x| < R_2$ 内。计算 $X$ 在碰到内球面 $|x| = R_1$ 之前碰到外球面 $|x| = R_2$ 的概率。

??? success "解答 (点击展开)"

    **第一步：利用 Itô 公式与鞅的性质导出 Dirichlet 边值问题**
    
    设 $\tau$ 为过程 $X$ 首次离开球壳区域 $D = \{x: R_1 < |x| < R_2\}$ 的退出时间。令 $u(x)$ 为从初始点 $x$ 出发，先击中外球面的概率，即 $u(x) = P^x(|X_\tau| = R_2)$。
    
    我们构造一个随机过程 $M_t = u(X_{t \wedge \tau})$。对 $u(X_t)$ 应用多维 Itô 公式展开：

    \[
    u(X_{t \wedge \tau}) = u(X_0) + \int_0^{t \wedge \tau} \nabla u(X_s) dW_s + \int_0^{t \wedge \tau} \frac{1}{2} \Delta u(X_s) ds
    \]

    为了使期望概率在停时前保持守恒，过程 $M_t = u(X_{t \wedge \tau})$ 必须是一个**鞅 (Martingale)**。
    
    鞅的核心性质要求其微分的漂移项（$ds$ 项）期望恒为 0。因此，被积函数必须满足：

    \[
    \frac{1}{2} \Delta u(x) = 0 \implies \Delta u(x) = 0 \quad (x \in D)
    \]

    对于边界条件：当过程击中外球面时，事件必然发生，概率为 1；击中内球面时，事件未发生，概率为 0。这就严谨地将概率问题转化为了如下的 Dirichlet 边值问题：

    \[
    \begin{cases}
    \Delta u(x) = 0, & R_1 < |x| < R_2 \\
    u(x) = 0, & |x| = R_1 \\
    u(x) = 1, & |x| = R_2
    \end{cases}
    \]

    **第二步：求解偏微分方程**

    由于区域是对称的，假设解具有径向形式 $u(x) = v(r)$，其中 $r = |x|$。在 $n \ge 3$ 时，Laplace 算子的基本解形式为：

    \[
    v(r) = A + B r^{2-n}
    \]

    代入边界条件：
    当 $r = R_1$ 时：

    \[
    A + B R_1^{2-n} = 0 \implies A = -B R_1^{2-n}
    \]

    当 $r = R_2$ 时：

    \[
    A + B R_2^{2-n} = 1
    \]

    将第一式代入第二式，解得常数 $B$ 和 $A$：

    \[
    B = \frac{1}{R_2^{2-n} - R_1^{2-n}}, \quad A = -\frac{R_1^{2-n}}{R_2^{2-n} - R_1^{2-n}}
    \]

    代回原函数并化简，得到最终的击中概率：

    \[
    u(x) = \frac{|x|^{2-n} - R_1^{2-n}}{R_2^{2-n} - R_1^{2-n}}
    \]

---
<br>

## Part II: 欧式期权定价实战

### 习题 3
**题目**
推导欧式看跌期权 (European Put Option) 的 Black-Scholes 公式。

??? success "解答 (点击展开)"

    **1. 建立 Black-Scholes 偏微分方程 (PDE)**
    设欧式看跌期权价格为 $P(S, t)$。基于无套利原理和几何布朗运动的假设，看跌期权同样满足 BS PDE：

    \[
    \frac{\partial P}{\partial t} + rS \frac{\partial P}{\partial S} + \frac{1}{2}\sigma^2 S^2 \frac{\partial^2 P}{\partial S^2} - rP = 0
    \]

    终端条件为：$P(S, T) = \max(K - S, 0) = (K - S)^+$。

    **2. 变量代换与方程简化**
    令 $x = \ln S$ (即 $S = e^x$)，并反转时间 $\tau = T - t$。设 $V(x, \tau) = P(e^x, T - \tau)$。代入方程后，为了消去一阶偏导和常数项，再作代换 $V(x, \tau) = u(x, \tau) e^{\alpha \tau + \beta x}$，解得系数：

    \[
    \beta = \frac{1}{2} - \frac{r}{\sigma^2}, \quad \alpha = -r - \frac{1}{2\sigma^2}\left(r - \frac{1}{2}\sigma^2\right)^2
    \]

    此时方程化为标准的热传导方程：

    \[
    \frac{\partial u}{\partial \tau} = \frac{1}{2}\sigma^2 \frac{\partial^2 u}{\partial x^2}
    \]

    **3. 利用高斯核 (Gaussian Kernel) 积分求解**
    转化后的初始条件为 $u(y, 0) = e^{-\beta y} \max(K - e^y, 0)$。热传导方程的解可以利用高斯核的卷积给出：

    \[
    u(x, \tau) = \frac{1}{\sigma \sqrt{2\pi \tau}} \int_{-\infty}^{\ln K} \exp\left(-\frac{(x-y)^2}{2\sigma^2 \tau}\right) e^{-\beta y} (K - e^y) dy
    \]

    (注意积分上限为 $\ln K$，因为当 $e^y > K$ 时收益为 0)
    
    我们将积分拆分为 $I_1$ (含 $K$) 和 $I_2$ (含 $e^y$) 两项。对指数部分进行“配方法 (Completing the square)”，并通过变量替换 $z = \frac{y - \mu}{\sigma}$ 将其化为标准正态分布累积函数 $N(\cdot)$。

    **4. 还原得到最终定价公式**
    最后利用 $P(S, t) = e^{\alpha \tau + \beta x} u(x, \tau)$ 还原变量，指数项完美相消，得到看跌期权的 Black-Scholes 公式：

    \[
    P(S, t) = K e^{-r(T-t)} N(-d_2) - S N(-d_1)
    \]

    其中：

    \[
    d_1 = \frac{\ln(S/K) + (r + \frac{1}{2}\sigma^2)(T-t)}{\sigma \sqrt{T-t}}
    \]

    \[
    d_2 = d_1 - \sigma \sqrt{T-t}
    \]

    推导完毕。