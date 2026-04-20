# 🌊 Stochastic Processes

本模块聚焦随机过程的相关内容。

## 一、 条件期望的代数性质 (Conditional Expectation)

!!! abstract "核心运算性质"
    条件期望 $E[X \mid \mathcal{F}]$ 本身是一个关于 $\mathcal{F}$ 可测的随机变量。
    
    1. **提出已知性 (Pulling out knowns)**：若 $Y$ 是 $\mathcal{F}$ 可测的，则 $E[XY \mid \mathcal{F}] = Y E[X \mid \mathcal{F}]$。
    2. **塔牌性质 (Tower Property)**：若 $\mathcal{F}_1 \subset \mathcal{F}_2$，则 $E[E[X \mid \mathcal{F}_2] \mid \mathcal{F}_1] = E[X \mid \mathcal{F}_1]$。
    3. **Jensen 不等式**：对于凸函数 $\phi$，有 $\phi(E[X \mid \mathcal{F}]) \le E[\phi(X) \mid \mathcal{F}]$。**该不等式常用于证明凸函数的鞅为下鞅**。

---

## 二、 鞅与停时 (Martingales and Stopping Times)

!!! info "鞅的定义 (Martingale)"
    适应于信息流 $\mathcal{F}_n$ 的可积序列 $X_n$，满足：

    \[ 
    E[X_{n+1} \mid \mathcal{F}_n] = X_n 
    \]

    （若为 $\ge$ 则为下鞅 Submartingale，若为 $\le$ 则为上鞅 Supermartingale）。

!!! info "停时 (Stopping Time)"
    一个取值为 $\mathbb{N} \cup \{\infty\}$ 的随机变量 $T$，若对任意 $n$，事件 $\{T \le n\}$ 完全由 $\mathcal{F}_n$ 决定（即发生与否仅依赖于当下的信息），则称 $T$ 为停时。
    * **示例**：首中时 $\inf\{n: S_n = 0\}$ 是停时。

??? success "可选停时定理 (Optional Stopping Theorem)"
    在一定条件下，鞅在停时 $T$ 处的期望等于其初始期望：$E[X_T] = E[X_0]$。满足以下**任一**条件即可成立：
    
    1. **时间有界**：停时 $T$ 几乎处处有界（即存在常数 $N$ 使 $P(T \le N) = 1$）。
    2. **过程一致有界**：对所有 $n \le T$，存在常数 $M$ 使得 $|X_n| \le M$ 几乎处处成立。
    3. **期望时间有限且增量有界**：$E[T] < \infty$，且鞅的单步增量有界（$|X_m - X_n| \le C$）。

---

## 三、 鞅的核心极限定理

!!! abstract "Doob 极大值不等式 (Doob's Maximal Inequality)"
    对于非负下鞅 $X_n$（或鞅的绝对值），其历史最大值的尾概率受末尾期望控制：

    \[ 
    \lambda P(\max_{1\le k\le n} X_k \ge \lambda) \le E[X_n I(\max_{1\le k\le n} X_k \ge \lambda)] \le E[X_n] 
    \]

!!! abstract "Doob 鞅收敛定理 (Doob's Convergence Theorem)"
    若 $X_n$ 为下鞅，且其正部期望一致有界（即 $\sup_n E[X_n^+] < \infty$），则存在一个可积的随机变量 $X$，使得：

    \[ 
    X_n \xrightarrow{a.s.} X 
    \]

    * 注意：仅保证几乎处处收敛，不保证 $L^1$ 收敛。若需 $E[X_n] \to E[X]$，必须外加**一致可积性 (Uniform Integrability)** 条件。

??? note "Doob-Meyer 分解定理"
    任何一个下鞅 $X_n$ 都可以被唯一分解为一个鞅 $M_n$ 和一个非减的可预测序列 (Predictable sequence) $A_n$ 的和：

    \[ 
    X_n = M_n + A_n 
    \]
    
    其中 $A_n$ 可预测意味着 $A_n$ 在 $n-1$ 时刻就已经完全已知（即 $\mathcal{F}_{n-1}$ 可测）。