# 🛰️ Research Frontier

这里是我通过脚本制作的论文速递版面，它会调用 GitHub 机器人，帮我在北京时间每天早上八点会自动扫描 [arXiv](https://arxiv.org/) 上的关于 Knockoff 与 Conformal Prediction 的统计学论文。


> 更新于: 2026-04-22

---
### Trustworthy Feature Importance Avoids Unrestricted Permutations 

- [ ] **分类**: Knockoff | **日期**: 2026-04-13
- **链接**: [PDF](http://arxiv.org/abs/2604.11253v1)

!!! note "AI 核心解读"

    本文创新点在于：通过条件模型依赖性与高斯变换的Knockoffs方法，从统计推导上规避了无限制置换带来的外推误差；同时提出受限ALE图设计，从理论性质上确保特征重要性评估的可靠性。

### Distribution-free screening of spatially variable genes in spatial transcriptomics 

- [ ] **分类**: Knockoff | **日期**: 2026-03-10
- **链接**: [PDF](http://arxiv.org/abs/2603.09061v1)

!!! note "AI 核心解读"

    本文提出了一种基于准似然比统计量MM-test的无分布空间可变基因筛选方法，该方法通过结合辅助空间距离信息构造检验统计量，并引入knockoff程序控制错误发现率；同时，理论证明了方法具有选择一致性、错误发现率控制及选择后聚类误差界等统计保证。

### Variable selection via knockoffs for clustered data 

- [ ] **分类**: Knockoff | **日期**: 2026-02-23
- **链接**: [PDF](http://arxiv.org/abs/2602.19398v1)

!!! note "AI 核心解读"

    本文创新点在于针对聚类数据提出了一种分层构造辅助变量的两步法：先将观测层预测变量分解为聚类均值与组内偏差两个正交分量，再分别在聚类层与观测层独立执行敲除筛选。该方法通过分层构造的敲除变量矩阵，严格保持了组内相关性结构，在控制错误发现率的同时提升了筛选功效。

### Improving the adjusted Benjamini--Hochberg method using e-values in knockoff-assisted variable selection 

- [ ] **分类**: Knockoff | **日期**: 2026-02-12
- **链接**: [PDF](http://arxiv.org/abs/2602.11610v1)

!!! note "AI 核心解读"

    本文创新点在于：将Sarkar和Tang（2022）的方法识别为未归一化e值加权Benjamini-Hochberg程序的特例，并引入有界p值到e值的校准器，实现更精细灵活的权重分配。在此基础上，构建了三种能控制错误发现率的程序，其模拟与实证分析在多种场景下均展现出优于原方法的性能。

### GRIP2: A Robust and Powerful Deep Knockoff Method for Feature Selection 

- [ ] **分类**: Knockoff | **日期**: 2026-01-30
- **链接**: [PDF](http://arxiv.org/abs/2602.00218v1)

!!! note "AI 核心解读"

    本文提出了一种基于二维正则化曲面积分的深度 knockoff 特征重要性统计量，通过块随机采样高效聚合不同正则化强度下的特征活动强度，其构造具有天然的反对称性，从而保证了有限样本下的错误发现率控制。该方法在高度相关和低信噪比场景下表现出更强的稳健性和检验功效。

### Conformal Robust Set Estimation 

- [ ] **分类**: Conformal | **日期**: 2026-04-20
- **链接**: [PDF](http://arxiv.org/abs/2604.18441v1)

!!! note "AI 核心解读"

    本文提出了一种基于半质量半径（即点到其第k近邻距离）的非共形评分函数，构建了具有边际覆盖保证的鲁棒共形预测集，并从理论上证明了该经验预测集以指数浓度速率概率收敛于基于距离测度泛函定义的总体鲁棒中心集。

### Online Conformal Prediction with Adversarial Semi-bandit Feedback via Regret Minimization 

- [ ] **分类**: Conformal | **日期**: 2026-04-20
- **链接**: [PDF](http://arxiv.org/abs/2604.17984v1)

!!! note "AI 核心解读"

    本文创新点在于：将在线共形预测重构为对抗性半赌博机问题，通过将候选预测集视为臂并利用已有对抗性赌博机算法，建立了学习器遗憾与长期覆盖率保证之间的显式理论联系。该方法在仅当真实标签落入预测集时才获得部分反馈的设定下，首次实现了无需分布假设的长期覆盖率保证，并控制了预测集规模。

### On a Probability Inequality for Order Statistics with Applications to Bootstrap, Conformal Prediction, and more 

- [ ] **分类**: Conformal | **日期**: 2026-04-16
- **链接**: [PDF](http://arxiv.org/abs/2604.15229v1)

!!! note "AI 核心解读"

    本文创新点在于：通过构建基于顺序统计量的概率不等式，在放宽独立同分布假设的条件下，推导出该不等式的近似成立形式，并严格证明其理论边界。该不等式为自助法、共形预测等方法的有效性提供了统一的理论框架，尤其在高维数据中实现了有限样本下的统计保证。

### Differentially Private Conformal Prediction 

- [ ] **分类**: Conformal | **日期**: 2026-04-16
- **链接**: [PDF](http://arxiv.org/abs/2604.14621v2)

!!! note "AI 核心解读"

    本文提出了一种不依赖数据分割的差分私有保形预测框架，通过利用差分隐私机制的稳定性，建立了与理想保形预测的直接理论关联，从而在保证隐私的同时继承了原有的有效性。在此基础上，作者构建了完整的私有训练与私有分位数校准流程，在相同隐私预算下，相比现有分割方法能获得更紧的预测集，并严格证明了其端到端隐私保障与覆盖性质。

### Conformal Prediction with Time-Series Data via Sequential Conformalized Density Regions 

- [ ] **分类**: Conformal | **日期**: 2026-04-08
- **链接**: [PDF](http://arxiv.org/abs/2604.07325v1)

!!! note "AI 核心解读"

    本文提出了一种针对时间序列数据的顺序保形化密度区域方法，通过分位数随机森林进行保形调整，在非交换数据下实现了渐近条件覆盖率的理论保证。该方法具有双重稳健性：只要预测密度模型设定正确或评分函数遵循正确阶数的非线性自回归模型，均可确保覆盖有效性。

