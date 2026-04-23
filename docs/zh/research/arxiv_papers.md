# 🛰️ Research Frontier

这里是我通过 GitHub 机器人脚本制作的论文速递版面，它北京时间每天早上八点会自动扫描 ArXiv 上的关于 Knockoff 与 Conformal Prediction 的统计学论文。


> 更新于: 2026-04-23

---
### Power of masking methods for adaptive testing in a multivariate normal means problem 

- [ ] **分类**: Knockoff | **日期**: 2026-01-12
- **链接**: [PDF](http://arxiv.org/abs/2601.07764v2)

!!! note "AI 核心解读"

    本文在多元正态均值模型中，通过渐近推导统一分析了三种掩蔽方法的检验功效：样本分割、零样本增广及理想样本内基准。主要创新在于严格证明增广法在匹配调参下优于分割法，且其功效最优的零样本数量为测试数量的可忽略比例，此时功效逼近理想基准；进一步推导出增广法的近似最优零样本量随测试量平方根增长，揭示了掩蔽策略在控制第一类错误与提升功效间的量化权衡。

### Improve Power of Knockoffs with Annotation Information of Covariates 

- [ ] **分类**: Knockoff | **日期**: 2026-01-05
- **链接**: [PDF](http://arxiv.org/abs/2601.02583v1)

!!! note "AI 核心解读"

    本文创新点在于：提出了一种基于注释信息的Knockoff方法（AnnoKn），通过将Knockoff过程与自适应Lasso回归结合，在统一的贝叶斯框架下整合功能注释信息，以提升变量选择的统计功效；同时，该方法进一步拓展至仅需汇总统计量的场景，在严格控制错误发现率的前提下，实现了比现有方法更强的因果变量检测能力。

### A General Stability Approach to False Discovery Rate Control 

- [ ] **分类**: Knockoff | **日期**: 2025-12-19
- **链接**: [PDF](http://arxiv.org/abs/2512.17401v1)

!!! note "AI 核心解读"

    本文提出一种基于多次运行基础FDR控制算法并聚合特征重要性排序的稳定性框架，通过构造具有理论保证的“稳定化松弛e值”并应用e-BH程序，在有限样本下严格控制FDR且渐近保持检验功效；同时证明了随着重复次数增加，所选特征集会收敛到确定性极限，从而在理论层面确保了方法的稳定性和可重复性。

### StaRQR-K: False Discovery Rate Controlled Regional Quantile Regression 

- [ ] **分类**: Knockoff | **日期**: 2025-11-26
- **链接**: [PDF](http://arxiv.org/abs/2511.21562v1)

!!! note "AI 核心解读"

    本文创新点在于：提出了一种结合区域分位数回归与模型X敲除技术的稳定化框架，通过构建基于截尾化处理的模型X敲除筛选器，实现了对超高维数据中区域分位数效应的错误发现率严格控制；同时设计了高效的区域分位数确定独立筛选流程，在保证理论上的错误发现率可控性的前提下，显著提升了检测尾部敏感效应与异质性关联的统计功效。

### When Features Beat Noise: A Feature Selection Technique Through Noise-Based Hypothesis Testing 

- [ ] **分类**: Knockoff | **日期**: 2025-11-25
- **链接**: [PDF](http://arxiv.org/abs/2511.20851v2)

!!! note "AI 核心解读"

    本文提出一种基于噪声比较的特征选择方法，其核心创新点在于：通过引入多个随机噪声特征，并构建以噪声特征最大重要性值为参照的非参数自助法假设检验框架，为特征筛选提供了统计严格的停止准则与显著性评估依据。该方法在统计推导上确保了算法设计的理论严谨性，通过模拟与真实数据验证，其在信号恢复能力上优于Boruta等对比方法，实现了理论性质与实用效能的统一。

### Conformal Margin Risk Minimization: An Envelope Framework for Robust Learning under Label Noise 

- [ ] **分类**: Conformal | **日期**: 2026-04-07
- **链接**: [PDF](http://arxiv.org/abs/2604.06468v2)

!!! note "AI 核心解读"

    本文提出了一种基于置信度边际的保形正则化框架，通过估计每批样本的保形分位数阈值，动态筛选高置信度样本并抑制可能误标的样本，从而提升任意分类损失在标签噪声下的鲁棒性。理论方面，该工作仅需边际分布满足温和的正则性条件，便推导出适用于任意标签噪声的学习误差界，其构造的评分函数无需先验噪声信息或训练流程修改，具有普适的稳定性。

### Weighted Bayesian Conformal Prediction 

- [ ] **分类**: Conformal | **日期**: 2026-04-07
- **链接**: [PDF](http://arxiv.org/abs/2604.06464v1)

!!! note "AI 核心解读"

    本文提出加权贝叶斯共形预测，通过将均匀狄利克雷先验替换为以有效样本量加权的狄利克雷分布，将贝叶斯求积框架推广至非独立同分布场景；理论证明了该加权参数是唯一能使频率学派与贝叶斯方差匹配的浓度参数，并推导出后验标准差以有效样本量平方根速率衰减，同时扩展了数据依条件随机支配性保证。

### Debiased Machine Learning for Conformal Prediction of Counterfactual Outcomes Under Runtime Confounding 

- [ ] **分类**: Conformal | **日期**: 2026-04-04
- **链接**: [PDF](http://arxiv.org/abs/2604.03772v1)

!!! note "AI 核心解读"

    本文提出了一种基于半参数效率理论的去偏机器学习框架，用于处理目标群体中仅观测到部分混杂变量（即运行时混杂）的反事实结果预测问题。其核心创新在于构造了具有双重稳健性质的评分函数，通过渐进线性估计量实现覆盖率的理论保证，并在有限样本下相比传统方法具有更快的收敛速度。

### Online Reasoning Calibration: Test-Time Training Enables Generalizable Conformal LLM Reasoning 

- [ ] **分类**: Conformal | **日期**: 2026-04-01
- **链接**: [PDF](http://arxiv.org/abs/2604.01170v1)

!!! note "AI 核心解读"

    本文提出一种基于元学习的在线校准模块，通过测试时训练动态调整置信度，利用保形预测理论为分布漂移下的推理过程提供风险可控的统计保证；所构造的评分函数融合了监督与自洽标签，在理论层面确保校准误差的同时，显著提升了多类推理任务中的采样效率与泛化性能。

### Conformal Prediction Assessment: A Framework for Conditional Coverage Evaluation and Selection 

- [ ] **分类**: Conformal | **日期**: 2026-03-28
- **链接**: [PDF](http://arxiv.org/abs/2603.27189v1)

!!! note "AI 核心解读"

    本文提出了一种将条件覆盖评估转化为监督学习任务的框架，通过训练可靠性估计器来预测实例级覆盖概率，并基于此构造了可分解为安全性（欠覆盖风险）与效率（过覆盖成本）的条件有效性指标。该框架为可靠性估计器建立了收敛速率，并证明了基于条件有效性指标的模型选择方法具有一致性。


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

