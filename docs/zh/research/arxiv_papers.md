# 🛰️ Research Frontier

这里是我通过GitHub机器人脚本制作的论文速递版面，它北京时间每天早上八点会自动扫描ArXiv上的关于Knockoff与Conformal Prediction的统计学论文。

> 更新于: 2026-04-21

---
### - [ ] Trustworthy Feature Importance Avoids Unrestricted Permutations 
- **分类**: Knockoff | **日期**: 2026-04-13
- **链接**: [PDF](http://arxiv.org/abs/2604.11253v1)
!!! note "AI 核心解读"
    该论文的创新点在于：  
    1. 揭示了传统特征重要性方法因无限制置换导致外推误差的根本缺陷，并提出三种新方法（条件模型依赖、高斯变换Knockoffs、受限ALE图设计）从理论上解决该问题。  
    2. 通过理论与数值实验证明，所提策略能有效减少或消除外推误差，从而提升特征重要性评估的可靠性。

### - [ ] xplainfi: Feature Importance and Statistical Inference for Machine Learning in R 
- **分类**: Knockoff | **日期**: 2026-03-16
- **链接**: [PDF](http://arxiv.org/abs/2603.15306v1)
!!! note "AI 核心解读"
    该论文的创新点在于：  
    1. 推出了首个在R语言中集成**条件特征重要性分析**与**统计推断框架**的综合性工具包，填补了现有方法在条件性采样与统计检验方面的空白。  
    2. 提供了**模块化的条件采样架构**（如对抗随机森林、基于Knockoff的采样器等），支持连续与混合数据的可解释性分析，并实现了多种统计推断方法（如方差校正置信区间）。

### - [ ] Distribution-free screening of spatially variable genes in spatial transcriptomics 
- **分类**: Knockoff | **日期**: 2026-03-10
- **链接**: [PDF](http://arxiv.org/abs/2603.09061v1)
!!! note "AI 核心解读"
    本论文的创新点在于：提出了一种基于新型准似然比统计量（MM检验）的无分布空间可变基因筛选方法，该方法结合Knockoff程序控制错误发现率，并首次系统适用于三维空间转录组数据；同时，该方法在理论层面建立了选择一致性、错误发现率控制及选择后聚类误差界的统计保证。

### - [ ] Variable selection via knockoffs for clustered data 
- **分类**: Knockoff | **日期**: 2026-02-23
- **链接**: [PDF](http://arxiv.org/abs/2602.19398v1)
!!! note "AI 核心解读"
    1. 针对聚类数据（如重复测量数据）中观测水平与群组水平变量混杂的问题，创新性地提出将原始预测变量分解为群组均值与组内偏差两部分，并分别在两个层次上独立进行变量选择。  
    2. 通过模拟研究验证了该方法在控制错误发现率与提高检验功效上的有效性，尤其证明分层处理的knockoff方法优于直接使用全数据集的传统方法（如Lasso）。

### - [ ] Improving the adjusted Benjamini--Hochberg method using e-values in knockoff-assisted variable selection 
- **分类**: Knockoff | **日期**: 2026-02-12
- **链接**: [PDF](http://arxiv.org/abs/2602.11610v1)
!!! note "AI 核心解读"
    1. 将Sarkar和Tang（2022）的knockoff辅助变量选择方法重新解释为一种未归一化的e值加权Benjamini-Hochberg程序，并在此基础上引入有界p值到e值的校准器，实现了更精细灵活的权重分配。  
    2. 提出了三种基于e值加权的FDR控制程序，在仿真和真实数据（如HIV-1耐药性分析）中均显示出更优的检验功效，尤其在低目标FDR、信号稀少且微弱的场景下性能提升显著。

### - [ ] Conformal Robust Set Estimation 
- **分类**: Conformal | **日期**: 2026-04-20
- **链接**: [PDF](http://arxiv.org/abs/2604.18441v1)
!!! note "AI 核心解读"
    本论文的创新点在于：提出了一种基于半质量半径（即点到其第k近邻距离）的非共形评分方法，构建了具有边际有效性的鲁棒共形预测集；该方法在任意样本量下均能保证覆盖度，且通过指数浓度界限严格量化了经验区域与总体稳健中心集之间的概率偏差。

### - [ ] Online Conformal Prediction with Adversarial Semi-bandit Feedback via Regret Minimization 
- **分类**: Conformal | **日期**: 2026-04-20
- **链接**: [PDF](http://arxiv.org/abs/2604.17984v1)
!!! note "AI 核心解读"
    该论文的创新点在于：  
    1. 首次将在线共形预测的**部分反馈**问题（仅当真实标签落入预测集时才被观测）建模为**对抗性半赌博机问题**，通过将候选预测集视为“臂”并利用对抗性赌博机算法实现长期覆盖率保证。  
    2. 通过**显式建立学习器遗憾与覆盖率之间的理论联系**，在无需分布假设的对抗性环境中实现了对错误覆盖率的控制，并在独立同分布与非独立同分布场景中验证了方法的有效性。


