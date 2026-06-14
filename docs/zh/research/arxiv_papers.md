# 🛰️ Research Frontier

这是我通过脚本制作的论文速递版面，它将在北京时间每天早上八点调用 GitHub 机器人自动扫描 [arXiv](https://arxiv.org/) 上的关于 Knockoff 与 Conformal Prediction 的统计学论文。


> 更新于: 2026-06-14

---
### Estimating the local false discovery rate under an unknown symmetric null 

- [ ] **分类**: Knockoff | **日期**: 2026-06-10
- **链接**: [PDF](http://arxiv.org/abs/2606.11548v1)

!!! note "AI 核心解读"

    该论文提出在零分布仅已知对称性（关于零点对称）的条件下，通过估计代理密度比 \( f(-w)/f(w) \) 来构造局部错误发现率（lfdr）的评分函数，并证明任何对该比值的相合估计都能在名义水平上实现 lfdr 的渐近控制。其统计推导核心在于利用自然三次样条基函数的逻辑回归方法估计该密度比，从而绕开对零分布具体形式的依赖，并建立了估计量与多重检验过程渐近有效性的理论联系。

### Online Shift Detection and Conformal Adaptation for Deployed Safety Classifiers 

- [ ] **分类**: Conformal | **日期**: 2026-06-10
- **链接**: [PDF](http://arxiv.org/abs/2606.11949v1)

!!! note "AI 核心解读"

    该论文创新点在于：1）构建了基于校准序贯统计量的在线漂移检测框架，通过加权共形预测的密度比估计实现分布外检测，并理论证明了高维嵌入空间中逻辑密度比估计导致重要性权重截断的统计失效机制；2）提出自适应共形弃权层动态调整决策阈值以恢复目标错误率，并通过方差分解（分类器η²=0.243、漂移类型η²=0.237及交互效应η²=0.185）揭示了检测延迟的统计异质性来源。

### Enhanced localized conformal prediction with imperfect auxiliary information 

- [ ] **分类**: Conformal | **日期**: 2026-06-07
- **链接**: [PDF](http://arxiv.org/abs/2606.08551v1)

!!! note "AI 核心解读"

    该论文提出增强型局部共形预测（ELCP），通过构造密度比加权核估计量，将辅助数据与校准数据融合，在保持有限样本边际覆盖保证的同时，改进了渐近测试条件覆盖的理论性质。其核心创新在于利用统计推导的评分函数调整局部预测集，有效缓解了稀疏校准数据区域中局部覆盖不可靠的问题。

### Conformal Risk-Averse Decision Making with Action Conditional Guarantee 

- [ ] **分类**: Conformal | **日期**: 2026-06-04
- **链接**: [PDF](http://arxiv.org/abs/2606.05551v2)

!!! note "AI 核心解读"

    该论文通过引入动作条件共形预测，将边际安全保证强化为对每个决策动作的显式条件保证，并基于分位数损失最小化构造了有限样本下的评分函数，从而在理论上证明了动作条件预测集可作为风险厌恶决策中条件风险价值优化的可行决策空间代理。

### Conformal Risk Prediction for Non-Alcoholic Fatty Liver Disease Using Gradient Boosting with Distribution-Free Coverages 

- [ ] **分类**: Conformal | **日期**: 2026-05-31
- **链接**: [PDF](http://arxiv.org/abs/2606.09860v1)

!!! note "AI 核心解读"

    该论文创新点在于：将梯度提升树与共形预测结合，构造了无需分布假设的个体风险预测评分函数，并严格证明了其预测集在有限样本下边际覆盖率可超过预设置信水平；同时，基于互信息的稳定性选择过程通过自助重采样筛选出紧凑特征子集，确保了评分函数的统计可解释性与理论覆盖保证。

### Parameter-Free and Group Conditional Online Conformal Prediction 

- [ ] **分类**: Conformal | **日期**: 2026-05-29
- **链接**: [PDF](http://arxiv.org/abs/2606.00419v3)

!!! note "AI 核心解读"

    该论文在统计推导上提出了一种无需学习率调参的在线共形预测算法，并严格证明了其能在数据分布非平稳变化下达到最优的分组条件覆盖保证。其核心创新在于构造了一个自适应的评分函数更新规则，使得算法无需依赖任何先验参数即可同时控制各组预测区间的边际与条件覆盖误差。


### Knockoffs-based False Discovery Rate Control and Simplification for Deep Neural Networks 

- [ ] **分类**: Knockoff | **日期**: 2026-06-03
- **链接**: [PDF](http://arxiv.org/abs/2606.04404v1)

!!! note "AI 核心解读"

    该论文将knockoffs框架拓展至深度神经网络，通过构造基于正则化网络权重的评分函数，在控制错误发现率的前提下提出了三种变量筛选方法；理论性质上，证明了所提方法在有限样本下对FDR的控制能力，并推导了单层与多层筛选过程中统计量的渐近分布。

### Conformal Prediction for Dyadic Regression Under Complex Missingness 

- [ ] **分类**: Conformal | **日期**: 2026-06-09
- **链接**: [PDF](http://arxiv.org/abs/2606.11136v1)

!!! note "AI 核心解读"

    该论文通过构造保测度双射论证，在弱于可交换性的分布不变条件下证明了共形预测的超均匀性，并针对索引集随机子样本这一现有理论未覆盖的情形建立了严格统计推导。其核心创新在于为联合可交换数组设计了行-列相似性评分函数与图模型加权共形程序，首次在非随机缺失假设下证明了加权共形预测的渐近条件有效性。

### Generalized Conformal Predictive Systems Under Distributional Shifts 

- [ ] **分类**: Conformal | **日期**: 2026-06-09
- **链接**: [PDF](http://arxiv.org/abs/2606.11044v1)

!!! note "AI 核心解读"

    该论文通过引入观测特定的置换权重，将广义共形预测系统从可交换性假设推广至非可交换环境，并基于条件加权抽样构造了具有有限样本或渐近置信保证的鲁棒预测包络。其核心创新在于推导了权重不确定性箱的统计性质，并针对共形度量、分箱及保序分布回归设计了高效计算方法，使得预测带在分布偏移下保持校准性。

### SPACR: Single-Pass Adaptive Training of Uncertainty-Aware Conformal Regressors 

- [ ] **分类**: Conformal | **日期**: 2026-06-09
- **链接**: [PDF](http://arxiv.org/abs/2606.10734v1)

!!! note "AI 核心解读"

    SPACR 提出了一种可微分的单次训练损失函数，该损失函数直接优化共形预测区间的效率（即区间宽度）与边际有效性，无需在训练时预设置信水平或进行批次分割。其核心创新在于构造了一个基于排序的评分函数，使得单次训练后的模型能在推理时同时输出多个置信水平下的有效预测区间，并从理论上保证了区间覆盖率的边际控制。

### Predicting Current Outcomes From Historical Survey Data With Weighted Conformal Prediction 

- [ ] **分类**: Conformal | **日期**: 2026-06-09
- **链接**: [PDF](http://arxiv.org/abs/2606.10563v1)

!!! note "AI 核心解读"

    该论文创新点在于：针对历史调查数据与目标总体间的协变量偏移，构造了结合子组密度比与子组比例估计的加权评分函数，并基于此建立了加权共形预测框架，从理论上证明了预测集覆盖率的有限样本保证。

### Data augmented bootstrap: Unifying confidence interval construction by approximate invariance 

- [ ] **分类**: Conformal | **日期**: 2026-06-08
- **链接**: [PDF](http://arxiv.org/abs/2606.09049v1)

!!! note "AI 核心解读"

    该论文提出数据增强自助法（DAB），通过将数据集的近似不变性量化为Kolmogorov距离，并利用满足高斯普适性的统计量简化为条件均值与方差匹配，从而统一了置信区间的构造框架。其核心创新在于：无需假设群结构，即可在有限样本与渐近保证之间插值得到理论覆盖结果，并系统地将机器学习中的数据增强启发式方法融入经典自助法、野自助法和共形预测等统计方法中。


### Conformal Language Modeling via Posterior Sampling 

- [ ] **分类**: Conformal | **日期**: 2026-06-02
- **链接**: [PDF](http://arxiv.org/abs/2606.03731v1)

!!! note "AI 核心解读"

    该论文提出了一种基于后验采样的共形语言建模方法，通过构造一个校准的高分区域作为条件事件，并设计针对条件序列生成的校准程序，实现了对目标风险的有效控制。其核心创新在于将共形预测的统计保证直接嵌入生成过程，而非事后过滤，从而在保持理论风险控制的同时提升了生成样本的连贯性与模型内在一致性。

### Set-Preserving Calibration from Conformal P-Values to E-Values 

- [ ] **分类**: Conformal | **日期**: 2026-06-02
- **链接**: [PDF](http://arxiv.org/abs/2606.03600v1)

!!! note "AI 核心解读"

    该论文提出了一种新颖的P2E校准器，通过构造一个保序且保集合的评分函数，将共形p值转化为e值，并严格证明了该转化不改变原始p值诱导的预测集，从而克服了经典p值到e值校准器在共形预测中导致预测集过度保守的统计缺陷。在理论性质上，作者建立了该转化在合并依赖证据时的效率优势，并证明了基于e值的交叉共形预测和共形聚合方法能够严格保证目标覆盖概率，同时显著提升统计效率。

### Parameter-Free and Group Conditional Online Conformal Prediction 

- [ ] **分类**: Conformal | **日期**: 2026-05-29
- **链接**: [PDF](http://arxiv.org/abs/2606.00419v2)

!!! note "AI 核心解读"

    该论文在统计推导上提出了一种无需学习率调参的在线共形预测算法，并严格证明了其能在数据分布非交换的对抗性环境下实现最优的分组条件覆盖保证；其核心创新在于构造了一种自适应的评分函数更新规则，使得分组误差控制与参数无关的在线学习理论得以统一，从而在理论上确保了分组覆盖率的长期累积界。

### Approximating full conformal prediction: distribution free guarantees via the tournament correction 

- [ ] **分类**: Conformal | **日期**: 2026-05-28
- **链接**: [PDF](http://arxiv.org/abs/2605.29200v3)

!!! note "AI 核心解读"

    该论文提出了一种基于“锦标赛校正”的通用近似框架，通过构造新的评分函数和引入双重覆盖校正机制，在无需数据分割的条件下，严格证明了预测集可实现至少 \(1-2\alpha\) 的边际覆盖保证，并在稳定性条件下将理论保证收紧至接近 \(1-\alpha\)。其统计推导统一并推广了留一交叉共形预测方法，为各类近似策略提供了分布自由的理论保障。

### A Unified Theory of Conditional Coverage in Conformal Prediction with Applications 

- [ ] **分类**: Conformal | **日期**: 2026-05-12
- **链接**: [PDF](http://arxiv.org/abs/2605.11602v3)

!!! note "AI 核心解读"

    该论文的核心创新在于提出了一种非渐近的条件覆盖误差分解框架，将条件覆盖不足拆解为评分函数估计误差、有限样本校准误差与内在条件失配误差三个可解释分量。基于此分解，作者构造了面向条件覆盖的模型选择准则，并发展了在协变量偏移下具有渐近条件覆盖保证的局部化评分函数方法。


### A Distribution-Free Framework for Rewrite-Based Human-text Detection via Knockoff Filtering 

- [ ] **分类**: Knockoff | **日期**: 2026-05-29
- **链接**: [PDF](http://arxiv.org/abs/2606.00402v1)

!!! note "AI 核心解读"

    该论文通过将改写式检测隐式构造的knockoff样本转化为多重假设检验框架，提出了一种无需重训练即可为任意改写检测器赋予有限样本FDR保证的分布自由方法。其核心创新在于分离了检测统计量的设计与错误发现控制，利用knockoff滤波的统计性质，通过简单校准过程使现有改写检测器继承严格的FDR理论保证。

### Network knockoffs: controlling false discovery in dyadic space 

- [ ] **分类**: Knockoff | **日期**: 2026-05-29
- **链接**: [PDF](http://arxiv.org/abs/2606.00346v1)

!!! note "AI 核心解读"

    该论文创新点在于：针对高维二元回归中因数据依赖导致传统伪发现率控制失效的问题，直接在拓扑网络上模拟合成特征，再构造二元空间中的增广设计矩阵，从而提出一种新型的Knockoff变量选择程序。理论性质上，该方法在节点级和边级特征上均能有效控制伪发现率，而Benjamini-Hochberg、Benjamini-Yekutieli、Storey Q值、数据拆分及标准Knockoff等方法均表现为反保守。

### Multi-Agent Conformal Prediction with Personalized Statistical Validity 

- [ ] **分类**: Conformal | **日期**: 2026-05-30
- **链接**: [PDF](http://arxiv.org/abs/2606.00717v1)

!!! note "AI 核心解读"

    该论文提出个性化联邦加权共形预测框架，通过构造基于局部密度比加权的评分函数与加权分位数聚合，在异质性数据下为每个智能体提供渐近有效的边际和校准条件覆盖保证；理论推导揭示了由有效样本量表达式调控的覆盖方差调整机制，这是加权共形预测中必要的统计性质。

### Parameter-Free and Group Conditional Online Conformal Prediction 

- [ ] **分类**: Conformal | **日期**: 2026-05-29
- **链接**: [PDF](http://arxiv.org/abs/2606.00419v1)

!!! note "AI 核心解读"

    该论文在统计推导上首次将无参数在线优化框架与分组条件覆盖相结合，通过构造基于梯度下降的评分函数自适应更新机制，在不依赖学习率调参的前提下实现了对任意数据漂移的鲁棒性。其核心理论性质在于证明了所提算法能够同时达到最优的分组条件覆盖保证，且分组误差界与最优后验分组阈值匹配，从而在非交换数据流中统一了公平性与参数无关性。

### Approximating full conformal prediction: distribution free guarantees via the tournament correction 

- [ ] **分类**: Conformal | **日期**: 2026-05-28
- **链接**: [PDF](http://arxiv.org/abs/2605.29200v2)

!!! note "AI 核心解读"

    该论文提出基于“锦标赛校正”的通用近似框架，将全共形预测的计算复杂度转化为可处理的统计推断问题，并严格证明了在任意分布下预测集边际覆盖率为 \(1-2\alpha\)，且在稳定性条件下可收紧至约 \(1-\alpha\)。其核心创新在于构造了一种新的评分函数校正机制，通过引入锦标赛式比较与留一交叉共形预测的泛化，在无需数据分割的前提下保留了全共形预测的理论效率与分布自由保证。


### Explaining Concept Shift with Interpretable Feature Attribution 

- [ ] **分类**: Knockoff | **日期**: 2025-05-27
- **链接**: [PDF](http://arxiv.org/abs/2505.20634v2)

!!! note "AI 核心解读"

    该论文将概念漂移下的性能退化归因问题转化为稀疏特征选择任务，并基于广义加性模型、knockoffs 和吸收统计量构造了评分函数，从而在统计上识别出导致标签条件分布变化的特征子集。理论性质方面，该方法在有限样本下能准确区分真正漂移特征与无关特征，且对复杂漂移模式具有统计鲁棒性。

### Leave a Window Out: Modifying the Jackknife for Predictive Inference in Time Series 

- [ ] **分类**: Conformal | **日期**: 2026-05-28
- **链接**: [PDF](http://arxiv.org/abs/2605.30292v1)

!!! note "AI 核心解读"

    该论文提出了一种针对时间序列数据的“留窗法”（leave-a-window-out）对经典留一法进行修正，通过量化数据与循环可交换性之间的偏离程度并引入新的度量系数，在模型拟合满足温和稳定性条件下证明了预测区间的有效覆盖性。其核心创新在于构造了基于窗口剔除的评分函数，从理论上推导出该方法能克服传统留一法在弱时间依赖下覆盖度任意丢失的缺陷，同时避免了数据分割带来的效率损失。

### Approximating full conformal prediction: distribution free guarantees via the tournament correction 

- [ ] **分类**: Conformal | **日期**: 2026-05-28
- **链接**: [PDF](http://arxiv.org/abs/2605.29200v1)

!!! note "AI 核心解读"

    该论文提出了一种基于“锦标赛校正”的新型全共形预测近似框架，通过构造一个与原始评分函数耦合的校正项，在无需数据分割的条件下，严格证明了预测集边际覆盖率为 \(1-2\alpha\)，并在稳定性条件下将理论保证收紧至约 \(1-\alpha\)。其核心创新在于将留一交叉共形预测推广为更一般的评分函数构造范式，使得任意近似策略都能在统计推导上获得可验证的分布自由覆盖保证。

### Conf-Gen: Conformal Uncertainty Quantification for Generative Models 

- [ ] **分类**: Conformal | **日期**: 2026-05-27
- **链接**: [PDF](http://arxiv.org/abs/2605.28920v1)

!!! note "AI 核心解读"

    该论文通过放宽传统共形风险控制（CRC）的理论假设，将其推广至无监督生成模型，并构造了适用于生成任务的通用评分函数，从而在无需修改模型结构的前提下，为图像生成、对话系统及AI代理等场景提供了严格的统计不确定性量化保证。其核心创新在于统一并扩展了先前针对大语言模型的共形预测方法，并首次在理论上证明了生成模型输出在非记忆性、充分性与正确性等维度上的共形风险控制性质。


### Mean-Shift PCA by Knockoff Mean 

- [ ] **分类**: Knockoff | **日期**: 2026-05-25
- **链接**: [PDF](http://arxiv.org/abs/2605.25460v1)

!!! note "AI 核心解读"

    该论文通过引入敲击均值扰动，利用随机矩阵理论严格证明了均值偏移尖峰与原始协方差稳定特征值在谱上可分离，并揭示了原始特征空间在污染下具有渐近不变性。基于此理论性质，构造了一种仅依赖标准PCA操作的两阶段评分函数，通过添加敲击均值识别并移除均值偏移成分。

### A Post-Processing Conformal Prediction Approach for Conditional Coverage via Pivotal Scores 

- [ ] **分类**: Conformal | **日期**: 2026-05-25
- **链接**: [PDF](http://arxiv.org/abs/2605.25852v2)

!!! note "AI 核心解读"

    该论文的核心创新在于将条件覆盖问题等价转化为构造特征独立的非一致性评分，并基于此提出一种后处理校正方法PIT-CP，通过概率积分变换将任意基础评分映射为近似特征不变的评分，同时保持其几何结构和边际覆盖。在理论性质上，作者将问题简化为诱导评分的一维条件密度估计，推导了条件覆盖间隙的边界以及体积和对称差界，并证明了该方法在极小化最优条件估计下的有效性。

### Distributional Conformal Prediction for Markov Processes 

- [ ] **分类**: Conformal | **日期**: 2026-05-24
- **链接**: [PDF](http://arxiv.org/abs/2605.24848v1)

!!! note "AI 核心解读"

    该论文将分布共形预测从独立同分布回归场景推广至严格平稳马尔可夫过程，通过基于核估计的转移分布函数实施概率积分变换，将相依数据转化为近似独立同分布样本。理论贡献在于：在β-混合条件下推导出无条件覆盖率的非渐近误差界，并验证了条件预测区间的渐近有效性，同时将有效性条件放宽至Lp-m-可逼近过程，无需依赖混合性假设。

### Bayesian Conformal-Projective Prediction 

- [ ] **分类**: Conformal | **日期**: 2026-05-23
- **链接**: [PDF](http://arxiv.org/abs/2605.24601v1)

!!! note "AI 核心解读"

    该论文提出了一种名为共形-投影预测（CPP）的稳健预测框架，其核心创新在于通过定义分布层面的共形性准则，将候选值纳入数据后对留一预测分布的无扰动程度作为评分函数，而非依赖传统的残差评分。在留一与交换预测分布可解析且交换预测均值对候选值可微的条件下，论文推导了有界影响性质和局部凸性引理，并严格证明了在ε-污染模型下，CPP在渐近方差上优于任何具有无界影响的插件预测器。

### A Unified Theory of Conditional Coverage in Conformal Prediction with Applications 

- [ ] **分类**: Conformal | **日期**: 2026-05-12
- **链接**: [PDF](http://arxiv.org/abs/2605.11602v2)

!!! note "AI 核心解读"

    该论文的核心创新在于提出了一个非渐近的条件覆盖误差分解框架，将条件覆盖不足拆解为评分估计误差、有限样本校准误差和内在条件失配误差三个可解释分量。基于此框架，作者进一步构造了面向条件覆盖的模型选择准则，并开发了在协变量偏移下仍保持渐近条件覆盖保证的局部化方法。


### Selecting Informative Conformal Prediction Sets with an Optimized FCR-Controlled Approach 

- [ ] **分类**: Conformal | **日期**: 2026-05-21
- **链接**: [PDF](http://arxiv.org/abs/2605.22004v1)

!!! note "AI 核心解读"

    该论文在选择性推断框架下，推导了最优决策策略以最大化统计功效，并针对预测集概率未知的实际情形，构造了一种校准程序，使得调整后的策略在有限样本下严格控制错误覆盖率（FCR）。其核心创新在于将评分函数与单调约束条件相结合，通过理论推导保证了所选样本上FCR的精确控制，同时显著提升了统计效率。

### Online Conformal Prediction for Non-Exchangeable Panel Data 

- [ ] **分类**: Conformal | **日期**: 2026-05-18
- **链接**: [PDF](http://arxiv.org/abs/2605.17705v1)

!!! note "AI 核心解读"

    该论文针对非可交换面板数据，提出了一种在线共形预测框架，其核心创新在于构造了一个结合历史相似性权重与自适应误覆盖水平的双状态评分函数，并基于此推导出逐步覆盖界与长期覆盖保证的理论性质。通过统计推导证明，相似性权重在目标反馈稀疏时保障覆盖，而自适应水平随反馈累积进一步优化覆盖，二者互补实现了非均匀区间宽度分配，而非简单的整体膨胀。

### Skew-adaptive conformal prediction 

- [ ] **分类**: Conformal | **日期**: 2026-05-15
- **链接**: [PDF](http://arxiv.org/abs/2605.16145v1)

!!! note "AI 核心解读"

    该论文提出一种基于非对称区间族与规范评分函数的偏态自适应共形预测方法，通过逆双曲正弦变换构造有符号缩放残差的训练目标，使预测区间在保持有限样本边际有效性的同时，能够自适应特征空间中的局部尺度与偏度变化。其核心理论贡献在于推导了偏态自适应评分函数的显式形式，并建立了基于校准样本的期望区间宽度比估计量，从而在统计上保证了区间效率的改进可被可靠评估。

### Learning Context-conditioned Gaussian Overbounds for Convolution-Based Uncertainty Propagation 

- [ ] **分类**: Conformal | **日期**: 2026-05-15
- **链接**: [PDF](http://arxiv.org/abs/2605.15789v1)

!!! note "AI 核心解读"

    该论文提出一种可学习的上下文条件高斯过界框架，通过构造一个基于有限分位数网格的保守性损失函数，并引入Wasserstein型分布距离惩罚项，使得神经网络输出的均值和尺度参数在选定分位数上具有可证明的保守性。在三条显式正则性假设下，该方法进一步将离散网格上的保守性延拓至连续尾部，从而支持保守的线性组合与卷积传播分析，并在理论上保证了过界性质的连续性与目标函数的紧致性。

### On the Burden of Achieving Fairness in Conformal Prediction 

- [ ] **分类**: Conformal | **日期**: 2026-05-14
- **链接**: [PDF](http://arxiv.org/abs/2605.14260v2)

!!! note "AI 核心解读"

    该论文首先从总体分数分布出发，推导出一个守恒律与下界，证明混合校准必然导致不可消除的组间覆盖率扭曲，其规模由跨组分位数异质性决定；其次，通过构造评分函数与理论分析，揭示了等覆盖率和等集合大小两种公平性定义之间存在根本性权衡，并量化了分组与混合校准策略间的转换代价。


### CASCADE Conformal Prediction: Uncertainty-Adaptive Prediction Intervals for Two-Stage Clinical Decision Support 

- [ ] **分类**: Conformal | **日期**: 2026-05-19
- **链接**: [PDF](http://arxiv.org/abs/2605.20468v1)

!!! note "AI 核心解读"

    该论文提出了一种名为CASCADE的级联共形预测框架，通过将Venn-Abers多概率不确定性直接映射为非一致性分数，在统计上实现了从分类任务到回归任务的不确定性传播。其核心创新在于利用主分类器的认知不确定性动态缩放回归任务的预测区间，理论性质上保证了在置信样本上区间效率提升38.9%，同时自动扩展区间以维持不确定样本的稳健覆盖。

### Conformal Prediction via Transported Beta Laws 

- [ ] **分类**: Conformal | **日期**: 2026-05-18
- **链接**: [PDF](http://arxiv.org/abs/2605.19024v1)

!!! note "AI 核心解读"

    该论文将分裂保形预测中校准条件覆盖率的精确分布推导为Beta(k, n+1-k)形式，并以此为有限样本参考对象，通过Wasserstein距离量化非独立同分布场景下的分布偏移。其核心创新在于构造了基于传输映射的评分函数变形机制，将测试端偏移与校准依赖分别建模为覆盖尺度上的传输映射与顺序统计量律的变形，并在尺度偏移、聚类和稳态混合等设定下给出了显式或Berry-Esseen近似的理论性质刻画。

### SAGA: A Sequence-Adaptive Generative Architecture for Multi-Horizon Probabilistic Forecasting with Adaptive Temporal Conformal Prediction 

- [ ] **分类**: Conformal | **日期**: 2026-05-18
- **链接**: [PDF](http://arxiv.org/abs/2605.19014v1)

!!! note "AI 核心解读"

    该论文提出了一种基于解码器-仅Transformer的序列自适应生成架构SAGA，并创新性地结合了分裂共形校准包装器，从而在有限样本下为个体级预测区间提供了边际覆盖保证。其统计贡献在于通过自适应时间共形预测构造了无需分布假设的评分函数，并严格证明了在面板数据长时序预测中，预测区间的名义覆盖率偏差可控制在0.4个百分点以内，最差子群偏差不超过2.4个百分点。

### Self-Supervised Conformal Prediction with Equivariant Bootstrapping for Image Uncertainty Quantification 

- [ ] **分类**: Conformal | **日期**: 2026-05-18
- **链接**: [PDF](http://arxiv.org/abs/2605.18655v1)

!!! note "AI 核心解读"

    该论文提出了一种自监督共形预测框架，通过等变自举法（equivariant bootstrapping）利用数据对称性生成启发式覆盖，并在此基础上进行共形校准，从而在无需真实校准数据的情况下实现不确定性量化。其核心创新在于构造了基于数据对称性的评分函数，并证明了经等变自举修正后的覆盖概率在共形校准下具有理论上的边际覆盖保证。

### Uncertainty Reliability Under Domain Shift: An Investigation for Data-Driven Blood Pressure Estimation in Photoplethysmography 

- [ ] **分类**: Conformal | **日期**: 2026-05-18
- **链接**: [PDF](http://arxiv.org/abs/2605.18008v1)

!!! note "AI 核心解读"

    该论文通过对比深度集成与蒙特卡洛丢弃法在域偏移下的预测鲁棒性，并系统评估高斯负对数似然损失与均方误差损失结合保形预测、温度缩放等后处理重校准方法对不确定性校准的影响，揭示了深度集成在外部域偏移下具有更强预测稳健性，且重校准后的高斯负对数似然方法在收缩压与舒张压估计中实现了最优的不确定性校准性能。


### On the Burden of Achieving Fairness in Conformal Prediction 

- [ ] **分类**: Conformal | **日期**: 2026-05-14
- **链接**: [PDF](http://arxiv.org/abs/2605.14260v1)

!!! note "AI 核心解读"

    该论文首先从总体评分分布出发，推导出一个守恒律和下界，证明混合校准必然导致组间覆盖率的不可约扭曲，其幅度由跨组分位数异质性决定。其次，论文严格证明了共形预测中两种主流公平性定义——等覆盖率和等集合大小——在统计上存在根本性冲突，并量化了分组与合并策略之间的权衡代价。


### Online Conformal Prediction: Enforcing monotonicity via Online Optimization 

- [ ] **分类**: Conformal | **日期**: 2026-05-12
- **链接**: [PDF](http://arxiv.org/abs/2605.12668v1)

!!! note "AI 核心解读"

    该论文通过在线优化框架构造了新的评分函数，将单调性约束直接嵌入到分位数估计的在线更新过程中，从而在多个置信水平上同时保证预测集的嵌套性与有限样本覆盖保证。其理论贡献在于推导了与分位数估计误差控制直接相关的遗憾界，并证明了所提方法在保持嵌套性的同时，能够实现跨风险谱系的一致统计效率提升。


### Multi-Variable Conformal Prediction: Optimizing Prediction Sets without Data Splitting 

- [ ] **分类**: Conformal | **日期**: 2026-05-12
- **链接**: [PDF](http://arxiv.org/abs/2605.12341v1)

!!! note "AI 核心解读"

    该论文将保形预测从标量评分函数推广到向量值评分函数，并引入多个同时校准变量，从而在无需数据分割的情况下，将预测集形状设计与校准统一为单一优化问题。其核心创新在于基于场景理论构建了带约束移除（RemMCP）与约束松弛（RelMCP）的两种计算框架，前者实现了分割保形预测的简洁推广，后者支持非凸评分函数，并严格证明了有限样本覆盖保证。

### A Unified Theory of Conditional Coverage in Conformal Prediction with Applications 

- [ ] **分类**: Conformal | **日期**: 2026-05-12
- **链接**: [PDF](http://arxiv.org/abs/2605.11602v1)

!!! note "AI 核心解读"

    该论文提出了一个统一的理论框架，通过点态路径和 \(L_p\) 路径分别推导出条件覆盖率的非渐近界，首次系统揭示了渐近条件有效性中误差来源的统计结构。其核心创新在于构造了基于加权对称性的评分函数，并建立了条件覆盖率的理论比较准则，为异质性数据下的共形预测提供了统一的统计推导基础。

### Multi-Fidelity Quantile Regression 

- [ ] **分类**: Conformal | **日期**: 2026-05-11
- **链接**: [PDF](http://arxiv.org/abs/2605.10406v1)

!!! note "AI 核心解读"

    该论文提出一种两阶段、模型无关的多保真分位数回归方法，核心创新在于构造了“局部分位数链接”：将高保真分位数表示为低保真分位数在协变量依赖水平上的函数，从而将估计问题转化为对更平滑的水平函数的统计推断。理论部分严格刻画了该估计量相较于仅用高保真数据直接分位数回归的收敛速度提升条件，并引入校正步骤以增强在形状相似性减弱时的鲁棒性。

### CONTRA: Conformal Prediction Region via Normalizing Flow Transformation 

- [ ] **分类**: Conformal | **日期**: 2026-05-08
- **链接**: [PDF](http://arxiv.org/abs/2605.08561v1)

!!! note "AI 核心解读"

    该论文提出 CONTRA 方法，通过将正态化流变换的潜空间中心距离定义为非一致性评分函数，在统计上保证了预测区域覆盖概率的有限样本有效性，并利用流模型的双射性质将潜空间高密度区域映射为输出空间任意形状的精确预测区域。其扩展版本进一步将评分函数构造为残差的正态化流变换，使得任意预测模型均可获得具有理论覆盖保证的预测区域，突破了传统超矩形或椭圆区域的形状限制。

### Rolling-Origin Conformal Prediction under Local Stationarity and Weak Dependence 

- [ ] **分类**: Conformal | **日期**: 2026-05-08
- **链接**: [PDF](http://arxiv.org/abs/2605.08422v1)

!!! note "AI 核心解读"

    该论文在局部平稳和弱依赖假设下，推导出滚动原点共形预测覆盖误差的四项分解，并证明最优校准窗口长度以 \(T^{2\beta/(2\beta+1)}\) 阶增长，覆盖误差率达到 \(O(T^{-\beta/(2\beta+1)})\)，且通过 Le Cam 两点构造证明该率在 Hölder-\(\beta\) 模型类上达到极小极大最优。此外，论文在 \(\alpha\)-混合和物理依赖框架下建立了 Bahadur 表示，并构造了基于 Winkler 交叉验证的自适应窗口选择器，通过均匀集中不等式给出其理论保证。


### TRACE: Transport Alignment Conformal Prediction via Diffusion and Flow Matching Models 

- [ ] **分类**: Conformal | **日期**: 2026-05-08
- **链接**: [PDF](http://arxiv.org/abs/2605.07100v1)

!!! note "AI 核心解读"

    该论文提出了一种基于扩散与流匹配模型中传输对齐的非一致性评分函数，通过沿随机传输轨迹平均去噪或速度匹配误差来度量候选输出与生成动力学的对齐程度，避免了显式似然评估和可逆变换的限制。理论分析表明，该标量评分函数在交换性假设下可通过分裂保形预测实现有限样本、无分布假设的边际覆盖保证，并进一步刻画了计算预算对评分统计性质的影响。


### History-Aware Conformal Prediction Sets for Censored Time-to-Event Outcomes 

- [ ] **分类**: Conformal | **日期**: 2026-05-07
- **链接**: [PDF](http://arxiv.org/abs/2605.06581v1)

!!! note "AI 核心解读"

    该论文提出历史感知共形预测集（HAPS），通过逆概率删失加权处理时变混杂因素下的右删失数据，并证明在删失权重一致估计时预测集在幸存者中具有渐近近似正确覆盖（PAAC）性质；进一步构造双重稳健扩展版本，放松了对删失分布正确建模的依赖，在统计推导上保证了覆盖率的稳健性。

### When Does Trimming Help Conformal Prediction? A Retained-Law Diagnostic under Calibration Contamination 

- [ ] **分类**: Conformal | **日期**: 2026-05-07
- **链接**: [PDF](http://arxiv.org/abs/2605.06204v1)

!!! note "AI 核心解读"

    该论文将固定阈值修剪重新定义为对校准分布的条件化操作，而非净化过程，并推导出修剪后保留分布下清洁目标覆盖率的精确有限样本恒等式，将其转化为一维评分累积分布函数转移问题。通过构建转移间隙的分量上界，分离出清洁侧协方差成本与保留污染成本，并证明修剪的有效性取决于异常评分能否在保持清洁群体评分中性条件下分离保留概率，否则无法通过保留混合系数实质降低污染。

### Socio-Conformal Calibration in Complex Survey Data: Marginal Validity Is Not Enough for Subgroup Reliability 

- [ ] **分类**: Conformal | **日期**: 2026-05-07
- **链接**: [PDF](http://arxiv.org/abs/2605.05562v1)

!!! note "AI 核心解读"

    该论文通过将Mondrian分位数回归与正则化收缩机制相结合，构造了一种在复杂调查数据中平衡子群覆盖率和预测集大小的评分函数，并严格推导了加权子群覆盖差距的渐近性质。其理论贡献在于证明了在存在校准单元碎片化和群体置信度不匹配时，边际有效性无法保证子群可靠性，且朴素分组校准反而会加剧公平性与效率之间的统计权衡。


### Impossibility of Distribution-Free Predictive Inference for Individual Treatment Effects 

- [ ] **分类**: Conformal | **日期**: 2026-05-06
- **链接**: [PDF](http://arxiv.org/abs/2605.05051v1)

!!! note "AI 核心解读"

    该论文通过建立有限样本与渐近不可能性定理，严格证明了在连续协变量存在时，任何分布自由的个体处理效应预测集若达到目标覆盖率，其期望长度必然为无穷大。其核心创新在于将个体处理效应推断的不可行性归约为条件独立性检验的统计难度，从而揭示了因果推断中缺失数据机制对评分函数构造的根本性限制。

### Classification-Powered Conformal Inference for Zero-inflated Outcomes 

- [ ] **分类**: Conformal | **日期**: 2026-05-05
- **链接**: [PDF](http://arxiv.org/abs/2605.04219v1)

!!! note "AI 核心解读"

    该论文提出一种结合分类步骤的共形推断框架，通过构造一个基于分类器输出的评分函数，将零膨胀结果拆解为零点与非零连续部分，从而在保持边际覆盖率的条件下，从理论上证明该方法能渐近达到区间长度的最优性。其核心创新在于利用分类信息修正共形预测的保守性，并证明无论分类或回归模型如何选择，所构造的预测集都能在交换性假设下实现目标覆盖与最小化区间长度的理论平衡。


### Graph Convolutional Support Vector Regression for Robust Spatiotemporal Forecasting of Urban Air Pollution 

- [ ] **分类**: Conformal | **日期**: 2026-05-05
- **链接**: [PDF](http://arxiv.org/abs/2605.03795v1)

!!! note "AI 核心解读"

    该论文的核心创新在于将图卷积网络与支持向量回归相结合，通过构造一个鲁棒的评分函数来同时捕捉站点间的空间依赖性与时间非线性动态，并在损失函数中引入对异常值的惩罚机制，从而在统计上提升了模型对污染尖峰和季节波动的稳健性。理论性质上，作者证明了所提方法在非平稳时空数据下的预测一致性，并利用共形预测为点预测提供了具有统计保证的校准预测区间。

### Conformalized Percentile Interval: Finite Sample Validity and Improved Conditional Performance 

- [ ] **分类**: Conformal | **日期**: 2026-05-04
- **链接**: [PDF](http://arxiv.org/abs/2605.03233v1)

!!! note "AI 核心解读"

    该论文提出一种在概率积分变换（PIT）空间中进行共形校准的方法，通过神经网络估计条件累积分布函数后，构造出由该估计决定最短长度的有限样本调整百分位区间。其核心创新在于：校准PIT值能利用其渐近特征独立性来缓解特征依赖的覆盖偏差，同时百分位校准对条件CDF的估计误差具有鲁棒性，并严格证明了有限样本边际覆盖和渐近条件覆盖的理论性质。


### Stable Localized Conformal Prediction via Transduction 

- [ ] **分类**: Conformal | **日期**: 2026-05-02
- **链接**: [PDF](http://arxiv.org/abs/2605.01452v1)

!!! note "AI 核心解读"

    该论文首次将共形预测集合大小的条件期望方差形式化为“集合稳定性”这一统计量，并基于此推导出稳定共形预测（StCP）方法，通过引入源任务标签数据与目标任务无标签数据构造了新的评分函数，在理论上严格证明了其边际覆盖率和稳定性性质。


### Optimal Spatio-Temporal Decoupling for Bayesian Conformal Prediction 

- [ ] **分类**: Conformal | **日期**: 2026-05-01
- **链接**: [PDF](http://arxiv.org/abs/2605.00432v1)

!!! note "AI 核心解读"

    该论文提出状态自适应贝叶斯共形预测（SA-BCP），通过空间核密度证据门控长期时间惯性，在统计上实现了时间适应性与结构稳定性的最优解耦。其核心创新在于严格证明了由证据阈值K控制的极小极大偏差-方差权衡，并构造了能同时消除ACI系统欠覆盖并缩减贝叶斯CP未校准区间膨胀的评分函数。

### SPLICE: Latent Diffusion over JEPA Embeddings for Conformal Time-Series Inpainting 

- [ ] **分类**: Conformal | **日期**: 2026-04-30
- **链接**: [PDF](http://arxiv.org/abs/2605.00126v1)

!!! note "AI 核心解读"

    该论文的核心创新在于将基于JEPA的潜变量扩散生成模型与自适应共形推断（ACI）相结合，首次为时间序列插补提供了具有有限样本覆盖保证的在线预测区间。其统计贡献在于：通过构造一个条件潜桥扩散过程，在64维潜空间中生成候选填补轨迹，并利用ACI动态调整共形分位数，从而在无需分布假设的情况下，严格校正了静态共形预测中高达7.5个百分点的覆盖不足问题。


### VLM Judges Can Rank but Cannot Score: Task-Dependent Uncertainty in Multimodal Evaluation 

- [ ] **分类**: Conformal | **日期**: 2026-04-28
- **链接**: [PDF](http://arxiv.org/abs/2604.25235v2)

!!! note "AI 核心解读"

    该论文首次将共形预测（conformal prediction）框架系统性地应用于视觉语言模型（VLM）作为自动评估者的场景，仅利用评分令牌的对数概率构造出无需重训练的校准预测区间，并严格证明了评估不确定性具有显著的任务依赖性。其核心创新在于揭示了“排序-评分解耦”这一统计失效模式：即便VLM在排序相关性指标上表现优异，其构造的区间仍可能过宽且无信息量，从而在理论上区分了相对排序能力与绝对评分可靠性的本质差异。


### VLM Judges Can Rank but Cannot Score: Task-Dependent Uncertainty in Multimodal Evaluation 

- [ ] **分类**: Conformal | **日期**: 2026-04-28
- **链接**: [PDF](http://arxiv.org/abs/2604.25235v1)

!!! note "AI 核心解读"

    该论文首次将共形预测框架系统性地应用于视觉语言模型（VLM）作为自动评分者的不确定性量化，仅利用评分令牌的对数概率构造出无需重训的校准预测区间。其核心理论发现是评分区间宽度与任务类型强相关，并揭示了“排序-评分解耦”这一统计失效模式：即模型在排序相关性指标上表现良好，但生成的绝对评分区间却过宽且无信息量。

### Tail allocation for conformal prediction intervals 

- [ ] **分类**: Conformal | **日期**: 2026-04-28
- **链接**: [PDF](http://arxiv.org/abs/2604.25202v1)

!!! note "AI 核心解读"

    该论文提出了一种新的尾部分配方法（TA-CQR），通过搜索分位数定义的核心来估计最优的单区间预测集，并利用非负加性分裂共形校准，在可交换性假设下严格保证了有限样本的边际覆盖。理论贡献包括：刻画了单区间预言机在单峰性下的最高密度几何解释与不连通高密度集的正连接成本，证明了所选分配与核心的局部恢复性，以及在端点密度条件下校准半径的渐近可忽略性，并给出了有限样本校准长度与预言机之间的显式不等式。

### Conflict Forecasting via Conformal Prediction for Markov Processes 

- [ ] **分类**: Conformal | **日期**: 2026-04-28
- **链接**: [PDF](http://arxiv.org/abs/2604.25139v1)

!!! note "AI 核心解读"

    该论文将保形预测框架拓展至离散状态马尔可夫过程，针对时序依赖数据构造了非交换性条件下的预测集，并推导了其有限样本覆盖保证。通过对比基于似然的预测策略，作者证明了所提评分函数在模型误设定下仍能维持有效的置信度校准，且无需依赖数据可交换性假设。


### Privacy-preserving Meta-analysis through Low-Rank Basis Hunting 

- [ ] **分类**: Conformal | **日期**: 2026-04-26
- **链接**: [PDF](http://arxiv.org/abs/2604.23847v1)

!!! note "AI 核心解读"

    该论文提出MetaHunt方法，通过将各研究的真实函数假设为共享低秩基函数的凸组合，并扩展函数型逐次投影算法实现去噪基函数搜寻，在温和正则条件下证明了基函数估计的一致性。进一步利用半参数或非参数模型刻画研究层面协变量与混合权重的关联，并基于共形预测构造预测区间，在可交换性与弱估计误差条件下证明了其渐近有效的边际覆盖性质。


### Confirmatory Biomarker Identification with k-FWER Control Using Derandomized Knockoffs with Cox Regression 

- [ ] **分类**: Knockoff | **日期**: 2025-04-04
- **链接**: [PDF](http://arxiv.org/abs/2504.03907v2)

!!! note "AI 核心解读"

    该论文创新点在于：针对Cox回归模型，提出一种去随机化knockoffs方法，通过聚合多次随机knockoffs实现的特征选择结果，构造了更稳定的评分函数，并在理论上严格证明了该方法能控制k族错误率（k-FWER）。其核心统计贡献是，在保持错误率控制的同时，通过去随机化策略显著降低了传统knockoffs因单次随机实现导致的选择不稳定性。

### Interpretable Feature Interaction via Statistical Self-supervised Learning on Tabular Data 

- [ ] **分类**: Knockoff | **日期**: 2025-03-23
- **链接**: [PDF](http://arxiv.org/abs/2503.18048v1)

!!! note "AI 核心解读"

    该论文提出Spofe方法，通过核主成分分析生成自监督信号，并将其蒸馏为稀疏多项式函数，从而在非线性特征交互建模中实现可解释性。其核心统计贡献在于：构建了基于多目标knockoff选择的特征筛选框架，并严格证明了错误发现率（FDR）的控制与误差界，为高维表格数据提供了兼具统计严谨性与可解释性的特征交互识别方法。

### One-at-a-time knockoffs: controlled false discovery rate with higher power 

- [ ] **分类**: Knockoff | **日期**: 2025-02-26
- **链接**: [PDF](http://arxiv.org/abs/2502.18750v1)

!!! note "AI 核心解读"

    该论文提出了一种逐列生成 knockoff 设计矩阵的新方法（OATK），通过仅替换原始设计矩阵的单一列来保持 Gram 矩阵，从而大幅简化了 Barber 与 Candès 的联合约束框架；在原始设计矩阵满足温和相关性假设的条件下，通过构造原始系数与 knockoff 系数的比较统计量，OATK 能够渐近控制错误发现率（FDR），并显著提升统计功效。

### Asymptotic FDR Control with Model-X Knockoffs: Is Moments Matching Sufficient? 

- [ ] **分类**: Knockoff | **日期**: 2025-02-09
- **链接**: [PDF](http://arxiv.org/abs/2502.05969v1)

!!! note "AI 核心解读"

    该论文通过将模型-X敲除框架中的分布可交换性条件替换为近似敲除统计量的三个新条件，建立了渐进FDR控制的理论统一框架；并严格证明了基于前两阶矩匹配的高斯敲除生成器在采用相应统计量时能实现渐进FDR控制，首次从理论上验证了该方法的有效性与稳健性。

### Can linear algebra create perfect knockoffs? 

- [ ] **分类**: Knockoff | **日期**: 2025-02-04
- **链接**: [PDF](http://arxiv.org/abs/2502.02148v1)

!!! note "AI 核心解读"

    该论文创新点在于：通过线性代数方法直接构造“完美”的Model-X knockoffs，并基于特征与knockoffs之间的平均绝对相关系数构建了新的质量度量指标，从而绕开了对复杂条件分布的精确推导；同时，作者提出了多种计算加速策略以降低优化算法的高昂计算成本，并严格证明了所构造knockoffs在统计性质上的伪完美性。

### Conformalized Super Learner 

- [ ] **分类**: Conformal | **日期**: 2026-04-24
- **链接**: [PDF](http://arxiv.org/abs/2604.22391v1)

!!! note "AI 核心解读"

    该论文创新性地将共形预测框架与超级学习器集成，通过加权多数投票机制融合各学习器的共形得分，从而在有限样本下为集成预测提供严格的覆盖保证。其核心理论贡献在于，在可交换性假设及异方差、稀疏性等复杂数据生成机制下，严格刻画了所构造预测区间的有限样本覆盖性质，无需依赖渐近论证或重采样过程。

### Cross-Domain Uncertainty Quantification for Selective Prediction: A Comprehensive Bound Ablation with Transfer-Informed Betting 

- [ ] **分类**: Conformal | **日期**: 2026-03-09
- **链接**: [PDF](http://arxiv.org/abs/2603.08907v1)

!!! note "AI 核心解读"

    该论文的核心创新在于提出了**迁移知情投注（Transfer-Informed Betting, TIB）**，通过用源域风险分布预热WSR财富过程，构造了一个在任意源-目标域散度下仍为有效鞅的评分函数，并证明了当域匹配时TIB严格优于标准WSR，且任何数据无关的预热策略都无法达到更优的收敛速率。此外，论文系统性地消融了九种有限样本界族，揭示了**投注置信序列与LTT单调检验及跨域迁移的三重组合**在统计推导上的理论优势，即通过消除联合界中的ln(K)惩罚项，在数据稀缺场景下实现了覆盖率的显著提升。

### Beyond Data Splitting: Full-Data Conformal Prediction by Differential Privacy 

- [ ] **分类**: Conformal | **日期**: 2026-03-08
- **链接**: [PDF](http://arxiv.org/abs/2603.07522v1)

!!! note "AI 核心解读"

    该论文提出了一种无需数据分割的全数据差分隐私保形预测框架，通过利用差分隐私引入的稳定性来约束样本内与样本外保形分数之间的差距，并构造了一个保守的私有分位数估计程序以防止覆盖不足。理论上，该工作证明了通用差分隐私保证可提供一个普适的覆盖下限，但无法恢复名义水平 \(1-\alpha\)，进而通过机制特定的稳定性分析实现了名义水平的渐近恢复。

### ConfHit: Conformal Generative Design with Oracle Free Guarantees 

- [ ] **分类**: Conformal | **日期**: 2026-03-07
- **链接**: [PDF](http://arxiv.org/abs/2603.07371v1)

!!! note "AI 核心解读"

    该论文提出了一种无需实验验证的分布自由框架，通过构建基于加权可交换性的多样本密度比加权共形p值，解决了生成模型在预算限制和无真实标签条件下的命中率统计保证问题。其核心创新在于设计了嵌套检验程序，在保持多重比较校正的同时，对生成批次进行认证与精炼，从而在有限样本下严格控制置信水平并保证覆盖率的有效性。

### CREDO: Epistemic-Aware Conformalized Credal Envelopes for Regression 

- [ ] **分类**: Conformal | **日期**: 2026-03-06
- **链接**: [PDF](http://arxiv.org/abs/2603.06826v1)

!!! note "AI 核心解读"

    该论文提出一种“先构造认知包络、再执行分裂共形校准”的两阶段框架，通过基于极端后验预测端点修剪的快速算法构建随局部证据强度自适应增宽的认知包络，并证明其能在无需额外假设下保证边际覆盖有效性，同时将预测区间宽度分解为偶然噪声、认知膨胀和分布自由校准松弛三个可解释分量。


### Gold after Randomized Sand: Model-X Split Knockoffs for Controlled Transformation Selection 

- [ ] **分类**: Knockoff | **日期**: 2025-07-02
- **链接**: [PDF](http://arxiv.org/abs/2507.01732v2)

!!! note "AI 核心解读"

    该论文提出Model-X Split Knockoffs方法，通过引入新颖的辅助随机化设计，解决了随机协变量设计与确定性线性变换之间的根本性矛盾，从而在随机设计下实现了对变换选择的有限样本FDR控制。其核心创新在于构造了一种新型评分函数，使得在已知或可准确估计的协变量分布下，无论响应变量的条件分布如何，都能保证与经典Model-X Knockoffs相同或更优的统计功效，并具备可证明的有限样本FDR控制理论性质。

### Knockoffs Inference under Privacy Constraints 

- [ ] **分类**: Knockoff | **日期**: 2025-06-11
- **链接**: [PDF](http://arxiv.org/abs/2506.09690v1)

!!! note "AI 核心解读"

    该论文在差分隐私约束下，为模型-X knockoff框架构造了经严格统计推导的噪声注入机制，并证明了该机制能同时保证精确的有限样本FDR控制与隐私保护。通过建立评分函数的渐近理论性质，作者给出了噪声不损害统计功效的充分条件，从而在理论上统一了隐私保护与变量选择的有效性。

### Knockoff-Guided Compressive Sensing: A Statistical Machine Learning Framework for Support-Assured Signal Recovery 

- [ ] **分类**: Knockoff | **日期**: 2025-05-30
- **链接**: [PDF](http://arxiv.org/abs/2505.24727v1)

!!! note "AI 核心解读"

    该论文的核心创新在于将Knockoff变量筛选框架与压缩感知结合，通过构造一种新型评分函数实现有限样本下对支撑集识别的严格FDR控制，并基于此推导出比传统ℓ1方法更弱的理论恢复条件，从而在保证信号重构精度的同时，为支撑恢复提供了统计可解释的误差控制。

### Explaining Concept Shift with Interpretable Feature Attribution 

- [ ] **分类**: Knockoff | **日期**: 2025-05-27
- **链接**: [PDF](http://arxiv.org/abs/2505.20634v1)

!!! note "AI 核心解读"

    该论文提出基于广义可加模型（GAM）的SGShift框架，通过构造条件标签分布的评分函数并引入knockoffs控制错误发现率，实现了对概念漂移特征的稀疏识别与统计推断。理论性质上，该方法在特征归因中融合了吸收项以修正模型欠拟合偏差，并证明了所提评分函数在特征选择中的渐近有效性。

### Controlling the false discovery rate in high-dimensional linear models using model-X knockoffs and $p$-values 

- [ ] **分类**: Knockoff | **日期**: 2025-05-22
- **链接**: [PDF](http://arxiv.org/abs/2505.16124v2)

!!! note "AI 核心解读"

    该论文的核心创新在于：将模型-X knockoff框架与去偏惩罚估计量结合，构造了两组天然配对的高维检验统计量及其p值，并证明了第一组统计量渐近相互独立，从而为Benjamini-Hochberg程序提供了理论依据；进一步通过两步法利用配对结构提升检验功效，在一般设计矩阵相关结构下严格证明了渐近FDR控制性，并刻画了功效增益的理论性质。

### Posterior inference via Hill's prediction model 

- [ ] **分类**: Conformal | **日期**: 2026-03-20
- **链接**: [PDF](http://arxiv.org/abs/2603.20071v1)

!!! note "AI 核心解读"

    该论文基于Hill的$A_n$预测模型，通过构造一步前向预测分布函数，在无需先验分布的前提下，直接生成完整数据集以实现任意统计量的后验推断。其核心创新在于利用共形预测框架下的等概率区间划分机制，推导出后验分布的闭合形式，并证明了该评分函数在有限样本下具有频率学派覆盖性质。

### Conformalized Robust Principal Component Analysis 

- [ ] **分类**: Conformal | **日期**: 2026-03-15
- **链接**: [PDF](http://arxiv.org/abs/2603.14233v1)

!!! note "AI 核心解读"

    该论文的创新点在于：将共形预测框架与稳健主成分分析相结合，构建了一种无需分布假设的评分函数，并引入加权校准机制以处理异质性观测概率，从而在有限样本下为低秩矩阵恢复提供了严格的覆盖概率理论保证。

### Efficient Federated Conformal Prediction with Group-Conditional Guarantees 

- [ ] **分类**: Conformal | **日期**: 2026-03-15
- **链接**: [PDF](http://arxiv.org/abs/2603.14198v2)

!!! note "AI 核心解读"

    该论文提出了一种基于组分层核集（group-stratified coresets）的联邦共形预测框架，通过构造可合并的加权评分摘要，在无需共享原始数据的前提下实现分布自由的组条件覆盖保证。其核心创新在于推导了跨客户端异构数据分布下的评分函数聚合理论，证明了所构造的加权评分集能够保持原始局部校准过程的覆盖性质，并严格控制了组条件覆盖误差的有限样本上界。

### Beyond Accuracy: Reliability and Uncertainty Estimation in Convolutional Neural Networks 

- [ ] **分类**: Conformal | **日期**: 2026-03-11
- **链接**: [PDF](http://arxiv.org/abs/2603.10731v1)

!!! note "AI 核心解读"

    该论文通过对比蒙特卡洛丢弃法（贝叶斯近似）与共形预测（非参数框架）在卷积神经网络中的不确定性量化表现，揭示了高精度模型（如H-CNN VGG16）往往伴随过度自信的校准缺陷，而共形预测凭借其构造的统计保证预测集，在理论上确保了覆盖率的有效性与无分布假设的稳健性。

### Conformal prediction for high-dimensional functional time series: Applications to subnational mortality 

- [ ] **分类**: Conformal | **日期**: 2026-03-11
- **链接**: [PDF](http://arxiv.org/abs/2603.10674v1)

!!! note "AI 核心解读"

    该论文将分裂与序列两种共形预测框架拓展至高维函数型时间序列，通过验证集校准经验覆盖概率以优化调参，或利用自回归过程序贯更新预测分位数，从而在无需模型假设下构造出具有有限样本有效性的预测区间。其理论贡献在于证明了所构造评分函数在非平稳依赖结构下的覆盖保真性，并给出了序列共形预测中分位数更新过程的渐近一致性。


### Diffusion-Driven High-Dimensional Variable Selection 

- [ ] **分类**: Knockoff | **日期**: 2025-08-19
- **链接**: [PDF](http://arxiv.org/abs/2508.13890v1)

!!! note "AI 核心解读"

    该论文提出一种基于扩散模型生成多重伪数据集的再抽样聚合框架，通过构造稳定性评分函数并证明其在温和假设下的选择相合性，实现了对高维强相关数据的稳定变量选择。其核心创新在于将扩散模型的生成能力与聚合推断理论结合，使得评分函数在理论上可校准，并拓展至图模型选择与统计推断，为高维变量选择提供了新的统计推导路径。

### Variable selection via knockoffs in missing data settings with categorical predictors 

- [ ] **分类**: Knockoff | **日期**: 2025-08-08
- **链接**: [PDF](http://arxiv.org/abs/2508.06138v1)

!!! note "AI 核心解读"

    该论文将多重插补与knockoffs框架结合，针对含缺失值的分类预测变量场景，构造了适用于无序类别变量的评分函数，并证明了在多重插补后各数据集上独立应用knockoffs筛选仍能控制错误发现率。其理论贡献在于将knockoffs的有限样本保证扩展至缺失数据与多水平随机效应并存的情形，无需对缺失机制或类别顺序施加额外假设。

### Differentially Private Model-X Knockoffs via Johnson-Lindenstrauss Transform 

- [ ] **分类**: Knockoff | **日期**: 2025-08-06
- **链接**: [PDF](http://arxiv.org/abs/2508.04800v1)

!!! note "AI 核心解读"

    该论文的核心创新在于：通过高斯Johnson-Lindenstrauss变换对数据knockoff矩阵进行私有化，利用近似等距性在满足(ε,δ)-差分隐私的同时保留协变量间的依赖结构，从而避免传统噪声注入破坏交换性条件；理论层面，作者提出了一种针对高维私有knockoff过程的“去偏技术”，并严格刻画了降维比、信噪比、隐私参数等因子对FDR控制与统计功效的渐近影响，建立了功效趋于1的充分条件。

### A powerful procedure that controls the false discovery rate with directional information 

- [ ] **分类**: Knockoff | **日期**: 2025-07-21
- **链接**: [PDF](http://arxiv.org/abs/2507.15631v1)

!!! note "AI 核心解读"

    该论文提出了 signed-knockoff 程序，通过构造一种结合符号信息的评分函数，在有限样本下严格控制错误发现率（FDR），并基于统计推导证明了其在利用方向性信息时相比传统 p 值方法具有更强的检验功效。其理论核心在于将方向性信息嵌入 knockoff 框架的统计量构造中，从而在不依赖渐近近似的前提下实现 FDR 的精确控制。

### Where to Intervene: Action Selection in Deep Reinforcement Learning 

- [ ] **分类**: Knockoff | **日期**: 2025-07-05
- **链接**: [PDF](http://arxiv.org/abs/2507.04187v1)

!!! note "AI 核心解读"

    该论文提出了一种基于数据驱动的无模型动作选择方法，通过构造knockoff采样机制在统计上控制错误发现率，从而在无需先验知识的情况下筛选出最小充分动作集。其核心理论贡献在于将变量选择中的假阳性控制框架与深度强化学习的在线训练过程相融合，并证明了所选动作集在保持策略最优性方面的统计性质。

### Elements of Conformal Prediction for Statisticians 

- [ ] **分类**: Conformal | **日期**: 2026-03-25
- **链接**: [PDF](http://arxiv.org/abs/2603.23923v1)

!!! note "AI 核心解读"

    该论文的核心创新在于：通过构造基于交换性假设的评分函数（如非一致性得分），将预测区间构建转化为分位数估计问题，并严格证明了在有限样本下无需分布假设即可实现精确的边际覆盖概率。其理论贡献在于揭示了评分函数与保序回归、核密度估计等经典统计工具的深层联系，并推导出条件覆盖率的渐近性质。

### Shape-Adaptive Conditional Calibration for Conformal Prediction via Minimax Optimization 

- [ ] **分类**: Conformal | **日期**: 2026-03-24
- **链接**: [PDF](http://arxiv.org/abs/2603.23374v1)

!!! note "AI 核心解读"

    该论文通过极小极大优化框架，将条件覆盖率的边际矩约束转化为对集值映射的灵活校准，从而突破了传统评分函数的结构限制，实现了形状自适应的预测集构造。理论方面，作者推导了非渐近的预言不等式，并证明在正则条件下覆盖误差的收敛速率达到最优阶，同时保证了基于校准阶段敏感属性（测试时不可观测）的条件推断有效性。

### A plug-and-play approach with fast uncertainty quantification for weak lensing mass mapping 

- [ ] **分类**: Conformal | **日期**: 2026-03-23
- **链接**: [PDF](http://arxiv.org/abs/2603.22006v1)

!!! note "AI 核心解读"

    该论文提出一种即插即用的弱透镜质量映射方法，通过交替执行带有精心设计数据保真项的梯度下降步骤与单一高斯白噪声去噪网络，实现了无需针对观测噪声协方差重新训练的快速点估计；同时构造基于矩网络的快速无采样不确定性量化方案，并引入保形预测理论保证误差条覆盖概率的严格理论性质。

### CRPS-Optimal Binning for Univariate Conformal Regression 

- [ ] **分类**: Conformal | **日期**: 2026-03-23
- **链接**: [PDF](http://arxiv.org/abs/2603.22000v2)

!!! note "AI 核心解读"

    该论文提出了一种基于协变量排序后分箱的非参数条件分布估计方法，通过最小化留一法连续排序概率评分（LOO-CRPS）的闭式代价函数，并利用动态规划在 \(O(n^2 K)\) 时间内全局最优地确定分箱边界。理论贡献在于证明了直接最小化样本内 LOO-CRPS 会导致过拟合，转而采用 K 折交叉验证选择箱数 \(K\)，从而得到具有明确最小值的 U 型准则；同时，基于 CRPS 作为非一致性评分构造的共形预测集，在有限样本下保证了边际覆盖率的精确控制。

### CoNBONet: Conformalized Neuroscience-inspired Bayesian Operator Network for Reliability Analysis 

- [ ] **分类**: Conformal | **日期**: 2026-03-23
- **链接**: [PDF](http://arxiv.org/abs/2603.21678v1)

!!! note "AI 核心解读"

    该论文创新点在于：将分治共形预测（split conformal prediction）与神经科学启发的脉冲神经元模型相结合，构造了一个具有理论保证的校准不确定性量化评分函数，从而在贝叶斯算子网络框架下实现了对非线性动力系统时变可靠性的高效推断。其理论性质体现在，通过共形预测的有限样本覆盖保证，所提方法能够在不依赖传统贝叶斯后验采样的前提下，提供具有统计严格性的失效概率置信区间。


### Choosing the nominal level post-hoc with knockoffs using e-values 

- [ ] **分类**: Knockoff | **日期**: 2025-11-14
- **链接**: [PDF](http://arxiv.org/abs/2511.11166v2)

!!! note "AI 核心解读"

    该论文通过引入 e 值框架，构造了一种可在数据观测后动态调整名义 FDR 水平的后验评分函数，突破了传统 knockoff 滤波器对名义水平的刚性依赖；其核心统计推导证明了该后验过程在任意数据实现下均能严格保持 FDR 控制，且在不增加任何统计代价的前提下，显著提升了低维稀疏信号场景下的选择效能与精度。

### DiffKnock: Diffusion-based Knockoff Statistics for Neural Networks Inference 

- [ ] **分类**: Knockoff | **日期**: 2025-10-01
- **链接**: [PDF](http://arxiv.org/abs/2510.01418v1)

!!! note "AI 核心解读"

    该论文提出基于扩散模型的DiffKnock框架，通过训练扩散模型生成满足有限样本FDR控制的knockoff变量，并利用神经网络构造的梯度与过滤统计量构建反对称特征重要性度量，从而在保留复杂特征依赖结构的同时实现非线性关联检测。其理论贡献在于将深度生成模型的灵活性嵌入knockoff统计框架，在有限样本下严格保证FDR控制，并通过构造的统计量在非线性高维场景中显著提升统计功效。

### Sparse minimum Redundancy Maximum Relevance for feature selection 

- [ ] **分类**: Knockoff | **日期**: 2025-08-26
- **链接**: [PDF](http://arxiv.org/abs/2508.18901v1)

!!! note "AI 核心解读"

    该论文将经典离散mRMR准则连续化，并引入非凸惩罚项构造稀疏化评分函数，通过零系数估计实现非活跃特征的精确识别；进一步结合knockoff滤波器设计多阶段筛选流程，在控制错误发现率（FDR）的理论框架下证明了特征剔除的统计一致性。

### Mapping beyond diseases: Controlled variable selection for secondary phenotypes using tilted knockoffs 

- [ ] **分类**: Knockoff | **日期**: 2025-08-25
- **链接**: [PDF](http://arxiv.org/abs/2508.18548v2)

!!! note "AI 核心解读"

    该论文通过引入“倾斜分布”对模型-X knockoff框架进行修正，将选择概率纳入评分函数的构造，从而在偏倚抽样下实现了对条件独立性假设的FDR控制；其核心创新在于证明了基于倾斜分布构造的knockoff变量能够校正collider偏差，并推导了该评分函数在有限样本下控制FDR的理论性质。

### Novel Knockoff Generation and Importance Measures with Heterogeneous Data via Conditional Residuals and Local Gradients 

- [ ] **分类**: Knockoff | **日期**: 2025-08-20
- **链接**: [PDF](http://arxiv.org/abs/2508.14882v1)

!!! note "AI 核心解读"

    该论文提出基于条件残差的 knockoff 生成框架，无需假设数据分布或类型同质性，通过残差构造实现异质数据下的伪变量生成；同时引入平均绝对局部导数（MALD）作为变量重要性度量，该评分函数兼容任意非线性模型，并在理论上保证了错误发现率的控制与统计功效的提升。

### A Kernel Nonconformity Score for Multivariate Conformal Prediction 

- [ ] **分类**: Conformal | **日期**: 2026-04-23
- **链接**: [PDF](http://arxiv.org/abs/2604.21595v1)

!!! note "AI 核心解读"

    该论文提出了一种多元核评分函数（MKS），通过将残差向量压缩为标量并显式保留其几何结构，其统计推导表明该评分可分解为各向异性最大均值差异（MMD），从而在核密度估计与协方差加权距离之间实现插值。理论性质上，论文证明了有限样本覆盖保证，并建立了基于核协方差算子有效秩的收敛速率，实现了与维度无关的自适应，而非依赖环境维度。

### Conformal Prediction Assessment: A Framework for Conditional Coverage Evaluation and Selection 

- [ ] **分类**: Conformal | **日期**: 2026-03-28
- **链接**: [PDF](http://arxiv.org/abs/2603.27189v2)

!!! note "AI 核心解读"

    该论文将条件覆盖率的评估重新构造为监督学习问题，通过训练一个可靠性估计器来预测实例层面的覆盖概率，并建立了该估计量的收敛速率。在此基础上，论文定义了条件有效性指数（CVI），将可靠性分解为安全性（欠覆盖风险）与效率（过覆盖代价），并证明了基于CVI的模型选择方法的一致性。

### Conformalized Signal Temporal Logic Inference under Covariate Shift 

- [ ] **分类**: Conformal | **日期**: 2026-03-28
- **链接**: [PDF](http://arxiv.org/abs/2603.27062v1)

!!! note "AI 核心解读"

    该论文的核心创新在于：针对训练与部署数据存在协变量偏移的场景，首次将加权保形预测框架与信号时序逻辑（STL）推理相结合，通过估计似然比构造了基于STL鲁棒性的加权评分函数。理论性质上，该方法在非交换性数据下仍能保证有限样本的边际覆盖有效性，且推导了加权保形集在分布偏移下的误差界。

### Contrastive Conformal Sets 

- [ ] **分类**: Conformal | **日期**: 2026-03-27
- **链接**: [PDF](http://arxiv.org/abs/2603.26261v1)

!!! note "AI 核心解读"

    该论文将共形预测扩展到对比学习场景，通过引入可学习的广义多范数约束构造最小体积覆盖集，并理论证明了体积最小化可作为负样本排除的代理指标。其评分函数基于集合几何的优化设计，在保证正样本无分布假设覆盖的同时，通过统计推导实现了负样本排除的最大化。

### Conformal Prediction for Nonparametric Instrumental Regression 

- [ ] **分类**: Conformal | **日期**: 2026-03-26
- **链接**: [PDF](http://arxiv.org/abs/2603.25509v1)

!!! note "AI 核心解读"

    该论文将条件覆盖保证重新表述为工具变量偏移类 \(\mathcal{F}\) 上的边际覆盖，从而将共形推断框架扩展到非参数工具变量回归中；其理论贡献在于证明了对于任意选定的偏移类，所构造的预测区间具有分布自由且有限样本的覆盖性质，且该评分函数可与任意NPIV估计量（如筛分2SLS或神经网络极小极大方法）结合使用。


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

