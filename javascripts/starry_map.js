document.addEventListener("DOMContentLoaded", function() {
    const container = document.getElementById('starry-map');
    if (!container) return;

    const width = container.clientWidth;
    const height = container.clientHeight;

    // 🌍 自动检测当前语言环境 (zh 或 en)
    const isEn = window.location.pathname.includes('/en/');

    // --- 🌌 星云配色系统 (Nebula Palette) ---
    const colors = {
        core: "#FFD700",       // 数学科学：恒星金
        analysis: "#00CED1",   // 分析学：深青色
        analysisSub: "#48D1CC",// 泛函分析：中青色
        analysisLeaf: "#E0FFFF",// 章节小行星：亮星白
        prob: "#4169E1",       // 概率与统计：皇家蓝
        probSub: "#87CEFA",    // 具体课程：天蓝色
        algebra: "#9370DB",    // 代数学：星云紫
        geometry: "#3CB371",   // 几何学：翡翠绿
        applied: "#FF7F50",    // 应用与计算：珊瑚橙
    };

    // --- 1. 定义整个宇宙的完整数据 ---
    const allNodes = [
        { id: "Math-Core", label: isEn ? "Mathematical Sciences" : "数学科学", radius: 50, color: colors.core, collapsed: false },

        { id: "Analysis", parentId: "Math-Core", label: isEn ? "Analysis" : "分析学", radius: 30, color: colors.analysis, collapsed: false },
        { id: "Algebra", parentId: "Math-Core", label: isEn ? "Algebra" : "代数学", radius: 30, color: colors.algebra, collapsed: false },
        { id: "Geometry", parentId: "Math-Core", label: isEn ? "Geometry" : "几何学", radius: 30, color: colors.geometry, collapsed: false },
        { id: "Prob-Stat", parentId: "Math-Core", label: isEn ? "Probability & Statistics" : "概率与统计", radius: 30, color: colors.prob, collapsed: false },
        { id: "Applied-Comp", parentId: "Math-Core", label: isEn ? "Applied Mathematics" : "应用与计算", radius: 30, color: colors.applied, collapsed: false },

        { id: "FA", parentId: "Analysis", label: isEn ? "Functional Analysis" : "泛函分析", radius: 15, color: colors.analysisSub, collapsed: true },
        { id: "AS", parentId: "Prob-Stat", label: isEn ? "Asymptotic Statistics" : "渐近统计", radius: 15, color: colors.probSub, collapsed: true },
        
        // ✨ LDP 节点更新：添加了跳转 URL，点击星星直接跳转
        { id: "LDP", parentId: "Prob-Stat", label: isEn ? "LDP Seminar" : "大偏差原理", radius: 15, color: colors.probSub, url: "Lectures/Course_Lectures/LDP/" },
        
        { id: "SDE", parentId: "Prob-Stat", label: isEn ? "Stochastic Diff Eq" : "随机微分方程", radius: 15, color: colors.probSub, collapsed: true },

        { id: "Top", parentId: "Geometry", label: isEn ? "General Topology" : "点集拓扑", radius: 15, color: colors.geometry, collapsed: true },

        { id: "RA", parentId: "Analysis", label: isEn ? "Real Analysis" : "实分析", radius: 15, color: colors.analysisSub, collapsed: true },

        { id: "SL", parentId: "Prob-Stat", label: isEn ? "Statistical Learning" : "统计学习", radius: 15, color: colors.probSub, collapsed: true },

        // 泛函的内部章节
        { id: "FA-Intro", label: isEn ? "0. Intro" : "课程简介", parentId: "FA", radius: 8, color: colors.analysisLeaf, url: "Functional_Analysis/0Intro/" },
        { id: "FA-Ch1", label: isEn ? "1. Metric Spaces" : "1.度量空间", parentId: "FA", radius: 8, color: colors.analysisLeaf, url: "Functional_Analysis/1distance_space/" },
        { id: "FA-Ch2", label: isEn ? "2. Linear Operators" : "2.线性算子", parentId: "FA", radius: 8, color: colors.analysisLeaf, url: "Functional_Analysis/2completeness/" },
        { id: "FA-Ch3", label: isEn ? "3. Compactness" : "3.紧性", parentId: "FA", radius: 8, color: colors.analysisLeaf, url: "Functional_Analysis/3compactness/" },
        { id: "FA-Ch4", label: isEn ? "4. Fixed Point Theorems" : "4.压缩映射定理", parentId: "FA", radius: 8, color: colors.analysisLeaf, url: "Functional_Analysis/4CMT/" },
        { id: "FA-Ch5", label: isEn ? "5. Normed Spaces" : "5.赋范线性空间", parentId: "FA", radius: 8, color: colors.analysisLeaf, url: "Functional_Analysis/5Linear/" },
        { id: "FA-Ch6", label: isEn ? "6. Banach Spaces" : "6.巴拿赫空间", parentId: "FA", radius: 8, color: colors.analysisLeaf, url: "Functional_Analysis/6Banach/" },
        { id: "FA-Ch7", label: isEn ? "7. Hilbert Spaces" : "7.希尔伯特空间", parentId: "FA", radius: 8, color: colors.analysisLeaf, url: "Functional_Analysis/7Hilbert/" },
        { id: "FA-Ch8", label: isEn ? "8. Hilbert Spaces Theory" : "8.希尔伯特空间理论", parentId: "FA", radius: 8, color: colors.analysisLeaf, url: "Functional_Analysis/8Hilbert/" },
        { id: "FA-Ch9", label: isEn ? "9. Linear Operators" : "9.线性算子", parentId: "FA", radius: 8, color: colors.analysisLeaf, url: "Functional_Analysis/9Linear/" },
        { id: "FA-Ch10", label: isEn ? "10. Open Mapping Theorem and Closed Graph Theorem" : "10.开映射与闭图定理", parentId: "FA", radius: 8, color: colors.analysisLeaf, url: "Functional_Analysis/10Open_Closed/" },
        { id: "FA-Ch11", label: isEn ? "11. Dual Spaces" : "11.对偶空间", parentId: "FA", radius: 8, color: colors.analysisLeaf, url: "Functional_Analysis/11dual/" },
        { id: "FA-Ch12", label: isEn ? "12. Hahn-Banach Theorem" : "12.Hahn-Banach定理", parentId: "FA", radius: 8, color: colors.analysisLeaf, url: "Functional_Analysis/12Hahn_Banach/" },
        // 渐近统计章节
        { id: "AS-Intro", label: isEn ? "0. Intro" : "课程简介", parentId: "AS", radius: 8, color: colors.analysisLeaf, url: "Asymptotic-statistics/0Intro/" },
        { id: "AS-Ch1", label: isEn ? "1. Convergence" : "1.随机收敛", parentId: "AS", radius: 8, color: colors.analysisLeaf, url: "Asymptotic-statistics/1Stochastic_convergence/" },
        { id: "AS-Ch2", label: isEn ? "2. Char Functions" : "2.特征函数", parentId: "AS", radius: 8, color: colors.analysisLeaf, url: "Asymptotic-statistics/2asymptotic_efficiency/" },
        { id: "AS-Ch3", label: isEn ? "3. CLT (I)" : "3.中心极限定理（一）", parentId: "AS", radius: 8, color: colors.analysisLeaf, url: "Asymptotic-statistics/3CLT/" },
        { id: "AS-Ch4", label: isEn ? "4. CLT (II)" : "4.中心极限定理（二）", parentId: "AS", radius: 8, color: colors.analysisLeaf, url: "Asymptotic-statistics/4CLT/" },
        { id: "AS-Ch5", label: isEn ? "5. Weak Dependence" : "5.弱相关数据", parentId: "AS", radius: 8, color: colors.analysisLeaf, url: "Asymptotic-statistics/5Weakly_Dep/" },
        { id: "AS-Ch6", label: isEn ? "6. Stationary Process" : "6.平稳过程", parentId: "AS", radius: 8, color: colors.analysisLeaf, url: "Asymptotic-statistics/6Stantionary_process/" },
        { id: "AS-Ch7", label: isEn ? "7. Delta Method" : "7.Delta方法", parentId: "AS", radius: 8, color: colors.analysisLeaf, url: "Asymptotic-statistics/7Delta_Method/" },
        { id: "AS-Ch8", label: isEn ? "8. Method of Moments" : "8.矩估计", parentId: "AS", radius: 8, color: colors.analysisLeaf, url: "Asymptotic-statistics/12ME/" },
        { id: "AS-Ch9", label: isEn ? "9. MLE I" : "9.极大似然估计（一）", parentId: "AS", radius: 8, color: colors.analysisLeaf, url: "Asymptotic-statistics/8MLE/" },
        { id: "AS-Ch10", label: isEn ? "10. MLE II" : "10.极大似然估计（二）", parentId: "AS", radius: 8, color: colors.analysisLeaf, url: "Asymptotic-statistics/9MLE/" },
        { id: "AS-Ch11", label: isEn ? "11. M-Estimation and Z-Estimation I" : "11.M估计和Z估计 (一)", parentId: "AS", radius: 8, color: colors.analysisLeaf, url: "Asymptotic-statistics/10M&Z1/" },
        { id: "AS-Ch12", label: isEn ? "12. M-Estimation and Z-Estimation II" : "12.M估计和Z估计 (二)", parentId: "AS", radius: 8, color: colors.analysisLeaf, url: "Asymptotic-statistics/11M&Z2/" },
        { id: "AS-Ch13", label: isEn ? "13.  U-Statistics I" : "13.U统计量(一)", parentId: "AS", radius: 8, color: colors.analysisLeaf, url: "Asymptotic-statistics/13U/" },
        { id: "AS-Ch14", label: isEn ? "14.  U-Statistics II" : "14.U统计量(二)", parentId: "AS", radius: 8, color: colors.analysisLeaf, url: "Asymptotic-statistics/14U/" },
        // 随机微分方程章节
        { id: "SDE-Intro", label: isEn ? "0. Intro" : "课程简介", parentId: "SDE", radius: 8, color: colors.analysisLeaf, url: "Stochastic-Differential-Equation/0Intro/" },
        { id: "SDE-Ch1", label: isEn ? "1. Cond Expectation" : "1.条件期望", parentId: "SDE", radius: 8, color: colors.analysisLeaf, url: "Stochastic-Differential-Equation/1conditional_expectation/" },
        { id: "SDE-Ch2", label: isEn ? "2. Brownian Motion" : "2.布朗运动", parentId: "SDE", radius: 8, color: colors.analysisLeaf, url: "Stochastic-Differential-Equation/2Brown_Motion/" },
        { id: "SDE-Ch3", label: isEn ? "3. Stoch Integral" : "3.随机积分", parentId: "SDE", radius: 8, color: colors.analysisLeaf, url: "Stochastic-Differential-Equation/3stochastic_integrate/" },
        { id: "SDE-Ch4", label: isEn ? "4. Ito Integral" : "4.伊藤积分", parentId: "SDE", radius: 8, color: colors.analysisLeaf, url: "Stochastic-Differential-Equation/4Ito_integrate/" },
        { id: "SDE-Ch5", label: isEn ? "5. Multivariate Ito" : "5.多元伊藤积分", parentId: "SDE", radius: 8, color: colors.analysisLeaf, url: "Stochastic-Differential-Equation/5SDEs/" },
        { id: "SDE-Ch6", label: isEn ? "6. SDEs" : "6.随机微分方程", parentId: "SDE", radius: 8, color: colors.analysisLeaf, url: "Stochastic-Differential-Equation/6SDEs/" },
        { id: "SDE-Ch7", label: isEn ? "7. Option Pricing" : "7.欧式期权定价", parentId: "SDE", radius: 8, color: colors.analysisLeaf, url: "Stochastic-Differential-Equation/7Option_Pricing/" },
        
        // 点集拓扑全部章节
        { id: "Top-Intro", label: isEn ? "0. Intro" : "课程简介", parentId: "Top", radius: 8, color: colors.geometry, url: "Topology/0Intro/" },
        { id: "Top-Ch1", label: isEn ? "1. Set Theory" : "1.集合论", parentId: "Top", radius: 8, color: colors.geometry, url: "Topology/1/" },
        { id: "Top-Ch2", label: isEn ? "2. Functions and Relations" : "2.函数与关系", parentId: "Top", radius: 8, color: colors.geometry, url: "Topology/2/" },
        { id: "Top-Ch3", label: isEn ? "3. Countable Set" : "3.可数集", parentId: "Top", radius: 8, color: colors.geometry, url: "Topology/3/" },
        { id: "Top-Ch4", label: isEn ? "4. Open Sets and Topology" : "4.开集与拓扑", parentId: "Top", radius: 8, color: colors.geometry, url: "Topology/4/" },
        { id: "Top-Ch5", label: isEn ? "5. Induced Topology" : "5.诱导拓扑", parentId: "Top", radius: 8, color: colors.geometry, url: "Topology/5/" },
        { id: "Top-Ch6", label: isEn ? "6.Continuity and Homeomorphisms" : "6.连续与同胚", parentId: "Top", radius: 8, color: colors.geometry, url: "Topology/6/" },
        { id: "Top-Ch7", label: isEn ? "7. Connectedness" : "7.连通性", parentId: "Top", radius: 8, color: colors.geometry, url: "Topology/7/" },
        { id: "Top-Ch8", label: isEn ? "8. Compactness" : "8.紧性", parentId: "Top", radius: 8, color: colors.geometry, url: "Topology/8/" },
        { id: "Top-Ch9", label: isEn ? "9. Compactness in Euclidean Spaces" : "9.欧氏空间紧性", parentId: "Top", radius: 8, color: colors.geometry, url: "Topology/9/" },
        { id: "Top-Ch10", label: isEn ? "10. Sequential Compactness" : "10.序列紧", parentId: "Top", radius: 8, color: colors.geometry, url: "Topology/10/" },
        { id: "Top-Ch11", label: isEn ? "11. One-Point Compactification" : "11.单点紧化", parentId: "Top", radius: 8, color: colors.geometry, url: "Topology/11/" },
        { id: "Top-Ch12", label: isEn ? "12. Countability" : "可数性", parentId: "Top", radius: 8, color: colors.geometry, url: "Topology/12/" },
        { id: "Top-Ch13", label: isEn ? "13. Separation Axioms" : "分离公理", parentId: "Top", radius: 8, color: colors.geometry, url: "Topology/13/" },
        //实分析全部章节
        { id: "RA-Intro", label: isEn ? "0. Intro" : "实分析课程概览", parentId: "RA", radius: 8, color: colors.analysisSub, url: "Real_analysis/0Intro/" },
        { id: "RA-Ch1", label: isEn ? "1. Sigma-Algebras" : "1.sigma-代数", parentId: "RA", radius: 8, color: colors.analysisSub, url: "Real_analysis/1_sigma_algebra/" },
        { id: "RA-Ch2", label: isEn ? "2. Outer Measure" : "2.外测度", parentId: "RA", radius: 8, color: colors.analysisSub, url: "Real_analysis/2_outer_measure/" },
        { id: "RA-Ch3", label: isEn ? "3. Lebesgue Measure" : "3.Lebesgue测度", parentId: "RA", radius: 8, color: colors.analysisSub, url: "Real_analysis/3_ls_measure/" },
        { id: "RA-Ch4", label: isEn ? "4. Measurable Functions I" : "4.可测函数（一）", parentId: "RA", radius: 8, color: colors.analysisSub, url: "Real_analysis/4_measurable_functions/" },
        { id: "RA-Ch5", label: isEn ? "5. Measurable Functions II" : "5.可测函数（二）", parentId: "RA", radius: 8, color: colors.analysisSub, url: "Real_analysis/5_limits_decomposition/" },
        { id: "RA-Ch6", label: isEn ? "6. Integration Construction" : "6.积分构造", parentId: "RA", radius: 8, color: colors.analysisSub, url: "Real_analysis/6_integration_construction/" },
        { id: "RA-Ch7", label: isEn ? "7. Convergence Theorems" : "7.收敛定理", parentId: "RA", radius: 8, color: colors.analysisSub, url: "Real_analysis/7_convergence_theorems/" },
        { id: "RA-Ch8", label: isEn ? "8. Fubini Theorems" : "8.Fubini定理", parentId: "RA", radius: 8, color: colors.analysisSub, url: "Real_analysis/8_fubini_theorems/" },
        { id: "RA-Ch9", label: isEn ? "9. Hahn-Jordan" : "9.Hahn-Jordan分解", parentId: "RA", radius: 8, color: colors.analysisSub, url: "Real_analysis/9_hahn_jordan/" },
        { id: "RA-Ch10", label: isEn ? "10. Radon-Nikodym" : "10.Radon-Nikodym定理", parentId: "RA", radius: 8, color: colors.analysisSub, url: "Real_analysis/10_radon_nikodym/" },
        { id: "RA-Ch11", label: isEn ? "11. FTC in R" : "11.微积分基本定理", parentId: "RA", radius: 8, color: colors.analysisSub, url: "Real_analysis/11_differentiation_ir/" },
        { id: "RA-Ch12", label: isEn ? "12. Lp Spaces" : "12.L^p空间", parentId: "RA", radius: 8, color: colors.analysisSub, url: "Real_analysis/12_lp_spaces/" },

        //统计学习
        { id: "SL-Intro", label: isEn ? "0. Intro" : "课程简介", parentId: "SL", radius: 8, color: colors.probSub, url: "Statistical_Learning/0Intro/" },
        { id: "SL-Ch1", label: isEn ? "1. KNN" : "1.KNN", parentId: "SL", radius: 8, color: colors.probSub, url: "Statistical_Learning/1KNN/" },
        { id: "SL-Ch2", label: isEn ? "2. Linear Regression" : "2.线性回归", parentId: "SL", radius: 8, color: colors.probSub, url: "Statistical_Learning/2Linear_Regression/" }, 
        { id: "SL-Ch3", label: isEn ? "3. Penalized Regression I" : "3.惩罚回归（一）", parentId: "SL", radius: 8, color: colors.probSub, url: "Statistical_Learning/3Penalized/" },
        { id: "SL-Ch4", label: isEn ? "4. Penalized Regression II" : "4.惩罚回归（二）", parentId: "SL", radius: 8, color: colors.probSub, url: "Statistical_Learning/4Penalized/" },
        { id: "SL-Ch5", label: isEn ? "5. Classification" : "5.分类问题", parentId: "SL", radius: 8, color: colors.probSub, url: "Statistical_Learning/5Classification/" },
        { id: "SL-Ch6", label: isEn ? "6. SVM" : "6.SVM", parentId: "SL", radius: 8, color: colors.probSub, url: "Statistical_Learning/6SVM/" },
        { id: "SL-Ch7", label: isEn ? "7. RKHS" : "7.核希尔伯特空间", parentId: "SL", radius: 8, color: colors.probSub, url: "Statistical_Learning/7RKHS/" },
        { id: "SL-Ch8", label: isEn ? "8. Kernel Methods" : "8.核方法", parentId: "SL", radius: 8, color: colors.probSub, url: "Statistical_Learning/8Kernel_Density/" },
        { id: "SL-Ch9", label: isEn ? "9. Random Forest" : "9.随机森林", parentId: "SL", radius: 8, color: colors.probSub, url: "Statistical_Learning/9Tree_Model/" },
        { id: "SL-Ch10", label: isEn ? "10. Causal Inference" : "10.因果推断", parentId: "SL", radius: 8, color: colors.probSub, url: "Statistical_Learning/10Causal_Inference/" },
    ];

    const allLinks = [
        { source: "Math-Core", target: "Analysis" },
        { source: "Math-Core", target: "Algebra" },
        { source: "Math-Core", target: "Geometry" },
        { source: "Math-Core", target: "Prob-Stat" },
        { source: "Math-Core", target: "Applied-Comp" },

        { source: "Analysis", target: "FA" },
        { source: "Prob-Stat", target: "AS" },
        { source: "Prob-Stat", target: "LDP" },
        { source: "Prob-Stat", target: "SDE" },
        { source: "Analysis", target: "SDE" },
        { source: "Geometry", target: "Top" },
        { source: "Analysis", target: "RA" },
        { source: "Prob-Stat", target: "SL" },

        { source: "FA", target: "FA-Intro" }, { source: "FA", target: "FA-Ch1" }, { source: "FA", target: "FA-Ch2" }, { source: "FA", target: "FA-Ch3" },
        { source: "FA", target: "FA-Ch4" }, { source: "FA", target: "FA-Ch5" }, { source: "FA", target: "FA-Ch6" }, { source: "FA", target: "FA-Ch7" },
        { source: "FA", target: "FA-Ch8" }, { source: "FA", target: "FA-Ch9" }, { source: "FA", target: "FA-Ch10" }, { source: "FA", target: "FA-Ch11" }, { source: "FA", target: "FA-Ch12" },
        { source: "AS", target: "AS-Intro" }, { source: "AS", target: "AS-Ch1" }, { source: "AS", target: "AS-Ch2" }, { source: "AS", target: "AS-Ch3" },
        { source: "AS", target: "AS-Ch4" }, { source: "AS", target: "AS-Ch5" }, { source: "AS", target: "AS-Ch6" }, { source: "AS", target: "AS-Ch7" },
        { source: "AS", target: "AS-Ch8" }, { source: "AS", target: "AS-Ch9" }, { source: "AS", target: "AS-Ch10" }, { source: "AS", target: "AS-Ch11" },
        { source: "AS", target: "AS-Ch12" }, { source: "AS", target: "AS-Ch13" }, { source: "AS", target: "AS-Ch14" },
        { source: "SDE", target: "SDE-Intro" }, { source: "SDE", target: "SDE-Ch1" }, { source: "SDE", target: "SDE-Ch2" }, { source: "SDE", target: "SDE-Ch3" },
        { source: "SDE", target: "SDE-Ch4" }, { source: "SDE", target: "SDE-Ch5" }, { source: "SDE", target: "SDE-Ch6" }, { source: "SDE", target: "SDE-Ch7" },
        { source: "Top", target: "Top-Intro" }, { source: "Top", target: "Top-Ch1" }, { source: "Top", target: "Top-Ch2" }, { source: "Top", target: "Top-Ch3" },
        { source: "Top", target: "Top-Ch4" }, { source: "Top", target: "Top-Ch5" }, { source: "Top", target: "Top-Ch6" }, { source: "Top", target: "Top-Ch7" },
        { source: "Top", target: "Top-Ch8" }, { source: "Top", target: "Top-Ch9" }, { source: "Top", target: "Top-Ch10" }, { source: "Top", target: "Top-Ch11" },
        { source: "Top", target: "Top-Ch12" }, { source: "Top", target: "Top-Ch13" },
        { source: "RA", target: "RA-Intro" }, { source: "RA", target: "RA-Ch1" }, { source: "RA", target: "RA-Ch2" }, { source: "RA", target: "RA-Ch3" },
        { source: "RA", target: "RA-Ch4" }, { source: "RA", target: "RA-Ch5" }, { source: "RA", target: "RA-Ch6" }, { source: "RA", target: "RA-Ch7" },
        { source: "RA", target: "RA-Ch8" }, { source: "RA", target: "RA-Ch9" }, { source: "RA", target: "RA-Ch10" }, { source: "RA", target: "RA-Ch11" },
        { source: "RA", target: "RA-Ch12" },
        { source: "SL", target: "SL-Intro" }, { source: "SL", target: "SL-Ch1" }, { source: "SL", target: "SL-Ch2" }, { source: "SL", target: "SL-Ch3" },
        { source: "SL", target: "SL-Ch4" }, { source: "SL", target: "SL-Ch5" }, { source: "SL", target: "SL-Ch6" }, { source: "SL", target: "SL-Ch7" },
        { source: "SL", target: "SL-Ch8" }, { source: "SL", target: "SL-Ch9" }, { source: "SL", target: "SL-Ch10" },

    ];

    // --- 2. 绘图逻辑 ---
    const svg = d3.select("#starry-map").append("svg")
        .attr("width", "100%")
        .attr("height", "100%")
        .attr("viewBox", [0, 0, width, height]);

    const linkGroup = svg.append("g").attr("stroke", "rgba(255,255,255,0.25)").attr("stroke-width", 2.2);
    const nodeGroup = svg.append("g");

    let simulation = d3.forceSimulation()
        .force("link", d3.forceLink().id(d => d.id).distance(d => {
            if (d.target.radius === 30) return 180;
            if (d.target.radius === 15) return 110;
            return 60;
        }))
        .force("charge", d3.forceManyBody().strength(-400))
        .force("center", d3.forceCenter(width / 2, height / 2))
        .force("collide", d3.forceCollide().radius(d => d.radius + 25));

    function update(sourceNode) {
        const visibleNodes = [];
        const isVisible = (node) => {
            if (!node.parentId) return true;
            const parent = allNodes.find(n => n.id === node.parentId);
            return parent && !parent.collapsed && isVisible(parent);
        };
        
        allNodes.forEach(n => { if (isVisible(n)) visibleNodes.push(n); });

        const visibleLinks = allLinks.filter(l => {
            const srcId = typeof l.source === "object" ? l.source.id : l.source;
            const tgtId = typeof l.target === "object" ? l.target.id : l.target;
            return visibleNodes.some(n => n.id === srcId) && visibleNodes.some(n => n.id === tgtId);
        });

        const link = linkGroup.selectAll("line").data(visibleLinks, d => {
            const src = typeof d.source === "object" ? d.source.id : d.source;
            const tgt = typeof d.target === "object" ? d.target.id : d.target;
            return src + "-" + tgt;
        });
        link.exit().remove();
        const linkEnter = link.enter().append("line");
        const linkMerge = linkEnter.merge(link);

        const node = nodeGroup.selectAll("g.node").data(visibleNodes, d => d.id);
        node.exit().transition().duration(300).attr("r", 0).remove();

        const nodeEnter = node.enter().append("g")
            .attr("class", "node")
            .style("cursor", "pointer")
            .attr("transform", d => {
                if (sourceNode) { d.x = sourceNode.x; d.y = sourceNode.y; }
                return `translate(${d.x || width/2},${d.y || height/2})`;
            })
            .call(d3.drag().on("start", dragstarted).on("drag", dragged).on("end", dragended))
            .on("click", (event, d) => {
                if (d.url) { 
                    window.location.href = d.url; 
                } 
                else { d.collapsed = !d.collapsed; update(d); }
            });

        nodeEnter.append("circle")
            .attr("r", d => d.radius)
            .attr("fill", d => d.color)
            .attr("stroke", d => d.collapsed === false ? "rgba(255,255,255,0.8)" : "none")
            .attr("stroke-width", 2)
            .style("filter", "drop-shadow(0 0 10px rgba(255,255,255,0.2))")
            .on("mouseenter", function() {
                d3.select(this).transition().duration(200).attr("r", d => d.radius * 1.25).style("filter", "drop-shadow(0 0 25px rgba(255,255,255,0.6))");
            })
            .on("mouseleave", function() {
                d3.select(this).transition().duration(200).attr("r", d => d.radius).style("filter", "drop-shadow(0 0 10px rgba(255,255,255,0.2))");
            });

        nodeEnter.append("text")
            .text(d => d.label || d.id)
            .attr("x", d => d.radius + 10)
            .attr("y", 4)
            .style("fill", "#ecf0f1")
            .style("font-size", d => d.radius < 10 ? "11px" : "14px")
            .style("font-family", "'JetBrains Mono', 'PingFang SC', sans-serif")
            .style("text-shadow", "0 2px 4px rgba(0,0,0,0.8)");

        const nodeMerge = nodeEnter.merge(node);
        simulation.nodes(visibleNodes).on("tick", () => {
            linkMerge.attr("x1", d => d.source.x).attr("y1", d => d.source.y).attr("x2", d => d.target.x).attr("y2", d => d.target.y);
            nodeMerge.attr("transform", d => `translate(${d.x},${d.y})`);
        });
        simulation.force("link").links(visibleLinks);
        simulation.alpha(0.5).restart(); 
    }

    // 拖拽逻辑
    function dragstarted(event, d) { if (!event.active) simulation.alphaTarget(0.3).restart(); d.fx = d.x; d.fy = d.y; }
    function dragged(event, d) { d.fx = event.x; d.fy = event.y; }
    function dragended(event, d) { if (!event.active) simulation.alphaTarget(0); d.fx = null; d.fy = null; }

    update();

    // 🌟 核心星球呼吸动画
    setTimeout(() => {
        const coreCircle = d3.select("#starry-map svg circle");
        if (!coreCircle.empty()) {
            (function pulse() {
                coreCircle.transition().duration(2500).ease(d3.easeSinInOut).attr("r", 58).style("filter", "drop-shadow(0 0 30px rgba(255,215,0,0.7))")
                    .transition().duration(2500).ease(d3.easeSinInOut).attr("r", 50).style("filter", "drop-shadow(0 0 15px rgba(255,215,0,0.3))").on("end", pulse);
            })();
        }
    }, 1000);
});