document.addEventListener("DOMContentLoaded", function() {
    const container = document.getElementById('starry-map');
    if (!container) return;

    const width = container.clientWidth;
    const height = container.clientHeight;

    // 🌍 自动检测当前语言环境
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

    // 🛠️ 路径修正逻辑：如果是在 en/ 或 zh/ 子目录下，跳转路径需要返回上一级
    const pathPrefix = "../";

    // --- 1. 定义整个宇宙的完整数据 ---
    // id 为固定标识符，label 为显示的名称
    const allNodes = [
        { id: "Math-Core", label: isEn ? "Mathematical Sciences" : "数学科学", radius: 50, color: colors.core, collapsed: false },

        { id: "Analysis", parentId: "Math-Core", label: isEn ? "Analysis" : "分析学", radius: 30, color: colors.analysis, collapsed: false },
        { id: "Algebra", parentId: "Math-Core", label: isEn ? "Algebra" : "代数学", radius: 30, color: colors.algebra, collapsed: true },
        { id: "Geometry", parentId: "Math-Core", label: isEn ? "Geometry" : "几何学", radius: 30, color: colors.geometry, collapsed: true },
        { id: "Prob-Stat", parentId: "Math-Core", label: isEn ? "Probability & Statistics" : "概率与统计", radius: 30, color: colors.prob, collapsed: false },
        { id: "Applied-Comp", parentId: "Math-Core", label: isEn ? "Applied Mathematics" : "应用与计算", radius: 30, color: colors.applied, collapsed: true },

        { id: "FA", parentId: "Analysis", label: isEn ? "Functional Analysis" : "泛函分析", radius: 15, color: colors.analysisSub, collapsed: true },
        { id: "AS", parentId: "Prob-Stat", label: isEn ? "Asymptotic Statistics" : "渐近统计", radius: 15, color: colors.probSub, collapsed: true },
        { id: "LDP", parentId: "Prob-Stat", label: isEn ? "LDP Seminar" : "大偏差原理", radius: 15, color: colors.probSub, collapsed: true },
        { id: "SDE", parentId: "Prob-Stat", label: isEn ? "Stochastic Diff Eq" : "随机微分方程", radius: 15, color: colors.probSub, collapsed: true },

        // 泛函的内部章节
        { id: "FA-Intro", label: isEn ? "0. Intro" : "课程简介", parentId: "FA", radius: 8, color: colors.analysisLeaf, url: pathPrefix + "Functional_Analysis/0Intro/" },
        { id: "FA-Ch1", label: isEn ? "1. Metric Spaces" : "1.度量空间", parentId: "FA", radius: 8, color: colors.analysisLeaf, url: pathPrefix + "Functional_Analysis/1distance_space/" },
        { id: "FA-Ch2", label: isEn ? "2. Linear Operators" : "2.线性算子", parentId: "FA", radius: 8, color: colors.analysisLeaf, url: pathPrefix + "Functional_Analysis/2completeness/" },
        { id: "FA-Ch3", label: isEn ? "3. Compactness" : "3.紧性", parentId: "FA", radius: 8, color: colors.analysisLeaf, url: pathPrefix + "Functional_Analysis/3compactness/" },

        // 渐近统计章节
        { id: "AS-Intro", label: isEn ? "0. Intro" : "课程简介", parentId: "AS", radius: 8, color: colors.analysisLeaf, url: pathPrefix + "Asymptotic-statistics/0Intro/" },
        { id: "AS-Ch1", label: isEn ? "1. Convergence" : "1.随机收敛", parentId: "AS", radius: 8, color: colors.analysisLeaf, url: pathPrefix + "Asymptotic-statistics/1Stochastic_convergence/" },
        { id: "AS-Ch2", label: isEn ? "2. Char Functions" : "2.特征函数", parentId: "AS", radius: 8, color: colors.analysisLeaf, url: pathPrefix + "Asymptotic-statistics/2asymptotic_efficiency/" },
        { id: "AS-Ch3", label: isEn ? "3. CLT (I)" : "3.中心极限定理（一）", parentId: "AS", radius: 8, color: colors.analysisLeaf, url: pathPrefix + "Asymptotic-statistics/3CLT/" },
        { id: "AS-Ch4", label: isEn ? "4. CLT (II)" : "4.中心极限定理（二）", parentId: "AS", radius: 8, color: colors.analysisLeaf, url: pathPrefix + "Asymptotic-statistics/4CLT/" },
        { id: "AS-Ch5", label: isEn ? "5. Weak Dependence" : "5.弱相关数据", parentId: "AS", radius: 8, color: colors.analysisLeaf, url: pathPrefix + "Asymptotic-statistics/5Weakly_Dep/" },
        { id: "AS-Ch6", label: isEn ? "6. Stationary Process" : "6.平稳过程", parentId: "AS", radius: 8, color: colors.analysisLeaf, url: pathPrefix + "Asymptotic-statistics/6Stantionary_process/" },
        { id: "AS-Ch7", label: isEn ? "7. Delta Method" : "7.Delta方法", parentId: "AS", radius: 8, color: colors.analysisLeaf, url: pathPrefix + "Asymptotic-statistics/7Delta_Method/" },
    
        // 随机微分方程章节
        { id: "SDE-Intro", label: isEn ? "0. Intro" : "课程简介", parentId: "SDE", radius: 8, color: colors.analysisLeaf, url: pathPrefix + "Stochastic-Differential-Equation/0Intro/" },
        { id: "SDE-Ch1", label: isEn ? "1. Cond Expectation" : "1.条件期望", parentId: "SDE", radius: 8, color: colors.analysisLeaf, url: pathPrefix + "Stochastic-Differential-Equation/1conditional_expectation/" },
        { id: "SDE-Ch2", label: isEn ? "2. Brownian Motion" : "2.布朗运动", parentId: "SDE", radius: 8, color: colors.analysisLeaf, url: pathPrefix + "Stochastic-Differential-Equation/2Brown_Motion/" },
        { id: "SDE-Ch3", label: isEn ? "3. Stoch Integral" : "3.随机积分", parentId: "SDE", radius: 8, color: colors.analysisLeaf, url: pathPrefix + "Stochastic-Differential-Equation/3stochastic_integrate/" },
        { id: "SDE-Ch4", label: isEn ? "4. Ito Integral" : "4.伊藤积分", parentId: "SDE", radius: 8, color: colors.analysisLeaf, url: pathPrefix + "Stochastic-Differential-Equation/4Ito_integrate/" },
        { id: "SDE-Ch5", label: isEn ? "5. Multivariate Ito" : "5.多元伊藤积分", parentId: "SDE", radius: 8, color: colors.analysisLeaf, url: pathPrefix + "Stochastic-Differential-Equation/5SDEs/" },
        { id: "SDE-Ch6", label: isEn ? "6. SDEs" : "6.随机微分方程", parentId: "SDE", radius: 8, color: colors.analysisLeaf, url: pathPrefix + "Stochastic-Differential-Equation/6SDEs/" },
        { id: "SDE-Ch7", label: isEn ? "7. Option Pricing" : "7.欧式期权定价", parentId: "SDE", radius: 8, color: colors.analysisLeaf, url: pathPrefix + "Stochastic-Differential-Equation/7Option_Pricing/" }
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

        { source: "FA", target: "FA-Intro" }, { source: "FA", target: "FA-Ch1" }, { source: "FA", target: "FA-Ch2" }, { source: "FA", target: "FA-Ch3" },
        { source: "AS", target: "AS-Intro" }, { source: "AS", target: "AS-Ch1" }, { source: "AS", target: "AS-Ch2" }, { source: "AS", target: "AS-Ch3" },
        { source: "AS", target: "AS-Ch4" }, { source: "AS", target: "AS-Ch5" }, { source: "AS", target: "AS-Ch6" }, { source: "AS", target: "AS-Ch7" },
        { source: "SDE", target: "SDE-Intro" }, { source: "SDE", target: "SDE-Ch1" }, { source: "SDE", target: "SDE-Ch2" }, { source: "SDE", target: "SDE-Ch3" },
        { source: "SDE", target: "SDE-Ch4" }, { source: "SDE", target: "SDE-Ch5" }, { source: "SDE", target: "SDE-Ch6" }, { source: "SDE", target: "SDE-Ch7" }
    ];

    // --- 2. 绘图逻辑 (已优化跳转和渲染) ---
    const svg = d3.select("#starry-map").append("svg")
        .attr("width", "100%")
        .attr("height", "100%")
        .attr("viewBox", [0, 0, width, height]);

    const linkGroup = svg.append("g").attr("stroke", "rgba(255,255,255,0.30)").attr("stroke-width", 2.5);
    const nodeGroup = svg.append("g");

    let simulation = d3.forceSimulation()
        .force("link", d3.forceLink().id(d => d.id).distance(d => {
            if (d.target.radius === 30) return 180;
            if (d.target.radius === 15) return 120;
            return 70;
        }))
        .force("charge", d3.forceManyBody().strength(-350))
        .force("center", d3.forceCenter(width / 2, height / 2))
        .force("collide", d3.forceCollide().radius(d => d.radius + 20));

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
                if (d.url) { window.location.href = d.url; } 
                else { d.collapsed = !d.collapsed; update(d); }
            });

        nodeEnter.append("circle")
            .attr("r", d => d.radius)
            .attr("fill", d => d.color)
            .attr("stroke", d => d.collapsed === false ? "rgba(255,255,255,0.8)" : "none")
            .attr("stroke-width", 2)
            .style("filter", "drop-shadow(0 0 10px rgba(255,255,255,0.2))")
            .on("mouseenter", function() {
                d3.select(this).transition().duration(200).attr("r", d => d.radius * 1.3).style("filter", "drop-shadow(0 0 25px rgba(255,255,255,0.6))");
            })
            .on("mouseleave", function() {
                d3.select(this).transition().duration(200).attr("r", d => d.radius).style("filter", "drop-shadow(0 0 10px rgba(255,255,255,0.2))");
            });

        nodeEnter.append("text")
            .text(d => d.label || d.id) // 优先使用 label
            .attr("x", d => d.radius + 8)
            .attr("y", 4)
            .style("fill", "#ecf0f1")
            .style("font-size", d => d.radius < 10 ? "11px" : "13px")
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

    // 🌟 核心星球呼吸光晕动画
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