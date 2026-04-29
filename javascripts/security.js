// 使用更稳健的路径获取方式
const siteRoot = window.location.origin + '/Notes/';

// 监听文档内容变化
const observer = new MutationObserver((mutations) => {
    const pageText = document.body.innerText;
    
    // 匹配密码错误的提示语
    const triggerText = "看来我们没有那么默契";
    
    if (pageText.includes(triggerText)) {
        observer.disconnect();
        
        // ⚽️ 图片替换指南：
        // 建议你在网上找一张搞笑的裁判看VAR的GIF（比如 Lahoz 或者经典的画框手势）
        // 把它保存到你的 mkdocs 项目的 docs/assets/ 目录下，命名为 var_check.gif
        // 然后把下面这行换成本地路径：const varImageUrl = siteRoot + "assets/var_check.gif";
        // 这里我先用一个网络图库的通用裁判GIF占位
        const varImageUrl = siteRoot + "assets/var_check.gif";
        // 强行覆盖整个页面的 HTML，营造转播屏幕的氛围
        document.body.innerHTML = `
            <div style="
                display: flex; 
                flex-direction: column; 
                align-items: center; 
                justify-content: center; 
                height: 100vh; 
                background-color: #111; /* 纯黑背景，更有氛围 */
                color: white; 
                font-family: 'Helvetica Neue', Helvetica, Arial, sans-serif; 
                text-align: center;
                margin: 0;
            ">
                <img src="${varImageUrl}" alt="VAR Check" style="
                    max-width: 300px; 
                    border-radius: 8px; 
                    margin-bottom: 25px; 
                    box-shadow: 0 0 20px rgba(255,255,255,0.1);
                ">
                <h1 style="
                    color: #f1c40f; /* 裁判黄牌的颜色 */
                    margin-bottom: 15px; 
                    font-weight: 900; 
                    letter-spacing: 2px;
                    font-size: 2.5rem;
                ">📺 VAR CHECK IN PROGRESS...</h1>
                
                <h2 style="
                    color: #e74c3c; /* 红牌的颜色 */
                    font-weight: normal;
                    font-size: 1.5rem;
                ">❌ OFFSIDE! 越位在先，进球无效。</h2>
            </div>
        `;
        
        // ⏱️ 延长至 3 秒
        // 因为有图片和梗，1.5秒太快了，对方还没看清怎么回事就被踢出去了
        // 留 3 秒钟让他们把图看清，笑出来，然后再平滑地踢回主页
        setTimeout(() => {
            window.location.href = siteRoot;
        }, 5000);
    }
});

// 开始监听整个 body 的变化
observer.observe(document.body, { 
    childList: true, 
    subtree: true 
});