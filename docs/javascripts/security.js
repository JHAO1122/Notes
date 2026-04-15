// 使用更稳健的路径获取方式
const siteRoot = window.location.origin + '/Notes/';

// 监听文档内容变化
const observer = new MutationObserver((mutations) => {
    const pageText = document.body.innerText;
    
    // 建议只匹配前面的核心文字，避开末尾可能变动的三个点
    const triggerText = "看来我们没有那么默契";
    
    if (pageText.includes(triggerText)) {
    observer.disconnect();
    
    // 修改页面文字，让他知道自己被驱逐了，但不弹窗
    document.body.innerHTML = "<h1 style='text-align:center;margin-top:20%;color:red;'>⚠️ 安全警报：非法访问，启动强制退出...</h1>";
    
    // 1.5秒后自动踢出，不需要点确定
    setTimeout(() => {
        window.location.href = siteRoot;
    }, 1500);
}
});

// 开始监听整个 body 的变化
observer.observe(document.body, { 
    childList: true, 
    subtree: true 
});