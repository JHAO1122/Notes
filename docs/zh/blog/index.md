<div style="width: 100%; height: 180px; overflow: hidden; border-radius: 12px; margin-bottom: 24px; box-shadow: 0 4px 10px rgba(0,0,0,0.1);" markdown>

![Mes Que Un Club](assets/camp_nou.jpg){: style="width: 100%; height: 180px; object-fit: cover; object-position: center 85%;" }

</div>

# ✍️ 个人随笔

## 关于我

> "我若等你，又何惧几个秋。"

我不仅会与概率测度和统计方法打交道，我还是一个狂热的足球迷。我有十一年的球龄（大概四年级起），但是现在由于事情繁多、球友陌生等因素已与绿茵场渐行渐远。我始终怀念初中和同班同学踢球的日子，那大抵是我人生目前为止最无忧无虑、没心没肺的时光。还记得我们班在班级联赛斩获冠军，我以五场比赛四粒进球位列全校射手榜第二名。

虽然现在踢球的机会逐渐变少，但就像一句很经典的话说的那样——一个男孩会在他十几岁时选择一只球队并从此坚守一生，我选择的是巴塞罗那。遗憾的是，我开始看足球的那一年巴萨拿到欧冠冠军，直到现在，巴萨都没有再拿过欧冠冠军。不过十年饮冰，难凉热血，我始终坚信那一天终会到来。

![FC Barcelona](assets/FCB.png)

我特意在个人笔记主页开设随笔一栏，一来是为了记录生活，给自己放空大脑的空间；二来是可以在此发发牢骚，面对巨大的学业、科研压力，不吐不快。

---

!!! quote "入口声明"
    这是我用来记录日常的树洞，学术与笔记相关的内容请移步别处。这里只谈热爱，只讲情怀，不推公式。

## 📜 记忆存档

* **2026-04-15** | [男人的一生只哭十一次](./posts/private/2026-4-15.md) 🔒
    > *关于 25-26 欧冠四分之一决赛、“下次再来”的循环。*

* **2026-04-23** | [何时收敛](./posts/private/2026-4-23.md) 🔒
    > *对于27fall申请的疲惫、一些对未来的想法以及期中总结*

* **2026-05-05** | [被烫伤的孩子仍然爱火](./posts/public/2026-5-5.md) 
    > *二十岁的青春就在那里，眼里都是十岁时的影子*

* **2026-06-13** | [暂时留白的情感总结](./posts/private/2026-6-13.md) 🔒
    > *考试周前不久分手了，暂时留白...*


<div style="background-color: rgba(67, 161, 213, 0.05); border-radius: 8px; padding: 20px; text-align: center; box-shadow: 0 2px 8px rgba(0,0,0,0.05); margin: 20px 0; border: 1px solid rgba(67, 161, 213, 0.2);">
  <h3 style="margin-top: 0; font-weight: 600;">⭐⭐⭐ 距离阿根廷捧起大力神杯已过去</h3>
  <div id="arg-countup" style="display: flex; justify-content: center; gap: 20px; font-size: 1.6rem; font-weight: bold; color: #43A1D5; margin-top: 15px;">
    <div><span id="arg-days">--</span> <span style="font-size: 0.9rem; font-weight: normal; color: #757575;">天</span></div>
    <div><span id="arg-hours">--</span> <span style="font-size: 0.9rem; font-weight: normal; color: #757575;">时</span></div>
    <div><span id="arg-minutes">--</span> <span style="font-size: 0.9rem; font-weight: normal; color: #757575;">分</span></div>
    <div><span id="arg-seconds">--</span> <span style="font-size: 0.9rem; font-weight: normal; color: #757575;">秒</span></div>
  </div>
  <div style="margin-top: 10px; font-size: 0.85rem; color: #D4AF37;">"Muchachos, ahora nos volvimos a ilusionar..."</div>
</div>

<script>
  // 锁定北京时间：2022年12月19日 02:00:00 (夺冠捧杯时刻)
  const argDate = new Date("2022-12-19T02:00:00+08:00").getTime();
  
  setInterval(function() {
    const now = new Date().getTime();
    const distance = now - argDate; // 注意：这里是现在减去过去
    
    document.getElementById("arg-days").innerText = Math.floor(distance / (1000 * 60 * 60 * 24));
    document.getElementById("arg-hours").innerText = Math.floor((distance % (1000 * 60 * 60 * 24)) / (1000 * 60 * 60)).toString().padStart(2, '0');
    document.getElementById("arg-minutes").innerText = Math.floor((distance % (1000 * 60 * 60)) / (1000 * 60)).toString().padStart(2, '0');
    document.getElementById("arg-seconds").innerText = Math.floor((distance % (1000 * 60)) / 1000).toString().padStart(2, '0');
  }, 1000);
</script>


<div style="background-color: rgba(128, 128, 128, 0.05); border-radius: 8px; padding: 20px; text-align: center; box-shadow: 0 2px 8px rgba(0,0,0,0.05); margin: 20px 0; border: 1px solid rgba(128,128,128,0.2);">
  <h3 style="margin-top: 0; font-weight: 600;">🏆 距离 2026 美加墨世界杯开幕还有</h3>
  <div id="world-cup-countdown" style="display: flex; justify-content: center; gap: 20px; font-size: 1.6rem; font-weight: bold; color: #D32F2F; margin-top: 15px;">
    <div><span id="wc-days">--</span> <span style="font-size: 0.9rem; font-weight: normal; color: #757575;">天</span></div>
    <div><span id="wc-hours">--</span> <span style="font-size: 0.9rem; font-weight: normal; color: #757575;">时</span></div>
    <div><span id="wc-minutes">--</span> <span style="font-size: 0.9rem; font-weight: normal; color: #757575;">分</span></div>
    <div><span id="wc-seconds">--</span> <span style="font-size: 0.9rem; font-weight: normal; color: #757575;">秒</span></div>
  </div>
</div>

<script>
  // 使用 ISO 8601 格式强制指定为北京时间 (UTC+8) 的 2026年6月12日 03:00:00
  const wcDate = new Date("2026-06-12T03:00:00+08:00").getTime();
  
  const wcTimer = setInterval(function() {
    const now = new Date().getTime();
    const distance = wcDate - now;
    
    if (distance < 0) {
      clearInterval(wcTimer);
      document.getElementById("world-cup-countdown").innerHTML = "⚽ 比赛已经开始啦！熬夜看球走起！";
      document.getElementById("world-cup-countdown").style.color = "#2E7D32"; // 变绿色
      return;
    }
    
    document.getElementById("wc-days").innerText = Math.floor(distance / (1000 * 60 * 60 * 24));
    document.getElementById("wc-hours").innerText = Math.floor((distance % (1000 * 60 * 60 * 24)) / (1000 * 60 * 60)).toString().padStart(2, '0');
    document.getElementById("wc-minutes").innerText = Math.floor((distance % (1000 * 60 * 60)) / (1000 * 60)).toString().padStart(2, '0');
    document.getElementById("wc-seconds").innerText = Math.floor((distance % (1000 * 60)) / 1000).toString().padStart(2, '0');
  }, 1000);
</script>