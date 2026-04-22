import arxiv
import openai
import os
import json
from datetime import datetime

# 配置与路径
DB_PATH = "scripts/processed_papers.json"
MD_PATH = "docs/zh/research/arxiv_papers.md"
DEEPSEEK_KEY = os.getenv("DEEPSEEK_API_KEY")

client_ds = openai.OpenAI(api_key=DEEPSEEK_KEY, base_url="https://api.deepseek.com")

def load_db():
    if os.path.exists(DB_PATH):
        with open(DB_PATH, "r", encoding="utf-8") as f:
            return json.load(f)
    return {}

def save_db(db):
    with open(DB_PATH, "w", encoding="utf-8") as f:
        json.dump(db, f, ensure_ascii=False, indent=4)

def get_ai_summary(title, abstract):
    # 强化 Prompt，要求 AI 站在数院统计系的视角
    prompt = (
        "你是一个深耕数理统计领域的专家。请用两句中文概括该论文创新点，"
        "重点突出其统计推导、构造的评分函数或理论性质，不要泛泛而谈应用场景。"
        "不要使用 Markdown 标题或复杂的列表嵌套。"
        f"\n标题: {title}\n摘要: {abstract}"
    )
    try:
        response = client_ds.chat.completions.create(
            model="deepseek-chat",
            messages=[{"role": "user", "content": prompt}],
            temperature=0.3
        )
        return response.choices[0].message.content
    except Exception as e:
        return f"AI 总结生成失败: {e}"

# 1. 运行抓取
db = load_db()

# --- 关键修改：加入统计学分类过滤 (stat.ME=方法论, stat.TH=理论, stat.ML=统计机器学习) ---
# 这样可以过滤掉大量的 cs.RO (机器人) 或 cs.AI (工程应用) 的水文
queries = {
    "Knockoff": "all:knockoff AND (cat:stat.ME OR cat:stat.TH OR cat:stat.ML)",
    "Conformal": 'all:"conformal prediction" AND (cat:stat.ME OR cat:stat.TH OR cat:stat.ML)'
}

# 配置一个“礼貌”的客户端
client = arxiv.Client(page_size=10, delay_seconds=3.0, num_retries=5)

new_entries = ""
for tag, q in queries.items():
    # 关键修改 1：把 max_results 调大（比如 30），扩大我们“淘金”的沙池
    search = arxiv.Search(query=q, max_results=30, sort_by=arxiv.SortCriterion.SubmittedDate)
    
    # 关键修改 2：加一个计数器，记录今天这个标签下抓到了几篇全新的
    new_found_count = 0 
    
    for result in client.results(search):
        paper_id = result.entry_id.split('/')[-1]
        
        if paper_id in db:
            continue # 见过的，直接跳过，去看下一篇
            
        print(f"New Academic Paper Found: {result.title}")
        summary = get_ai_summary(result.title, result.summary)
        
        # 确保多行 summary 也能完美缩进
        indented_summary = summary.replace("\n", "\n    ")
        
        db[paper_id] = {
            "title": result.title,
            "date": str(result.published.date()),
            "tag": tag,
            "status": "todo"
        }
        
        new_entries += f"### {result.title} \n\n" 
        new_entries += f"- [ ] **分类**: {tag} | **日期**: {result.published.date()}\n"
        new_entries += f"- **链接**: [PDF]({result.entry_id})\n\n"
        new_entries += f'!!! note "AI 核心解读"\n\n'
        new_entries += f"    {indented_summary}\n\n"
        
        # 关键修改 3：每找到一篇全新的，计数器 +1。凑满 5 篇就收工
        new_found_count += 1
        if new_found_count >= 5:
            break
# 2. 更新文件
if new_entries:
    old_content = ""
    if os.path.exists(MD_PATH):
        with open(MD_PATH, "r", encoding="utf-8") as f:
            content = f.read()
            if "---" in content:
                old_content = content.split("---")[-1]
            else:
                old_content = content

    # 1. 定义你专属的抬头模板
    intro_text = """这里是我通过 GitHub 机器人脚本制作的论文速递版面，它北京时间每天早上八点会自动扫描 ArXiv 上的关于 Knockoff 与 Conformal Prediction 的统计学论文。
"""

    # 2. 构造完整的 header
    header = f"# 🛰️ Research Frontier\n\n{intro_text}\n\n> 更新于: {datetime.now().strftime('%Y-%m-%d')}\n\n---\n"
    
    with open(MD_PATH, "w", encoding="utf-8") as f:
        f.write(header + new_entries + old_content)
    
    save_db(db)