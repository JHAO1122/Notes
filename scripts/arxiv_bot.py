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
    # (同之前的 AI 摘要函数)
    prompt = f"你是统计学专家，请用两句中文概括该论文创新点：\n标题: {title}\n摘要: {abstract}"
    response = client_ds.chat.completions.create(
        model="deepseek-chat",
        messages=[{"role": "user", "content": prompt}],
        temperature=0.3
    )
    return response.choices[0].message.content

# 1. 运行抓取
db = load_db()
queries = {
    "Knockoff": "all:knockoff",
    "Conformal": 'all:"conformal prediction"'
}

new_entries = ""
for tag, q in queries.items():
    search = arxiv.Search(query=q, max_results=5, sort_by=arxiv.SortCriterion.SubmittedDate)
    for result in arxiv.Client().results(search):
        paper_id = result.entry_id.split('/')[-1]
        
        # 去重逻辑
        if paper_id in db:
            continue
            
        print(f"New Paper Found: {result.title}")
        summary = get_ai_summary(result.title, result.summary)
        
        # 记录到数据库
        db[paper_id] = {
            "title": result.title,
            "date": str(result.published.date()),
            "tag": tag,
            "status": "todo" # 初始状态
        }
        
        # 构造 Markdown 条目 (使用 Task List 语法)
        new_entries += f"### - [ ] {result.title} \n"
        new_entries += f"- **分类**: {tag} | **日期**: {result.published.date()}\n"
        new_entries += f"- **链接**: [PDF]({result.entry_id})\n"
        new_entries += f"!!! note \"AI 核心解读\"\n    {summary}\n\n"

# 2. 更新文件（追加到页面顶部，保留旧内容）
if new_entries:
    old_content = ""
    if os.path.exists(MD_PATH):
        with open(MD_PATH, "r", encoding="utf-8") as f:
            old_content = f.read().split("---")[-1] # 保留分割线后的历史记录

    header = f"# 🛰️ Research Frontier\n\n> 更新于: {datetime.now().strftime('%Y-%m-%d')}\n\n---\n"
    with open(MD_PATH, "w", encoding="utf-8") as f:
        f.write(header + new_entries + old_content)
    
    save_db(db)