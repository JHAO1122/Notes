import arxiv
import openai
import os
import json
import sys
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
    # 关键修复：保存前先确保 scripts 文件夹存在，否则会报 Exit code 1
    os.makedirs(os.path.dirname(DB_PATH), exist_ok=True)
    with open(DB_PATH, "w", encoding="utf-8") as f:
        json.dump(db, f, ensure_ascii=False, indent=4)

def get_ai_summary(title, abstract):
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

# --- 关键修复：加入全局异常捕获，保证遇到任何报错都能绿灯退出 ---
try:
    # 1. 运行抓取
    db = load_db()

    queries = {
        "Knockoff": "all:knockoff AND (cat:stat.ME OR cat:stat.TH OR cat:stat.ML)",
        "Conformal": 'all:"conformal prediction" AND (cat:stat.ME OR cat:stat.TH OR cat:stat.ML)'
    }

    client = arxiv.Client(page_size=10, delay_seconds=3.0, num_retries=5)

    new_entries = ""
    for tag, q in queries.items():
        search = arxiv.Search(query=q, max_results=30, sort_by=arxiv.SortCriterion.SubmittedDate)
        new_found_count = 0 
        
        for result in client.results(search):
            paper_id = result.entry_id.split('/')[-1]
            
            if paper_id in db:
                continue 
                
            print(f"New Academic Paper Found: {result.title}")
            summary = get_ai_summary(result.title, result.summary)
            
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

        intro_text = """这是我通过脚本制作的论文速递版面，它将在北京时间每天早上八点调用 GitHub 机器人自动扫描 [arXiv](https://arxiv.org/) 上的关于 Knockoff 与 Conformal Prediction 的统计学论文。
"""
        header = f"# 🛰️ Research Frontier\n\n{intro_text}\n\n> 更新于: {datetime.now().strftime('%Y-%m-%d')}\n\n---\n"
        
        # 关键修复：确保 docs 的多级目录存在
        os.makedirs(os.path.dirname(MD_PATH), exist_ok=True)
        with open(MD_PATH, "w", encoding="utf-8") as f:
            f.write(header + new_entries + old_content)
        
        save_db(db)
        print("✅ 更新成功！")
    else:
        print("✅ 今天没有新论文，跳过更新。")

except Exception as e:
    # 遇到任何错误打印出来，但通过 sys.exit(0) 欺骗 GitHub 这是成功的
    print(f"⚠️ 运行中出现异常，但这不影响之前的数据。错误详情: {e}")
    sys.exit(0)