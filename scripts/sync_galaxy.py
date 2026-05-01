import os
import re
import json
from pathlib import Path
from dotenv import load_dotenv
from openai import OpenAI

# 1. 配置加载
load_dotenv(Path(__file__).parent / ".env")
client = OpenAI(
    api_key=os.getenv("DEEPSEEK_API_KEY"),
    base_url="https://api.deepseek.com"
)

# 📂 路径配置 (已根据你的 docs/javascripts 结构修正)
PROJECT_ROOT = Path(__file__).parent.parent
DOCS_ZH = PROJECT_ROOT / "docs" / "zh"
JS_FILE = PROJECT_ROOT / "docs" / "javascripts" / "starry_map.js"

def scan_docs():
    """扫描 zh 文件夹下的所有课程和章节"""
    structure = {}
    # 排除不需要扫描的文件夹
    exclude = ['javascripts', 'overrides', 'assets', 'Barca']
    
    for folder in os.listdir(DOCS_ZH):
        folder_path = DOCS_ZH / folder
        if folder_path.is_dir() and folder not in exclude:
            # 这里的扫描逻辑会抓取文件夹下的所有 .md 文件
            files = [f.stem for f in folder_path.glob("*.md") if f.stem != "index"]
            structure[folder] = files
    return structure

def get_ai_suggestions(structure):
    """请求 AI 生成星图代码片段，并强制保留 1899 节点"""
    prompt = f"""
    我的 MkDocs 目录结构如下：
    {json.dumps(structure, ensure_ascii=False)}

    任务：为 D3.js 星空图生成 `allNodes` 和 `allLinks`。
    
    约束条件：
    1. 根节点 ID: "Math-Core"。
    2. 一级学科 ID: "Analysis", "Algebra", "Geometry", "Prob-Stat", "Applied-Comp"。
    3. ！！！绝对保护节点！！！：必须包含以下节点，严禁删除：
       {{ id: "Barca-1899", parentId: "Math-Core", label: "1899", radius: 10, color: colors.barca, url: "Barca/" }}
    4. 文件夹为 Level 2 节点（radius: 15, color: colors.analysisSub/probSub），连向一级学科。
    5. 文件为 Level 3 节点（radius: 8, color: colors.analysisLeaf），连向所属文件夹。
    6. label 必须严格使用三元表达式：isEn ? "English Name" : "中文名称"。
    7. URL 格式：文件夹名/文件名/ (Intro文件映射为 0Intro/)。

    请直接返回 JSON，不要包含 Markdown 代码块标签。
    """
    
    response = client.chat.completions.create(
        model="deepseek-chat",
        messages=[{"role": "user", "content": prompt}]
    )
    raw_res = response.choices[0].message.content.strip().replace("```json", "").replace("```", "")
    return json.loads(raw_res)

def update_js_file(nodes, links):
    """将数据写回 JS，并处理三元表达式的特殊语法"""
    with open(JS_FILE, 'r', encoding='utf-8') as f:
        content = f.read()

    def dict_to_js_obj(d):
        fields = []
        for k, v in d.items():
            if k == 'label':
                # 核心修正：label 的值直接作为代码写入，不加引号
                fields.append(f'{k}: {v}')
            elif k == 'color' and 'colors.' in str(v):
                # 颜色变量也不加引号
                fields.append(f'{k}: {v}')
            elif isinstance(v, str):
                fields.append(f'{k}: "{v}"')
            else:
                fields.append(f'{k}: {str(v).lower() if isinstance(v, bool) else v}')
        return "        { " + ", ".join(fields) + " },"

    nodes_js = "\n".join([dict_to_js_obj(n) for n in nodes])
    links_js = "\n".join([f'        {{ source: "{l["source"]}", target: "{l["target"]}" }},' for l in links])

    # 正则替换数组内容
    content = re.sub(r'const allNodes = \[.*?\];', f'const allNodes = [\n{nodes_js}\n    ];', content, flags=re.DOTALL)
    content = re.sub(r'const allLinks = \[.*?\];', f'const allLinks = [\n{links_js}\n    ];', content, flags=re.DOTALL)

    with open(JS_FILE, 'w', encoding='utf-8') as f:
        f.write(content)
    print(f"✨ 同步完成！已更新至: {JS_FILE}")

if __name__ == "__main__":
    try:
        print("🔍 扫描目录中...")
        struct = scan_docs()
        print("🤖 请求 AI 脑暴连线中...")
        data = get_ai_suggestions(struct)
        print("💾 注入 JS 文件...")
        update_js_file(data['nodes'], data['links'])
    except Exception as e:
        print(f"❌ 出错了: {e}")