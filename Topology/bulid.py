import os
import re
import subprocess

# 你的文件列表
files = ["1.md", "2.md", "3.md", "4.md", "5.md", "6.md", "7.md"]
temp_file = "temp_combined.md"

with open(temp_file, "w", encoding="utf-8") as outfile:
    for fname in files:
        if os.path.exists(fname):
            with open(fname, "r", encoding="utf-8") as infile:
                content = infile.read()
                
                # 1. 处理 MkDocs 的 Admonition (!!! info "标题")
                # 将其转换为加粗标题
                content = re.sub(r'!!! (?:info|success|warning|error|note|abstract) "(.*?)"', r'### \1', content)
                
                # 2. 处理缩进问题：删除行首的 4 个空格或 1 个制表符
                # 这是最关键的一步，防止 Pandoc 把正文当成代码块
                lines = content.split('\n')
                fixed_lines = []
                for line in lines:
                    if line.startswith('    '):
                        fixed_lines.append(line[4:])
                    elif line.startswith('\t'):
                        fixed_lines.append(line[1:])
                    else:
                        fixed_lines.append(line)
                content = '\n'.join(fixed_lines)
                
                outfile.write(content)
                outfile.write("\n\n\\newpage\n\n") # 自动添加换页符

# 调用 Pandoc 生成 PDF
cmd = [
    "pandoc", temp_file, "-o", "Topology_Review.pdf",
    "--pdf-engine=xelatex",
    "--toc",
    "--number-sections",
    "-V", "mainfont=Noto Serif CJK SC",
    "-V", "CJKmainfont=Noto Serif CJK SC",
    "-V", "geometry:margin=1in"
]

print("正在生成 PDF...")
subprocess.run(cmd)
# os.remove(temp_file) # 可选：生成后删除临时文件
print("生成成功：Topology_Review.pdf")