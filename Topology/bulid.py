import os
import re
import subprocess

# 1. 更新你的文件列表
# 按照你希望在 PDF 中出现的顺序排列。如果你有 8.md，记得也加上。
files = [
    "1.md", "2.md", "3.md", "4.md", "5.md", "6.md", "7.md",
    "sol_12.md", "sol_34.md", "sol_56.md"
]
temp_file = "temp_combined.md"

with open(temp_file, "w", encoding="utf-8") as outfile:
    for fname in files:
        if os.path.exists(fname):
            with open(fname, "r", encoding="utf-8") as infile:
                content = infile.read()
                
                # 2. 关键修改：同时兼容 !!! 和 ??? 两种 MkDocs 语法
                # 将它们统一转换为标准的 ### 三级标题
                content = re.sub(r'(?:!!!|\?\?\?) (?:info|success|warning|error|note|abstract) "(.*?)"', r'### \1', content)
                
                # 3. 处理缩进问题：删除行首的 4 个空格或 1 个制表符
                # 这保证了解答框里的内容不会被 Pandoc 误认为成代码块
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
                outfile.write("\n\n\\newpage\n\n") # 自动添加换页符，保证每一章/每个作业从新的一页开始

# 调用 Pandoc 生成 PDF
# 这里我把输出文件名改成了 Topology_Full_Review.pdf，以免覆盖你之前的纯知识点版本
cmd = [
    "pandoc", temp_file, "-o", "Topology_Full_Review.pdf",
    "--pdf-engine=xelatex",
    "--toc",
    "--number-sections",
    "-V", "mainfont=Noto Serif CJK SC",
    "-V", "CJKmainfont=Noto Serif CJK SC",
    "-V", "geometry:margin=1in"
]

print("正在整合文档并生成 PDF...")
subprocess.run(cmd)

# 如果你不想保留中间生成的 temp_combined.md，可以取消下面这行的注释
# os.remove(temp_file) 

print("生成成功：Topology_Full_Review.pdf")