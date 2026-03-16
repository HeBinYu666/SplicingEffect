import os
import re

# 我们要精准捕获的三个目标转录本 ID
target_tx_ids = ["ENST00000375213", "ENST00000338983", "ENST00000409176"]

# 你的原始大文件路径
gtf_in = "/data5/lab/hebinyu/splicing_impact/ensembl_v110_gtf/Homo_sapiens.GRCh38.110.gtf"
fasta_in = "/data5/lab/hebinyu/splicing_impact/ensembl_v110_fasta/Homo_sapiens.GRCh38.pep.all.fa"
csv_in = "/data5/lab/hebinyu/splicing_impact/pfam/all_mapped_domains.csv"

# 你要在 GitHub 仓库里存放测试数据的路径
out_dir = "/data5/lab/hebinyu/1_Github/splicing/example_data"

# 如果文件夹不存在，自动创建
os.makedirs(out_dir, exist_ok=True)

gtf_out = os.path.join(out_dir, "mini(1).gtf")
fasta_out = os.path.join(out_dir, "mini(1).fasta")
csv_out = os.path.join(out_dir, "mini_domain(1).csv")

# =======================================================
# 核心修复：准备一个空集合，用来动态存放抓取到的 ENSP 蛋白 ID
target_prot_ids = set()
# =======================================================

print("🚀 开始提取迷你测试集...")

# ================= 1. 提取 GTF，并顺手牵羊抓取 ENSP =================
print("正在扫描 GTF 大文件 (同时自动解析 ENSP 蛋白 ID)...")
with open(gtf_in, 'r') as fin, open(gtf_out, 'w') as fout:
    for line in fin:
        if line.startswith('#'):
            fout.write(line)
            continue
        
        # 如果这行内容包含我们的目标 ENST ID
        if any(tx in line for tx in target_tx_ids):
            fout.write(line)
            
            # 【关键修复】立刻用正则把这行的 protein_id (ENSP) 挖出来，加入集合！
            match_prot = re.search(r'protein_id "([^"]+)"', line)
            if match_prot:
                target_prot_ids.add(match_prot.group(1))

print(f"✅ GTF 提取完成 -> {gtf_out}")
print(f"🎯 成功从 GTF 中捕获对应的 Protein IDs: {target_prot_ids}")

# ================= 2. 提取 FASTA (双重保险匹配) =================
print("正在扫描 FASTA 大文件...")
with open(fasta_in, 'r') as fin, open(fasta_out, 'w') as fout:
    keep = False
    for line in fin:
        if line.startswith('>'):
            # 遇到新的序列头，检查是不是我们要的 ENST 或者刚才抓到的 ENSP
            keep = any(tx in line for tx in target_tx_ids) or any(prot in line for prot in target_prot_ids)
        
        if keep:
            fout.write(line)
print(f"✅ FASTA 提取完成 -> {fasta_out}")

# ================= 3. 提取 Domain CSV (使用 ENSP 匹配) =================
print("正在扫描 Domain CSV 大文件...")
with open(csv_in, 'r') as fin, open(csv_out, 'w') as fout:
    header = fin.readline()
    fout.write(header)
    
    for line in fin:
        # 【关键修复】这里必须用刚才抓到的 ENSP 去匹配，绝对不能用 ENST！
        if any(prot in line for prot in target_prot_ids):
            fout.write(line)
print(f"✅ Domain CSV 提取完成 -> {csv_out}")

print("🎉 全部搞定！这下数据逻辑彻底完美了！")