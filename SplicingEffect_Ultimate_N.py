import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
import os
import re
import subprocess
import argparse  # 新增：用来接收命令行参数的模块
from Bio.SeqUtils.ProtParam import ProteinAnalysis

# =====================================================================
# 1. 命令行接口解析 (完美实现 Single Command)
# =====================================================================
parser = argparse.ArgumentParser(description="SplicingEffect: Automated Isoform Functional Analysis Pipeline")
parser.add_argument('-i', '--ids', nargs='+', required=True, help="Input Multi-Transcript IDs separated by space")
parser.add_argument('-g', '--gtf', default="example_data/mini.gtf", help="Path to GTF file")
parser.add_argument('-f', '--fasta', default="example_data/mini.fasta", help="Path to FASTA file")
parser.add_argument('-d', '--domain', default="example_data/mini_domain.csv", help="Path to Domain CSV file")

args = parser.parse_args()

# 把用户输入的参数赋值给全局变量
GTF_FILE = args.gtf
FASTA_FILE = args.fasta
DOMAIN_CSV = args.domain
INPUT_TXS = args.ids

domain_colors = {"PK_Tyr_Ser-Thr": "#e76f51", "SAM_2": "#5bc0be"}

# =====================================================================
# 1.5 智能调度 R 脚本自动化模块 (论文核心承诺)
# =====================================================================
R_SCRIPT_PATH = "/data5/lab/hebinyu/splicing_impact/map_all_domains.R" # 你的真实R脚本路径

if os.path.exists(DOMAIN_CSV):
    print(f"✅ 检测到本地已存在结构域坐标库，直接极速加载...")
else:
    print(f"🔄 未检测到结构域库，正在后台静默启动 R 脚本进行绝对坐标映射，请稍候...")
    try:
        # 这一行就是 Python 自动调 R 的核心！
        subprocess.run(["Rscript", R_SCRIPT_PATH], check=True)
        print("✅ R 脚本执行完毕，成功生成全新的绝对坐标库！")
    except subprocess.CalledProcessError as e:
        print(f"❌ R 脚本执行失败，请检查 R 环境或网络连接: {e}")
        # 如果 R 失败了，生成一个空的 DataFrame 防止后面画图崩溃
        pd.DataFrame(columns=['Protein_ID', 'Domain_Name', 'Genomic_Start', 'Genomic_End']).to_csv(DOMAIN_CSV, index=False)
        
# =====================================================================
# 2. 超高维特征扩充模块 (Biopython + propy3 稳定版)
# =====================================================================
def get_high_dim_propy3_features(seq):
    """调用 propy3 稳定模块，瞬间提取高阶机器学习王牌特征"""
    if not seq or len(seq) < 35: 
        return {} 
    feats = {}
    try:
        from propy.PyPro import GetProDes
        des = GetProDes(seq)
        feats.update(des.GetPAAC())             # 伪氨基酸组成 (约50维)
        feats.update(des.GetAPAAC())            # 双亲性伪氨基酸 (约80维)
        feats.update(des.GetMoreauBrotoAuto())  # MoreauBroto 拓扑自相关
        feats.update(des.GetMoranAuto())        # Moran 拓扑自相关
        feats.update(des.GetGearyAuto())        # Geary 拓扑自相关
    except Exception as e:
        print(f"⚠️ 提示: 提取 PAAC/Auto 特征时遇到个别异常序列，已安全跳过: {e}")
    return feats

def calculate_ctd(seq):
    """保留 147 维基础 CTD 特征"""
    if not seq or len(seq) == 0: return {}
    properties = {
        'Hydrophobicity': [{'R','K','E','D','Q','N'}, {'G','A','S','T','P','H','Y'}, {'C','V','L','I','M','F','W'}],
        'NormVDWVolume': [{'G','A','S','T','P','D','C'}, {'N','V','E','Q','I','L'}, {'M','H','K','F','R','Y','W'}],
        'Polarity': [{'L','I','F','W','C','M','V','Y'}, {'P','A','T','G','S'}, {'H','Q','R','K','N','E','D'}],
        'Polarizability': [{'G','A','S','D','T'}, {'C','P','N','V','E','Q','I','L'}, {'K','M','H','F','R','Y','W'}],
        'Charge': [{'K','R'}, {'A','N','C','Q','G','I','L','M','F','P','S','T','W','Y','V','H'}, {'D','E'}],
        'SecondaryStruct': [{'E','A','L','M','Q','K','R','H'}, {'V','I','Y','C','W','F','T'}, {'G','N','P','S','D'}],
        'SolventAccess': [{'A','L','F','C','G','I','V','W'}, {'R','K','Q','E','N','D'}, {'M','S','P','T','H','Y'}]
    }
    ctd = {}
    seq_len = len(seq)
    for prop_name, groups in properties.items():
        counts = [sum(1 for aa in seq if aa in g) for g in groups]
        for i in range(3): ctd[f"CTD_C_{prop_name}_G{i+1}"] = round(counts[i] / seq_len, 4)
        t_counts = [0, 0, 0]
        for i in range(seq_len - 1):
            g1 = next((j for j, g in enumerate(groups) if seq[i] in g), -1)
            g2 = next((j for j, g in enumerate(groups) if seq[i+1] in g), -1)
            if g1 != -1 and g2 != -1 and g1 != g2:
                if (g1==0 and g2==1) or (g1==1 and g2==0): t_counts[0] += 1
                elif (g1==0 and g2==2) or (g1==2 and g2==0): t_counts[1] += 1
                elif (g1==1 and g2==2) or (g1==2 and g2==1): t_counts[2] += 1
        ctd[f"CTD_T_{prop_name}_12"] = round(t_counts[0] / (seq_len-1), 4) if seq_len>1 else 0
        ctd[f"CTD_T_{prop_name}_13"] = round(t_counts[1] / (seq_len-1), 4) if seq_len>1 else 0
        ctd[f"CTD_T_{prop_name}_23"] = round(t_counts[2] / (seq_len-1), 4) if seq_len>1 else 0
        for i, g in enumerate(groups):
            pos_list = [idx+1 for idx, aa in enumerate(seq) if aa in g]
            for percent in [1, 25, 50, 75, 100]:
                feat_name = f"CTD_D_{prop_name}_G{i+1}_{percent}p"
                if len(pos_list) == 0: ctd[feat_name] = 0
                else:
                    target_idx = max(0, int((percent/100) * len(pos_list)) - 1)
                    if percent == 1: target_idx = 0
                    if percent == 100: target_idx = len(pos_list) - 1
                    ctd[feat_name] = round(pos_list[target_idx] / seq_len, 4)
    return ctd

def get_all_features(seq):
    """整合基础特征与超高维稳定特征"""
    if not seq or len(seq) == 0: return {"MW": 0}
    seq = seq.replace("X", "").replace("U", "").replace("Z", "").replace("B", "")
    pa = ProteinAnalysis(seq)
    feats = {"MW": round(pa.molecular_weight(), 2), 
             "pI": round(pa.isoelectric_point(), 4), 
             "GRAVY": round(pa.gravy(), 4), 
             "Instability_Index": round(pa.instability_index(), 4)}
    extinction = pa.molar_extinction_coefficient()
    feats["Extinction_Coefficient_Reduced"] = extinction[0]
    feats["Extinction_Coefficient_Oxidized"] = extinction[1]
    sec = pa.secondary_structure_fraction()
    feats["SecStruct_Helix"], feats["SecStruct_Turn"], feats["SecStruct_Sheet"] = round(sec[0],4), round(sec[1],4), round(sec[2],4)
    aac = pa.amino_acids_percent
    for aa, pct in aac.items(): feats[f"AAC_{aa}"] = round(pct, 4)
    mole_pct_A, mole_pct_V, mole_pct_I, mole_pct_L = aac.get('A', 0)*100, aac.get('V', 0)*100, aac.get('I', 0)*100, aac.get('L', 0)*100
    feats["Aliphatic_Index"] = round(mole_pct_A + 2.9 * mole_pct_V + 3.9 * (mole_pct_I + mole_pct_L), 4)
    
    feats.update(calculate_ctd(seq))
    feats.update(get_high_dim_propy3_features(seq))
    return feats

# =====================================================================
# 3. 解析 GTF 获取基础结构
# =====================================================================
tx_info = {tx: {'exons': [], 'cds': [], 'strand': '+', 'is_canonical': False, 
                'protein_id': '', 'gene_id': '', 'biotype': 'unknown',
                'span_start': float('inf'), 'span_end': 0} for tx in INPUT_TXS}

if os.path.exists(GTF_FILE):
    with open(GTF_FILE, 'r') as f:
        for line in f:
            if line.startswith('#'): continue
            parts = line.split('\t')
            if len(parts) < 9: continue
            info_str = parts[8]
            for tx in INPUT_TXS:
                if tx in info_str:
                    tx_info[tx]['strand'] = parts[6]
                    if 'Ensembl_canonical' in info_str or 'MANE_Select' in info_str: tx_info[tx]['is_canonical'] = True
                    match_prot = re.search(r'protein_id "([^"]+)"', info_str)
                    if match_prot: tx_info[tx]['protein_id'] = match_prot.group(1)
                    match_gene = re.search(r'gene_id "([^"]+)"', info_str)
                    if match_gene: tx_info[tx]['gene_id'] = match_gene.group(1)
                    match_biotype = re.search(r'transcript_biotype "([^"]+)"', info_str)
                    if not match_biotype: match_biotype = re.search(r'transcript_type "([^"]+)"', info_str)
                    if match_biotype: tx_info[tx]['biotype'] = match_biotype.group(1)
                    
                    start, end = int(parts[3]), int(parts[4])
                    if parts[2] == 'exon':
                        tx_info[tx]['span_start'] = min(tx_info[tx]['span_start'], start)
                        tx_info[tx]['span_end'] = max(tx_info[tx]['span_end'], end)
                        tx_info[tx]['exons'].append({'start': start, 'end': end})
                    elif parts[2] == 'CDS':
                        tx_info[tx]['cds'].append({'start': start, 'end': end})

can_tx = next((tx for tx in INPUT_TXS if tx_info[tx]['is_canonical']), None)
if not can_tx:
    can_tx = max(INPUT_TXS, key=lambda x: tx_info[x]['span_end'] - tx_info[x]['span_start'])
alt_txs = [tx for tx in INPUT_TXS if tx != can_tx]

# =====================================================================
# 4. 提取 FASTA 并构建矩阵
# =====================================================================
target_prots = [tx_info[tx]['protein_id'] for tx in INPUT_TXS if tx_info[tx]['protein_id']]
sequences = {pid: "" for pid in target_prots}
current_id = None
if os.path.exists(FASTA_FILE):
    with open(FASTA_FILE, 'r') as f:
        for line in f:
            if line.startswith('>'):
                current_id = next((pid for pid in target_prots if pid in line), None)
            elif current_id: sequences[current_id] += line.strip()

can_prot_id = tx_info[can_tx]['protein_id']
feat_c = get_all_features(sequences.get(can_prot_id, ""))

df_domains = pd.read_csv(DOMAIN_CSV) if os.path.exists(DOMAIN_CSV) else pd.DataFrame(columns=['Protein_ID', 'Domain_Name'])
d_c_list = df_domains[df_domains['Protein_ID'].str.contains(can_prot_id, na=False)]['Domain_Name'].tolist() if can_prot_id else []
sp_c = tx_info[can_tx]['span_end'] - tx_info[can_tx]['span_start'] + 1

rows = []
for alt_tx in alt_txs:
    alt_prot_id = tx_info[alt_tx]['protein_id']
    feat_a = get_all_features(sequences.get(alt_prot_id, ""))
    
    d_a_list = df_domains[df_domains['Protein_ID'].str.contains(alt_prot_id, na=False)]['Domain_Name'].tolist() if alt_prot_id else []
    d_shared, d_lost, d_gained = list(set(d_c_list) & set(d_a_list)), list(set(d_c_list) - set(d_a_list)), list(set(d_a_list) - set(d_c_list))
    sp_a = tx_info[alt_tx]['span_end'] - tx_info[alt_tx]['span_start'] + 1
    
    row = {
        "Gene_ID": tx_info[can_tx]['gene_id'],
        "Canonical_Transcript_ID": can_tx,
        "Alternative_Transcript_ID": alt_tx,
        "Canonical_Protein_ID": can_prot_id,
        "Alternative_Protein_ID": alt_prot_id,
        "Coding_Potential_Canonical": tx_info[can_tx]['biotype'],
        "Coding_Potential_Alternative": tx_info[alt_tx]['biotype'],
        "Protein_Length_Canonical": len(sequences.get(can_prot_id, "")),
        "Protein_Length_Alternative": len(sequences.get(alt_prot_id, "")),
        "Protein_Length_Difference": len(sequences.get(alt_prot_id, "")) - len(sequences.get(can_prot_id, "")),
        "Molecular_Weight_Canonical_Da": feat_c.get("MW", 0),
        "Molecular_Weight_Alternative_Da": feat_a.get("MW", 0),
        "Molecular_Weight_Difference_Da": round(feat_a.get("MW", 0) - feat_c.get("MW", 0), 2),
        "Exon_Count_Canonical": len(tx_info[can_tx]['exons']),
        "Exon_Count_Alternative": len(tx_info[alt_tx]['exons']),
        "Exon_Count_Difference": len(tx_info[alt_tx]['exons']) - len(tx_info[can_tx]['exons']),
        "Genomic_Span_Canonical_bp": sp_c if sp_c != float('inf') else 0,
        "Genomic_Span_Alternative_bp": sp_a if sp_a != float('inf') else 0,
        "Genomic_Span_Difference_bp": (sp_a - sp_c) if sp_c != float('inf') else 0,
        "Domains_List_Canonical": ";".join(d_c_list) if d_c_list else "None",
        "Domains_List_Alternative": ";".join(d_a_list) if d_a_list else "None",
        "Domains_List_Shared": ";".join(d_shared) if d_shared else "None",
        "Domains_List_Lost": ";".join(d_lost) if d_lost else "None",
        "Domains_List_Gained": ";".join(d_gained) if d_gained else "None",
        "Domain_Count_Lost": len(d_lost),
        "Domain_Count_Gained": len(d_gained),
        "Domain_Count_Total_Change": len(d_lost) + len(d_gained)
    }
    
    feat_c_copy = feat_c.copy()
    feat_c_copy.pop("MW", None); feat_a.pop("MW", None)
    for k in feat_c_copy.keys():
        row[f"Can_{k}"] = feat_c_copy[k]
        row[f"Alt_{k}"] = feat_a.get(k, 0)
        try:
            row[f"Delta_{k}"] = round(feat_a.get(k, 0) - feat_c_copy[k], 4)
        except:
            row[f"Delta_{k}"] = 0
            
    rows.append(row)

pd.DataFrame(rows).to_csv("SplicingEffect_Ultimate_Features(N).csv", index=False)
print(f"成功生成纯净满血版 CSV 表格！当前特征总列数: {len(rows[0])} 列。")

# =====================================================================
# 5. 动态多轨可视化模块
# =====================================================================
y_spacing = 1.5
tx_pairs = [{"tx": can_tx, "prot": can_prot_id, "label": f"Canonical:\n{can_tx}\n(Strand: {tx_info[can_tx]['strand']})", "y_pos": len(alt_txs) * y_spacing + 1}]
for i, alt_tx in enumerate(alt_txs):
    tx_pairs.append({"tx": alt_tx, "prot": tx_info[alt_tx]['protein_id'], "label": f"Alternative:\n{alt_tx}\n(Strand: {tx_info[alt_tx]['strand']})", "y_pos": (len(alt_txs) - i - 1) * y_spacing + 1})

exons_dict = {tx: sorted(tx_info[tx]['exons'], key=lambda x: x['start']) for tx in INPUT_TXS}
cds_dict = {tx: sorted(tx_info[tx]['cds'], key=lambda x: x['start']) for tx in INPUT_TXS}

boundaries = set()
for ex_list in exons_dict.values():
    for e in ex_list: boundaries.update([e['start'], e['end']])
for c_list in cds_dict.values():
    for c in c_list: boundaries.update([c['start'], c['end']])
boundaries = sorted(list(boundaries))

segments = []
if boundaries:
    for i in range(len(boundaries)-1):
        s, e = boundaries[i], boundaries[i+1]
        is_cds = any(c['start'] <= s and e <= c['end'] for c_list in cds_dict.values() for c in c_list)
        is_exon = any(ex['start'] <= s and e <= ex['end'] for ex_list in exons_dict.values() for ex in ex_list)
        segments.append({'start': s, 'end': e, 'is_cds': is_cds, 'is_exon': is_exon})

MAX_INTRON_LEN = 200
MAX_UTR_LEN = 150
mapping = {boundaries[0]: 0} if boundaries else {}
current_vis_pos = 0

for seg in segments:
    length = seg['end'] - seg['start']
    if not seg['is_exon']: vis_len = min(length, MAX_INTRON_LEN)
    elif not seg['is_cds']: vis_len = min(length, MAX_UTR_LEN)
    else: vis_len = length         
    current_vis_pos += vis_len
    mapping[seg['end']] = current_vis_pos

def map_pos(pos):
    if not boundaries: return 0
    if pos <= boundaries[0]: return mapping[boundaries[0]]
    if pos >= boundaries[-1]: return mapping[boundaries[-1]]
    for i in range(len(boundaries)-1):
        if boundaries[i] <= pos <= boundaries[i+1]:
            s, e = boundaries[i], boundaries[i+1]
            if s == e: return mapping[s]
            ratio = (pos - s) / (e - s)
            return mapping[s] + ratio * (mapping[e] - mapping[s])
    return 0

fig, ax = plt.subplots(figsize=(16, max(5, len(INPUT_TXS) * 1.5))) 
H_utr, H_cds, H_domain = 0.15, 0.35, 0.45 

for tx in tx_pairs:
    y = tx["y_pos"]
    exons = exons_dict[tx["tx"]]
    cds_list = cds_dict[tx["tx"]]
    if not exons: continue
    
    vis_gene_start, vis_gene_end = map_pos(exons[0]['start']), map_pos(exons[-1]['end'])
    ax.plot([vis_gene_start, vis_gene_end], [y, y], color='dimgray', linewidth=1.2, zorder=1)
    ax.text(vis_gene_start - 300, y, tx["label"], ha='right', va='center', fontweight='bold', color='#333333', fontsize=10)

    for i in range(len(exons) - 1):
        intron_start = map_pos(exons[i]['end'])
        intron_end = map_pos(exons[i+1]['start'])
        intron_len = intron_end - intron_start
        if intron_len >= 40: 
            if intron_len < 120: ax.plot(intron_start + intron_len/2, y, marker='>', markersize=4, color='#888888', linestyle='None', zorder=2) 
            else:
                ax.plot(intron_start + intron_len/3, y, marker='>', markersize=4, color='#888888', linestyle='None', zorder=2) 
                ax.plot(intron_start + 2*intron_len/3, y, marker='>', markersize=4, color='#888888', linestyle='None', zorder=2) 

    for exon in exons:
        v_start, v_end = map_pos(exon['start']), map_pos(exon['end'])
        ax.add_patch(patches.Rectangle((v_start, y - H_utr/2), v_end - v_start, H_utr, linewidth=0.8, edgecolor='gray', facecolor='#e0e0e0', zorder=3)) 

    for cds in cds_list:
        v_start, v_end = map_pos(cds['start']), map_pos(cds['end'])
        ax.add_patch(patches.Rectangle((v_start, y - H_cds/2), v_end - v_start, H_cds, linewidth=1, edgecolor='dimgray', facecolor='#a0a0a0', zorder=4)) 
        
    tx_domains = df_domains[df_domains['Protein_ID'].str.contains(tx["prot"]) if tx["prot"] else False]
    for _, row_domain in tx_domains.iterrows():
        d_start, d_end = row_domain['Genomic_Start'], row_domain['Genomic_End'] 
        color = domain_colors.get(row_domain['Domain_Name'], "orange") 
        for cds in cds_list:
            overlap_start = max(cds['start'], d_start) 
            overlap_end = min(cds['end'], d_end) 
            if overlap_start < overlap_end:
                v_start, v_end = map_pos(overlap_start), map_pos(overlap_end) 
                ax.add_patch(patches.Rectangle((v_start, y - H_domain/2), v_end - v_start, H_domain, linewidth=0.8, edgecolor='black', facecolor=color, alpha=1.0, zorder=5)) 

legend_elements = [patches.Patch(facecolor='#e0e0e0', edgecolor='gray', label='UTR'), patches.Patch(facecolor='#a0a0a0', edgecolor='dimgray', label='CDS')]
for d_name, color in domain_colors.items(): legend_elements.append(patches.Patch(facecolor=color, edgecolor='dimgray', label=d_name)) 

ax.legend(handles=legend_elements, loc='upper right', bbox_to_anchor=(1.15, 1.25), frameon=False) 
ax.autoscale_view() 
ax.set_ylim(0, len(INPUT_TXS) * y_spacing + 1)
ax.axis('off')
plt.title("Figure 1. Isoform Structure & Domain Mapping(N)", fontsize=16, fontweight='bold', loc='left', pad=20)
plt.tight_layout()
plt.savefig("SplicingEffect_Ultimate_Output(N).png", dpi=300, bbox_inches='tight')
print("多序列全自动流水线执行完毕！图像与特征矩阵已完美生成！")