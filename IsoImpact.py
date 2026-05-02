import os
import re
import argparse
import gzip

# =====================================================================
# 1. Command-line interface
# =====================================================================
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
DEFAULT_HUMAN_DOMAIN_CSV = os.path.join(BASE_DIR, "data", "human_domain.csv")
DEFAULT_MOUSE_DOMAIN_CSV = os.path.join(BASE_DIR, "data", "mouse_domain.csv")
REQUIRED_DOMAIN_COLUMNS = [
    "Protein_ID",
    "Domain_ID",
    "Domain_Name",
    "Genomic_Start",
    "Genomic_End",
]


def open_text(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")

parser = argparse.ArgumentParser(
    description="IsoImpact: automated isoform functional impact analysis"
)

group = parser.add_mutually_exclusive_group(required=True)
group.add_argument(
    "-i", "--isoform", nargs="+",
    help="Compare at least two transcript IDs, for example: -i ENST1 ENST2"
)
group.add_argument(
    "-q", "--query-exon", nargs=3, metavar=("TRANSCRIPT_ID", "START", "END"),
    help="Query domains overlapping one exon, for example: -q ENST1 1000 2000"
)

parser.add_argument("-g", "--gtf", help="Path to the matching GTF annotation file")
parser.add_argument("-f", "--fasta", help="Path to the matching protein FASTA file")
parser.add_argument(
    "-d", "--domain", default=None,
    help="Optional custom domain-coordinate CSV. Standard human/mouse runs use built-in files."
)
parser.add_argument(
    "-o", "--outdir", default=".",
    help="Directory for output files in isoform comparison mode"
)
parser.add_argument(
    "--prefix", default="IsoImpact",
    help="Prefix for output files in isoform comparison mode"
)

args = parser.parse_args()

GTF_FILE = args.gtf
FASTA_FILE = args.fasta
OUT_DIR = args.outdir
OUT_PREFIX = args.prefix

try:
    import pandas as pd
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.patches as patches
    from matplotlib.gridspec import GridSpec
    import numpy as np
    from Bio.SeqUtils.ProtParam import ProteinAnalysis
except ImportError as exc:
    parser.error(
        "Missing Python dependency. Install the required packages with "
        "`pip install -r requirements.txt`. Original error: " + str(exc)
    )

# =====================================================================
# 2. Built-in human/mouse domain routing, with optional custom override
# =====================================================================
target_id_for_routing = args.query_exon[0] if args.query_exon else args.isoform[0]

if args.domain:
    DOMAIN_CSV = args.domain
    print(f"[IsoImpact] Using user-provided domain-coordinate file: {DOMAIN_CSV}")
elif target_id_for_routing.startswith("ENSMUS"):
    DOMAIN_CSV = DEFAULT_MOUSE_DOMAIN_CSV
    print(f"[IsoImpact] Mouse Ensembl ID detected; using built-in mouse domain database: {DOMAIN_CSV}")
elif target_id_for_routing.startswith("ENS"):
    DOMAIN_CSV = DEFAULT_HUMAN_DOMAIN_CSV
    print(f"[IsoImpact] Human Ensembl ID detected; using built-in human domain database: {DOMAIN_CSV}")
else:
    parser.error(
        "Could not infer a built-in human/mouse domain database from the transcript ID. "
        "For novel or custom isoforms, provide a compatible domain CSV with -d/--domain."
    )

if not GTF_FILE:
    parser.error("A matching GTF annotation file is required with -g/--gtf.")

if args.isoform and not FASTA_FILE:
    parser.error("Isoform comparison mode requires a matching protein FASTA file with -f/--fasta.")

if not os.path.exists(GTF_FILE):
    parser.error(f"GTF file not found: {GTF_FILE}")

if args.isoform and not os.path.exists(FASTA_FILE):
    parser.error(f"FASTA file not found: {FASTA_FILE}")

if not os.path.exists(DOMAIN_CSV):
    parser.error(f"Domain-coordinate CSV not found: {DOMAIN_CSV}")

df_domains = pd.read_csv(DOMAIN_CSV)
missing_domain_cols = [c for c in REQUIRED_DOMAIN_COLUMNS if c not in df_domains.columns]
if missing_domain_cols:
    parser.error(
        "Domain-coordinate CSV is missing required column(s): "
        + ", ".join(missing_domain_cols)
    )

df_domains["Protein_ID"] = df_domains["Protein_ID"].astype(str)

# =====================================================================
# 3. Exon-domain query mode (-q)
# =====================================================================
if args.query_exon:
    target_enst = args.query_exon[0]
    try:
        exon_start = int(args.query_exon[1])
        exon_end = int(args.query_exon[2])
    except ValueError:
        print("Error: exon start and end coordinates must be integers.")
        exit(1)
        
    print(f"[IsoImpact] Querying domains overlapping {target_enst}:{exon_start}-{exon_end}...")
    
    prot_id = None
    if os.path.exists(GTF_FILE):
        with open_text(GTF_FILE) as f:
            for line in f:
                if line.startswith('#'): continue
                if target_enst in line:
                    match_prot = re.search(r'protein_id "([^"]+)"', line)
                    if match_prot:
                        prot_id = match_prot.group(1)
                        break
    
    if not prot_id:
        print(f"No protein_id was found for {target_enst} in the GTF file.")
        exit(0)

    target_doms = df_domains[df_domains['Protein_ID'].str.contains(prot_id, na=False, regex=False)]
    
    if target_doms.empty:
        print(f"No annotated domain was found for protein {prot_id}.")
        exit(0)
        
    overlapping_domains = []
    for _, row in target_doms.iterrows():
        overlap_start = max(exon_start, row['Genomic_Start'])
        overlap_end = min(exon_end, row['Genomic_End'])
        
        if overlap_start <= overlap_end:
            overlap_length = overlap_end - overlap_start + 1
            dom_name = row.get('Domain_Name', 'Unknown')
            dom_id = row.get('Domain_ID', '')
            if pd.notna(dom_id) and str(dom_id).strip():
                display_name = f"{dom_name} ({dom_id})"
            else:
                display_name = dom_name
                
            overlapping_domains.append(f"{display_name} [overlap: {overlap_length} bp]")
            
    if overlapping_domains:
        print("Overlapping domain(s): " + " | ".join(overlapping_domains))
    else:
        print("No domain overlaps this exon interval.")
        
    # Query mode ends here; FASTA is not required.
    exit(0)

# =====================================================================
# 4. Isoform comparison and visualization mode (-i)
# =====================================================================
if args.isoform:
    INPUT_TXS = args.isoform
    if len(INPUT_TXS) < 2:
        print("Error: isoform comparison mode requires at least two transcript IDs.")
        exit(1)

os.makedirs(OUT_DIR, exist_ok=True)

# =====================================================================
# 4. Protein sequence feature extraction
# =====================================================================
def get_high_dim_propy3_features(seq):
    if not seq or len(seq) < 35: return {}
    feats = {}
    try:
        from propy.PyPro import GetProDes
        des = GetProDes(seq)
        try: feats.update(des.GetPAAC())
        except: pass
        try: feats.update(des.GetAPAAC())
        except: pass
        try: feats.update(des.GetMoreauBrotoAuto())
        except: pass
        try: feats.update(des.GetMoranAuto())
        except: pass
        try: feats.update(des.GetGearyAuto())
        except: pass
    except Exception:
        pass
    return feats

def calculate_ctd(seq):
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
    if not seq or len(seq) == 0: return {"MW": 0}
    seq = seq.replace("X", "").replace("U", "").replace("Z", "").replace("B", "")
    if not seq:
        return {"MW": 0}
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
# 5. Parse the GTF annotation
# =====================================================================
tx_info = {tx: {'exons': [], 'cds': [], 'strand': '+', 'is_canonical': False, 
                'protein_id': '', 'gene_id': '', 'biotype': 'unknown',
                'span_start': float('inf'), 'span_end': 0} for tx in INPUT_TXS}

if os.path.exists(GTF_FILE):
    with open_text(GTF_FILE) as f:
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

missing_txs = [tx for tx in INPUT_TXS if tx_info[tx]['span_start'] == float('inf')]
if missing_txs:
    parser.error(
        "The following transcript ID(s) were not found in the GTF exon records: "
        + ", ".join(missing_txs)
    )

missing_proteins = [tx for tx in INPUT_TXS if not tx_info[tx]['protein_id']]
if missing_proteins:
    parser.error(
        "The following transcript ID(s) do not have protein_id annotations in the GTF: "
        + ", ".join(missing_proteins)
    )

can_tx = next((tx for tx in INPUT_TXS if tx_info[tx]['is_canonical']), None)
if not can_tx:
    can_tx = max(INPUT_TXS, key=lambda x: tx_info[x]['span_end'] - tx_info[x]['span_start'])
alt_txs = [tx for tx in INPUT_TXS if tx != can_tx]

# =====================================================================
# 6. Extract protein FASTA sequences and build the feature matrix
# =====================================================================
target_prots = [tx_info[tx]['protein_id'] for tx in INPUT_TXS if tx_info[tx]['protein_id']]
sequences = {pid: "" for pid in target_prots}
current_id = None
if os.path.exists(FASTA_FILE):
    with open_text(FASTA_FILE) as f:
        for line in f:
            if line.startswith('>'):
                current_id = next((pid for pid in target_prots if pid in line), None)
            elif current_id: sequences[current_id] += line.strip()

missing_sequences = [pid for pid, seq in sequences.items() if not seq]
if missing_sequences:
    parser.error(
        "The following protein sequence(s) were not found in the FASTA file: "
        + ", ".join(missing_sequences)
    )

can_prot_id = tx_info[can_tx]['protein_id']
feat_c = get_all_features(sequences.get(can_prot_id, ""))

d_c_list = df_domains[df_domains['Protein_ID'].str.contains(can_prot_id, na=False, regex=False)]['Domain_Name'].tolist() if can_prot_id and not df_domains.empty else []
sp_c = tx_info[can_tx]['span_end'] - tx_info[can_tx]['span_start'] + 1

rows = []
for alt_tx in alt_txs:
    alt_prot_id = tx_info[alt_tx]['protein_id']
    feat_a = get_all_features(sequences.get(alt_prot_id, ""))
    
    d_a_list = df_domains[df_domains['Protein_ID'].str.contains(alt_prot_id, na=False, regex=False)]['Domain_Name'].tolist() if alt_prot_id and not df_domains.empty else []
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

df_all_features = pd.DataFrame(rows)
if df_all_features.empty:
    print("Error: no isoform comparison rows were generated. Please check transcript IDs, GTF, FASTA, and domain files.")
    exit(1)

delta_cols = [c for c in df_all_features.columns if c.startswith('Delta_')]
if not delta_cols:
    parser.error(
        "No protein feature delta columns were generated. Please check the protein FASTA sequences."
    )

csv_filename = os.path.join(OUT_DIR, f"{OUT_PREFIX}_features.csv")
df_all_features.to_csv(csv_filename, index=False)
print(f"[IsoImpact] Feature matrix written to {csv_filename} ({len(rows[0])} columns).")

# =====================================================================
# 7. Dynamic visualization layout
# =====================================================================
print("[IsoImpact] Rendering domain map and ranked feature-difference plot...")

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

# Prepare plot data and domain colors.
all_txs_to_plot = [can_tx] + alt_txs
current_used_domains = []
for tx in all_txs_to_plot:
    prot_id = tx_info[tx]['protein_id']
    if prot_id and not df_domains.empty:
        tx_doms = df_domains[df_domains['Protein_ID'].str.contains(prot_id, na=False, regex=False)]
        current_used_domains.extend(tx_doms['Domain_Name'].dropna().tolist())

unique_domains = sorted(list(set(current_used_domains)))
DOMAIN_PALETTE = [
    '#F4A460',
    '#5DCAA5',
    '#9B8DC4',
    '#E8735A',
    '#87CEEB',
    '#FFB347',
    '#77DD77',
    '#DDA0DD',
    '#48D1CC',
    '#F08080',
    '#AEC6CF',
    '#C8A96A',
    '#6FB8A3',
    '#E0A080',
    '#A8C8A0',
    '#C4A0C4',
    '#80C4C4',
    '#D4B080',
    '#90B8D0',
    '#C8A090',
]
dynamic_domain_colors = {dom: DOMAIN_PALETTE[i % len(DOMAIN_PALETTE)] 
                         for i, dom in enumerate(unique_domains)}

# Use an A4-style figure height that scales with the number of isoforms.
fig_height = max(12, 4 + len(INPUT_TXS) * 2.0)
fig = plt.figure(figsize=(20, fig_height), dpi=300, constrained_layout=True)
gs = GridSpec(2, 1, figure=fig, height_ratios=[len(INPUT_TXS)*0.8, 3.0])

# ==================== [A] Genomic structure and domain map ====================
ax_a = fig.add_subplot(gs[0, 0])
ax_a.set_title('A', loc='left', fontsize=22, fontweight='bold', pad=10)

H_utr, H_cds, H_domain = 0.20, 0.40, 0.55

# Use one shared left label anchor to keep transcript labels aligned.
global_min_x = float('inf')
for tx in all_txs_to_plot:
    if exons_dict[tx]: global_min_x = min(global_min_x, map_pos(exons_dict[tx][0]['start']))
text_align_x = global_min_x - 400 if global_min_x != float('inf') else -400

y_positions = np.arange(len(all_txs_to_plot), 0, -1) * 2.0 
track_labels = [f"Canonical:\n{can_tx}"] + [f"Alternative:\n{tx}" for tx in alt_txs]

for tx, y, label in zip(all_txs_to_plot, y_positions, track_labels):
    exons = exons_dict[tx]
    cds_list = cds_dict[tx]
    if not exons: continue
    
    vis_gene_start, vis_gene_end = map_pos(exons[0]['start']), map_pos(exons[-1]['end'])
    ax_a.plot([vis_gene_start, vis_gene_end], [y, y], color='dimgray', linewidth=1.5, zorder=1)
    
    ax_a.text(text_align_x, y, label, ha='right', va='center', fontweight='bold', color='#222222', fontsize=20)

    for i in range(len(exons) - 1):
        intron_start = map_pos(exons[i]['end'])
        intron_end = map_pos(exons[i+1]['start'])
        intron_len = intron_end - intron_start
        if intron_len >= 40: 
            if intron_len < 120: ax_a.plot(intron_start + intron_len/2, y, marker='>', markersize=6, color='#888888', linestyle='None', zorder=2) 
            else:
                ax_a.plot(intron_start + intron_len/3, y, marker='>', markersize=6, color='#888888', linestyle='None', zorder=2) 
                ax_a.plot(intron_start + 2*intron_len/3, y, marker='>', markersize=6, color='#888888', linestyle='None', zorder=2) 

    for exon in exons:
        v_start, v_end = map_pos(exon['start']), map_pos(exon['end'])
        ax_a.add_patch(patches.Rectangle((v_start, y - H_utr/2), v_end - v_start, H_utr, linewidth=1.0, edgecolor='gray', facecolor='#e0e0e0', zorder=3)) 

    for cds in cds_list:
        v_start, v_end = map_pos(cds['start']), map_pos(cds['end'])
        ax_a.add_patch(patches.Rectangle((v_start, y - H_cds/2), v_end - v_start, H_cds, linewidth=1.2, edgecolor='#2C5F8A' , facecolor='#5B8DB8', zorder=4)) 
        
    if not df_domains.empty and tx_info[tx]['protein_id']:
        tx_domains = df_domains[df_domains['Protein_ID'].str.contains(tx_info[tx]['protein_id'], na=False, regex=False)]
        for _, row_domain in tx_domains.iterrows():
            d_start, d_end = row_domain['Genomic_Start'], row_domain['Genomic_End'] 
            color = dynamic_domain_colors.get(row_domain['Domain_Name'], "orange") 
            for cds in cds_list:
                overlap_start = max(cds['start'], d_start) 
                overlap_end = min(cds['end'], d_end) 
                if overlap_start < overlap_end:
                    v_start, v_end = map_pos(overlap_start), map_pos(overlap_end) 
                    ax_a.add_patch(patches.Rectangle((v_start, y - H_domain/2), v_end - v_start, H_domain, linewidth=1.0, edgecolor='black', facecolor=color, alpha=1.0, zorder=5)) 

ax_a.set_ylim(0.5, len(all_txs_to_plot) * 2.0 + 1.0) 
ax_a.axis('off')

legend_elements = [patches.Patch(facecolor='#e0e0e0', edgecolor='#999999', label='UTR'), 
                   patches.Patch(facecolor='#5B8DB8', edgecolor='#2C5F8A', label='CDS')]
for d_name, color in dynamic_domain_colors.items(): 
    legend_elements.append(patches.Patch(facecolor=color, edgecolor='dimgray', label=d_name)) 
ax_a.legend(handles=legend_elements, loc='center left', bbox_to_anchor=(1.01, 0.5), frameon=False, 
            prop={'size': 20, 'weight': 'bold'}, title="Regions & Domains", title_fontproperties={'size': 22, 'weight': 'bold'})

# ==================== [B] Ranked signed feature differences ====================
ax_new = fig.add_subplot(gs[1, :])

if not df_all_features.empty:
    mean_delta = df_all_features[delta_cols].mean(axis=0)
    sorted_delta_cols = mean_delta.sort_values(ascending=False).index.tolist()
    sorted_feat_names = [c.replace('Delta_', '') for c in sorted_delta_cols]
    n_feats = len(sorted_feat_names)

    isoform_styles = [
        {'color': '#0072B2'},
        {'color': '#E69F00'},
        {'color': '#009E73'},
        {'color': '#CC79A7'},
        {'color': '#56B4E9'},
        {'color': '#F0E442'},
        {'color': '#D55E00'},
    ]
    all_y_vals = []
    isoform_plot_data = []
    for i, (_, row) in enumerate(df_all_features.iterrows()):
        alt_id = row['Alternative_Transcript_ID']
        y_vals = [float(row[f'Delta_{feat}']) if f'Delta_{feat}' in row.index else 0.0
                  for feat in sorted_feat_names]
        isoform_plot_data.append((alt_id, y_vals, isoform_styles[i % len(isoform_styles)]))
        all_y_vals.extend(y_vals)

    y_arr = np.array(all_y_vals)
    y_lo = np.percentile(y_arr, 1)
    y_hi = np.percentile(y_arr, 99)
    y_margin = (y_hi - y_lo) * 0.15
    if y_margin == 0:
        y_margin = max(abs(y_hi) * 0.15, 1.0)
    ax_new.set_ylim(y_lo - y_margin, y_hi + y_margin)

    n_iso = len(isoform_plot_data)
    offsets = np.linspace(-0.3, 0.3, n_iso) if n_iso > 1 else [0]

    for idx, (alt_id, y_vals, style) in enumerate(isoform_plot_data):
        # Each alternative isoform is ranked independently, as described in the paper.
        y_sorted = sorted(y_vals, reverse=True)
        x_indices = list(range(len(y_sorted)))
        
        ax_new.vlines(
            x_indices, 0, y_sorted,
            color=style['color'],
            alpha=0.55,
            linewidth=0.6,
            label=alt_id,
            zorder=3
        )

    ax_new.axhline(0, color='black', linewidth=1.5, zorder=2)
    ax_new.set_xlabel(
        f'Feature Rank',
        fontsize=20, fontweight='bold'
    )
    ax_new.set_ylabel('Difference  (Alt - Canonical)', fontsize=20, fontweight='bold')
    ax_new.tick_params(axis='both', labelsize=18)
    plt.setp(ax_new.get_xticklabels(), fontweight='bold')
    plt.setp(ax_new.get_yticklabels(), fontweight='bold')
    ax_new.spines['top'].set_visible(False)
    ax_new.spines['right'].set_visible(False)
    # Label one positive and one negative extreme for each alternative isoform.
    y_vis_hi = ax_new.get_ylim()[1]
    y_vis_lo = ax_new.get_ylim()[0]

    def truncate_name(name, max_len=22):
        return name if len(name) <= max_len else name[:max_len] + '...'

    for iso_idx, (alt_id, y_vals, style) in enumerate(isoform_plot_data):
        paired_iso = sorted(zip(sorted_delta_cols, y_vals), key=lambda x: x[1], reverse=True)
        feat_names_iso = [c.replace('Delta_', '') for c, _ in paired_iso]
        y_iso = [v for _, v in paired_iso]

        for rank in range(len(y_iso)):
            val = y_iso[rank]
            if val <= 0:
                break
            ax_new.annotate(
                feat_names_iso[rank],
                xy=(rank, min(val, y_vis_hi * 0.98)),
                xytext=(rank + 40, y_vis_hi * (0.92 - iso_idx * 0.18)),
                fontsize=11, fontweight='bold', color=style['color'],
                ha='left', va='center',
                arrowprops=dict(arrowstyle='->', color=style['color'], lw=0.9),
                clip_on=False
            )
            break

        for rank in range(len(y_iso)):
            real_x = len(y_iso) - 1 - rank
            val = y_iso[real_x]
            if val >= 0:
                break
            ax_new.annotate(
                feat_names_iso[real_x],
                xy=(real_x, max(val, y_vis_lo * 0.98)),
                xytext=(real_x - 200, y_vis_lo * (0.92 - iso_idx * 0.18)),
                fontsize=11, fontweight='bold', color=style['color'],
                ha='left', va='center',
                arrowprops=dict(arrowstyle='->', color=style['color'], lw=0.9),
                clip_on=False
            )
            break
    legend_elements_new = [
        patches.Patch(facecolor=style['color'], edgecolor='dimgray', label=alt_id)
        for alt_id, _, style in isoform_plot_data
    ]
    ax_new.legend(
        handles=legend_elements_new,
        loc='center left', bbox_to_anchor=(1.01, 0.5), frameon=False,
        prop={'size': 18, 'weight': 'bold'},
        title='Alternative Isoforms',
        title_fontproperties={'size': 20, 'weight': 'bold'}
    )
output_filename = os.path.join(OUT_DIR, f"{OUT_PREFIX}_figure.png")
plt.savefig(output_filename, dpi=300, bbox_inches='tight', facecolor='white')

print("[IsoImpact] Analysis complete.")
print(f"[IsoImpact] Figure written to {output_filename}")
