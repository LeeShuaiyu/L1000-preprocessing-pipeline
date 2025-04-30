
* **SMILES → 重复次数分布**（解释 “ threshold = 3 仅剩 7 000 条” 的原因）  
* **新参数**：`MIN_REPLICA = 1`（保留全部 10 670 个唯一 SMILES）

---

# LINCS-L1000 数据预处理与 SMILES 聚合管线  
*(阶段性总结 v2  2025-04-30)*

---

## 1   原始数据概览

| 文件 | 维度 / 行数 | 关键字段 | 说明 |
|------|-------------|----------|------|
| **level5_beta_trt_cp_n720216x12328.gctx** | 12328 genes × 720 216 columns | `cid` | 处理后表达 z-scores |
| **instinfo_beta.txt** | 3 026 460 rows | `sample_id`, `pert_id`, `pert_type`, `qc_pass` | 每列的实验元数据 |
| **compoundinfo_beta.txt** | 19 811 rows | `pert_id`, `canonical_smiles` | 化合物字典 |
| **geneinfo_beta.txt** | 12 328 rows | `gene_id`, `gene_symbol`, `feature_space` | 行注释 (`landmark` = 978) |

---

## 2   cid 解析与匹配策略

| cid 形式 | 例子 | 解析键 | 占比 |
|----------|------|--------|-----:|
| `plate:well` (1 colon) | `AICHI001_BJAB_24H:A07` | **plate-well → pert_id** 哈希 | 462 523 |
| `plate:pert_id:dose` (2 colons) | `AML001_CD34_24H:BRD-A03772856:10` | `parts[1]` = `pert_id` | 256 893 |
| `plate:pert_id:dose:time` (3 colons) | `ABY001_A375_XH:BRD-A61304759:0.625:24` | 同上 | 800 |

### 最终正则（Mac / Win 均通过）

```python
clean_plate = re.compile(
    r'(_X\d+)?'          # 重复编号
    r'(_F\d+B\d+)?'      # 子板块
    r'(_DUO.*)?$'        # 双孔标签
).sub
plate_core = clean_plate('', raw_plate)      # e.g. AICHI001_BJAB_24H
key = f'{plate_core}:{well}'
```

---

## 3   匹配结果

| 指标 | 数量 | 占 720 216 (%) |
|------|------|---------------|
| 列 → `pert_id` 成功 | **699 937** | 97.2 |
| 其中 `pert_type == 'trt_cp'` | **465 354** | 64.6 |
| 且有 `canonical_smiles` | **465 354** | 64.6 |
| **唯一 SMILES** | **10 670** | — |

> 说明：19 811 个化合物 ID 中约一半映射到同一 canonical SMILES——盐酸盐、批次前缀等已折叠。

---

## 4   重复次数分布（trt_cp + SMILES）

| 重复次数 (`expr_cnt`) | SMILES 数 | 累积占比 |
|----------------------|----------:|----------|
| **1 次** | 2 958 | 27.7 % |
| **2 次** | 742 | 34.7 % |
| **3–5 次** | 2 892 | 61.8 % |
| **≥ 6 次** | 4 078 | 100 % |

> *原因*：LINCS Phase II 中大量化合物只测 **1–2 条剂量-时间组合**，导致 **threshold = 3** 时仅保留 7 012 条 SMILES。  
> 科研场景往往需要 **覆盖性优先**：因此本版本将  
> ```python
> MIN_REPLICA = 1
> ```  
> 并在模型中使用 **样本数加权** 或 **校正方差** 来抵消重复不足的影响。

---

## 5   Pipeline 流程（适配 16 GB RAM）

| 步骤 (Cell) | 任务 | 关键点 & 进度监控 |
|-------------|------|------------------|
| **0** | 路径 & 参数 | `CHUNK_SIZE = 2_000`, `float32` |
| **1** | 导入依赖 | — |
| **2** | *plate-well → pert_id* 哈希、`trt_cp_set`, `pid2smiles` | `print` 键数 |
| **3** | 读取 geneinfo | 行注释、landmark mask |
| **4** | **Pass-1** 扫描列元数据，映射到 SMILES | `tqdm` 进度；输出可用列/SMILES |
| **5** | 初始化 **memmap** 累加矩阵 | 峰值内存 ≤ 4 GB |
| **6** | **Pass-2** 分块 2 000 列累加 | `tqdm` “accumulate” |
| **7** | 计算均值，过滤重复 < MIN_REPLICA | 统计保留 SMILES 数 |
| **8** | 导出 **12328 × SMILES** CSV | 每行 `gene_id, gene_symbol, feature_space, …` |
| **9** | 导出 **978 Landmark** CSV | — |

---

## 6   输出文件

| 文件 | 尺寸 | 描述 |
|------|------|------|
| `expr_by_smiles_12328.csv` | ~ 1.1 GB | 12 328 genes × 10 670 SMILES |
| `expr_by_smiles_978.csv` | ~ 90 MB | 978 genes × 10 670 SMILES |
| `unmatched_cid_summary.csv` | 4 行 | R1–R3 统计 |
| `unmatched_cid_detail_examples.csv` | 90 行 | 每类 30 例示范 |

---

## 7   结论与下一步

1. **数据完整度**：65 % 列映射成功且保留全部 10 670 SMILES；足以训练化合物-转录组模型。  
2. **阈值决策**：  
   * `MIN_REPLICA = 1` → 覆盖 10 670 SMILES  
   * `MIN_REPLICA = 3` → 仅 ~7 000 SMILES（高重复、低覆盖）  
   根据实验需求切换即可。  
3. **后续工作**：  
   - PCA/UMAP 可视化 SMILES-表达空间  
   - 构建 GNN → 表达回归模型；样本权重 = 重复次数  
   - 交叉比对 `siginfo_beta.txt` 的 exemplar signature 质量标志

