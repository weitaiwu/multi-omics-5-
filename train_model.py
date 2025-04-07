import scanpy as sc
import pandas as pd
import numpy as np
import scipy.sparse as sp
from celltypist import train
import gc
import os
import tempfile  # 新增导入

os.environ['OPENBLAS_NUM_THREADS'] = '1'
from anndata import concat

def load_and_preprocess(csv_path, chunksize=20000):
    """分块加载并二次验证标准化"""
    with open(csv_path, 'r') as f:
        genes = f.readline().strip().split(',')[1:]
    
    adata = None
    # 修正点：添加 enumerate 获取迭代索引
    for i, chunk in enumerate(pd.read_csv(csv_path, chunksize=chunksize, index_col=0)):
        # 转换为32位稀疏矩阵
        X = sp.csr_matrix(chunk.values.astype(np.float32))
        
        # 分块CPM标准化（精确到每个细胞）
        sc.pp.normalize_total(sc.AnnData(X), target_sum=1e4, inplace=True)
        X.data = np.log1p(X.data)
        
        # 构建分块对象
        chunk_adata = sc.AnnData(
            X, 
            obs=pd.DataFrame(index=chunk.index.astype(str)),
            var=pd.DataFrame(index=genes)
        )
        
        # 合并数据（修正合并方向）
        adata = chunk_adata if adata is None else concat([adata, chunk_adata], axis=0)
        
        # 立即释放内存
        del X, chunk, chunk_adata
        gc.collect()
        print(f"✅ 已处理第 {i+1} 块 | 总细胞数: {adata.n_obs}")
    
    # 移除全局二次标准化（避免重复处理）
    # 验证关键指标
    sample_idx = np.random.choice(adata.n_obs)
    cell_total = adata.X[sample_idx].sum()
    print(f"随机细胞总计数（应≈1e4）: {cell_total:.0f}")
    print(f"最大值（应<15）: {adata.X.max():.2f}")
    
    return adata

# 加载数据（强制使用float32）
adata = load_and_preprocess("matrix.csv", chunksize=20000)
adata.X = adata.X.astype(np.float32)

# 加载元数据（按需加载）
cell_metadata = pd.read_csv(
    "metadata.csv", 
    index_col="sample_name", 
    usecols=['sample_name', 'cell_type_designation_label']
)
cell_metadata.index = cell_metadata.index.astype(str)

# 整合元数据（惰性索引）
common_cells = adata.obs_names[adata.obs_names.isin(cell_metadata.index)]
adata = adata[common_cells].copy()
adata.obs = cell_metadata.loc[common_cells]

# 训练前压缩数据
sc.pp.filter_genes(adata, min_cells=1)
adata.var_names_make_unique()
adata = adata[:, ~adata.var_names.duplicated()]

# 强制垃圾回收
gc.collect()

# 训练参数优化（极限内存模式）
model = train(
    adata,
    labels=adata.obs['cell_type_designation_label'],
    use_SGD=True,
    mini_batch=True,
    batch_size=200,
    epochs=30,
    n_jobs=1,
    feature_selection=False,
    top_genes=None,
    balance_cell_type=False,
    pretrain=False
)

# 分块验证（使用临时文件）
accuracy = []
for i in range(0, adata.n_obs, 5000):
    end = min(i+5000, adata.n_obs)
    chunk = adata[i:end].copy()
    chunk.X = chunk.X.astype(np.float32)
    
    with tempfile.NamedTemporaryFile(suffix=".h5ad") as tmp:
        chunk.write(tmp.name)
        pred = model.predict(sc.read_h5ad(tmp.name))
    
    acc = (pred.predicted_labels['majority_voting'] == chunk.obs['cell_type_designation_label']).mean()
    accuracy.append(acc)
    print(f"块 {i}-{end} 准确率: {acc:.2%}")
    del chunk, pred
    gc.collect()

print(f"平均准确率: {np.mean(accuracy):.2%}")
model.write("optimized_model.pkl", save_adata=False)
