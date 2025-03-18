# -*- coding: utf-8 -*-
import h5py
import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import SGDClassifier
from sklearn.metrics import classification_report
import joblib
import time

def incremental_fit(model, scaler, X_chunk, y_chunk, all_classes, is_first_chunk=False):
    """增量拟合数据和模型"""
    # 每次分块都更新缩放器
    scaler.partial_fit(X_chunk)  # 修改点：移除is_first_chunk判断
    
    X_scaled = scaler.transform(X_chunk)
    
    if is_first_chunk:
        model.partial_fit(X_scaled, y_chunk, classes=all_classes)
    else:
        model.partial_fit(X_scaled, y_chunk)
    
    return model, scaler

def load_and_align_data(hdf5_path, metadata_path):
    """数据加载与对齐核心函数"""
    # 加载并清洗元数据
    metadata = pd.read_csv(metadata_path)
    metadata["sample_name"] = metadata["sample_name"].str.strip()
    
    # 加载并清洗HDF5样本
    with h5py.File(hdf5_path, "r") as f:
        hdf5_samples = np.char.strip(f["/data/samples"][:].astype(str))
    
    # 创建有效样本掩码
    valid_mask = np.isin(hdf5_samples, metadata["sample_name"])
    valid_sample_ids = hdf5_samples[valid_mask]
    
    # 生成数据对齐报告
    mismatch_count = len(hdf5_samples) - len(valid_sample_ids)
    print(f"\n=== 数据对齐报告 ===")
    print(f"原始HDF5样本数: {len(hdf5_samples):,}")
    print(f"有效样本数:     {len(valid_sample_ids):,}")
    print(f"缺失样本数:     {mismatch_count} ({mismatch_count/len(hdf5_samples):.2%})")
    
    # 过滤并重新对齐元数据
    metadata = metadata.set_index("sample_name").loc[valid_sample_ids].reset_index()
    return metadata, valid_sample_ids, valid_mask

def main():
    # 配置参数
    CHUNK_SIZE = 10000           # 根据内存调整
    TEST_SIZE = 0.2
    RANDOM_STATE = 42
    MIN_SAMPLES_PER_CLASS = 2   # 新增：最小样本阈值
    HDF5_PATH = "./expression_matrix.hdf5"
    METADATA_PATH = "./metadata.csv"
    
    # 1. 数据加载与对齐
    start_time = time.time()
    metadata, valid_samples, valid_mask = load_and_align_data(HDF5_PATH, METADATA_PATH)
    y = metadata["cell_type_designation_label"].values
    
    # 类别过滤系统
    print("\n=== 类别分布分析 ===")
    unique_classes, class_counts = np.unique(y, return_counts=True)
    print(f"原始类别数: {len(unique_classes)}")
    print(f"最小类别样本数: {class_counts.min()}")
    
    # 过滤样本数不足的类别
    valid_classes_mask = class_counts >= MIN_SAMPLES_PER_CLASS
    valid_classes = unique_classes[valid_classes_mask]
    class_filter = np.isin(y, valid_classes)
    
    # 应用过滤
    y_filtered = y[class_filter]
    valid_samples = valid_samples[class_filter]
    metadata = metadata[class_filter]
    original_indices = np.where(valid_mask)[0][class_filter]  # 修正原始索引
    assert len(np.unique(y_filtered)) == len(valid_classes), "过滤后类别数异常"
    
    print(f"\n=== 过滤后数据 ===")
    print(f"保留类别数: {len(valid_classes)}")
    print(f"有效样本数: {len(y_filtered):,}")
    print(f"过滤样本数: {len(y) - len(y_filtered)}")
    

    
    print(f"\n=== 过滤后数据 ===")
    print(f"保留类别数: {len(valid_classes)}")
    print(f"有效样本数: {len(y_filtered):,}")
    print(f"过滤样本数: {len(y) - len(y_filtered)}")
    
    # 2. 数据划分（保持分层抽样）
    indices = np.arange(len(valid_samples))
    try:
        train_idx, test_idx = train_test_split(
            indices, 
            test_size=TEST_SIZE, 
            random_state=RANDOM_STATE,
            stratify=y_filtered
        )
    except ValueError as e:
        print(f"\n分层抽样失败: {e}")
        print("正在使用常规随机划分...")
        train_idx, test_idx = train_test_split(
            indices, 
            test_size=TEST_SIZE, 
            random_state=RANDOM_STATE
        )
    train_idx.sort()
    
    print(f"\n=== 训练配置 ===")
    print(f"总样本数:      {len(valid_samples):,}")
    print(f"训练集样本数:  {len(train_idx):,}")
    print(f"测试集样本数:  {len(test_idx):,}")
    print(f"分块大小:      {CHUNK_SIZE}")
    print(f"预计训练批次:  {len(train_idx)//CHUNK_SIZE + 1}")
    
    # 3. 模型初始化
    scaler = StandardScaler()
    model = SGDClassifier(
        loss='log_loss',
        penalty='l2',
        learning_rate='adaptive',
        eta0=0.01,
        max_iter=1000,
        n_jobs=-1,
        random_state=RANDOM_STATE,
        tol=1e-4
    )
    
    # 4. 训练流程
    print("\n=== 开始训练 ===")
    with h5py.File(HDF5_PATH, "r") as f:
        dset = f["/data/counts"]
        
        # 分块训练
        for batch_num, i in enumerate(range(0, len(train_idx), CHUNK_SIZE), 1):
            batch_start = time.time()
            
            # 获取当前分块的映射索引
            chunk_train_idx = train_idx[i:i+CHUNK_SIZE]
            hdf5_idx = original_indices[chunk_train_idx]
            
            # 读取数据（基因×细胞 → 细胞×基因）
            X_chunk = dset[:, hdf5_idx].T.astype(np.float32)
            y_chunk = y_filtered[chunk_train_idx]
            
            # 处理缺失值
            if np.isnan(X_chunk).any():
                X_chunk = np.nan_to_num(X_chunk)
            
           # 增量训练
            is_first_chunk = (batch_num == 1)
            # !!! 传递valid_classes到训练函数
            model, scaler = incremental_fit(
                model, scaler, 
                X_chunk, y_chunk,
                all_classes=valid_classes,  # !!! 关键修改
                is_first_chunk=is_first_chunk
            )
            
            # 进度跟踪
            batch_time = time.time() - batch_start
            processed = min(i+CHUNK_SIZE, len(train_idx))
            print(f"批次 {batch_num:03d} | 已处理 {processed:,}/{len(train_idx):,} "
                  f"({processed/len(train_idx):.1%}) | 耗时 {batch_time:.1f}s")

    # 5. 保存模型
    pipeline = Pipeline([('scaler', scaler), ('classifier', model)])
    model_name = f"models/scClassifier_{time.strftime('%Y%m%d-%H%M')}.pkl"
    joblib.dump(pipeline, model_name)
    print(f"\n模型已保存至: {model_name}")

    # 6. 验证评估
    def evaluate_in_chunks():
        y_true, y_pred = [], []
        with h5py.File(HDF5_PATH, "r") as f:
            dset = f["/data/counts"]
            original_test_idx = original_indices[test_idx]
            
            for i in range(0, len(original_test_idx), CHUNK_SIZE):
                chunk_idx = original_test_idx[i:i+CHUNK_SIZE]
                X_chunk = dset[:, chunk_idx].T.astype(np.float32)
                X_scaled = scaler.transform(X_chunk)
                y_true.extend(y_filtered[test_idx[i:i+CHUNK_SIZE]])
                y_pred.extend(model.predict(X_scaled))
        return y_true, y_pred
    
    print("\n=== 开始验证 ===")
    y_true, y_pred = evaluate_in_chunks()
    print("\n分类报告:")
    print(classification_report(y_true, y_pred, target_names=np.unique(y_true)))
    print(f"总耗时: {time.time()-start_time:.1f}秒")

if __name__ == "__main__":
    main()
