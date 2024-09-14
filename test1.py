import pandas as pd
import numpy as np
from scipy.optimize import curve_fit
import statsmodels.api as sm

# 读取Excel文件
file_path = 'D:/新建文件夹/photoperoid relative/1-E.xlsx'
data = pd.read_excel(file_path)

# 时间点（对应3, 7, 11, 15, 19, 23小时）
time_points = np.array([3, 7, 11, 15, 19, 23])

# 将时间点转换为弧度
time_in_radians = 2 * np.pi * time_points / 24

# 定义 Cosinor 模型函数
def cosinor_model(time, mesor, amplitude, acrophase):
    return mesor + amplitude * np.cos(time - acrophase)

# 保存每个基因的节律分析结果
results = []

# 为每个基因进行分析
for i, row in data.iterrows():
    gene_id = row['ID']

    # 提取每个时间点的重复样本表达量的平均值
    gene_expression = np.array([
        row[['1-3-E1_FPKM', '1-3-E2_FPKM', '1-3-E3_FPKM', '1-3-E4_FPKM']].mean(),
        row[['1-7-E1_FPKM', '1-7-E2_FPKM', '1-7-E3_FPKM', '1-7-E4_FPKM']].mean(),
        row[['1-11-E1_FPKM', '1-11-E2_FPKM', '1-11-E3_FPKM', '1-11-E4_FPKM']].mean(),
        row[['1-15-E1_FPKM', '1-15-E2_FPKM', '1-15-E3_FPKM', '1-15-E4_FPKM']].mean(),
        row[['1-19-E1_FPKM', '1-19-E2_FPKM', '1-19-E3_FPKM', '1-19-E4_FPKM']].mean(),
        row[['1-23-E1_FPKM', '1-23-E2_FPKM', '1-23-E3_FPKM', '1-23-E4_FPKM']].mean()
    ])

    # 检查基因表达的标准差是否足够大，避免分析那些波动极小的基因
    if np.std(gene_expression) < 1e-3:
        results.append({
            'Gene_ID': gene_id,
            'Mesor': np.nan,
            'Amplitude': np.nan,
            'Acrophase': np.nan,
            'P_value': np.nan,
            'Rhythmic': False
        })
        continue

    # 初始参数猜测
    initial_guess = [
        gene_expression.mean(),  # mesor
        (gene_expression.max() - gene_expression.min()) / 2,  # amplitude
        0  # acrophase
    ]

    try:
        # 拟合 Cosinor 模型
        params, params_covariance = curve_fit(
            cosinor_model,
            time_in_radians,
            gene_expression,
            p0=initial_guess
        )

        # 提取参数
        mesor, amplitude, acrophase = params

        # 构建用于线性回归的余弦分量和正弦分量
        cos_component = np.cos(time_in_radians - acrophase)
        sin_component = np.sin(time_in_radians - acrophase)

        # 构建线性模型
        model = sm.OLS(gene_expression, sm.add_constant(np.column_stack([cos_component, sin_component])))
        results_model = model.fit()

        # 提取 F-检验的p值
        p_value = results_model.f_pvalue

        # 判断是否有显著的节律性，通常 p 值小于 0.05 表示显著
        rhythmic = p_value < 0.05

        # 将结果添加到列表
        results.append({
            'Gene_ID': gene_id,
            'Mesor': mesor,
            'Amplitude': amplitude,
            'Acrophase': acrophase,
            'P_value': p_value,
            'Rhythmic': rhythmic
        })

    except RuntimeError:
        # 如果拟合失败，记录结果为NaN
        results.append({
            'Gene_ID': gene_id,
            'Mesor': np.nan,
            'Amplitude': np.nan,
            'Acrophase': np.nan,
            'P_value': np.nan,
            'Rhythmic': False
        })

# 将结果保存为DataFrame
results_df = pd.DataFrame(results)

# 将结果保存为Excel文件
output_path = 'D:/新建文件夹/photoperoid relative/1-E-rhythm.xlsx'
results_df.to_excel(output_path, index=False)

print(f"结果已保存至 {output_path}")

