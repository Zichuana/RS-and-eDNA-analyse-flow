#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
alpha 多样性 4 指标 1 行 4 图
- 每个图用不同可爱色系
- 横坐标文本自动 45° 倾斜，不重叠
- 保存为横向长图，适合论文插图
"""
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

#  ------------------------------------------------------------------
df = pd.read_excel("alpha.xlsx", index_col=0)
df.columns = df.columns.str.strip()
if 'Group' in df.columns:
    df = df.rename(columns={'Group': 'group'})
df['group'] = df['group'].astype(str)

# ----------------------------------------------------------
metrics = ['simpson_diversity', 'shannon_entropy', 'pielou_evenness', 'faith_pd']
rmb_palettes = {
    # 'observed_features': ['#7B1FA2', '#9C27B0', '#BA68C8', '#E1BEE7'],
    # 100 元红
    'simpson_diversity': ['#C83C3C', '#D85A5A', '#E37373', '#ED8B8B'],
    # 50 元绿
    'shannon_entropy':   ['#2E7D32', '#43A047', '#66BB6A', '#81C784'],
    # 20 元棕
    'pielou_evenness':   ['#8D6E63', '#A1887F', '#BCAAA4', '#D7CCC8'],
    # 10 元蓝
    'faith_pd':          ['#1565C0', '#1976D2', '#42A5F5', '#64B5F6']
}


# ----------------------------------------------------------
fig, axes = plt.subplots(1, 4, figsize=(10, 4), sharey=False)   # 横向长图
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size']   = 12

for ax, m in zip(axes, metrics):
    # 散点
    sns.stripplot(
        data=df, x='group', y=m,
        color='black', edgecolor='gray', size=3, alpha=0.7, jitter=True, ax=ax,
        order=["2500-3000", "3000-3500", "3500-4000", "4000-4500", "4500-5000"]
    )
    # 箱子
    sns.boxplot(
        data=df, x='group', y=m,
        hue='group', palette=rmb_palettes[m],
        width=0.5, flierprops={"marker": ""}, legend=False, ax=ax,
        order=["2500-3000", "3000-3500", "3500-4000", "4000-4500", "4500-5000"]
    )
    # 标题 & 轴标签
    ax.set_title(m.replace('_', ' ').title(), fontsize=12, weight='bold')
    ax.set_xlabel('Group', fontsize=12)
    ax.set_ylabel(m.replace('_', ' ').title(), fontsize=12)
    # 长标签倾斜
    ax.tick_params(axis='x', rotation=45)
    # 去掉上边右边边框
    sns.despine(ax=ax, trim=True)

# ----------------------------------------------------------
plt.tight_layout()
plt.savefig('alpha_metrics_cute_row.png', dpi=300)
plt.show()

