"""Plot ECM spatial dynamics from the ROI-aware mixed-effects pipeline output.

Run this script after HistologyData_Pipeline_MixedEffects.py.
"""

import os

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


PROJECT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS_DIR = os.path.join(PROJECT_DIR, 'outputs', 'ecm')
INPUT_PATH = os.path.join(RESULTS_DIR, '0_SectionRegion_Summary.csv')
OUTPUT_PATH = os.path.join(RESULTS_DIR, 'Comprehensive_3Factor_Spatial_Dynamics.png')

if not os.path.exists(INPUT_PATH):
    raise FileNotFoundError(
        '0_SectionRegion_Summary.csv was not found. Run '
        'HistologyData_Pipeline_MixedEffects.py first.'
    )

print(f'Loading section-level ECM summaries from: {INPUT_PATH}')
data = pd.read_csv(INPUT_PATH)

required_columns = {'Metric', 'Sample_ID', 'region', 'Coverage'}
missing_columns = required_columns - set(data.columns)
if missing_columns:
    raise ValueError(f'Missing required columns: {sorted(missing_columns)}')

region_display_map = {
    'Wound Skin': 'Wound Centre',
    'Wound Centre': 'Wound Centre',
    'Proximal Edge': 'Proximal Edge',
    'Distant Edge': 'Distant Edge',
    'Normal Skin': 'Normal Skin',
}
region_order = ['Wound Centre', 'Proximal Edge', 'Distant Edge', 'Normal Skin']

plot_data = data.copy()
plot_data['Region'] = plot_data['region'].map(region_display_map)
plot_data = plot_data.dropna(subset=['Region', 'Coverage'])
plot_data['Region'] = pd.Categorical(
    plot_data['Region'], categories=region_order, ordered=True
)

metric_styles = {
    'Necrosis': {'marker': 's', 'color': '#d62728'},
    'Mucin': {'marker': 'o', 'color': '#2ca02c'},
    'Collagen': {'marker': '^', 'color': '#1f77b4'},
}

sns.set_theme(style='ticks')
plt.rcParams['font.sans-serif'] = 'Arial'
plt.rcParams['axes.unicode_minus'] = False

plt.figure(figsize=(10, 6))
for metric, style in metric_styles.items():
    metric_data = plot_data[plot_data['Metric'] == metric]
    if metric_data.empty:
        print(f'Warning: no data found for {metric}.')
        continue

    sns.lineplot(
        data=metric_data,
        x='Region',
        y='Coverage',
        marker=style['marker'],
        color=style['color'],
        label=metric,
        linewidth=2.5,
        markersize=10,
        errorbar=('ci', 95),
    )

# No figure title or x-axis title: use the thesis figure caption instead.
plt.ylabel('Tissue Coverage (%)', fontsize=12)
plt.xlabel('')
plt.ylim(bottom=0)
plt.grid(axis='y', linestyle='--', alpha=0.7)
sns.despine()
plt.legend(title='Tissue Component', loc='upper right')
plt.tight_layout()
plt.savefig(OUTPUT_PATH, dpi=300, bbox_inches='tight')
plt.close()

print(f'Spatial dynamics figure saved to: {OUTPUT_PATH}')
