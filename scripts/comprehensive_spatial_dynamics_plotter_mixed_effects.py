"""Plot a descriptive regional ECM profile from section-region summaries.

Run this script after HistologyData_Pipeline_MixedEffects.py.

This figure is descriptive only. It displays individual section-region values
and mean +/- SD; it does not fit a continuous spatial model or calculate
correlations. Formal inference is performed by the mixed-effects models in
HistologyData_Pipeline_MixedEffects.py.
"""

import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


# GitHub repository layout: this script is stored in the scripts directory.
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
    'Necrosis': {'marker': 's', 'color': '#d62728', 'offset': -0.16},
    'Mucin': {'marker': 'o', 'color': '#2ca02c', 'offset': 0.00},
    'Collagen': {'marker': '^', 'color': '#1f77b4', 'offset': 0.16},
}

sns.set_theme(style='ticks')
plt.rcParams['font.sans-serif'] = 'Arial'
plt.rcParams['axes.unicode_minus'] = False

figure, axis = plt.subplots(figsize=(9, 6))
region_positions = {region: index for index, region in enumerate(region_order)}
rng = np.random.default_rng(20260330)

for metric, style in metric_styles.items():
    metric_data = plot_data[plot_data['Metric'] == metric]
    if metric_data.empty:
        print(f'Warning: no data found for {metric}.')
        continue

    # Individual points are section-region observations: n = 6 for each wound
    # region and n = 4 for Normal Skin after averaging its three ROIs.
    x_values = np.array(
        [region_positions[str(region)] for region in metric_data['Region']]
    )
    x_values = x_values + style['offset'] + rng.uniform(-0.035, 0.035, len(metric_data))
    axis.scatter(
        x_values,
        metric_data['Coverage'],
        marker=style['marker'],
        s=38,
        color=style['color'],
        edgecolor='black',
        linewidth=0.45,
        alpha=0.62,
        zorder=2,
    )

    summary = (
        metric_data.groupby('Region', observed=True)['Coverage']
        .agg(['mean', 'std'])
        .reindex(region_order)
    )
    positions = np.array([region_positions[region] for region in region_order])
    positions = positions + style['offset']
    axis.plot(
        positions,
        summary['mean'],
        marker=style['marker'],
        color=style['color'],
        label=metric,
        linewidth=2.0,
        markersize=7,
        zorder=4,
    )
    axis.errorbar(
        positions,
        summary['mean'],
        yerr=summary['std'],
        fmt='none',
        ecolor=style['color'],
        elinewidth=1.35,
        capsize=4,
        capthick=1.35,
        zorder=3,
    )

# Normal Skin is an independent reference group, not the next measured point
# along a continuous wound-distance scale.
axis.axvline(2.5, color='0.55', linestyle='--', linewidth=1.0, zorder=1)
axis.set_xticks(range(len(region_order)), region_order)
axis.set_ylabel('Tissue Coverage (%)', fontsize=12)
axis.set_xlabel('')
axis.set_ylim(bottom=0)
axis.grid(axis='y', linestyle='--', alpha=0.45)
sns.despine()
axis.legend(title='Tissue Component', loc='upper right')
figure.tight_layout()
figure.savefig(OUTPUT_PATH, dpi=300, bbox_inches='tight')
plt.close(figure)

print(f'Regional ECM profile figure saved to: {OUTPUT_PATH}')
