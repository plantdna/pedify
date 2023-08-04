import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

from .color_transform import hex_to_rgb


def genome_image(data, output, color_map):
    cmp = {
        "": np.array([255, 255, 255]),
        "centromere": np.array([0, 0, 0])
    }
    legend_elements = []
    for key, color in color_map.items():
        legend_elements.append(
            Patch(facecolor=color, edgecolor=color, label=key)
        )
        cmp[key] = np.array(hex_to_rgb(color))

    fig, ax = plt.subplots(figsize=(20, 10))
    heatmap_data = np.array(data)
    data_3d = np.ndarray(shape=(heatmap_data.shape[0], heatmap_data.shape[1], 3), dtype=int)

    for i in range(0, heatmap_data.shape[0]):
        for j in range(0, heatmap_data.shape[1]):
            data_3d[i][j] = cmp[heatmap_data[i][j]]

    ax.imshow(data_3d, interpolation=None)
    ax.set_aspect('auto')

    ax.set_xticks([])
    ax.set_yticks(np.arange(heatmap_data.shape[0]))

    ax.set_yticklabels(['Chr{}'.format(i) for i in range(1, 11)], fontsize=36)
    ax.set_yticks(np.arange(heatmap_data.shape[0]+1)-.5, minor=True)
    ax.grid(which="minor", color="w", linestyle='-', linewidth=20)
    ax.tick_params(which="minor", bottom=False, left=False)

    ax.legend(
        handles=legend_elements,
        bbox_to_anchor=(1, 0),
        loc=3,
        prop={"size": 32}
    )

    plt.tight_layout()
    plt.savefig('{}/pedigree.png'.format(output), format='png', dpi=300)
    plt.close()
