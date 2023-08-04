from pandas import DataFrame
from numpy import array, zeros, arange
import matplotlib.pyplot as plt
from itertools import chain
from matplotlib.patches import Patch


def statistics_image(data, output, color_map):
    fig = plt.figure(figsize=(12, 6))
    grid = plt.GridSpec(2, 3, hspace=0.2, wspace=0.4)

    ax1 = fig.add_subplot(grid[:2, :2])
    ax2 = fig.add_subplot(grid[0, 2])
    ax3 = fig.add_subplot(grid[1, 2])

    chrs = [f"chr{i+1}" for i in range(len(data))]
    chain_data = list(chain.from_iterable(data))
    
    pie_data = []
    pie_labels = []

    plot_data = {}

    data_sets = list(set(chain_data))

    if "" in data_sets:
        data_sets.remove("")
    if "centromere" in data_sets:
        data_sets.remove("centromere")

    ancestor_counts = {}
    for key in data_sets:
        ancestor_counts[key] = array([row.count(key) for row in data])
        pie_data.append(chain_data.count(key))
        pie_labels.append(key)

        plot_data[key] = []
        for index, values in enumerate(data):
            plot_data[key].append(values.count(key))


    statistics_df = DataFrame({"Chr": chrs, **ancestor_counts})
    statistics_df.to_csv(f"{output}/statistics.txt", sep="\t", index=False, encoding="utf-8_sig")

    width = 0.8

    bottom = zeros(len(data))

    for ancestor, ancestor_count in ancestor_counts.items():
        p = ax1.bar(
            chrs, 
            ancestor_count, 
            width, 
            label=ancestor, 
            bottom=bottom,
            color=color_map[ancestor]
        )
        bottom += ancestor_count

        ax1.bar_label(p, label_type='center')

    # 饼图
    ax2.pie(
        pie_data, 
        labels=pie_labels, 
        autopct='%1.1f%%', 
        shadow=True,
        colors=[color_map[key] for key in pie_labels]
    )

    # 折线图
    for key, value in plot_data.items():
        ax3.plot(
            arange(1, len(data)+1), 
            value, 
            'o-', 
            color=color_map[key], 
            label=key
        )
        ax3.set_xticks(arange(1, len(data)+1))

    # legend

    legend_elements = []

    for key, color in color_map.items():
        legend_elements.append(
            Patch(facecolor=color, edgecolor=color, label=key)
        )

    ax2.legend(
        handles=legend_elements,
        bbox_to_anchor=(1, 0, 0.5, 1),
        loc="center left",
        prop={"size": 12}
    )

    plt.savefig('{}/statistics.png'.format(output), format='png', dpi=300)
    plt.close()