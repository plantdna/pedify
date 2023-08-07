import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
sns.set_theme(style="whitegrid")


def get_change_count(data):
    count = 0
    record_value = ""
    for index, v in enumerate(data):
        if index == 0:
            record_value = v
        else:
            if v != record_value:
                record_value = v
                count+=1
    return count


def exchange_report(dir_path, name):
    data = {}
    for chr_id in range(1, 11):
        chr_df = pd.read_csv(f"{dir_path}/{name}/htp_denoise/chr{chr_id}.csv")
        chr_counts = []
        for index, row in chr_df.iterrows():
            row_values = row.values.tolist()[1:]
            chr_counts.append(get_change_count(row_values))
        
        data[f"chr{chr_id}"] = chr_counts

    df = pd.DataFrame(data)

    df.to_csv(f"{dir_path}/{name}.txt", index=False, encoding="utf-8_sig", sep="\t")

    # Set up the matplotlib figure
    f, ax = plt.subplots(figsize=(11, 6))

    # Draw a violinplot with a narrower bandwidth than the default
    sns.violinplot(data=df, palette="Set3", bw=.2, cut=2, linewidth=1)

    sns.despine(left=True, bottom=True)

    plt.savefig(f"{dir_path}/{name}.png", dpi=300, format="png")
