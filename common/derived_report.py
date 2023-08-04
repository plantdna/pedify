
import pandas as pd
import numpy as np
import traceback
import matplotlib.pyplot as plt
import os
from .color_transform import hex_to_rgb

config = {}


null_cols = 20
null_rows = 1

fig, ax = plt.subplots(figsize=(100, 40), dpi=300)


def draw_image(data, ylabels):
    heatmap_data = np.array(data)
    color_map = {
        0: np.array([255, 255, 255]),
        1: np.array(hex_to_rgb(config["param"]["colors"]["ancestor"])),
        2: np.array(hex_to_rgb(config["param"]["colors"]["derived"])),
    }

    data_3d = np.ndarray(shape=(heatmap_data.shape[0], heatmap_data.shape[1], 3), dtype=int)

    for i in range(0, heatmap_data.shape[0]):
        for j in range(0, heatmap_data.shape[1]):
            data_3d[i][j] = color_map[heatmap_data[i][j]]

    ax.imshow(data_3d, interpolation=None)
    ax.set_aspect('auto')

    ax.set_xticks([])
    ax.set_yticks(range(0, len(data), 1))

    # 设置ylabels
    ax.set_yticklabels(ylabels, fontdict={"fontsize": 48, "fontweight": 800} )

    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['right'].set_visible(False)




def heatmap(df):
    rates = {config["param"]["ancestorName"]: 1}
    cols = df.columns[4:-1]
    for col in cols:
        rates[col] = round(df[col].tolist().count(config["param"]["ancestorName"])/len(df), 4)

    image_df = pd.DataFrame()

    for chr_id in range(1, 11):
        chr_df = df[df['Chr_id']==chr_id].copy()
        chr_df_t = chr_df.T
        chr_df_t.reset_index(inplace=True)

        htps = chr_df_t.iloc[0].tolist()[1:]

        chr_df_t.drop([1, 2, 3, len(chr_df_t)-1], axis=0, inplace=True)

        chr_df_t.iloc[0] = [config["param"]["ancestorName"]] * chr_df_t.shape[1]


        chr_df_t["rate"] = chr_df_t["index"].apply(lambda x: rates[x])

        chr_df_t.sort_values(by=['rate'], ascending=False, inplace=True)
        chr_df_t.drop('rate', axis=1, inplace=True)

        chr_df_t.columns=["code"]+htps

        for index, row in chr_df_t.iterrows():
            for htp in htps:
                if row[htp] == config["param"]["ancestorName"]:
                    chr_df_t.at[index, htp] = 1
                else:
                    chr_df_t.at[index, htp] = 2


        for htp in htps:
            chr_df_t[htp] = chr_df_t[htp].astype(int)

        if chr_id == 1:
            image_df = chr_df_t
        else:
            image_df = pd.merge(image_df, chr_df_t, on="code")

        # 染色体之间插入空列
        if chr_id != 10:
            null_col_df = chr_df_t[["code"]].copy()
            for nc in range(null_cols):
                null_col_df[f"chr{chr_id}{nc}"] = 0

            image_df = pd.merge(image_df, null_col_df, on="code")

    image_data = []
    ylabels = []

    for index, row in image_df.iterrows():
        row_data = row.values.tolist()[1:]
        image_data.append(row_data)
        image_data.append(row_data)
        image_data.append(row_data)

        ylabels.append("")
        ylabels.append(row["code"])
        ylabels.append("")

        for i in range(null_rows):
            image_data.append([0]*len(row_data))
            ylabels.append("")

    draw_image(image_data, ylabels)


    plt.savefig(f'{config["output"]}/DPSA.png', format='png')
    plt.close()




def derived_report(params):
    try:
        global config
        config = params

        ndf = pd.DataFrame()

        filepath = [f"{params['output']}/{fp}" for fp in os.listdir(params["output"]) if fp[-3:] == "txt"]

        for index, fp in enumerate(filepath):
            data = []
            df = pd.read_csv(fp, sep=params["sep"])

            df_chr_group = df.groupby('chr_id')
            for chr_id, chr_df in df_chr_group:
                if params["HTP"]:
                    chr_htp_group = chr_df.groupby("HTP")
                else:
                    chr_htp_group = chr_df.groupby("locus")

                for htp, htp_df in chr_htp_group:
                    if index == 0:
                        data.append([
                            htp,
                            chr_id,
                            htp_df['position'].tolist()[0],
                            htp_df['position'].tolist()[-1],
                            htp_df['HTP_Source'].tolist()[0] if params["HTP"] else htp_df['source'].tolist()[0],
                        ])
                    else:
                        data.append([
                            htp,
                            htp_df['HTP_Source'].tolist()[0] if params["HTP"] else htp_df['source'].tolist()[0],
                        ])
            if index == 0:
                sdf = pd.DataFrame(data, columns=['Locus', 'Chr_id', 'Start', 'Stop', fp.split('/')[-1].split('.')[0]])
                ndf = sdf
            else:
                sdf = pd.DataFrame(data, columns=['Locus', fp.split('/')[-1].split('.')[0]])
                ndf = pd.merge(ndf, sdf, on="Locus")

        ndf["Density"] = ndf.apply(lambda x: round(x.values.tolist()[4:].count(params["param"]["ancestorName"]) / len(x.values.tolist()[4:]), 3), axis=1)
        ndf.sort_values(by=["Chr_id","Start"], ascending=[True, True], inplace=True)
        ndf.to_csv(f"{params['output']}/statistics.txt", index=False, encoding="utf-8_sig", sep=params["sep"])
        heatmap(ndf)

    except Exception as e:
        traceback.print_exc(e)
