import traceback
import os

from .genome import genome_image
from .statistics import statistics_image
from .load_json import load_json


def generate_report(params, df):
    try:
        genome_data = load_json(
            os.path.abspath(os.path.join(os.getcwd(), f"./dataset/{params['genome']}.json"))
        )

        # 生成系谱在染色体上的分布图
        centromere_index = {}
        chr_groups = df.groupby('chr_id')
        for chr_id, chr_df in chr_groups:
            distances = []
            chr_htp_groups = chr_df.groupby("HTP")
            chr_centromere = genome_data[f"chr{chr_id}"]["centromere"]
            for htp_name, htp_df in chr_htp_groups:
                htp_range = [htp_df["position"].tolist()[0], htp_df["position"].tolist()[-1]]
                distances.append(abs(htp_range[0] - chr_centromere[0]))

            centromere_index[str(chr_id)] = distances.index(min(distances))
        
        pedigree_data = []
        for chr_id, chr_df in chr_groups:
            chr_data = []
            chr_index = 0
            chr_htp_groups = chr_df.groupby("HTP")
            htp_index = 0
            for htp_name, htp_df in chr_htp_groups:
                if htp_index == centromere_index[str(chr_id)]:
                    chr_data+=["centromere"] * 5
                    chr_data.append(htp_df["HTP_Source"].tolist()[0])
                else:
                    chr_data.append(htp_df["HTP_Source"].tolist()[0])
                htp_index += 1

            pedigree_data.append(chr_data)

            chr_index += 1

        # 格式化pedigree_data
        max_cols = 0
        for row in pedigree_data:
            if len(row) > max_cols:
                max_cols = len(row)

        for index, row in enumerate(pedigree_data):
            pedigree_data[index] = row + [""]* (max_cols-len(row))

        genome_image(
            pedigree_data,
            params["output"],
            params["colors"]
        )

        # 生成统计图
        statistics_image(
            pedigree_data,
            params["output"],
            params["colors"]
        )

    except Exception as e:
        traceback.print_exc(e)