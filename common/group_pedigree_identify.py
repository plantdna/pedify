import traceback
import multiprocessing
import os

import pandas as pd
from .snp_compare import snp_compare
from .makedirs import makedirs
from .BRE import background_reply_evaluation
from .gpi_report import exchange_report

config = {}


def write_result(df):
    try:
        if os.path.exists(f'{config["output"]}/pedigree_identify.txt'):
            df.to_csv(f'{config["output"]}/pedigree_identify.txt', index=False, encoding="utf-8_sig", mode="a+", header=None, sep=config["sep"])
        else:
            df.to_csv(f'{config["output"]}/pedigree_identify.txt', index=False, encoding="utf-8_sig", sep=config["sep"])

    except Exception as e:
        traceback.print_exc(e)
    

def chunk_compare(df, pdf, useHTP, chunk_t, output):
    try:
        results = []
        samples = df.columns.tolist()[2:]
        for sample in samples:
            item_df = pd.merge(
                df[["locus","gt",sample]],
                pdf[["locus", "chr_id", "position", "HTP"]],
                on="locus",
                how="left"
            )

            item_df.dropna(axis=0, subset=["gt", sample], inplace=True)
            item_df["CR"] = item_df.apply(lambda x: snp_compare(x["gt"], x[sample]), axis=1)

            res_list = item_df["CR"].tolist()
            snp_gs = round(res_list.count(1) / (len(res_list) - res_list.count(0)), 3) if (len(res_list) - res_list.count(0)) != 0 else 0

            if snp_gs >= 0.9:
                parent_df = chunk_t[["locus", sample]].copy()
                parent_df.columns = ["locus", "gt"]
                makedirs(f'{output}/{sample}')
                parent_df.to_csv(f'{output}/{sample}/parent.txt', index=False, encoding="utf-8_sig", sep="\t")


            if useHTP:
                htp_gss = []
                item_htp_groups = item_df.groupby("HTP")
                for htp, group_df in item_htp_groups:
                    htp_res = group_df["CR"].tolist()
                    htp_gs = round(htp_res.count(1) / (len(group_df) - htp_res.count(0)), 3) if (len(group_df) - htp_res.count(0)) != 0 else 0
                    htp_gss.append(htp_gs)

                results.append([
                    sample,
                    res_list.count(0),
                    res_list.count(1),
                    res_list.count(2),
                    res_list.count(3),
                    snp_gs,
                    round(sum(htp_gss)/len(htp_gss), 3)
                ])
            else:
                results.append([
                    sample,
                    res_list.count(0),
                    res_list.count(1),
                    res_list.count(2),
                    res_list.count(3),
                    snp_gs,
                ])
        

        result_df = pd.DataFrame(
            results, 
            columns=[
                'CODE', 
                'MISS_SNP_NUMBER', 
                'EQUAL_SNP_NUMBER', 
                'SIMILARITY_SNP_NUMBER',
                'DIFF_SNP_NUMBER',
                'SNP_GENETIC_SIMILARITY',
                'HTP_GENETIC_SIMILARITY',
            ] if useHTP else [
                'CODE', 
                'MISS_SNP_NUMBER', 
                'EQUAL_SNP_NUMBER', 
                'SIMILARITY_SNP_NUMBER',
                'DIFF_SNP_NUMBER',
                'SNP_GENETIC_SIMILARITY'
            ]
        )
        return result_df

    except Exception as e:
        traceback.print_exc(e)


def whole_library_search(df, pdf):
    try:
        pool = multiprocessing.Pool(config["cores"])
        file_reader = pd.read_csv(config["param"]["dataset"], sep=config["sep"], iterator=True) 
        loop= True
        while loop:
            try:
                chunk = file_reader.get_chunk(50)
                chunk_t = chunk.T.copy()
                chunk_t.reset_index(inplace=True)
                chunk_t.columns = chunk_t.iloc[0].tolist()
                chunk_t.drop(0, axis=0, inplace=True)
                chunk_t.rename(columns={"Name":"locus"}, inplace=True)

                ndf = pd.merge(df[["locus", "gt"]], chunk_t, on="locus", how="left")

                pool.apply_async(
                    func=chunk_compare,
                    args=(
                        ndf,
                        pdf,
                        config["HTP"],
                        chunk_t,
                        config["output"],
                    ), 
                    callback=write_result
                )

            except StopIteration:
                loop = False
        
        pool.close()
        pool.join()

    except Exception as e:
        traceback.print_exc(e)

def find_monomorphic_locus(df):
    locus = {}
    for index, row in df.iterrows():
        row_gts = row.values.tolist()[1:]

        counts = pd.value_counts(row_gts).to_dict()

        count_keys = list(counts.keys())
        count_values = list(counts.values())

        if count_values[0] / len(row_gts) >= 0.8:
            locus[row['locus']] = count_keys[0]

    return locus

def format_chr_data(folder, pdf, gdf):
    makedirs(f"{config['output']}/{folder}/chr_data")
    parent_df = pd.read_csv(f"{config['output']}/{folder}/parent.txt", sep="\t")

    parent_position_df = pd.merge(pdf[["locus", "chr_id", "position"]], parent_df, on="locus", how="right")

    group_df = pd.merge(parent_position_df, gdf, on="locus", how="right")

    group_df.dropna(subset=["chr_id"], inplace=True)

    group_df.sort_values(by=["chr_id","position"], ascending=[True, True], inplace=True)

    group_df["chr_id"] = group_df["chr_id"].astype(int)

    for chr_id in range(1, 11):
        chr_df = group_df[group_df["chr_id"]==chr_id].copy()
        chr_df.drop(["chr_id","position"], inplace=True, axis=1)
        chr_df_t = chr_df.T
        chr_df_t.reset_index(inplace=True)
        chr_df_t.columns = chr_df_t.iloc[0]
        chr_df_t.drop(0, inplace=True, axis=0)
        chr_df_t.rename(columns={"locus": "call_code"}, inplace=True)
        chr_df_t.to_csv(f"{config['output']}/{folder}/chr_data/chr{chr_id}.csv", index=False, encoding="utf-8_sig")


def group_pedigree_identify(params):
    try:
        global config
        config = params

        pdf = pd.read_csv(config["locus_file"], sep=config["sep"], low_memory=False)
        gt_df = pd.read_csv(config["param"]["target"], sep=config["sep"], low_memory=False)

        monomorphic_data = find_monomorphic_locus(gt_df)
        monomorphic_df = pd.DataFrame({"locus": list(monomorphic_data.keys()), "gt": list(monomorphic_data.values())})

        # 去数据库比对
        whole_library_search(monomorphic_df, pdf)

        parent_folders = [folder for folder in os.listdir(config["output"]) if os.path.isdir(f"{config['output']}/{folder}")]
        if len(parent_folders) <= 0:
            raise ValueError("Can't find parents")

        # 删除单态标记
        gdf = gt_df[~gt_df["locus"].isin(list(monomorphic_data.keys()))]

        for folder in parent_folders:
            format_chr_data(folder, pdf, gdf)

            background_reply_evaluation(
                f"{config['output']}/{folder}",
                config["locus_file"]
            )

            exchange_report(config['output'], folder)

    except Exception as e:
        traceback.print_exc(e)
