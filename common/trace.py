import traceback
import multiprocessing
import os

from pandas import read_csv, merge, DataFrame
from .snp_compare import snp_compare

output = ""

def write_result(df):
    try:
        if os.path.exists(f'{output}/trace.txt'):
            df.to_csv(f'{output}/trace.txt', index=False, encoding="utf-8_sig", mode="a+", header=None, sep="\t")
        else:
            df.to_csv(f'{output}/trace.txt', index=False, encoding="utf-8_sig", sep="\t")

    except Exception as e:
        traceback.print_exc(e)
    

def chunk_compare(df, pdf):
    try:
        results = []
        samples = df.columns.tolist()[2:]
        for sample in samples:
            item_df = merge(
                df[["probeset_id","gt",sample]],
                pdf[["probeset_id", "Chr_id", "Start", "HTP"]],
                on="probeset_id",
                how="left"
            )

            item_df["CR"] = item_df.apply(lambda  x: snp_compare(x["gt"], x[sample]), axis=1)

            res_list = item_df["CR"].tolist()
            snp_gs = round(res_list.count(1) / (len(res_list) - res_list.count(0)), 3) if (len(res_list) - res_list.count(0)) != 0 else 0

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
        
        result_df = DataFrame(
            results, 
            columns=[
                'CODE', 
                'MISS_SNP_NUMBER', 
                'EQUAL_SNP_NUMBER', 
                'SIMILARITY_SNP_NUMBER',
                'DIFF_SNP_NUMBER',
                'SNP_GENETIC_SIMILARITY',
                'HTP_GENETIC_SIMILARITY',
            ]
        )
        return result_df

    except Exception as e:
        traceback.print_exc(e)


def whole_library_search(cpus, df, library_path, pdf):
    try:
        pool = multiprocessing.Pool(cpus)
        file_reader = read_csv(library_path, sep='\t', iterator=True) 
        loop= True
        while loop:
            try:
                chunk = file_reader.get_chunk(50)
                chunk_t = chunk.T.copy()
                chunk_t.reset_index(inplace=True)
                chunk_t.columns = chunk_t.iloc[0].tolist()
                chunk_t.drop(0, axis=0, inplace=True)
                chunk_t.rename(columns={"call_code":"probeset_id"}, inplace=True)

                ndf = merge(df[["probeset_id", "gt"]], chunk_t, on="probeset_id", how="left")

                pool.apply_async(
                    func=chunk_compare,
                    args=(
                        ndf,
                        pdf,
                    ), 
                    callback=write_result
                )

            except StopIteration:
                loop = False
        
        pool.close()
        pool.join()

    except Exception as e:
        traceback.print_exc(e)



def get_inbredX(target_df, ancestors, miss_string, pdf):
    nuclear_snps = pdf["probeset_id"].tolist()

    merge_df = target_df[["probeset_id", "gt"]].copy()
    inbredx_snps = []

    for fp in ancestors:
        ancestor_df = read_csv(fp)
        merge_df = merge(merge_df, ancestor_df[["probeset_id", "gt"]])

    for index, row in merge_df.iterrows():
        gts = row.values.tolist()
        if gts[0] not in miss_string and gts[1] not in gts[2:] and row["probeset_id"] in nuclear_snps:
            inbredx_snps.append(row["probeset_id"])

    inbredx_df = target_df[target_df["probeset_id"].isin(inbredx_snps)][["probeset_id", "gt"]].copy()

    return inbredx_df



def trace(params:dict, logging):
    try:
        
        global output
        output = params["output"]
        if os.path.exists(f'{output}/trace.txt'):
            os.system(f'rm {output}/trace.txt')

        target_df = read_csv(params["target"])
        pdf = read_csv(params["probeset"])

        if len(params["ancestors"]) == 0:
            logging.info("系谱完全未知, 执行全库比对程序")
            # 完全未知, 执行全库比对程序
            whole_library_search(
                cpus=params["cpu"],
                df=target_df,
                library_path=params["dataset_path"],
                pdf=pdf
            )
        else:
            logging.info("部分未知，执行InbredX比对程序")
            # 部分未知，执行InbredX比对程序
            inbredX_df = get_inbredX(target_df, params["ancestors"], params["miss_string"], pdf)
            whole_library_search(
                params["cpu"],
                inbredX_df,
                params["dataset_path"],
                pdf
            )

    except Exception as e:
        traceback.print_exc(e)