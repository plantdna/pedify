import traceback
import multiprocessing
import subprocess
import math

from pandas import read_csv, merge, DataFrame, concat
from snp_compare import snp_compare


def write_result(df):
    df.to_csv('../files/test.csv', index=False, encoding="utf-8_sig", mode="a+", header=None)
    

def chunk_compare(df, pdf):
    try:
        results = []
        samples = df.columns.tolist()[2:]
        print(samples)
        for sample in samples:
            item_df = merge(
                pdf[["probeset_id", "Chr_id", "Start", "HTP"]],
                df[["probeset_id","gt",sample]],
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
        print(result_df)
        return result_df

    except Exception as e:
        traceback.print_exc(e)


if __name__ == "__main__":
    pool = multiprocessing.Pool(5)

    df = read_csv("../files/test/u8112.csv")
    pdf = read_csv("../dataset/NUCLEAR_SNP_PROBESET.csv")

    library_path = "../dataset/Inbred_Lines.txt"

    print("数据加载完毕，准备进入多进程")

    file_reader = read_csv(library_path, sep='\t', iterator=True, nrows=1000) 
    loop= True
    while loop:
        try:
            print('======>')
            chunk = file_reader.get_chunk(20)
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
