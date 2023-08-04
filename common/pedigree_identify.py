import traceback
import multiprocessing
import os

from pandas import read_csv, merge, DataFrame
from .snp_compare import snp_compare

config = {}

def write_result(df):
    try:
        if os.path.exists(f'{config["output"]}/pedigree_identify.txt'):
            df.to_csv(f'{config["output"]}/pedigree_identify.txt', index=False, encoding="utf-8_sig", mode="a+", header=None, sep=config["sep"])
        else:
            df.to_csv(f'{config["output"]}/pedigree_identify.txt', index=False, encoding="utf-8_sig", sep=config["sep"])

    except Exception as e:
        traceback.print_exc(e)
    

def chunk_compare(df, pdf, useHTP):
    try:
        results = []
        samples = df.columns.tolist()[2:]
        for sample in samples:
            item_df = merge(
                df[["locus","gt",sample]],
                pdf[["locus", "chr_id", "position", "HTP"]],
                on="locus",
                how="left"
            )

            item_df.dropna(axis=0, subset=["gt", sample], inplace=True)
            item_df["CR"] = item_df.apply(lambda x: snp_compare(x["gt"], x[sample]), axis=1)

            res_list = item_df["CR"].tolist()
            snp_gs = round(res_list.count(1) / (len(res_list) - res_list.count(0)), 3) if (len(res_list) - res_list.count(0)) != 0 else 0

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
        file_reader = read_csv(config["param"]["dataset"], sep=config["sep"], iterator=True) 
        loop= True
        while loop:
            try:
                chunk = file_reader.get_chunk(50)
                chunk_t = chunk.T.copy()
                chunk_t.reset_index(inplace=True)
                chunk_t.columns = chunk_t.iloc[0].tolist()
                chunk_t.drop(0, axis=0, inplace=True)
                chunk_t.rename(columns={"Name":"locus"}, inplace=True)

                ndf = merge(df[["locus", "gt"]], chunk_t, on="locus", how="left")

                pool.apply_async(
                    func=chunk_compare,
                    args=(
                        ndf,
                        pdf,
                        config["HTP"],
                    ), 
                    callback=write_result
                )

            except StopIteration:
                loop = False
        
        pool.close()
        pool.join()

    except Exception as e:
        traceback.print_exc(e)



def find_inbredX(target_df, ancestors, miss_string, sep, pdf):

    nuclear_snps = pdf["locus"].tolist()

    merge_df = target_df[["locus", "gt"]].copy()
    inbredx_snps = []

    for fp in ancestors:
        ancestor_df = read_csv(fp, sep=sep)
        merge_df = merge(merge_df, ancestor_df[["locus", "gt"]], on="locus")


    for index, row in merge_df.iterrows():
        gts = row.values.tolist()[1:]
        if gts[0] != miss_string and gts[0] not in gts[1:] and row["locus"] in nuclear_snps:
            inbredx_snps.append(row["locus"])

    inbredx_df = target_df[target_df["locus"].isin(inbredx_snps)][["locus", "gt"]].copy()

    return inbredx_df


def find_LCF(df, pdf):
    ndf = merge(df, pdf[["locus", "chr_id", "position"]], on="locus")
    ndf.sort_values(by=["chr_id", "position"], ascending=[True, True], inplace=True)

    max_lcfs = []

    for chr_id in range(1, 11):
        chr_df = ndf[ndf["chr_id"] == chr_id]
        chr_max_lcfs = []
        lcfs = []
        for index, row in chr_df.iterrows():
            if index > 0:
                if row["position"] - ndf.at[index-1, "position"] <= 5000000:
                    lcfs.append(row["locus"])
                else:
                    if len(chr_max_lcfs) < len(lcfs):
                        chr_max_lcfs = lcfs
                        lcfs = []

        max_lcfs+=chr_max_lcfs

    return df[df["locus"].isin(max_lcfs)]


def pedigree_identify(params):
    try:
        global config
        config = params

        if os.path.exists(f'{config["output"]}/pedigree_identify.txt'):
            os.system(f'rm {config["output"]}/pedigree_identify.txt')

        target_df = read_csv(params["param"]["target"], sep=params["sep"], low_memory=False)
        pdf = read_csv(params["locus_file"], sep=params["sep"], low_memory=False)

        if len(params["param"]["ancestors"]) == 0:
            whole_library_search(
                df=target_df,
                pdf=pdf
            )
        else:
            inbredX_df = find_inbredX(target_df, params["param"]["ancestors"], params["miss"], params["sep"], pdf)
            if "algorithm" in params["param"] and params["param"]["algorithm"] == "LPI":
                lcf_inbredX_df = find_LCF(inbredX_df, pdf)

                whole_library_search(
                    df=lcf_inbredX_df,
                    pdf=pdf
                )
            else:
                whole_library_search(
                    df=inbredX_df,
                    pdf=pdf
                )

        df = read_csv(f'{config["output"]}/pedigree_identify.txt',  sep=config["sep"])
        if config["HTP"]:
            df.sort_values(by="HTP_GENETIC_SIMILARITY", ascending=False, inplace=True)
        else:
            df.sort_values(by="SNP_GENETIC_SIMILARITY", ascending=False, inplace=True)

        df.to_csv(f'{config["output"]}/pedigree_identify_all.txt', index=False, encoding="utf-8_sig", sep=config["sep"])

        if config["HTP"]:
            df.sort_values(by="HTP_GENETIC_SIMILARITY", ascending=False, inplace=True)
            if len(df[df["HTP_GENETIC_SIMILARITY"]>0.9])>10:
                df[df["HTP_GENETIC_SIMILARITY"]>0.9].head(10).to_csv(f'{config["output"]}/pedigree_identify.txt', index=False, encoding="utf-8_sig", sep=config["sep"])

            else:
                df.head(10).to_csv(f'{config["output"]}/pedigree_identify.txt', index=False, encoding="utf-8_sig", sep=config["sep"])

        else:
            df.sort_values(by="SNP_GENETIC_SIMILARITY", ascending=False, inplace=True)
            if len(df[df["SNP_GENETIC_SIMILARITY"]>0.9])>10:
                df[df["SNP_GENETIC_SIMILARITY"]>0.9].head(10).to_csv(f'{config["output"]}/pedigree_identify.txt', index=False, encoding="utf-8_sig", sep=config["sep"])

            else:
                df.head(10).to_csv(f'{config["output"]}/pedigree_identify.txt', index=False, encoding="utf-8_sig", sep=config["sep"])


    except Exception as e:
        traceback.print_exc(e)