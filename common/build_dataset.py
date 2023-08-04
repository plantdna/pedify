# build dataset
import os
import pandas as pd
import traceback

def assess_locus(data):
    maf = 0
    pic = 0
    alleles = []
    for v in data:
        if '/' in v:
            gt = v.split('/')
            alleles+=gt

    allele = list(set(alleles))

    if len(allele) == 2:
        af = [alleles.count(allele[0])/len(alleles), alleles.count(allele[1])/len(alleles)]
        maf = round(min(af), 3)
        pic = round(1-(pow(maf, 2) + pow(1-maf, 2)), 3)
    else:
        # 单态
        maf = 0
        pic = 0

    return (maf,pic)
    

def build_dataset(param):
    try:
        locus_file = param["locus_file"]
        sep = param["sep"]
        foldepath = param["param"]["genotypes"]
        pic_limit = param["param"]["pic_limit"]
        maf_limit = param["param"]["maf_limit"]
        output = param["output"]

        locus_df = pd.read_csv(locus_file, sep=sep)
        common_locus = locus_df["locus"].to_list()

        txtfiles = [filename for filename in os.listdir(foldepath) if filename[-3:] == "txt"]
        database = pd.DataFrame()
        for tf in txtfiles:
            name = tf.split('.')[0]
            df = pd.read_csv(f"{foldepath}/{tf}", sep=sep)
            df.rename(columns={"gt":name}, inplace=True)
            if len(database) > 0:
                database = pd.merge(database, df[df["locus"].isin(common_locus)], on="locus", how="left")
            else:
                database = df[df["locus"].isin(common_locus)]

        database.set_index("locus", inplace=True)
        dt = database.T
        dt.reset_index(inplace=True)
        dt.rename(columns={"index":"Name"}, inplace=True)

        #  get maf and pic of locus
        locus = dt.columns.tolist()[1:]
        delete_locus = []
        for loci in locus:
            maf, pic = assess_locus(dt[loci].tolist())
            if maf <= maf_limit or pic <= pic_limit:
                delete_locus.append(loci)

        # filter locus
        dt.drop(delete_locus, axis=1, inplace=True)
        dt.to_csv(f"{output}/dataset.txt", index=False, encoding="utf-8_sig", sep="\t")

        print(f"Remove substandard loci: {len(delete_locus)}")

    except Exception as e:
        traceback.print_exc(e)
