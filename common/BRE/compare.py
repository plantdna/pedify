import pandas as pd
import os

from .base import append_write_csv_by_df


def compare_genotype(mother, children):
    '''
    群体的基因型数据与受体比较:
    群体: 却失        结果: 0
    群体: AA 受体: AA 结果: 1
    群体: AA 受体: AB 结果: 2
    群体: AA 受体: BB 结果: 3
    '''
    if not pd.isna(children) and len(children)>=2:
        if children in ['NoCall', '---']:
            return 0
        else:
            mother_allele = mother.split('/')
            children_allele = children.split('/')
            if mother == children:
                return 1
            elif children_allele[0] in mother_allele or children_allele[-1] in mother_allele:
                return 2
            else:
                return 3
    else:
        return 0

def compare(input_path, output_path, code):
    '''
    与受体数据进行比对
    '''
    if os.path.exists(output_path):
        os.system('rm -rf {}'.format(output_path))

    for chr_id in range(1, 11):
        chr_df = pd.read_csv('{}/chr{}.csv'.format(input_path, chr_id))
        chr_probesets = chr_df.columns.values.tolist()[1:]
        chr_df.set_index('call_code', inplace=True)
        group_chr_df_copy = chr_df.copy()

        for index, row in chr_df.iterrows():
            for probeset in chr_probesets:
                group_chr_df_copy.at[index, probeset] = compare_genotype(chr_df.loc[code, probeset], row[probeset])
        group_chr_df_copy.reset_index(inplace=True)

        append_write_csv_by_df(group_chr_df_copy, output_path, 'chr{}'.format(chr_id))

