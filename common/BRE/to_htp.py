import pandas as pd
import os

from .base import snp2htp

def to_htp(input_path, output_path, htp_df):
    if os.path.exists(output_path):
        os.system('rm -rf {}'.format(output_path))

    for chr_id in range(1, 11):
        if os.path.exists('{}/chr{}.csv'.format(input_path, chr_id)):
            chr_htp_df = htp_df[htp_df['chr_id'].isin([chr_id, str(chr_id)])]
            snp2htp('{}/chr{}.csv'.format(input_path, chr_id), output_path, chr_id, chr_htp_df)

