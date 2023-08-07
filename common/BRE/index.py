
import pandas as pd

from .compare import compare
from .denoise import denoise
from .to_htp import to_htp

def background_reply_evaluation(output, locus_file):
    chr_dir = f"{output}/chr_data"
    compare_dir = f"{output}/compare"
    denoise_snp_dir = f"{output}/snp_denoise"
    chr_htp_dir = f"{output}/htp"
    probeset_info_df = pd.read_csv(locus_file, sep="\t", low_memory=False)
    denoise_htp_dir = f"{output}/htp_denoise"

    probeset_info_df.set_index('locus', inplace=True)

    compare(
        chr_dir, 
        compare_dir,
        'gt'
    )
    denoise(compare_dir, denoise_snp_dir)
    to_htp(denoise_snp_dir, chr_htp_dir, probeset_info_df)
    denoise(chr_htp_dir, denoise_htp_dir)
