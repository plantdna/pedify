import os
from .base import smoothing


def denoise(input_path, output_path):
    '''
    对比对后的数据进行降噪处理
    '''
    if os.path.exists(output_path):
        os.system('rm -rf {}'.format(output_path))

    for chr_id in range(1, 11):
        if os.path.exists('{}/chr{}.csv'.format(input_path, chr_id)):
            smoothing('{}/chr{}.csv'.format(input_path, chr_id), output_path, chr_id)

