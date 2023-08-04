import pandas as pd
from itertools import chain
import math
import json
from .snp_compare import snp_compare


def get_top_key_value(data):
    value_counts = pd.value_counts(data).to_dict()
    value_keys = list(value_counts.keys())
    return (value_keys[0], value_counts[value_keys[0]])

def denoise(window, data):
    """根据滑动窗口进行降噪

    Args:
        window (int): 窗口大小
        data (list): 需要降噪的数据
    """
    # * 先将缺失数据填充
    c_data = data.copy()
    
    if '' in c_data and c_data.count('') == len(c_data):
        return c_data
    else:
        for i in range(len(c_data)):
            if c_data[i] == '':
                if i == 0:
                    j = i
                    while c_data[j] == '':
                        j+=1
                    c_data[i] = c_data[j]
                else:
                    j = i
                    while c_data[j] == '':
                        j-=1
                    c_data[i] = c_data[j]

        new_data = []
        if len(c_data) < window:
            (top_key, top_value) = get_top_key_value(c_data)
            new_data+=[top_key]*len(c_data)
        else:
            slices = math.ceil(len(c_data)/window)
            front_key = None
            for i in range(slices):
                slice_data = c_data[i*window:(i+1)*window]
                if len(slice_data) < window:
                    new_data += [front_key]*len(slice_data)
                else:
                    end_data = c_data[(i+1)*window:(i+2)*window] if i<slices-1 else c_data[i*window:]

                    (top_key, top_value) = get_top_key_value(slice_data)
                    top_ratio = round(top_value/len(slice_data), 2)

                    (end_key, end_value) = get_top_key_value(end_data)
                    end_ratio = round(end_value/len(end_data), 2)

                    if front_key:
                        if front_key == end_key:
                            new_data += [front_key]*len(slice_data)
                        elif front_key == top_key:
                            new_data += [front_key]*len(slice_data)
                        elif top_ratio >= 0.5:
                            new_data += [top_key]*len(slice_data)
                            front_key = top_key
                        else:
                            new_data += [front_key]*len(slice_data)

                    else:
                        front_key = top_key
                        new_data += [top_key]*len(slice_data)
        return new_data



def format_source(c, n):
    '''
    格式化溯源结果
    '''
    if c == 1:
        return n
    else:
        return ''


def get_certain_source(data):
    '''
    确定唯一的来源
    '''
    s = list(set([v for v in data if v != '']))
    if len(s) == 1:
        return s[0]
    else:
        return ''



def smooth(previous, next, data):
    '''
    HTP内部平滑，检测内部的异常小片段
    '''
    data_counts = pd.value_counts(data).to_dict()
    keys = list(data_counts.keys())
    
    if len(data) <= 10:
        # * 小HTP
        return [previous]*len(data)
    else:
        # * 大HTP
        if len(keys) > 1:
            # * HTP内部发生了交换，根据前后检测是否是真正的交换
            if keys[0] == previous:
                return [previous]*len(data)
            elif keys[0] == next:
                return [next]*len(data)
            else:
                return [previous]*len(data)
        else:
            # * HTP内部不发生交换
            return data



def reconstruction(params):
    '''
    系谱重建
    '''
    target_df = pd.read_csv(params['target'])

    probeset_df = pd.read_csv(params['probeset'], low_memory=False)

    chr_df = probeset_df[['probeset_id', 'Chr_id', 'Start', 'HTP']].copy()

    # * 第一步进行数据比对
    for key, value in params['ancestors'].items():
        item_df = pd.read_csv(value)
        cdf = pd.merge(target_df, item_df, on='probeset_id')
        cdf[key] = cdf.apply(lambda x: format_source(
            snp_compare(x['gt_x'], x['gt_y']), key), axis=1)
        chr_df = pd.merge(
            chr_df, cdf[['probeset_id', key]], on='probeset_id', how='left')


    # * 确定唯一的来源
    chr_df['source'] = chr_df.apply(
        lambda x: get_certain_source(x.values.tolist()[-len(params['ancestors']):]), axis=1)


    # * SNP未知来源预测
    for index, row in chr_df.iterrows():
        if row['source'] == '':
            # * 只处理未知来源
            if index == 0:
                # * 第一个
                chr_df.at[index, 'source'] = chr_df.iloc[index+1]['source']
            elif index == len(chr_df)-1:
                # * 最后一个
                chr_df.at[index, 'source'] = chr_df.iloc[index-1]['source']
            else:
                sources = row.values.tolist()[4:-1:]
                if sources.count('') == len(params['ancestors']):
                    # * 未知来源
                    chr_df.at[index, 'source'] = chr_df.iloc[index-1]['source']
                elif len(sources) == len(params['ancestors']):
                    # * 共同来源
                    chr_df.at[index, 'source'] = chr_df.iloc[index-1]['source']
                elif chr_df.iloc[index-1]['source'] == chr_df.iloc[index+1]['source'] and chr_df.iloc[index-1]['source'] != '':
                    # * 前后相同确定中间
                    chr_df.at[index, 'source'] = chr_df.iloc[index-1]['source']
                elif chr_df.iloc[index-1]['source'] in sources:
                    # * 前后不同，前在其中
                    chr_df.at[index, 'source'] = chr_df.iloc[index-1]['source']
                else:
                    # * 前后不同，前不在其中
                    htp_df = chr_df[chr_df['HTP'] == row['HTP']]
                    range_data = [v for v in list(chain.from_iterable(
                        htp_df[list(params['ancestors'].keys())].values.tolist())) if v != '']
                    if len(range_data) > 0:
                        data_count = pd.value_counts(range_data).to_dict()
                        data_keys = list(data_count.keys())
                        top_key = data_keys[0]
                        chr_df.at[index, 'source'] = top_key

    # * HTP内部降噪
    htp_denoise = []
    # * HTP内部平滑
    htp_smooth = []

    htps = list(set(chr_df['HTP'].tolist()))
    htps.sort()

    # * HTP降噪
    for index, value in enumerate(htps):
        htp_df = chr_df[chr_df['HTP']==value]

        hdd = denoise(10, htp_df['source'].tolist())
        htp_denoise += hdd

    chr_df['HTP_Denoise'] = htp_denoise


    # * HTP内部平滑
    for index, value in enumerate(htps):
        htp_df = chr_df[chr_df['HTP']==value]

        if index == 0:
            next_htp_df = chr_df[chr_df['HTP']==htps[1]]
            data_counts = next_htp_df['HTP_Denoise'].value_counts().to_dict()
            sd = smooth(list(data_counts.keys())[0], list(data_counts.keys())[0], htp_df['HTP_Denoise'].tolist())

        elif index == len(htps)-1:
            previous_htp_df = chr_df[chr_df['HTP']==htps[len(htps)-2]]
            data_counts = previous_htp_df['HTP_Denoise'].value_counts().to_dict()
            sd = smooth(list(data_counts.keys())[0], list(data_counts.keys())[0], htp_df['HTP_Denoise'].tolist())

        else:
            previous_htp_df = chr_df[chr_df['HTP']==htps[index-1]]
            next_htp_df = chr_df[chr_df['HTP']==htps[index+1]]
            provious_counts = previous_htp_df['HTP_Denoise'].value_counts().to_dict()
            next_counts = next_htp_df['HTP_Denoise'].value_counts().to_dict()
            sd = smooth(list(provious_counts.keys())[0], list(next_counts.keys())[0], htp_df['HTP_Denoise'].tolist())

        htp_smooth += sd

    chr_df['HTP_Smooth'] = htp_smooth


    # * 分染色体进行 HTP间降噪
    htp_source = []

    for chr_id in range(1, 11):
        chr_htp_values = []
        cdf = chr_df[chr_df['Chr_id']==chr_id].copy()
        chr_htp_groups = cdf.groupby('HTP')

        for chr_htp_name, chr_htp_df in chr_htp_groups:
            chr_htp_values.append(list(chr_htp_df['HTP_Smooth'].value_counts().to_dict().keys())[0])

        chr_htp_value_denoise = denoise(10, chr_htp_values)
        htp_source+=chr_htp_value_denoise

    htp_smooth_2 = []
    for index, value in enumerate(htps):
        htp_df = chr_df[chr_df['HTP']==value]
        htp_smooth_2+=[htp_source[index]]*len(htp_df)

                

    chr_df['HTP_Source'] = htp_smooth_2

    return chr_df

