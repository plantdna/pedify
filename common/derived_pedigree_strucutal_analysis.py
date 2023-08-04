import pandas as pd
from itertools import chain
import math
import traceback
from .snp_compare import snp_compare


config = {}

def get_top_key_value(data):
    value_counts = pd.value_counts(data).to_dict()
    value_keys = list(value_counts.keys())
    return (value_keys[0], value_counts[value_keys[0]])


def denoise(window, data):
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
    if c == 1:
        return n
    else:
        return ''


def get_certain_source(data):
    s = list(set([v for v in data if v != '']))
    if len(s) == 1:
        return s[0]
    elif len(s) == 0:
        return 'Unknow'
    else:
        return ''



def smooth(previous, next, data):
    data_counts = pd.value_counts(data).to_dict()
    keys = list(data_counts.keys())
    
    if len(data) <= config["window_size"]:
        return [previous]*len(data)
    else:
        if len(keys) > 1:
            if keys[0] == previous:
                return [previous]*len(data)
            elif keys[0] == next:
                return [next]*len(data)
            else:
                return [previous]*len(data)
        else:
            return data


def derived_pedigree_strucutal_analysis(params):
    try:
        global config
        config = params

        target_df = pd.read_csv(params["param"]['ancestor'], sep=params["sep"], low_memory=False)
        locus_df = pd.read_csv(params['locus_file'], sep=params["sep"], low_memory=False)

        for key, file_path in params["param"]["derivedlines"].items():
            df = locus_df[['locus', 'chr_id', 'position', 'HTP']].copy()
            item_df = pd.read_csv(file_path, sep=params["sep"], low_memory=False)
            cdf = pd.merge(item_df, target_df, on='locus')

            cdf[params["param"]["ancestorName"]] = cdf.apply(lambda x: format_source(
                snp_compare(x['gt_x'], x['gt_y']), params["param"]["ancestorName"]), axis=1)
            df = pd.merge(
                df, cdf[['locus', params["param"]["ancestorName"]]], on='locus', how='left')

            df.sort_values(by=["chr_id", "position"], ascending=[True, True], inplace=True)

            df['source'] = df.apply(
                lambda x: get_certain_source(x.values.tolist()[4:]), axis=1)

            for index, row in df.iterrows():
                if row['source'] == '':
                    if index == 0:
                        df.at[index, 'source'] = df.iloc[index+1]['source'] 
                    elif index == len(df)-1:
                        df.at[index, 'source'] = df.iloc[index-1]['source']
                    else:
                        sources = row.values.tolist()[4:-1:]
                        if sources.count('') == len(params["param"]['ancestors']):
                            df.at[index, 'source'] = df.iloc[index-1]['source']
                        elif len(sources) == len(params["param"]['ancestors']):
                            df.at[index, 'source'] = df.iloc[index-1]['source']
                        elif df.iloc[index-1]['source'] == df.iloc[index+1]['source'] and df.iloc[index-1]['source'] != '':
                            df.at[index, 'source'] = df.iloc[index-1]['source']
                        elif df.iloc[index-1]['source'] in sources:
                            df.at[index, 'source'] = df.iloc[index-1]['source']

                        if params["HTP"]:
                            htp_df = df[df['HTP'] == row['HTP']]
                            range_data = [v for v in list(chain.from_iterable(
                                htp_df[list(params["param"]['ancestors'].keys())].values.tolist())) if v != '']
                            if len(range_data) > 0:
                                data_count = pd.value_counts(range_data).to_dict()
                                data_keys = list(data_count.keys())
                                top_key = data_keys[0]
                                df.at[index, 'source'] = top_key


            if params["HTP"]:
                htp_denoise = []
                htp_smooth = []

                htps = list(set(df['HTP'].tolist()))
                htps.sort()

                for index, value in enumerate(htps):
                    htp_df = df[df['HTP']==value]

                    hdd = denoise(params["window_size"], htp_df['source'].tolist())
                    htp_denoise += hdd

                df['HTP_Denoise'] = htp_denoise

                for index, value in enumerate(htps):
                    htp_df = df[df['HTP']==value]

                    if index == 0:
                        next_htp_df = df[df['HTP']==htps[1]]
                        data_counts = next_htp_df['HTP_Denoise'].value_counts().to_dict()
                        sd = smooth(list(data_counts.keys())[0], list(data_counts.keys())[0], htp_df['HTP_Denoise'].tolist())

                    elif index == len(htps)-1:
                        previous_htp_df = df[df['HTP']==htps[len(htps)-2]]
                        data_counts = previous_htp_df['HTP_Denoise'].value_counts().to_dict()
                        sd = smooth(list(data_counts.keys())[0], list(data_counts.keys())[0], htp_df['HTP_Denoise'].tolist())

                    else:
                        previous_htp_df = df[df['HTP']==htps[index-1]]
                        next_htp_df = df[df['HTP']==htps[index+1]]
                        provious_counts = previous_htp_df['HTP_Denoise'].value_counts().to_dict()
                        next_counts = next_htp_df['HTP_Denoise'].value_counts().to_dict()
                        sd = smooth(list(provious_counts.keys())[0], list(next_counts.keys())[0], htp_df['HTP_Denoise'].tolist())

                    htp_smooth += sd

                df['HTP_Smooth'] = htp_smooth


                htp_source = []

                for chr_id in range(1, 11):
                    chr_htp_values = []
                    cdf = df[df['chr_id']==chr_id].copy()
                    chr_htp_groups = cdf.groupby('HTP')

                    for chr_htp_name, chr_htp_df in chr_htp_groups:
                        chr_htp_values.append(list(chr_htp_df['HTP_Smooth'].value_counts().to_dict().keys())[0])

                    chr_htp_value_denoise = denoise(10, chr_htp_values)
                    htp_source+=chr_htp_value_denoise

                htp_smooth_2 = []
                for index, value in enumerate(htps):
                    htp_df = df[df['HTP']==value]
                    htp_smooth_2+=[htp_source[index]]*len(htp_df)


                df['HTP_Source'] = htp_smooth_2

                df.to_csv(
                    f"{params['output']}/{key}.txt",
                    index=False,
                    encoding="utf-8_sig",
                    sep=params["sep"]
                )

    except Exception as e:
        traceback.print_exc(e)

