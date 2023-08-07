import pandas as pd
import os
import json
import traceback


def load_json(file_path:str):
    '''
    load json file
    '''
    f = open(file_path,'r',encoding='utf-8')
    json_dict = json.load(f)
    f.close()
    return json_dict


def format_opts(opts:list):
    '''
    format options
    '''
    opt_dict = {}
    for item in opts:
        opt_dict[item[0]] = item[1]
    return opt_dict


def append_write_csv_by_df(df, dir_path:str, file_name:str):
    '''
    以追加的方式将DataFrame写成CSV文件
    '''
    try:
        file_path='{}/{}.csv'.format(dir_path, file_name)
        if os.path.exists(dir_path):
            if os.path.exists(file_path):
                df.to_csv(file_path, mode='a',index=False, header=False, encoding='utf-8_sig')
            else:
                df.to_csv(file_path, index=False, encoding='utf-8_sig')            
        else:
            os.makedirs(dir_path)
            df.to_csv(file_path, index=False, encoding='utf-8_sig')

        return
    except Exception as e:
        traceback.format_exc(e)


def smoothing(file_path, output_dir, chr_id):
    '''
    比对结果进行平滑
    '''
    try:
        chr_df = pd.read_csv(file_path)    
        new_data = []
        step_length=5
        
        for index, row in chr_df.iterrows():
            new_row = [row['call_code']]
            row_list = row.tolist()[1:]
            row_list_len = len(row_list)
            
            previous_value = 0
            next_value = 0
        
            for i in range(0, row_list_len, step_length):
                block_values = row_list[i:i+step_length]
                block_length = len(block_values)
                block_counts = pd.value_counts(block_values).to_dict()
                block_keys = list(block_counts.keys())
                top_key = block_keys[0]

                next_block_values = row_list[i+step_length: i+step_length*2]
                if len(next_block_values) > 0:
                    next_block_counts = pd.value_counts(next_block_values).to_dict()
                    next_value = list(next_block_counts.keys())[0]

                if top_key == 0:
                    new_row+=[previous_value]*block_length
                elif previous_value == next_value:
                    new_row+=[previous_value]*block_length
                elif block_counts[top_key]/block_length >= 0.5:
                    previous_value = top_key
                    new_row+=[top_key]*block_length
                else:
                    new_row+=[previous_value]*block_length

            new_data.append(new_row)
        new_df = pd.DataFrame(new_data, columns=chr_df.columns.tolist())
        append_write_csv_by_df(new_df, output_dir, 'chr{}'.format(chr_id))
    except Exception as e:
        traceback.format_exc(e)


def snp2htp(file_path:str, output_dir:str, chr_id, chr_htp_df):
    '''
    将平滑后的SNP数据转成HTP
    '''
    try:
        chr_df = pd.read_csv(file_path)
        chr_probesets = chr_df.columns.values.tolist()[1:]
        chr_htp_groups = chr_htp_df.groupby(by="HTP")

        new_htp_df = chr_df[['call_code']].copy()

        for htp_name, htp_df in chr_htp_groups:
            htp_probesets = htp_df.index.tolist()
            common_htp_probesets = [probeset_id for probeset_id in htp_probesets if probeset_id in chr_probesets]
            htp_col_values = []

            if len(common_htp_probesets) > 0:
                chr_htp_df = chr_df[common_htp_probesets]
                for index, row in chr_htp_df.iterrows():
                    row_values = row.tolist()
                    if 1 in row_values and row_values.count(1)/len(row_values) >= 0.6:
                        htp_col_values.append(1)
                    elif 3 in row_values and row_values.count(3)/len(row_values) >= 0.7:
                        htp_col_values.append(3)
                    else:
                        htp_col_values.append(2)
            
                col_df = pd.DataFrame({htp_name : htp_col_values})
                new_htp_df = pd.concat([new_htp_df, col_df], axis=1)

        append_write_csv_by_df(new_htp_df, output_dir, 'chr{}'.format(chr_id))

    except Exception as e:
        print(e)
        traceback.format_exc(e)


def color_excel(file_path:str, output_dir:str, file_name:str, skip_cols:int):
    '''
    生成带有背景色的Excel
    '''
    try:
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        df = pd.read_csv(file_path)
        columns = df.columns.values.tolist()

        writer = pd.ExcelWriter('{}/{}.xlsx'.format(output_dir, file_name), engine='xlsxwriter')

        df.to_excel(writer, sheet_name='Sheet1', header=True, index=False) 

        workbook  = writer.book
        worksheet = writer.sheets['Sheet1']

        # green
        format_green = workbook.add_format({'bg_color': '#15a43f'})
        # yellow
        format_yellow = workbook.add_format({'bg_color': '#feb406'})
        # red
        format_red = workbook.add_format({'bg_color': '#fb0006'})
        # black
        format_black = workbook.add_format({'bg_color': '#000000', 'font_color': '#ffffff'})
        
        
        # 给数据着色
        worksheet.conditional_format(0, skip_cols, len(df), len(columns),
                                    {'type':'cell',
                                    'criteria': '=',
                                    'value':  1,
                                    'format': format_green})
        # 给数据着色
        worksheet.conditional_format(0, skip_cols, len(df), len(columns),
                                    {'type':'cell',
                                    'criteria': '=',
                                    'value':  2,
                                    'format': format_yellow})
        # 给数据着色
        worksheet.conditional_format(0, skip_cols, len(df), len(columns),
                                    {'type':'cell',
                                    'criteria': '=',
                                    'value':  3,
                                    'format': format_red})
        
        worksheet.set_column(skip_cols, len(columns), 0.1) 
        worksheet.set_column(0, skip_cols, 15) 

        writer.close()

    except Exception as e:
        traceback.format_exc(e)


def all_chrs_htp_color_excel(worksheet, workbook, skip_cols, sheet_df):
    try:
        sheet_columns = sheet_df.columns.tolist()

        # green
        format_green = workbook.add_format({'bg_color': '#15a43f'})
        # yellow
        format_yellow = workbook.add_format({'bg_color': '#feb406'})
        # red
        format_red = workbook.add_format({'bg_color': '#fb0006'})
        
        
        # 给数据着色
        worksheet.conditional_format(0, skip_cols, len(sheet_df), len(sheet_columns),
                                    {'type':'cell',
                                    'criteria': '=',
                                    'value':  1,
                                    'format': format_green})
        # 给数据着色
        worksheet.conditional_format(0, skip_cols, len(sheet_df), len(sheet_columns),
                                    {'type':'cell',
                                    'criteria': '=',
                                    'value':  2,
                                    'format': format_yellow})
        # 给数据着色
        worksheet.conditional_format(0, skip_cols, len(sheet_df), len(sheet_columns),
                                    {'type':'cell',
                                    'criteria': '=',
                                    'value':  3,
                                    'format': format_red})
        
        for index, col in enumerate(sheet_columns):
            if index >= skip_cols:
                if 'call_code' in col:
                    worksheet.set_column(index, index, 15)
                else:
                    worksheet.set_column(index, index, 0.1)
            else:
                worksheet.set_column(0, skip_cols, 15)

        return worksheet
    except Exception as e:
        traceback.format_exc(e)

def statistics(snp_dir, htp_dir, output_dir):
    '''
    统计SNP与HTP标记在个染色体上的回复情况
    '''
    try:
        statistics_obj = {
            'call_code': [],
            'SNP_Num':[],
            'SNP_Miss_Num':[],
            'SNP_Rate':[],
            'HTP_Num':[],
            'HTP_Rate':[],
        }

        # Sheet2
        sheet2_data = [
            ['总标记数', 67984, ''], 
            ['总HTP数', 6164, ''], 
            ['跟踪的标记数', 0, ''], 
            ['跟踪的HTP数', 0, ''], 
            ['Chr','标记数','HTP数'],
        ]

        all_chrs_htp = pd.DataFrame()

        for chr_id in range(1, 11):
            statistics_obj['chr{}_SNP_Num'.format(chr_id)] = []    
            statistics_obj['chr{}_SNP_Rate'.format(chr_id)] = []    
            statistics_obj['chr{}_HTP_Num'.format(chr_id)] = []
            statistics_obj['chr{}_HTP_Rate'.format(chr_id)] = []

        for chr_id in range(1, 11):
            chr_snp_df = pd.read_csv('{}/chr{}.csv'.format(snp_dir, chr_id))
            chr_snp_df.rename(columns={'call_code': '{}_call_code'.format(chr_id)}, inplace=True)
            chr_htp_df = pd.read_csv('{}/chr{}.csv'.format(htp_dir, chr_id))
            chr_htp_df.rename(columns={'call_code': '{}_call_code'.format(chr_id)}, inplace=True)
            
            if len(all_chrs_htp) > 0:
                all_chrs_htp = pd.concat([all_chrs_htp, chr_htp_df], axis=1)
            else:
                all_chrs_htp = chr_htp_df

            sheet2_data[2][1]+=chr_snp_df.shape[1]-1
            sheet2_data[3][1]+=chr_htp_df.shape[1]-1
            sheet2_data.append(['chr{}'.format(chr_id), chr_snp_df.shape[1]-1, chr_htp_df.shape[1]-1])

            if len(statistics_obj['call_code']) == 0:
                statistics_obj['call_code'] = chr_snp_df['{}_call_code'.format(chr_id)].tolist()
                statistics_obj['SNP_Num'] = [0]*chr_snp_df.shape[0]
                statistics_obj['SNP_Miss_Num'] = [0]*chr_snp_df.shape[0]
                statistics_obj['SNP_Rate'] = [0]*chr_snp_df.shape[0]
                statistics_obj['HTP_Num'] = [0]*chr_snp_df.shape[0]
                statistics_obj['HTP_Rate'] = [0]*chr_snp_df.shape[0]
                
            df_len = len(chr_snp_df)
            
            for i in range(df_len):
                sam_snp_data = chr_snp_df.loc[i].tolist()[1:]        
                sam_htp_data = chr_htp_df.loc[i].tolist()[1:]
                
                snp_num = sam_snp_data.count(2)+sam_snp_data.count(3)
                htp_num = sam_htp_data.count(2)+sam_htp_data.count(3)
                
                statistics_obj['SNP_Num'][i]+=snp_num
                statistics_obj['SNP_Miss_Num'][i]+=sam_snp_data.count(0)
                statistics_obj['HTP_Num'][i]+=htp_num

                statistics_obj['chr{}_SNP_Num'.format(chr_id)].append(snp_num) 
                statistics_obj['chr{}_SNP_Rate'.format(chr_id)].append(round(snp_num/len(sam_snp_data), 3))
                statistics_obj['chr{}_HTP_Num'.format(chr_id)].append(htp_num)
                statistics_obj['chr{}_HTP_Rate'.format(chr_id)].append(round(htp_num/len(sam_htp_data),3))

        statistics_df = pd.DataFrame(statistics_obj)
        statistics_df['SNP_Rate'] = statistics_df[['SNP_Num']].apply(lambda x: round(x['SNP_Num']/sheet2_data[2][1], 3), axis=1)
        statistics_df['HTP_Rate'] = statistics_df[['HTP_Num']].apply(lambda x: round(x['HTP_Num']/sheet2_data[3][1], 3), axis=1)
        sheet2_df = pd.DataFrame(sheet2_data)
        
        statistics_df.to_csv('{}/数据统计.csv'.format(output_dir), index=False, encoding='utf-8_sig')
        sheet2_df.iloc[4:].to_csv('{}/跟踪标记统计.csv'.format(output_dir), index=False, header=None, encoding='utf-8_sig')
        # 对Sheet1进行拼接
        sheet1_df = pd.concat([statistics_df, all_chrs_htp], axis=1)
        # 对文件进行着色
        skip_cols = 47
        writer = pd.ExcelWriter('{}/{}.xlsx'.format(output_dir, '汇总表'), engine='xlsxwriter')
        sheet1_df.to_excel(writer, sheet_name='Sheet1', header=True, index=False)
        sheet2_df.to_excel(writer, sheet_name='Sheet2', header=False, index=False)

        workbook  = writer.book
        worksheet = writer.sheets['Sheet1']

        all_chrs_htp_color_excel(worksheet, workbook, skip_cols, sheet1_df)
                 
        writer.close()

    except Exception as e:
        traceback.format_exc(e)


def selected_sample(step7_dir, except_codes):
    try:
        statistics_df = pd.read_excel('{}/统计表.xlsx'.format(step7_dir), sheet_name='Sheet1')
        all_htp_num = statistics_df.shape[1] - 45 -10

        statistics_df.rename(
            columns={
                'SNP_Num':'Chr_HTP_Num',
                'SNP_Rate':'Chr_HTP_Rate',
            }, 
            inplace=True
        )
        statistics_df[['Chr_HTP_Num', 'Chr_HTP_Rate']] = 0

        # 全基因组挑选
        sort_by_htp_rate = statistics_df.sort_values(by='HTP_Rate')
        selected_df = sort_by_htp_rate[~sort_by_htp_rate['call_code'].isin(except_codes)].head(5).copy()
        selected_df.reset_index(inplace=True)

        for index, row in selected_df.iterrows():
            selected_df.loc[index, 'index'] = 'ALL_{}'.format(index+1)

        # 分染色体挑选
        for chr_id in range(1, 11):
            chr_sort_df = statistics_df.sort_values(by=['chr{}_HTP_Rate'.format(chr_id), "HTP_Rate"], ascending=[False, True])
            chr_selected_df = chr_sort_df[~chr_sort_df['call_code'].isin(except_codes)].head(5).copy()
            chr_selected_df.reset_index(inplace=True)
            
            for index, row in chr_selected_df.iterrows():
                chr_selected_df.loc[index, 'index'] = 'Chr{}_{}'.format(chr_id, index+1)
                chr_selected_df.loc[index, 'Chr_HTP_Num'] = row['chr{}_HTP_Num'.format(chr_id)]
                chr_selected_df.loc[index, 'Chr_HTP_Rate'] = row['chr{}_HTP_Rate'.format(chr_id)]
                chr_selected_df.loc[index, 'HTP_Num'] = row['HTP_Num'] - row['chr{}_HTP_Num'.format(chr_id)]
                chr_selected_df.loc[index, 'HTP_Rate'] = round((row['HTP_Num'] - row['chr{}_HTP_Num'.format(chr_id)])/all_htp_num, 3)

            selected_df = selected_df.append(chr_selected_df, ignore_index=True)

        drop_cols = [col for col in statistics_df.columns.tolist()[:46] if 'SNP' in col or 'chr' in col]
        selected_df.drop(drop_cols, inplace=True, axis=1)

        selected_df.rename(
            columns={
                'HTP_Num':'Other_HTP_Num',
                'HTP_Rate':'Other_HTP_Rate',
            }, 
            inplace=True
        )

        # 对文件进行着色
        skip_cols = 6
        writer = pd.ExcelWriter('{}/{}.xlsx'.format(step7_dir, '筛选结果'), engine='xlsxwriter')
        selected_df.to_excel(writer, sheet_name='Sheet1', header=True, index=False)

        workbook  = writer.book
        worksheet = writer.sheets['Sheet1']

        all_chrs_htp_color_excel(worksheet, workbook, skip_cols, selected_df)
                 
        writer.close()

    except Exception as e:
        traceback.format_exc(e)
