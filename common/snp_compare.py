def list_compare(lista:list, listb:list) -> int:
    '''
    计算两个列表相同元素个数
    '''
    num = 0
    compare = lista
    contrast = listb
    
    if len(lista) > len(listb):
        compare = listb
        contrast = lista
        
    for v in compare:
        if v in contrast:
            num+=1
    return num


def snp_compare(str1:str, str2:str) -> int:
    '''
    基因型比对
    0: 缺失
    1: 完全相同
    3: 完全不同
    2: 部分相同
    '''
    alleles_a = [v for v in str1.split('/') if v != ''] if '/' in str1 else []
    alleles_b = [v for v in str2.split('/') if v != ''] if '/' in str2 else []
    
    if len(alleles_a) == 2 and len(alleles_b) == 2:
        set_a = list(set(alleles_a))
        set_b = list(set(alleles_b))
        # 完全相同
        if len(set_a) == len(set_b) and len(set_a) == list_compare(set_a, set_b):
            return 1
        # 完全不同
        if list_compare(set_a, set_b) == 0:
            return 3
        # 部分相同
        if list_compare(set_a, set_b) == 1 and (len(set_a) == 1 or len(set_b) == 1):
            return 2
    else:
        return 0
