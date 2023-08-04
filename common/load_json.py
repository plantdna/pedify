# 读取json文件
import json
import os
import traceback

def load_json(filepath:str) -> dict:
    data = {}
    try:
        if os.path.exists(filepath):
            with open(filepath, 'r') as f:
                data = json.loads(f.read())
        else:
            raise Exception(f"this file path '{filepath}' is not exists")
        return data
    except Exception as e:
        traceback.print_exc(e)
