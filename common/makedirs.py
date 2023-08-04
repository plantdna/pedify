import os
import traceback

def makedirs(dirpath:str) -> dict:
    try:
        if os.path.exists(dirpath):
            os.system(f"rm -rf {dirpath}")

        os.makedirs(dirpath)
    except Exception as e:
        traceback.print_exc(e)