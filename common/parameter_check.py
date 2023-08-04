# 参数校验
import os
import traceback
from .load_json import load_json


def raise_error(key):
    print(f"Exception: Please set the correct parameter: '{key}'!")
    raise Exception(f"Please set the correct parameter: '{key}'!")


def parameter_check(filepath):
    try:
        config = load_json(filepath)
        
        # initial default values
        if "HTP" not in config:
            config["HTP"] = False
        if "miss" not in config:
            config["miss"] = "---"
        if "cores" not in config:
            config["cores"] = 1
        if "sep" not in config:
            config["sep"] = "\t"

        if "mode" in config and config["mode"] in ["BD", "PI", "PR", "DPR"]:
            if "locus_file" in config and os.path.exists(config["locus_file"]):
                pass
            else:
                raise_error("locus_file")

            if "output" in config:
                if not os.path.exists(config["output"]):
                    pass
                else:
                    v = input(f"'{config['output']}' directory already exists,  Do you want ot delete it ? [Y/N]: ")
                    if v in ["Y", "N"]:
                        if v == "Y":
                            os.system(f"rm -rf {config['output']}")
                        else:
                            raise_error("output directory already exists, please modify the parameter or delete the directory.")
                    else:
                        raise_error("output directory already exists, please modify the parameter or delete the directory.")
            else:
                raise_error("output")
    
            if "param" in config:
                param = config["param"]
            else:
                raise_error("param")

            if config["mode"] == "BD":
                if "maf_limit" not in config["param"]:
                    config["param"]["maf_limit"] = 0.05
                if "pic_limit" not in config["param"]:
                    config["param"]["pic_limit"] = 0.02

                if "genotypes" in param and os.path.exists(param["genotypes"]):
                    pass
                else:
                    raise_error("param.genotypes")

            if config["mode"] == "PI":
                if "target" in param and os.path.exists(param["target"]):
                    pass
                else:
                    raise_error("param.target")

                if "ancestors" in param and isinstance(param["ancestors"], list):
                    for fp in param["ancestors"]:
                        if os.path.exists(fp):
                            pass
                        else:
                            raise_error(f"param.ancestors: {fp} is not exists")
                else:
                    raise_error("param.ancestors")

            if config["mode"] == "PR":
                if "target" in param and os.path.exists(param["target"]):
                    pass
                else:
                    raise_error("param.target")

                if "ancestors" in param and isinstance(param["ancestors"], dict):
                    for key, value in param["ancestors"].items():
                        if os.path.exists(value):
                            pass
                        else:
                            raise_error(f"param.ancestors: {value} is not exists of {key}")

                if "colors" in param and isinstance(param["colors"], dict):
                    ancestor_keys = list(param["ancestors"].keys())
                    colors_keys = list(param["colors"].keys())

                    for a_k in ancestor_keys:
                        if a_k not in colors_keys:
                            raise_error(f"param.colors: {a_k} has not color value set")

                else:
                    raise_error("param.colors")

            if config["mode"] == "DPR":
                if "ancestor" in param and os.path.exists(param["ancestor"]):
                    pass
                else:
                    raise_error("param.ancestor")

                if "derivedlines" in param and isinstance(param["derivedlines"], dict):
                    for key, value in param["derivedlines"].items():
                        if os.path.exists(value):
                            pass
                        else:
                            raise_error(f"param.derivedlines: {value} is not exists of {key}")
                else:
                    raise_error("param.derivedlines")

            return config
        else:
          raise_error("mode")
    except Exception as e:
        traceback.print_exc(e)