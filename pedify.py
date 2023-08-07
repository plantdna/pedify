import traceback
import time
import logging
import sys
import argparse

from common import makedirs, parameter_check, build_dataset, pedigree_identify, pedigree_strucutal_analysis, generate_report, derived_pedigree_strucutal_analysis, derived_report, group_pedigree_identify


def pedify():
    try:
        parser = argparse.ArgumentParser(
            description="Pedifys Toolkit: Algorithm And Commend line program For Crop Pedigree Identify And Structural Analysis",
            epilog='Example: python3 pedify.py -c "./config.json"'
        )
        parser.add_argument(
            '-c',
            '--config',
            help="Input config file path",
            required=True,
            type=str
        )
        args = parser.parse_args()

        log_file_path = f"./{int(time.time())}.log"
        logging.basicConfig(filename=log_file_path, filemode='a', level=logging.DEBUG,
                            format='%(asctime)s %(name)s %(levelname)s %(message)s')
        logging.info('load config file')

        config_data = parameter_check(args.config)

        logging.info('check parameters')

        makedirs(config_data['output'])

        if config_data["mode"] == "BD":
            build_dataset(config_data)

        if config_data["mode"] == "PI":
            pedigree_identify(config_data)

        if config_data["mode"] == "GPI":
            group_pedigree_identify(config_data)

        if config_data["mode"] == "PR":
            df = pedigree_strucutal_analysis(config_data)
            if config_data["HTP"]:
                generate_report({
                    "genome": config_data["species"],
                    "output": config_data["output"],
                    "colors": config_data["param"]["colors"]
                }, df)

        if config_data["mode"] == "DPR":
            derived_pedigree_strucutal_analysis(config_data)
            derived_report(config_data)

    except Exception as e:
        exc_type, exc_value, exc_traceback_obj = sys.exc_info()
        print(exc_value)
        traceback.print_exception(
            exc_type,
            exc_value,
            exc_traceback_obj,
            file=open(
                log_file_path,
                'a'
            )
        )


if __name__ == "__main__":
    pedify()
