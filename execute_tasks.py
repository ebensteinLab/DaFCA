import argparse
import os
import subprocess
from multiprocessing.pool import Pool

import yaml


def call_proc(cmd):
    with open(f"outputs/{cmd[1]}_out.txt", "wb") as out, open(f"outputs/{cmd[1]}_err.txt", "wb") as err:
        print(f"{cmd[1]} return code:{subprocess.check_call(cmd[0], shell=True, stderr=err, stdout=out)}")



def generate_cmd(cmd_dict):
    key = next(iter(cmd_dict))
    cmd_dict = cmd_dict[key]
    cmd = ""
    for key,val in cmd_dict.items():
        if type(val)==list:
            val = ','.join(map(str, val))
        cmd += f"-{key} {val} "
    return cmd


def execute_script(python_alias, max_paralism, path_to_execution_script, is_path_is_global, executions):
    pool = Pool(int(max_paralism))
    if is_path_is_global:
        cmds = [(f"{python_alias} {path_to_execution_script} {generate_cmd(ex)}", next(iter(ex))) for ex in executions]
    else:
        cwd = os.path.dirname(os.path.realpath(__file__))
        cmds = [(f"{python_alias} {cwd}\\{path_to_execution_script} {generate_cmd(ex)}", next(iter(ex))) for ex in
                executions]
    pool.map(call_proc, cmds)
    pool.close()
    pool.join()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--execution_file_path', required=True, help='path to execution file')
    args = parser.parse_args()
    execution_dict = yaml.load(open(args.execution_file_path), Loader=yaml.FullLoader)
    execute_script(**execution_dict)
