""" Gernetate documentation for config """
from pathlib import Path
import re
import subprocess as sp
from typing import List
from jinja2 import Environment, FileSystemLoader

# pylint: disable=W1401
DOCKER_REPO_OPENMS = 'blcdsdockerregistry/openms-thirdparty:2.8.0'
MODULE_COMMANDS = {
        'PeptideIndexer': {
            'path': 'modules/index_peptides.nf',
            'image': 'OpenMS'
        },
        'FalseDiscoveryRate': {
            'path': 'modules/score_FDR_OpenMS.nf',
            'image': 'OpenMS'
        },
        'ConsensusID': {
            'path': 'modules/aggregate_consensus_id.nf',
            'image': 'OpenMS'
        },
        'IDFilter': {
            'path': 'modules/filter_PSM.nf',
            'image': 'OpenMS',
            'prefix': 'filter_psm'
        }
    }

def get_command_args(path:Path):
    """ Read NF module source code and extract the arguments needed. """
    with open(path, 'rt') as handle:
        scopes = ['PARAMETER_MAPPER', 'FLAG_MAPPER', 'ADDITIONAL_MAPPER']
        args = {}
        cur = None
        content = ''
        for line in handle:
            line = line.rstrip()
            if not cur:
                for scope in scopes:
                    if line.lstrip().startswith(scope):
                        cur = scope
                        break
            if not cur:
                continue
            content += line
            if line.endswith(']'):
                content = content.split('=')[1].strip().strip('[]').split(',')
                if content != [':']:
                    for item in content:
                        pattern = re.compile(r"'(.+)':\s*'(.+)'")
                        match = pattern.search(item)
                        val = match.group(1)
                        key = match.group(2)
                        args[key] = val
                if cur == 'FLAG_MAPPER':
                    return args
                cur = None
                content = ''
        return args

def wrap_text(val:str, width:int):
    """ Wrap text to certain width """
    val = val.replace(':', ': ')
    left = val.split(' ')
    right = []
    while len(' '.join(left)) > width:
        right.append(left.pop(-1))
    right.reverse()
    return ' '.join(left).replace(': ', ':'), ' '.join(right).replace(': ', ':')

def get_command_data_openms(command:str, args:List[str], prefix:str=None):
    """ Get command data for OpenMS """
    cmd = f"docker run --rm {DOCKER_REPO_OPENMS} {command} --helphelp"
    print(cmd)
    res = sp.run(cmd, shell=True, capture_output=True, check=False)
    message = res.stderr.decode('utf-8')
    arg_desc = []
    for line in message.split('\n'):
        line = line.strip()
        if not line.startswith('-'):
            continue
        pattern1 = re.compile('^(-\S+) <(.+?)>\*?\s+(.+)$')
        pattern2 = re.compile('^(-\S+) (\[.+\])\s+(.+)$')
        pattern3 = re.compile('^(-\S+)\s+(.+)$')
        match = pattern1.search(line)
        if match:
            arg = match.group(1)
            _type = match.group(2)
            desc = match.group(3)
        else:
            match = pattern2.search(line)
            if match:
                arg = match.group(1)
                _type = match.group(2)
                desc = match.group(3)
            else:
                match = pattern3.search(line)
                if not match:
                    raise ValueError
                arg = match.group(1)
                _type = 'bool'
                desc = match.group(2)
        if not arg in args:
            continue
        if _type == 'decoy_string':
            _type = 'str'
        param = args[arg]
        if prefix:
            param = param.replace(prefix, '<>')
        arg_desc.append([param, f"({arg})", _type, desc])
    return arg_desc

def get_command_desc(command:str, args:List[str], image, prefix:str=None):
    """ Get the paramter name, type, and description for the specified OpenMS
    command/program """
    # pylint: disable=R0914
    if image == 'OpenMS':
        arg_desc = get_command_data_openms(command, args, prefix)
    else:
        raise ValueError(f"Image {image} not supported")
    width_param = max(len(x[0]) for x in arg_desc)
    width_arg = max(len(x[1]) for x in arg_desc)
    width_type = max(len(x[2]) for x in arg_desc)

    width_arg = max([width_param, width_arg])

    width_desc = 130 - 13 - width_arg - width_type

    widths = [width_arg + 3, width_type + 3, width_desc + 3]

    table:List[List[List[str]]] = []
    for param in arg_desc:
        table.append([])
        left, param[3] = wrap_text(param[3], width_desc)
        table[-1].append([
            param[0].ljust(width_arg), param[2].ljust(width_type),
            left.ljust(width_desc)
        ])
        left, param[3] = wrap_text(param[3], width_desc)
        table[-1].append([
            param[1].ljust(width_arg), ' ' * width_type, left.ljust(width_desc)
        ])
        while len(param[3]) > width_desc:
            left, param[3] = wrap_text(param[3], width_desc)
            table[-1].append([
                ' ' * width_arg, ' ' * width_type, left.ljust(width_desc)
            ])
        table[-1].append([
            ' ' * width_arg, ' ' * width_type, param[3].ljust(width_desc)
        ])

    return widths, table


def get_params_data(commands:dict):
    """ get params data """
    data = {}
    for command in commands.keys():
        path = commands[command]['path']
        image = commands[command]['image']
        if 'prefix' in commands[command]:
            prefix = commands[command]['prefix']
        else:
            prefix = None
        args = get_command_args(path)
        widths, params = get_command_desc(command, args, image, prefix)
        data[command] = {'widths': widths, 'params': params}
    return data

def render_params_rst():
    """ render template """
    env = Environment(loader=FileSystemLoader('docs/partials'))
    template = env.get_template('_params.rst')
    data = get_params_data(MODULE_COMMANDS)
    print(data.keys())
    rendered = template.render(data=data)
    with open('docs/params.rst', "w") as handle:
        handle.write(rendered)

if __name__ == '__main__':
    render_params_rst()
