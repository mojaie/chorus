#
# (C) 2014-2017 Seiji Matsuoka
# Licensed under the MIT License (MIT)
# http://opensource.org/licenses/MIT
#

import re


def substitute(dict_, source):
    """ Perform re.sub with the patterns in the given dict
    Args:
      dict_: {pattern: repl}
      source: str
    """
    d_esc = (re.escape(k) for k in dict_.keys())
    pattern = re.compile('|'.join(d_esc))
    return pattern.sub(lambda x: dict_[x.group()], source)


def decode(byte_):
    try:
        res = byte_.decode('utf-8')
    except UnicodeDecodeError:
        try:
            res = byte_.decode("cp1252")
        except UnicodeDecodeError:
            try:
                res = byte_.decode("shift-jis")
            except UnicodeDecodeError:
                raise ValueError("Unsupported codec")
    return res


def decode_file(path):
    with open(path, 'rb') as f:
        for line in f:
            yield decode(line)
