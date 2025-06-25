import re

    # LINE_RE = re.compile(
    # r'^(?P<rec>[A-Z]{3}\d+)\s+' # REC_RESIDUE must be 3 capital letters followed by number (i.e. AAA353)
    # r'(?P<lig>[A-Z]{3}\d+)\s+' # LIG_RESIDUE must be 3 capital letters followed by number (i.e. AAA353)
    # r'(?P<dmin>(?:\d+(?:\.\d*)?|\.\d+))\s+' # dmin must be floating point number (2.0)
    # r'(?P<dmax>(?:\d+(?:\.\d*)?|\.\d+))\s*$')

def check_re(pattern, input):
    if re.fullmatch(pattern, input):
        print('Valid')

check_re(r'^[A-Z]{3}\d+$', 'HIS3')
check_re(r'^(?:0|0?\.\d+|[1-9]\d*(?:\.\d*)?)$', '11.0')