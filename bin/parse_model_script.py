#!/usr/bin/env python
import re
from sys import argv


with open(argv[1]) as infile:
    for line in infile:
        if line.startswith('altd'):
            break
match_result = re.match(r'altd\s+<-\sc\((\d+).*\)', line)
if match_result is not None:
    print(match_result.groups()[0], end='')
