import re
import sys

inputpath = sys.argv[1]
fout = None

for line in open(inputpath):
    if re.match('^[A-Za-z0-9_-]+.[a-z0-9]+:$', line.rstrip()):
        if fout:
            fout.close()
        fout = open(line.rstrip(':\n'), 'w')
    else:
        fout.write(line)
fout.close()
