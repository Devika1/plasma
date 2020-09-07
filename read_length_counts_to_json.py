#!/usr/bin/env python
import sys, json
_, f_out = sys.argv
cnt_dict = {}
for l in sys.stdin:
    try:
        length = int(l.strip())
    except ValueError:
        continue
    if length in cnt_dict:
        cnt_dict[length] += 1
    else:
        cnt_dict[length] = 1
with open(f_out, "w") as fp:
    json.dump([cnt_dict], fp)
