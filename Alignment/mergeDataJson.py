#!/share/public/software/python/3.6.3/bin/python3.6
import os, sys
import json

jl = []
for ffj in sys.argv[1:]:
    with open(ffj , 'r') as fh:
        jl.append(json.loads(fh.read()))

stat = {}

n = 0
name = []
rate = {}
ratename = []
for js in jl:
    for cc in ['filtering_result' , 'adapter_cutting']:
        for key in js[cc].keys():
            if isinstance(js[cc][key] , dict):
                continue
            if n == 0:
                name.append(key)
            if key in stat:
                stat[key] += js[cc][key]
            else:
                stat[key] = js[cc][key]

    for cc in ['after_filtering']:
        for key in js['summary'][cc].keys():
            if isinstance(js['summary'][cc][key] , dict):
                continue
            elif js['summary'][cc][key] > 0 and js['summary'][cc][key] < 1:
                if n == 0:
                    ratename.append(key)
                if key in rate:
                    rate[key] += js['summary'][cc][key] * js['summary'][cc]['total_bases']
                else:
                    rate[key] = js['summary'][cc][key] * js['summary'][cc]['total_bases']
            else:
                if n == 0:
                    name.append(key)
                if key in stat:
                    stat[key] += js['summary'][cc][key]
                else:
                    stat[key] = js['summary'][cc][key]
    n += 1

for key in name:
    print("%s\t%s" % (key , stat[key]))

for key in ratename:
    print("%s\t%s" % (key , rate[key]/stat['total_bases']))





