#!/share/public/software/python/3.6.3/bin/python3.6
import os, sys

def mkdir(dir):
    if not os.path.exists(dir):
        os.makedirs(dir)
    dir = os.path.abspath(dir)
    if dir[-1] != '/':
        dir += '/'
    return dir


txt = ''
with open(os.path.dirname(os.path.realpath(__file__))+"/wgsa.setting" , 'r') as fh:
    txt = fh.readlines()
txt = "".join(txt)
dr = mkdir(os.path.dirname(sys.argv[2]))
mkdir(os.path.dirname(sys.argv[2])+'/tmp')
apt = {'in':sys.argv[1],'out':sys.argv[2],'dir':dr}
out = txt.format(**apt)

print(out, end='')



