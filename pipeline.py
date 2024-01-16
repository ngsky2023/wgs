import sys, os, subprocess, re, shutil

class makeJob:
    def __init__(self, beforeJob='', outdir='./', config=''):
        self.order = beforeJob
        self.outdir = self.mkdir(outdir)
        self.soft = {}
        pipeDir = os.path.dirname(os.path.realpath(__file__))
        self.soft['pipeDir'] = pipeDir
        self.software(pipeDir+"/soft.list")
        if config != '':
            self.software(config)
        self.cmd = {}
        self.jopt = {}
        self.readCMD(pipeDir+'/cmd.txt')
        self.jobFile = sys.stdout

    def printsched(self, name, cmd, cpu=1, mem=1, queue='all.q', order=''):
        P, orderl = '', ''
        if 'mem' in queue:
            i = queue.find('mem')
            if i+4 == len(queue) or not queue[i+4].isdigit():
                P = queue[i:i+4]
            elif queue[i+4].isdigit():
                P = queue[i:i+5]
        elif 'joyce' in queue:
            P = 'joyce'
        if P.isalnum():
            P = ' -P %s' % (P)
        order = self.isorder(order, name)
        sjm = """job_begin
    name %s
    sched_options -V -cwd -l p=%s,vf=%sg -q %s%s
    cmd_begin
%s
    cmd_end
job_end
%s
""" % (name, cpu, mem, queue, P, self.mcmd(cmd), order)
        print(sjm, file=self.jobFile)

    def printlogin(self, name, cmd, order=''):
        order = self.isorder(order, name)
        sjm = """job_begin
    name %s
    host localhost
    cmd_begin
%s
    cmd_end
job_end
%s
""" % (name, self.mcmd(cmd), order)
        print(sjm, file=self.jobFile)

    def makeSched(self, name, cmd, cpu=0, mem=1, queue='all.q', order='', wd='./'):
        wd = self.mkdir(wd)
        shellfile = wd + name + '.sh'
        mark = wd + name + '.mark'
        with open(shellfile, 'w') as mfh:
            print("cd %s && rm -f ERROR %s && echo START $(date)" % (wd, mark), file=mfh)
            print(cmd, file=mfh)
            print("excode=$?", file=mfh)
            print("[ $excode -ne 0 ] && echo $(date) > ERROR && exit 1", file=mfh)
            print("[ $excode -eq 0 ] && echo END   $(date) && echo $(date) > "+mark, file=mfh)
#        mcc = """%s : %s
#\techo %s START $(date) && cd %s && qsub -cwd -V -l vf=%sG,p=%s -q %s -sync y %s && ls > /dev/null && %s/sleep.pl -i %s && echo %s END $(date) || exit 1
#""" % (mark, order, name, wd, mem, cpu, queue, shellfile, self.soft['pipeDir'], mark, name)
        mcc = """%s : %s
\techo %s START `date` && %s/qsub.pl -i %s -l vf=%sG,p=%s -q %s -s %s && echo %s END `date` || exit 1
""" % (mark, order, name, self.soft['pipeDir'], mark, mem, cpu, queue, shellfile, name)

        print(mcc, file=self.jobFile)
        return mark

    def makeLocal(self, name, cmd, order='', wd='./'):
        wd = self.mkdir(wd)
        shellfile = wd + name + '.sh'
        mark = wd + name + '.mark'
        with open(shellfile, 'w') as mfh:
            print("cd %s && rm -f ERROR %s && echo START $(date) &> %s.log" % (wd, mark, shellfile), file=mfh)
            print(cmd, file=mfh)
            print("excode=$?", file=mfh)
            print("[ $excode -ne 0 ] && echo $(date) > ERROR && exit 1", file=mfh)
            print("[ $excode -eq 0 ] && echo END   $(date) &>> "+shellfile+".log && echo $(date) > "+mark, file=mfh)
        mcc = """%s : %s
\techo %s START `date` && cd %s && sh %s && ls > /dev/null && [ -e %s ] && echo %s END `date` || exit 1
""" % (mark, order, name, wd, shellfile, mark, name)
        print(mcc, file=self.jobFile)
        return mark

    def software(self, file):
        softf = open(file , 'r')
        for line in softf:
            sp = line.strip(" \n").split("\t")
            if len(sp) != 2:
                continue
            sf, path = sp
            self.soft[sf] = path
        softf.close()
    
    def isorder(self, beforeJobs, job):
        order = ''
        if len(beforeJobs) > 0:
            for beforeJob in beforeJobs.split(','):
                order += "order %s after %s\n" % (job, beforeJob)
        return order
    
    def mkdir(self, dir):
        if not os.path.exists(dir):
            os.makedirs(dir)
        dir = os.path.abspath(dir)
        if dir[-1] != '/':
            dir += '/'
        return dir

    def mcmd(self, cmd):
        if isinstance(cmd, list):
            return " &&\n".join(map(self.addtab(), cmd))
        else:
            return self.addtab(cmd.strip(" \t\n&"))
    
    def addtab(self, cmd):
        #return " &&\n".join(map(lambda x: "\t\t"+x.strip(" \t\n&"), cmd.split("\n")))
        return " &&\n".join(map(lambda x: x.strip(" \t\n&"), cmd.split("\n")))
    
    def qsubt(self, sh, wd='', l='vf=1g,p=1', n=1, name='job'):
        n = int(n)
        if wd == '':
            wd = self.outdir
        i = 0
        if sh.strip().find(' ') != -1:
            cmd = sh.strip("\t\n ")
            sh = wd+'/'+name+'.sh'
            with open(sh, 'w') as cf:
                cf.write(cmd)
        shf = open(sh, 'r')
        for line in shf:
            if line.strip('\n\t ') != '':
                i += 1
        shf.close()
        m = i//n
        if i%n != 0:
            m += 1
        tc = m
        maxJob = 100
        if tc > maxJob:
            tc = maxJob
        CMD = ['qsub', '-N', name, '-wd', wd, '-t', '1-'+str(m), '-tc', str(tc), '-l', l, '-sync', 'y', self.soft['runtask'], sh, str(n)]
        print(' '.join(CMD))
        jobo = subprocess.Popen(CMD, stdout=subprocess.PIPE)
        jobo.wait()
        jout = jobo.stdout.read().decode()
        print(jout)
        jco = jout.split('\n')
        rem = re.match(r'Your job-array (\d+)', jco[0])
        pid = [1] * (m+1)
        ex = 0
        if rem:
            pid[0] = rem.group(1)
            for jo in jco[1:]:
                if m > 1:
                    rem = re.match(r'Job %s.(\d+) exited with exit code (\d+).'%(pid[0]), jo)
                    if rem:
                        pid[int(rem.group(1))] = rem.group(2)
                elif m == 1:
                    rem = re.match(r'Job %s exited with exit code (\d+).'%(pid[0]), jo)
                    if rem:
                        pid[1] = rem.group(1)
            logdir = self.mkdir(wd+'/log')
            for i in range(1, (m+1)):
                lof = "%s/%s.o%s.%s"%(wd, name, pid[0], i)
                lef = "%s/%s.e%s.%s"%(wd, name, pid[0], i)
                shutil.move(lof, logdir)
                shutil.move(lef, logdir)
                if pid[i] != '0':
                    ex = int(pid[i])
                    lcmd = ''
                    ck = 0
                    with open(lof, 'r') as lofh:
                        for jl in lofh:
                            if jl == '[CMD start]':
                                ck = 1
                            if jl == '[CMD end]':
                                ck = 0
                            elif ck == 1:
                                lcmd += jl
                    print("job %s.%s failed. please redo the command: \n%s\n" % (pid[0], i, lcmd))
        else:
            sys.exit(2)
        if ex != 0:
            sys.exit(ex)
        elif jobo.returncode != 0:
            sys.exit(jobo.returncode)

    def readCMD(self, txt):
        f = open(txt, 'r')
        name = ''
        for line in f:
            if line[0:2] == '#!':
                b = line[2:].strip(' \n').split("\t")
                name = b[0]
                self.cmd[name] = ''
                if len(b) > 1:
                    self.jopt[name] = b[1:]
            elif line[0] == '#' or line.isspace():
                continue
            else:
                self.cmd[name] += line
        f.close()

    def prefix(self, r1, r2):
        prx = ''
        r1 = os.path.basename(r1)
        r2 = os.path.basename(r2)
        for i in range(len(r1)):
            if r1[i] == r2[i]:
                prx += r1[i]
            else:
                break
        prx = prx.strip(' \t\n_.')
        if len(prx) > 0:
            return prx
        else:
            return r1+r2

    def mergeMarkL(self , omark):
        end = ''
        for emk in omark:
            for k in emk.keys():
                if isinstance(emk[k], dict):
                    for ch in emk[k].keys():
                        end += ' ' + emk[k][ch]
                else:
                    end += ' ' + emk[k]
        return end

    def ocmd(self, name, opt, wd='', k=1):
        if wd == '':
            wd = self.outdir
        adpp = ''
        for key in self.soft.keys():
            if key not in opt:
                opt[key] = self.soft[key]
        out = "%s" % (self.cmd[name].format(
            **opt
        ))
        out = out.strip("\n \t")
        if name not in self.jopt or self.jopt[name][0] != 'for':
            out = self.mcmd(out)
        out = adpp + out
        return out
    
    def jcmd(self, name, opt, sample, order='', wd='./'):
        if wd == '':
            wd = self.outdir
        if isinstance(opt, dict):
            cmd = self.ocmd(name, opt, wd=wd)
            if self.jopt[name][0] == 'sched':
                if len(self.jopt[name]) > 3:
                    if 'queue' in self.soft:
                        self.jopt[name][3] = self.soft['queue']
                    return self.makeSched(name+'_'+sample, cmd, cpu=self.jopt[name][1], mem=self.jopt[name][2], queue=self.jopt[name][3], order=order, wd=wd)
                else:
                    return self.makeSched(name+'_'+sample, cmd, cpu=self.jopt[name][1], mem=self.jopt[name][2], order=order, wd=wd)
            elif self.jopt[name][0] == 'login':
                return self.makeLocal(name+'_'+sample, cmd, order=order, wd=wd)
            else:
                print('No use jcmd 1\n', file=sys.stderr);
        elif isinstance(opt, list):
            if self.jopt[name][0] == 'for':
                os.chdir(wd)
                with open(name+'_'+sample+'_pipelineFor.sh', 'w') as forf:
                    for optl in opt:
                        cmd = self.ocmd(name, optl, wd=wd)
                        print(cmd, file=forf)
                subprocess.call(['sh', name+'_'+sample+'_pipelineFor.sh'])
                self.qsubt(wd+'/'+self.jopt[name][3], wd=wd, l=self.jopt[name][1], n=self.jopt[name][2], name=name+'_'+sample)
            else:
                print('No use jcmd 2\n', file=sys.stderr);
        else:
            print('No use jcmd 3\n', file=sys.stderr);
          



