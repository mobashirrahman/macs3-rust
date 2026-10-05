import json,os,subprocess,glob,sys,re
root=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
rows=[]
for xls in sorted(glob.glob(f'{root}/tests/golden/*/*/default/*_peaks.xls')):
    fx='/'.join(xls.split('/tests/golden/')[1].split('/')[:2])
    d=f'{root}/tests/fixtures/{fx}'
    if not os.path.exists(f'{d}/treat.bedpe'): continue
    cfg=f'{root}/tests/golden/{fx}/default/command.json'
    argv=json.load(open(cfg))['command']
    keep=[];i=0
    while i<len(argv):
        a=argv[i]
        if a=='-f': keep+=['--format',argv[i+1]]; i+=2; continue
        if a in ('--call-summits',): keep.append(a); i+=1; continue
        if a in ('-g','--gsize','--extsize','--slocal','--llocal','--qvalue'):
            keep+=[{'-g':'--gsize'}.get(a,a),argv[i+1]]; i+=2; continue
        i+=1
    env=dict(os.environ); env['CALLPEAK_LAMBDA']='1'
    r=subprocess.run(['cargo','run','-q','--bin','macs-callpeak-e2e','--',d]+keep+['--xls',xls],
                     capture_output=True,text=True,env=env,cwd=root)
    lam=[l.split('\t') for l in r.stderr.splitlines() if l.startswith('LAM\t')]
    if not lam: continue
    # upstream rows
    up=[]
    for l in open(xls):
        if l.startswith('#') or not l.strip(): continue
        f=l.split('\t')
        if len(f)<9: continue
        try: up.append((f[0],int(f[4]),float(f[5]),float(f[7])))
        except ValueError: continue
    for (c,su,pu,fo) in up:
        m=[x for x in lam if x[1]==c and int(x[2])==su]
        if not m: continue
        mine=float(m[0][4])
        rows.append((fx,fo,mine,(mine-fo)/fo*100))
rows.sort(key=lambda r:-abs(r[3]))
n=len(rows)
# F104: report the WHOLE distribution, not just the worst tail. Printing only the
# top 14 of 72 made "12 of 14 within 0.01%" true while leaving the other 58
# unmeasured -- an instrument that cannot distinguish "almost all good" from "a
# dozen bad and the rest unknown".
import math
buckets=[(0.005,'<= 0.005%'),(0.01,'<= 0.01%'),(0.05,'<= 0.05%'),(0.25,'<= 0.25%'),(1.0,'<= 1%'),(1e9,'> 1%')]
print('compared summits: %d'%n)
prev=0
for lim,lbl in buckets:
    c=sum(1 for r in rows if abs(r[3])<=lim)
    print('  %-12s %4d  (%.1f%% of total)'%(lbl,c,100.0*c/n))
    prev=c
print('median |dev| %.6f%%   p90 %.6f%%'%(sorted(abs(r[3]) for r in rows)[n//2],
                                          sorted(abs(r[3]) for r in rows)[int(n*0.9)]))
print('outside 0.01%%: %d'%(n-prev if False else sum(1 for r in rows if abs(r[3])>0.01)))
for fx,fo,mine,pct in rows[:10]:
    print('  %-40s up=%.5f mine=%.5f %+.4f%%'%(fx,fo,mine,pct))
