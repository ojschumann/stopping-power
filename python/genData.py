import os

tmpl = """---Stopping/Range Input Data (Number-format: Period = Decimal Point)
---Output File Name
"{fname}"
---Ion(Z), Ion Mass(u)
 {ion_Z}   {ion_A}
---Target Data: (Solid=0,Gas=1), Density(g/cm3), Compound Corr.
 {gas}    {rho}    0
---Number of Target Elements
 1
---Target Elements: (Z), Target name, Stoich, Target Mass(u)
 {target_Z}   "{target_name}"               1.0             {target_A}
---Output Stopping Units (1-8)
 6
---Ion Energy : E-Min(keV), E-Max(keV)
0 0
0.01
0.1
1
10
100
1000
10000
100000
1000000
2000000
4000000
5000000
0
"""


E=[]
for l in open("SCOEF03.dat", "rt"):
  sym = l[1:3].strip()
  name = l[6:18].strip()
  d = [float(v) for v in l[20:].split()]

  Z = int(d[0])
  A = d[2]
  rho = d[4]
  
  E.append((sym, name, Z, A, rho))


for sym1, name1, Z1, A1, rho1 in E:
  for sym2, name2, Z2, A2, rho2 in E:
    for gas in (0,1):
      fname = f'{sym1}_in_{sym2}_g{gas}.dat'
      if not os.path.exists(fname):
        s = tmpl.format(fname=fname, ion_Z=Z1, ion_A=A1, target_Z=Z2, target_name=name2, target_A=A2, rho=rho2, gas=gas)
        open("SR.IN", "wt").write(s)
        os.system('unix2dos SR.IN')
        os.system('wine SRModule.exe')



