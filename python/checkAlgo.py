from stopping_power import StoppingPower
import glob
import re




sym2za = {}
name2za = {}
for l in open("SCOEF03.dat", "rt"):
  sym = l[1:3].strip()
  name = l[6:18].strip()
  d = [float(v) for v in l[20:].split()]

  Z = int(d[0])
  A = d[2]

  sym2za[sym] = (Z,A)
  name2za[name] = (Z,A)



from collections import defaultdict
D=defaultdict(int)

for fn in glob.glob('data/*.dat'):
  print (fn)
  f = open(fn, "rt")

  M1 = re.match('Ion = ([A-Z][a-z]*)\s*, Mass = (.*)', f.readline())
  M2 = re.match('Tgt = (\s?[A-Z][a-z]?)\s+\(([A-Za-z]+)\), Density = ([0-9\.]+).*', f.readline())
  f.readline()
  f.readline()

  ion = M1.groups()[0].strip()
  ion_mass = float(M1.groups()[1])
  target = M2.groups()[0].strip()
  isGas = M2.groups()[1] == "Gas"
  rho = float(M2.groups()[2])
  ion_Z = name2za[ion][0]
  target_Z = sym2za[target][0]
  target_mass = sym2za[target][1]


  S=StoppingPower(ion_Z, ion_mass, target_Z, target_mass, rho, isGas)


  for l in f:
    l = l.replace(',', '.')
    E, Sele, Snuc = [float(v) for v in l.split()]

    SnucC = S.getSnucTbl(E)
    SeleC = S.getSele(E)

    dNuc = abs(Snuc-SnucC)/Snuc
    dEle = abs(Sele-SeleC)/Sele

    if True:
      print (f'{E:.3e} | {Snuc:.3e} {SnucC:.3e} {dNuc:.3e} | {Sele:.3e} {SeleC:.3e} {dEle:.3e}')
      if dNuc > 5.2e-4:
        raise Exception("Snuc difference to big")
      if dEle > 5.2e-4:
        raise Exception("Sele difference to big")
    else:
      if dNuc > 5e-4 or dEle > 5e-4:
        def f(x):
           if x > 5e-4:
              return "*"
           else:
              return " "
        print (f'{fn:15} {E:.3e} | {Snuc:.3e} {SnucC:.3e} {dNuc:.3e} {f(dNuc)} | {Sele:.3e} {SeleC:.3e} {dEle:.3e} {f(dEle)}')
        D[target_Z] += 1

for n,z in sorted([(D[z],z) for z in D]):
  print(n,z)