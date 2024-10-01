from numpy import *


Na = 6.0223e23 # at/mol


# Coeff
# 10-19: Parameter für SElec H
# 20-29: Parameter für SElec He
# 30-39: Parameter für SElec Li
# 40-52: Polynom in log(E) in SElec Other für Ion
# 53-65: Polynom in log(E) in SElec Other für Target
# 66-68: Parameter für calc_high -log(P) (isGas / ion_Z=1)
# 69-79: Polynom in log(E) in calc_high
# 80-87: Polynom in log(E) in calc_high für Gas-Correction
# 88-92: High energy correction in SElec Other (Ion)

class StoppingPower:
  SCOEF = None
  SNUC = None

  def __init__(self, ion_Z, ion_mass, target_Z, target_mass, rho, isGas=False):
    if self.SCOEF is None:
       self.loadCoeff()

    self._ion_Z = ion_Z
    self._ion_mass = ion_mass
    self._target_Z = target_Z
    self._target_mass = target_mass
    self._isGas = isGas


    # Constants for the calculation of S_nuclear
    self._zf = (ion_Z**0.23+target_Z**0.23)
    self._reduced_energy_factor  = 32.53         * target_mass / (ion_Z*target_Z*(ion_mass+target_mass)*self._zf)
    self._reduced_energy_factor2 = 1e3*0.4685 / 14.4 * target_mass / (ion_Z*target_Z*(ion_mass+target_mass)*self._zf)

    self._sNucTblFactor = 8.462 * ion_Z*target_Z*ion_mass/(ion_mass+target_mass)/self._zf * (Na * 1e-21 / target_mass)
    self._sNucFactor = 4 * 3.1415926 * 14.4*0.4685 * 0.6023 / 0.60223   * 60.222 * ion_Z*target_Z*ion_mass/target_mass/(ion_mass+target_mass)/self._zf


    # Constants for the calculation of S_electronic
    ID = tuple(sorted((ion_Z, target_Z)))
    if ID in self.SNUC:
      self._sNuc = self.SNUC[ID]
    else:
      self._sNuc = (1.1383, 0.01321, 0.21226, 0.19593)

    if ion_Z > 3:
      self._ion_poly = self.SCOEF[ion_Z-1, 40:53]
      self._ion_v5c = self.SCOEF[ion_Z-1, 88:93]
      self._target_poly = self.SCOEF[target_Z-1, 53:66]
    self._he_delta = -log(self.SCOEF[target_Z-1, (67 if ion_Z==1 else 68) if isGas else 66])
    self._he_poly = self.SCOEF[target_Z-1, 69:80]
    self._he_solid_poly = self.SCOEF[target_Z-1, 80:88]

    m = 10*ion_Z if ion_Z <=3 else 20
    self._A = self.SCOEF[target_Z-1, m:m+10]
    self._C = 7.0 if ion_Z == 1 else 4.0

    # correction coefficient
    if target_Z in (2, 10, 18, 36, 54, 86): # noble gas
       self._correction = 1
    else:
       gaseous = target_Z in (1, 2, 7, 8, 9, 10, 17, 18, 35, 36, 53, 54, 85, 86) # gaseous
       if gaseous == self._isGas:
          self._correction = 1
       else:
          self._correction = 0.92
          if not gaseous:
             self._correction = 1 / self._correction


  def loadCoeff(self):
     StoppingPower.SCOEF = array([[float(v) for v in l[20:].split()] for l in open("SCOEF03.dat", "rt")])
     StoppingPower.SNUC = { (int(a), int(b)): (float(c), float(d), float(e), float(f)) for a,b,c,d,e,f in [l.split() for l in open("SNUC03.dat", "rt")]}

  def getSnucTbl(self, E):
    eps = E * self._reduced_energy_factor
    if eps > 30:
      sn = 0.5 * log(eps) / eps
    else:
      a, b, c, d = self._sNuc
      sn = log(1+a*eps)/(2*(eps+b*eps**c + d*sqrt(eps)))

    return self._sNucTblFactor * sn

  def getSnuc(self, E):
      eps = E *  self._reduced_energy_factor2
      sn = log(1+1.1383*eps)/(2*(eps+0.01321*eps**0.21226 + 0.19594*sqrt(eps)))
      return sn * self._sNucFactor


  def getSele(self, E):

    C = 1 + (self._correction - 1) / (1 + exp(1.48 * (sqrt(E/self._ion_mass/25) - self._C)))

    if self._ion_Z <=3:
      return  60.222 / self._target_mass * 10 *C * self.calc_S_elec_H(E/self._ion_mass)
    else:
      return  60.222 / self._target_mass * 10 *C * self.calc_S_elec_other(E/self._ion_mass)


  def calc_S_elec_H(self, E):

    Emin = 2.0
    Emax = 1000.0 if self._ion_Z == 1 else 2000

    if E > 1.2*Emax:
      return self.calc_high_energy(E)

    eps = min(max(E, Emin), Emax)
    A = self._A

    v1 = A[0]*eps**A[1] + A[2]*eps**A[3] + A[8]*sqrt(eps)
    v2 = A[4] / eps**A[5] * (log(A[6]/eps + A[7]*eps) + A[9]*eps**-0.2)

    Selec =  v1*v2 / (v1 + v2)

    if E < Emin:
        Selec *= sqrt(E/Emin)
    elif E > Emax:
        S = self.calc_high_energy(1.2*Emax)
        Selec += (S - Selec) / (1.2*Emax - Emax) * (E - Emax)

    return Selec



  def calc_S_elec_other(self, E):
    Emin = 2.0
    Emax = 3e4
    eps = min(max(E, Emin), Emax)

    x = log(eps)
    p1 = dot(self._ion_poly, x**arange(13))
    p2 = dot(self._target_poly, x**arange(13))

    val = 0.25 * self._ion_Z**2 * p1 * p2

    save_Z = self._ion_Z
    self._ion_Z = 2

    val *= self.calc_S_elec_H(eps)

    if E < Emin:
       val *= sqrt(E/Emin)
    elif E > Emax:
       v5c = self._ion_v5c
       def f(e):
          x = log10(e) - v5c[1]
          m = 3 if x > 0 else 2
          y = v5c[0] + (1 - v5c[0]) / (1 + exp( -(x) / v5c[m]))
          return (y + v5c[4] * exp(-abs(x-1.3))) * self.calc_high_energy(e)

       val *= f(E) / f(Emax)

    self._ion_Z = save_Z

    return val



  def calc_high_energy(self, E_over_u):
    x = (1 + E_over_u / 931494)**-2
    var_40 = log(1022000 * (1-x) / x) + x - 1
    var_D0 = 1 - x

    x = log(E_over_u / 1000)
    if self._isGas:
      var_A0 = 1.35 / E_over_u**0.4
    else:
      var_A0 = 1.5 / (E_over_u**0.4) + 45000 / (self._target_Z * E_over_u**1.6)
      var_40 -= dot(self._he_solid_poly, x**arange(8))
    var_40 -= dot(self._he_poly, x**arange(11))
    var_40 += 1e-3 * self._ion_Z * E_over_u * var_A0  / (1e-3 * E_over_u + var_A0)
    var_40 += self._he_delta

    y2 = (self._ion_Z/137.036)**2 / var_D0
    var_40 -= y2 * (1.20206 - y2*(1.042-0.8549*y2+0.343*y2**2))

    return  self._ion_Z**2 * self._target_Z * var_40  / var_D0 / 1961.16
