#!/usr/bin/env python3
import numpy as np
import numba

class RefractiveIndex:
    def __init__(self):
        self.initConstants()
        self.materialList = ['PPKTP','KTP']
        self.material = 'PPKTP'

        PPKTPIndices = [['kato'], ['koenig', 'kato'], ['fradkin', 'kato', 'kato2']]
        KTPIndices = [['kato'], ['koenig', 'kato'], ['fradkin', 'kato', 'kato2']]
        self.AvailableIndices = []
        self.AvailableIndices.append(['PPKTP', PPKTPIndices])
        self.AvailableIndices.append(['KTP', KTPIndices])

    def setMaterial(self,material):
        if material in self.materialList:
            self.material = material
        else:
            print('ERROR: material unknown')

    def getMaterial(self):
        return self.material

    def getMaterialList(self):
        return self.materialList

    def getAvailableRefractiveIndices(self,material):
        success=False
        for i in range(0,len(self.AvailableIndices)):
            if (self.AvailableIndices[i][0]==material):
                index=i
                success=True
                break
        if success==True:
            return self.AvailableIndices[index][1]
        else:
            print('ERROR: Material unknown.')
            return -1

    def getIDX(self,crystaltype,IDXtype):
        global ridz, ridy, ridx
        success = 0
        if (len(IDXtype)!=3):
            print('Error: need 3 entries when requesting refractive indices (for nx, ny and nz')
            return -1
        if (crystaltype=='PPKTP'):
            if (IDXtype[0]=='kato'):
                ridx = self.PPKTP_nx_kato
                success = success + 1
            else:
                print('Error: no refractive index type ' , IDXtype[0] , 'known')

            if (IDXtype[1]=='könig'):
                ridy=self.PPKTP_ny_koenig
                success = success + 1
            elif (IDXtype[1]=='kato'):
                ridy=self.PPKTP_ny_kato
                success = success + 1
            else:
                print('Error: no refractive index type ' , IDXtype[1] , 'known')

            if (IDXtype[2]=='fradkin'):
                ridz = self.PPKTP_nz_fradkin
                success = success + 1
            elif (IDXtype[2] == 'kato'):
                ridz = self.PPKTP_nz_kato
                success = success + 1
            elif (IDXtype[2]=='kato2'):
                ridz=self.PPKTP_nz2_kato
                success = success + 1
            else:
                print('Error: no refractive index type ' , IDXtype[2] , 'known')

            if (success==3):
                return [ridx, ridy, ridz]
            else:
                print('Error occured.')
                return -1

        else:
            print('Error: Crystal type unknown')
            return -1

    def getSingleIDX(self,material,pol,paper):
        if (material == 'PPKTP' or material == 'KTP'):
            if (pol == 'X'):
                if (paper == 'kato'):
                    return self.PPKTP_nx_kato
                else:
                    print('Error: Paper unknown')
                    return -1
            elif (pol == 'Y'):
                if (paper == 'kato'):
                    return self.PPKTP_ny_kato
                elif (paper == 'koenig'):
                    return self.PPKTP_ny_koenig
                else:
                    print('Error: Paper unknown')
                    return -1
            elif (pol == 'Z'):
                if (paper == 'kato'):
                    return self.PPKTP_nz_kato
                if (paper == 'kato2'):
                    return self.PPKTP_nz2_kato
                elif (paper == 'fradkin'):
                    return self.PPKTP_nz_fradkin
                else:
                    print('Error: Paper unknown')
                    return -1
            else:
                print('Error: Polarizaion unknown')
                return -1
        else:
            print('Error: Material unknown')
            return -1

    def initConstants(self):
        # speed of light in µm/s
        c = 299792458000000

        # refractive indices fit parameters
        # emanueli & arie 2003
        self.a1z = np.array([9.9587 * 10 ** (-6), 9.9228 * 10 ** (-6), -8.9603 * 10 ** (-6), 4.1010 * 10 ** (-6)])
        self.a2z = np.array([-1.1882 * 10 ** (-8), 10.459 * 10 ** (-8), -9.8136 * 10 ** (-8), 3.1481 * 10 ** (-8)])
        self.a1y = np.array([6.2897 * 10 ** (-6), 6.3061 * 10 ** (-6), -6.0629 * 10 ** (-6), 2.6486 * 10 ** (-6)])
        self.a2y = np.array([-0.14445 * 10 ** (-8), 2.2244 * 10 ** (-8), -3.5770 * 10 ** (-8), 1.3470 * 10 ** (-8)])

        # refractive indices fit parameters
        # kato & takaoka 2002
        # thermal
        self.ax = np.array([0.1627, 0.8416, -0.5353, 0.1717])
        self.ay = np.array([0.5425, 0.5154, -0.4063, 0.1997])
        self.az = np.array([-0.1897, 3.6677, -2.9220, 0.9221])
        self.az2 = np.array([-0.5523, 3.3920, -1.7101, 0.3424])
        # sellmeier
        self.fnx = np.array([3.2910, 0.04140, -0.03978, 9.35522, -31.45571])
        self.fny = np.array([3.45018, 0.04341, -0.04597, 16.98825, -39.43799])
        self.fnz = np.array([4.59423, 0.06206, -0.04763, 110.80672, -86.12171])

        # Thermal expansion coefficients of KTP
        # (Emanueli 2003)
        self.TXCa = 6.7 * 10 ** (-6)
        self.TXCb = 11 * 10 ** (-9)
        self.TXrefT = 25

    # refractive indices
    def PPKTP_ny_koenig(self, lin, t):
        # input in meter, equations for µm
        return self.PPKTP_ny_koenig_numba(lin, t, self.a1y, self.a2y)
    @staticmethod
    @numba.njit(parallel=True)
    def PPKTP_ny_koenig_numba(lin, t, a1y, a2y):
        l = lin * 10 ** 6
        ll = l**2
        linv = 1/l
        # könig, wong 2002
        # emanueli 2003
        f1 = a1y[0] + linv * (a1y[1] + linv * (a1y[2] + (a1y[3] *linv)))
        f2 = a2y[0] + linv * (a2y[1] + linv * (a2y[2] + (a2y[3] *linv)))
        f3 = np.sqrt(np.abs(-0.0138408 * ll + 0.922683 / (1 - 0.0467695 / ll) + 2.0993))
        return (t - 25) * ( f1 +  (t - 25) * f2 )  + f3

    # refractive indices
    def PPKTP_nz_fradkin(self, lin, t):
        return self.PPKTP_nz_fradkin_numba(lin, t, self.a1z, self.a2z)
    @staticmethod
    @numba.njit(parallel=True)
    def PPKTP_nz_fradkin_numba(lin, t, a1z, a2z):
        # input in meter, equations for µm
        l = lin * 10 ** 6
        ll = l**2
        linv = 1/l
        # fradkin, arie, skilar, rosenman 1999
        # emanueli 2003
        f1 = a1z[0] + linv * (a1z[1] + linv * (a1z[2] + (a1z[3] *linv)))
        f2 = a2z[0] + linv * (a2z[1] + linv * (a2z[2] + (a2z[3] *linv)))
        f3 = np.sqrt(np.abs(-0.00968956 * ll + 1.18431 / (1 - 0.0514852 / ll) + 0.6603 / (1 - 100.005 / ll) + 2.12725))
        return (t - 25) * (f1 + (t - 25)  * f2 ) + f3

    # refractive indices
    def PPKTP_nx_kato(self, lin, t):
        # input in meter, equations for µm
        return self.PPKTP_nx_kato_numba(lin, t, self.fnx, self.ax)
    @staticmethod
    @numba.njit(parallel=True)
    def PPKTP_nx_kato_numba(lin, t, fnx, ax):
        l = lin * 10 ** 6
        linv=1/l
        ll= l**2
        # kato & takaoka 2002
        # note that this is at 20°C
        nx = np.sqrt(fnx[0] + fnx[1] / (ll + fnx[2]) + fnx[3] / (ll + fnx[4]))
        # change of nx per °C
        dnxdT = (ax[0] + linv * (ax[1] + linv*(ax[2] + linv*ax[3]))) * 10 ** (-5)
        dT = t - 20
        return (nx + dT * dnxdT)

    def PPKTP_ny_kato(self, lin, t):
        return self.PPKTP_ny_kato_numba(lin, t, self.fny, self.ay)
    @staticmethod
    @numba.njit(parallel=True)
    def PPKTP_ny_kato_numba(lin, t, fny, ay):
        # input in meter, equations for µm
        l = lin * 10 ** 6
        linv=1/l
        ll = l**2
        # kato & takaoka 2002
        # note that this is at 20°C
        ny = np.sqrt(fny[0] + fny[1] / (ll + fny[2]) + fny[3] / (ll + fny[4]))
        # change of nx per °C
        dnydT = (ay[0] + linv * (ay[1] + linv*(ay[2] + linv*ay[3] ))) * 10 ** (-5)
        dT = t - 20
        return (ny + dT * dnydT)

    def PPKTP_nz_kato(self, lin, t):
        return self.PPKTP_nz_kato_numba(lin, t, self.fnz, self.az)
    @staticmethod
    @numba.njit(parallel=True)
    def PPKTP_nz_kato_numba(lin, t, fnz, az):
        #input in meter, equations for µm
        l=lin*10**6
        linv=1/l
        ll = l**2
        # kato & takaoka 2002
        # note that this is at 20°C
        nz = np.sqrt(fnz[0] + fnz[1] / (ll + fnz[2]) + fnz[3] / (ll + fnz[4]))
        # change of nx per °C
        dnzdT = (az[0] + linv*(az[1] + linv*(az[2] + linv*az[3] ))) * 10 ** (-5)
        dT = t - 20
        return (nz + dT * dnzdT)

    def PPKTP_nz2_kato(self, lin, t):
        return PPKTP_nz2_kato_numba(lin, t, self.fnz, self.az2)
    @staticmethod
    @numba.njit(parallel=True)
    def PPKTP_nz2_kato_numba(lin, t, fnz, az2):
        #input in meter, equations for µm
        l=lin*10**6
        ll=l**2
        # kato & takaoka 2002
        # note that this is at 20°C
        nz = np.sqrt(fnz[0] + fnz[1] / (ll + fnz[2]) + fnz[3] / (ll + fnz[4]))
        # change of nx per °C
        dnzdT = (az2[0]/l + az2[1] + az2[2] * l + az2[3] * ll) * 10 ** (-5)
        dT = t - 20
        return (nz + dT * dnzdT)

    # thermal expansion factor
    def thermexpfactor(self, T):
        dT=T - self.TXrefT
        return (1 + self.TXCa * dT + self.TXCb * dT**2)
