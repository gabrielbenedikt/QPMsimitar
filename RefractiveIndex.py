#!/usr/bin/python3
import numpy
class RefractiveIndex():
    def __init__(self):
        self.initConstants()
        self.materialList = ['PPKTP','KTP']
        self.material = 'PPKTP'

        PPKTPIndices = [['kato'], ['koenig', 'kato'], ['fradkin', 'kato']]
        KTPIndices = [['kato'], ['koenig', 'kato'], ['fradkin', 'kato']]
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
            elif (IDXtype[2]=='kato'):
                ridz=self.PPKTP_nz_kato
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

    def initConstants(self):
        # speed of light in µm/s
        c = 299792458000000

        # refractive indices fit parameters
        # emanueli & arie 2003
        self.a1z = [9.9587 * 10 ** (-6), 9.9228 * 10 ** (-6), -8.9603 * 10 ** (-6), 4.1010 * 10 ** (-6)]
        self.a2z = [-1.1882 * 10 ** (-8), 10.459 * 10 ** (-8),-9.8136 * 10 ** (-8), 3.1481 * 10 ** (-8)]
        self.a1y = [6.2897 * 10 ** (-6), 6.3061 * 10 ** (-6), -6.0629 * 10 ** (-6), 2.6486 * 10 ** (-6)]
        self.a2y = [-0.14445 * 10 ** (-8), 2.2244 * 10 ** (-8), -3.5770 * 10 ** (-8), 1.3470 * 10 ** (-8)]

        # refractive indices fit parameters
        # kato & takaoka 2002
        # thermal
        self.ax = [0.1627, 0.8416, -0.5353, 0.1717]
        self.ay = [0.5425, 0.5154, -0.4063, 0.1997]
        self.az = [-0.1897, 3.6677, -2.9220, 0.9221]
        # sellmeier
        self.fnx = [3.2910, 0.04140, -0.03978, 9.35522, -31.45571]
        self.fny = [3.45018, 0.04341, -0.04597, 16.98825, -39.43799]
        self.fnz = [4.59423, 0.06206, -0.04763, 110.80672, -86.12171]

        # Thermal expansion coefficients of KTP
        # (Emanueli 2003)
        self.TXCa = 6.7 * 10 ** (-6)
        self.TXCb = 11 * 10 ** (-9)
        self.TXrefT = 25

    # refractive indices
    def PPKTP_ny_koenig(self, l, t):
        # könig, wong 2002
        # emanueli 2003
        f1 = self.a1y[0] + (1 / l) * (self.a1y[1] + (1 / l) * (self.a1y[2] + (self.a1y[3] / l)))
        f2 = self.a2y[0] + (1 / l) * (self.a2y[1] + (1 / l) * (self.a2y[2] + (self.a2y[3] / l)))
        f3 = numpy.sqrt(numpy.abs(-0.0138408 * l * l + 0.922683 / (1 - 0.0467695 / (l * l)) + 2.0993))
        return (t - 25) * ( f1 +  (t - 25) * f2 )  + f3

    # refractive indices
    def PPKTP_nz_fradkin(self, l, t):
        # fradkin, arie, skilar, rosenman 1999
        # emanueli 2003
        f1 = self.a1z[0] + (1 / l) * (self.a1z[1] + (1 / l) * (self.a1z[2] + (self.a1z[3] / l)))
        f2 = self.a2z[0] + (1 / l) * (self.a2z[1] + (1 / l) * (self.a2z[2] + (self.a2z[3] / l)))
        f3 = numpy.sqrt(numpy.abs(
            -0.00968956 * l * l + 1.18431 / (1 - 0.0514852 / (l * l)) + 0.6603 / (1 - 100.005 / (l * l)) + 2.12725))
        return (t - 25) * (f1 + (t - 25)  * f2 ) + f3

    # refractive indices
    def PPKTP_nx_kato(self, l, t):
        # kato & takaoka 2002
        # note that this is at 20°C
        nx = numpy.sqrt(self.fnx[0] + self.fnx[1] / (l * l + self.fnx[2]) + self.fnx[3] / (l * l + self.fnx[4]))
        # change of nx per °C
        dnxdT = (self.ax[0] + self.ax[1] / l + self.ax[2] / (l * l) + self.ax[3] / (l * l * l)) * 10 ** (-5)
        dT = t - 20
        return (nx + dT * dnxdT)

    def PPKTP_ny_kato(self, l, t):
        # kato & takaoka 2002
        # note that this is at 20°C
        ny = numpy.sqrt(self.fny[0] + self.fny[1] / (l * l + self.fny[2]) + self.fny[3] / (l * l + self.fny[4]))
        # change of nx per °C
        dnydT = (self.ay[0] + self.ay[1] / l + self.ay[2] / (l * l) + self.ay[3] / (l * l * l)) * 10 ** (-5)
        dT = t - 20
        return (ny + dT * dnydT)

    def PPKTP_nz_kato(self, l, t):
        # kato & takaoka 2002
        # note that this is at 20°C
        nz = numpy.sqrt(self.fnz[0] + self.fnz[1] / (l * l + self.fnz[2]) + self.fnz[3] / (l * l + self.fnz[4]))
        # change of nx per °C
        dnzdT = (self.az[0] + self.az[1] / l + self.az[2] / (l * l) + self.az[3] / (l * l * l)) * 10 ** (-5)
        dT = t - 20
        return (nz + dT * dnzdT)

    # thermal expansion factor
    def thermexpfactor(self, T):
        return (1 + self.TXCa * (T - self.TXrefT) + self.TXCb * (T - self.TXrefT) * (T - self.TXrefT))