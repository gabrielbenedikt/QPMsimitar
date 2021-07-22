import numpy

class Constants:
    def __init__(self):
        self.pi = numpy.pi
        # speed of light in m/s
        self.c = 299792458

        # factors to convert pulsewidth between autocorrelator measurement and pulsewidth for various pulseshapes
        self.taucfgauss = 1 / numpy.sqrt(2)  # gaussian shape
        self.taucfsech = 0.647  # sech2 shape
        self.taucflorentzian = 1/2.0

        #time bandwidth products
        self.tbwpgauss = 0.441
        self.tbwpsech = ((2 * numpy.log(1 + numpy.sqrt(2))) / (self.pi)) ** 2
        self.tbwpsinc = 1
    
    def usetaucf(self, useit):
        if useit:
            self.taucfgauss=0.1 / numpy.sqrt(2)
            self.taucfsech=0.647
        else:
            self.taucfgauss = 1
            self.taucfsech = 1
