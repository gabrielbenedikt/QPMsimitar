import numpy

class Constants:
    def __init__(self):
        self.pi = numpy.pi
        # speed of light in m/s
        self.c = 299792458

        # factors to convert pulsewidth between autocorrelator measurement and real pulsewidth
        self.taucfgauss = 1 / numpy.sqrt(2)  # gaussian shape
        self.taucfsech = 0.647  # sech2 shape

        self.tbwpgauss = 0.441
        self.tbwpsech = ((2 * numpy.log(1 + numpy.sqrt(2))) / (self.pi)) ** 2
    
    def usetaucf(self, useit):
        if useit:
            self.taucfgauss=0.1 / numpy.sqrt(2)
            self.taucfsech=0.647
        else:
            self.taucfgauss = 1
            self.taucfsech = 1
