#!/usr/bin/env python3

import ruamel.yaml as yaml
from pathlib import Path

class Settings:
    def __init__(self):
        pass

    def saveSettings(self):
        with open('config.yaml','w') as stream:
            try:
                print('Writing config file')
                tmpconfig = yaml.YAML(typ='safe', pure=True)
                tmpconfig.dump(self.config, stream)
            except yaml.YAMLError as exc:
                print(exc)

    def standardSettings(self):
        print('Initializing standard settings.')
        self.config=[]
        #self.config.append(["key",value])
        self.config.append(["Crystal Material", "PPKTP"])
        self.config.append(["Crystal Poling Period From", 46.0*10**(-6)])
        self.config.append(["Crystal Poling Period Single", 46.25*10**(-6)])
        self.config.append(["Crystal Poling Period To", 46.5*10**(-6)])
        self.config.append(["Crystal Refractive Index X", 'kato'])
        self.config.append(["Crystal Refractive Index Y", 'kato'])
        self.config.append(["Crystal Refractive Index Z", 'kato'])

        self.config.append(["QPM Order", -1])

        self.config.append(["Pump wavelength from", 760*10**(-9)])
        self.config.append(["Pump wavelength single", 775*10**(-9)])
        self.config.append(["Pump wavelength to", 790*10**(-9)])

        self.config.append(["SI wavelength from", 1500*10**(-9)])
        self.config.append(["SI wavelength single", 1550*10**(-9)])
        self.config.append(["SI wavelength to", 1600*10**(-9)])

        self.config.append(["Pump pulsewidth from", 1*10**(-12)])
        self.config.append(["Pump pulsewidth single", 2.9*10**(-12)])
        self.config.append(["Pump pulsewidth to", 5*10**(-12)])
        self.config.append(["Pump pulse shape", 'Sech^2'])
        self.config.append(["Pump pulsewidth apply deconvolution factor", True])

        self.config.append(["Crystal Temperature from", 20])
        self.config.append(["Crystal Temperature single", 35])
        self.config.append(["Crystal Temperature to", 50])

        self.config.append(["Crystal Length from", 20*10**(-3)])
        self.config.append(["Crystal Length single", 30*10**(-3)])
        self.config.append(["Crystal Length to", 40*10**(-3)])

        self.config.append(["JSI wavelength range", 8*10**(-9)])
        self.config.append(["JSI resolution", 50])

        self.config.append(["Purity wavelength resolution", 50])
        self.config.append(["Purity wavelength range", 8*10**(-9)])
        self.config.append(["Purity tau resolution", 20])

        self.config.append(["SI filter Idler Type", "None"])
        self.config.append(["SI filter Signal Type", "None"])
        self.config.append(["SI filter Idler center wavelength", 1550*10**(-9)])
        self.config.append(["SI filter Signal center wavelength", 1550*10**(-9)])
        self.config.append(["SI filter Idler FWHM", 3*10**(-9)])
        self.config.append(["SI filter Signal FWHM", 3*10**(-9)])

        self.config.append(["HOM interference plot resolution", 20])
        self.config.append(["HOM interference plot phase", 0])
        self.config.append(["HOM interference plot range", 6*10**(-12)])
        
        self.config.append(['FWHM plot resolution',10])
        self.config.append(['FWHM decimal precision',3])

    def loadSettings(self):
        if Path('config.yaml').is_file():
            with open('config.yaml') as stream:
                try:
                    print('Reading config file.')
                    tmpconfig = yaml.YAML(typ='safe', pure=True)
                    tmpconfig = tmpconfig.load(stream)
                    for i in range(0,len(tmpconfig)):
                        [key,val] = tmpconfig[i]
                        res = self.find(self.config,key)
                        if (res==-1):
                            pass
                        else:
                            self.set(key,val)
                    print(self.config)
                except yaml.YAMLError as exc:
                    print(exc)

    def set(self, key, val):
        idx,unused=self.find(self.config,key)
        self.config[idx][1]=val
        #print('set config value: ', self.config[idx])

    def get(self, key):
        idx,unused=self.find(self.config,key)
        return self.config[idx][1]

    def find(self, l, elem):
        for row, i in enumerate(l):
            try:
                column = i.index(elem)
            except ValueError:
                continue
            return row, column
        return -1
