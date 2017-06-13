#!/usr/bin/python3


class Crystal:
    def __init(self):
        from enum import Enum
        self.AvailableCrystals = Enum('PPKTP')

    def setCrystalType(self, type):
        success=False
        for crystal in self.AvailableCrystals:
            if (crystal.name==type):
                self.CrystalType = type
                success=True
                return 0
        if (success==False):
            print('Crystal type unknown')
            return -1

    def setPolingPeriod(self,period):
        self.PolingPeriod = period
        return 0


    def getPolingPeriod(self):
        assert isinstance(self.PolingPeriod, object)
        return self.PolingPeriod

    def getCrystalType(self):
        assert isinstance(self.CrystalType, object)
        return self.CrystalType
