#!/usr/bin/python3
from GUI import GUI
from RefractiveIndex import RefractiveIndex
from Constants import Constants
from Settings import Settings
from PyQt5.QtWidgets import QApplication
import sys

class MultiPhotonAnalysis():
    def __init__(self):
        self.config = Settings()
        self.config.standardSettings()
        self.config.loadSettings()
        constants=Constants()
        print(constants.pi)
        self.gui=GUI(self.config)

    def showGUI(self):
        self.gui.showWindow()

if __name__ == '__main__':
    app = QApplication(sys.argv)
    MPA = MultiPhotonAnalysis()
    app.exec_()
    MPA.config.saveSettings()
    #app = QApplication(sys.argv)
    #MPA=MultiPhotonAnalysis()

    #config = Settings()
    #config.standardSettings()
    #config.loadSettings()
    #cpp=config.get("Crystal Poling Period")[1]
    #config.set("Crystal Poling Period", cpp+1)

    #refIDX = RefractiveIndex()
    #[nx, ny, nz] = refIDX.getIDX('PPKTP',['kato', 'könig', 'fradkin'])
    #g = GUI()
    #MP = MultiPhotonAnalysis()
    #MP.showGUI()

    #app.exec_()

    #config.saveSettings()