#!/usr/bin/python3
from GUI import GUI
from RefractiveIndex import RefractiveIndex
from PMC import PMC
from Constants import Constants
from Settings import Settings
from PyQt5.QtWidgets import QApplication
import sys

class MultiPhotonAnalysis:
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
