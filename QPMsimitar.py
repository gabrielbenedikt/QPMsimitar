#!/usr/bin/env python3

has_qt6 = False
try:
    from PyQt6.QtWidgets import QApplication
    has_qt6 = True
except ModuleNotFoundError:
    print("Qt6 not found. Usingg Qt5 fallback")
    from PyQt5.QtWidgets import QApplication
from GUI import GUI
from RefractiveIndex import RefractiveIndex
from PMC import PMC
from Constants import Constants
from Settings import Settings
import sys

class QPMsimitar:
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
    MPA = QPMsimitar()
    if has_qt6:
        app.exec()
    else:
        app.exec_()
    MPA.config.saveSettings()
