#!/usr/bin/python3
from GUI import GUI
from RefractiveIndex import RefractiveIndex
from PyQt5.QtWidgets import QApplication
import sys

class MultiPhotonAnalysis():
    def __init__(self):
        self.gui=GUI()

    def showGUI(self):
        self.gui.showWindow()

if __name__ == '__main__':
    app = QApplication(sys.argv)
    MPA=MultiPhotonAnalysis()

    #refIDX = RefractiveIndex()
    #[nx, ny, nz] = refIDX.getIDX('PPKTP',['kato', 'könig', 'fradkin'])
    #g = GUI()
    #MP = MultiPhotonAnalysis()
    #MP.showGUI()

    app.exec_()