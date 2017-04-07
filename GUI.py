#!/usr/bin/python3
import sys
from PyQt5.QtCore import QObject, QVariant
from PyQt5.QtWidgets import (QApplication, QWidget, QToolTip)
from PyQt5.QtWidgets import (QPushButton, QSpinBox, QLabel, QDoubleSpinBox, QGroupBox, QComboBox, QCheckBox, QMainWindow)
from PyQt5.QtWidgets import (QGridLayout, QVBoxLayout, QHBoxLayout)
from PyQt5.QtGui import QFont
from RefractiveIndex import RefractiveIndex
#from Multiphoton import MultiPhotonAnalysis
from Settings import Settings

class GUI(QMainWindow):
    def __init__(self, config, parent=None):
        super(GUI, self).__init__(parent)

        self.mainWidget = QWidget(self)
        self.setCentralWidget(self.mainWidget)

        self.config=config

        self.initUI()

    def initUI(self):
        QToolTip.setFont(QFont('SansSerif', 10))
        self.setToolTip('Tooltip!')

        self.setWindowTitle('Multiphoton Analysis GUI')

        self.ui_layout = QGridLayout()

        self.getProperties()

        self.ui_layoutCrystal = self.initLayoutCrystal()
        self.ui_layoutPump = self.initLayoutPump()
        self.ui_layoutTemperature = self.initLayoutTemperature()

        #Group: set crystal properties
        self.ui_layout.addLayout(self.ui_layoutCrystal,     1, 1)
        #Group: set pump properties
        self.ui_layout.addLayout(self.ui_layoutPump,        2, 1)
        #Group: set temperature
        self.ui_layout.addLayout(self.ui_layoutTemperature, 3, 1)

        self.centralWidget().setLayout(self.ui_layout)
        self.resize(self.sizeHint())

        self.show()

    def initLayoutTemperature(self):
        self.ui_layoutTemp = QGridLayout()
        self.ui_TempGroupBox = QGroupBox()
        self.ui_layoutTempGroupBox = QGridLayout()

        self.ui_fromLabel = QLabel()
        self.ui_toLabel = QLabel()
        self.ui_singleLabel = QLabel()
        self.ui_fromLabel.setText('from')
        self.ui_toLabel.setText('to')
        self.ui_singleLabel.setText('single')

        self.ui_TsingleSB = QDoubleSpinBox()
        self.ui_TfromSB = QDoubleSpinBox()
        self.ui_TtoSB = QDoubleSpinBox()
        self.ui_TLabel = QLabel()
        self.ui_TLabel.setText('Temperature [°C]')

        self.ui_layoutTempGroupBox.addWidget(self.ui_fromLabel,     1, 2)
        self.ui_layoutTempGroupBox.addWidget(self.ui_singleLabel,   1, 3)
        self.ui_layoutTempGroupBox.addWidget(self.ui_toLabel,       1, 4)
        self.ui_layoutTempGroupBox.addWidget(self.ui_TLabel,        2, 1)
        self.ui_layoutTempGroupBox.addWidget(self.ui_TfromSB,       2, 2)
        self.ui_layoutTempGroupBox.addWidget(self.ui_TsingleSB,     2, 3)
        self.ui_layoutTempGroupBox.addWidget(self.ui_TtoSB,         2, 4)

        self.ui_TempGroupBox.setLayout(self.ui_layoutTempGroupBox)
        self.ui_layoutTemp.addWidget(self.ui_TempGroupBox)

        return self.ui_layoutTemp

    def initLayoutPump(self):
        self.ui_layoutPump = QGridLayout()
        self.ui_PumpGroupBox = QGroupBox()
        self.ui_layoutPumpGroupBox = QVBoxLayout()
        self.ui_layoutUpper = QGridLayout()
        self.ui_layoutLower = QGridLayout()

        self.ui_fromLabel = QLabel()
        self.ui_toLabel = QLabel()
        self.ui_singleLabel = QLabel()
        self.ui_fromLabel.setText('from')
        self.ui_toLabel.setText('to')
        self.ui_singleLabel.setText('single')

        self.ui_pumpwlsingleSB = QDoubleSpinBox()
        self.ui_pumpwlfromSB = QDoubleSpinBox()
        self.ui_pumpwltoSB = QDoubleSpinBox()
        self.ui_pumpwlLabel = QLabel()
        self.ui_pumpwlLabel.setText('Wavelength [nm]')

        self.ui_pulsewidthsingleSB = QDoubleSpinBox()
        self.ui_pulsewidthfromSB = QDoubleSpinBox()
        self.ui_pulsewidthtoSB = QDoubleSpinBox()
        self.ui_pulsewidthLabel = QLabel()
        self.ui_pulsewidthLabel.setText('Pulsewidth [ps]')

        self.ui_pumpShapeCB = QComboBox()
        self.ui_pumpShapeCB.addItem('Gaussian')
        self.ui_pumpShapeCB.addItem('Sech^2')
        self.ui_pumpShapeLabel = QLabel()
        self.ui_pumpShapeLabel.setText('Pump pulse shape')

        self.ui_pumpShapeCorrectionFactorCheckBox = QCheckBox()
        self.ui_pumpShapeCorrectionFactorLabel = QLabel()
        self.ui_pumpShapeCorrectionFactorLabel.setText('Apply deconvolution factor?')

        self.ui_layoutUpper.addWidget(self.ui_fromLabel,             1, 2)
        self.ui_layoutUpper.addWidget(self.ui_singleLabel,           1, 3)
        self.ui_layoutUpper.addWidget(self.ui_toLabel,               1, 4)
        self.ui_layoutUpper.addWidget(self.ui_pumpwlLabel,           2, 1)
        self.ui_layoutUpper.addWidget(self.ui_pumpwlfromSB,          2, 2)
        self.ui_layoutUpper.addWidget(self.ui_pumpwlsingleSB,        2, 3)
        self.ui_layoutUpper.addWidget(self.ui_pumpwltoSB,            2, 4)
        self.ui_layoutUpper.addWidget(self.ui_pulsewidthLabel,       3, 1)
        self.ui_layoutUpper.addWidget(self.ui_pulsewidthfromSB,      3, 2)
        self.ui_layoutUpper.addWidget(self.ui_pulsewidthsingleSB,    3, 3)
        self.ui_layoutUpper.addWidget(self.ui_pulsewidthtoSB,        3, 4)

        self.ui_layoutLower.addWidget(self.ui_pumpShapeLabel,                       1, 1)
        self.ui_layoutLower.addWidget(self.ui_pumpShapeCB,                          1, 2)
        self.ui_layoutLower.addWidget(self.ui_pumpShapeCorrectionFactorLabel,       2, 1)
        self.ui_layoutLower.addWidget(self.ui_pumpShapeCorrectionFactorCheckBox,    2, 2)

        self.ui_layoutPumpGroupBox.addLayout(self.ui_layoutUpper)
        self.ui_layoutPumpGroupBox.addLayout(self.ui_layoutLower)

        self.ui_PumpGroupBox.setLayout(self.ui_layoutPumpGroupBox)
        self.ui_layoutPump.addWidget(self.ui_PumpGroupBox)

        return self.ui_layoutPump

    def initLayoutCrystal(self):
        # crystal section
        self.ui_layoutCrystal = QGridLayout()
        self.ui_CrystalGroupBox = QGroupBox()
        self.ui_layoutCrystalGroupBox = QGridLayout()

        self.ui_CrystalPolingPeriodSpinBox = QDoubleSpinBox()
        self.ui_CrystalPolingPeriodLabel = QLabel()
        self.ui_CrystalMaterialComboBox = QComboBox()
        self.ui_CrystalMaterialLabel = QLabel()
        self.ui_CrystalNXComboBox = QComboBox()
        self.ui_CrystalNYComboBox = QComboBox()
        self.ui_CrystalNZComboBox = QComboBox()
        self.ui_CrystalNXLabel = QLabel()
        self.ui_CrystalNYLabel = QLabel()
        self.ui_CrystalNZLabel = QLabel()

        self.ui_CrystalPolingPeriodLabel.setText('Poling Period [µm]')
        self.ui_CrystalMaterialLabel.setText('Material')
        self.ui_CrystalNXLabel.setText('Refractive index (x)')
        self.ui_CrystalNYLabel.setText('Refractive index (y)')
        self.ui_CrystalNZLabel.setText('Refractive index (z)')
        self.ui_CrystalGroupBox.setTitle('Crystal')

        self.ui_layoutCrystalGroupBox.addWidget(self.ui_CrystalPolingPeriodSpinBox, 1, 1)
        self.ui_layoutCrystalGroupBox.addWidget(self.ui_CrystalPolingPeriodLabel,   1, 2)
        self.ui_layoutCrystalGroupBox.addWidget(self.ui_CrystalMaterialComboBox,    2, 1)
        self.ui_layoutCrystalGroupBox.addWidget(self.ui_CrystalMaterialLabel,       2, 2)
        self.ui_layoutCrystalGroupBox.addWidget(self.ui_CrystalNXComboBox,          3, 1)
        self.ui_layoutCrystalGroupBox.addWidget(self.ui_CrystalNXLabel,             3, 2)
        self.ui_layoutCrystalGroupBox.addWidget(self.ui_CrystalNYComboBox,          4, 1)
        self.ui_layoutCrystalGroupBox.addWidget(self.ui_CrystalNYLabel,             4, 2)
        self.ui_layoutCrystalGroupBox.addWidget(self.ui_CrystalNZComboBox,          5, 1)
        self.ui_layoutCrystalGroupBox.addWidget(self.ui_CrystalNZLabel,             5, 2)

        for material in self.CrystalMaterials:
            self.ui_CrystalMaterialComboBox.addItem(material)

        for refidx in self.CurrentAvailableRefractiveIndices[0]:
            self.ui_CrystalNXComboBox.addItem(refidx)
        for refidx in self.CurrentAvailableRefractiveIndices[1]:
            self.ui_CrystalNYComboBox.addItem(refidx)
        for refidx in self.CurrentAvailableRefractiveIndices[2]:
            self.ui_CrystalNZComboBox.addItem(refidx)

        #TODO: Set selected item
        #self.ui_CrystalNXComboBox.setCurrentText()

        self.ui_CrystalGroupBox.setLayout(self.ui_layoutCrystalGroupBox)

        self.ui_layoutCrystal.addWidget(self.ui_CrystalGroupBox)

        return self.ui_layoutCrystal

    def initLayoutTest(self):
        self.ui_layoutTest = QGridLayout()
        self.ui_Testbtn1 = QPushButton()
        self.ui_Testbtn1.setText('Button 1')
        self.ui_Testbtn2 = QPushButton()
        self.ui_Testbtn2.setText('Button 2')
        self.ui_TestLabel1 = QLabel()
        self.ui_TestLabel1.setText('Label 1')
        self.ui_TestLabel2 = QLabel()
        self.ui_TestLabel2.setText('Label 2')

        self.ui_layoutTest.addWidget(self.ui_Testbtn1,      1, 2)
        self.ui_layoutTest.addWidget(self.ui_Testbtn2,      2, 2)
        self.ui_layoutTest.addWidget(self.ui_TestLabel1,    1, 1)
        self.ui_layoutTest.addWidget(self.ui_TestLabel2,    2, 1)

        return self.ui_layoutTest

    def getProperties(self):
        self.CrystalMaterials = RefractiveIndex().materialList
        self.currentCrystalMaterial = self.CrystalMaterials[0]#fallback
        lastMaterial = self.config.get('Crystal Material')
        if lastMaterial in self.CrystalMaterials:
            #print('Found Material')
            self.currentCrystalMaterial == lastMaterial

        self.CurrentAvailableRefractiveIndices = RefractiveIndex().getAvailableRefractiveIndices(self.currentCrystalMaterial)

        lastRefractiveIndexX = self.config.get('Crystal Refractive Index X')
        lastRefractiveIndexY = self.config.get('Crystal Refractive Index Y')
        lastRefractiveIndexZ = self.config.get('Crystal Refractive Index Z')

        print (self.CurrentAvailableRefractiveIndices)
        #self.findChild()
        #self.CrystalRefractiveIndices = RefractiveIndex().getAvailableRefractiveIndices()

    def showWindow(self):
        g = GUI()
        sys.exit(app.exec_())


if __name__ == '__main__':
    app = QApplication(sys.argv)
    g = GUI()
    sys.exit(app.exec_())