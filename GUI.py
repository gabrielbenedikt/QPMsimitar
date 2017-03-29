#!/usr/bin/python3
import sys
from PyQt5.QtCore import QObject, QVariant
from PyQt5.QtWidgets import (QApplication, QWidget, QToolTip)
from PyQt5.QtWidgets import (QPushButton, QSpinBox, QLabel, QDoubleSpinBox, QGroupBox, QComboBox, QCheckBox, QMainWindow)
from PyQt5.QtWidgets import (QGridLayout, QVBoxLayout, QHBoxLayout)
from PyQt5.QtGui import QFont
from RefractiveIndex import RefractiveIndex

class GUI(QMainWindow):
    def __init__(self, parent=None):
        super(GUI, self).__init__(parent)

        self.mainWidget = QWidget(self)
        self.setCentralWidget(self.mainWidget)

        self.initUI()

    def initUI(self):
        QToolTip.setFont(QFont('SansSerif', 10))
        self.setToolTip('Tooltip!')

        self.setWindowTitle('Multiphoton Analysis GUI')

        layout = QGridLayout()

        self.getProperties()

        self.layoutCrystal = self.initLayoutCrystal()
        self.layoutTest = self.initLayoutTest()
        self.layoutPump = self.initLayoutPump()
        self.layoutTemperature = self.initLayoutTemperature()

        layout.addLayout(self.layoutCrystal,     1, 1)
        layout.addLayout(self.layoutTest,        1, 2)
        layout.addLayout(self.layoutPump,        2, 1)
        layout.addLayout(self.layoutTemperature, 3, 1)

        self.centralWidget().setLayout(layout)
        self.resize(self.sizeHint())


        self.show()


    def initLayoutTemperature(self):
        layoutTemp = QGridLayout()
        TempGroupBox = QGroupBox()
        layoutTempGroupBox = QGridLayout()

        fromLabel = QLabel()
        toLabel = QLabel()
        singleLabel = QLabel()
        fromLabel.setText('from')
        toLabel.setText('to')
        singleLabel.setText('single')

        TsingleSB = QDoubleSpinBox()
        TfromSB = QDoubleSpinBox()
        TtoSB = QDoubleSpinBox()
        TLabel = QLabel()
        TLabel.setText('Temperature [°C]')

        layoutTempGroupBox.addWidget(fromLabel,     1, 2)
        layoutTempGroupBox.addWidget(singleLabel,   1, 3)
        layoutTempGroupBox.addWidget(toLabel,       1, 4)
        layoutTempGroupBox.addWidget(TLabel,        2, 1)
        layoutTempGroupBox.addWidget(TfromSB,       2, 2)
        layoutTempGroupBox.addWidget(TsingleSB,     2, 3)
        layoutTempGroupBox.addWidget(TtoSB,         2, 4)

        TempGroupBox.setLayout(layoutTempGroupBox)
        layoutTemp.addWidget(TempGroupBox)

        return layoutTemp

    def initLayoutPump(self):
        layoutPump = QGridLayout()
        PumpGroupBox = QGroupBox()
        layoutPumpGroupBox = QVBoxLayout()
        layoutUpper = QGridLayout()
        layoutLower = QGridLayout()

        fromLabel = QLabel()
        toLabel = QLabel()
        singleLabel = QLabel()
        fromLabel.setText('from')
        toLabel.setText('to')
        singleLabel.setText('single')

        pumpwlsingleSB = QDoubleSpinBox()
        pumpwlfromSB = QDoubleSpinBox()
        pumpwltoSB = QDoubleSpinBox()
        pumpwlLabel = QLabel()
        pumpwlLabel.setText('Wavelength [nm]')

        pulsewidthsingleSB = QDoubleSpinBox()
        pulsewidthfromSB = QDoubleSpinBox()
        pulsewidthtoSB = QDoubleSpinBox()
        pulsewidthLabel = QLabel()
        pulsewidthLabel.setText('Pulsewidth [ps]')

        pumpShapeCB = QComboBox()
        pumpShapeCB.addItem('Gaussian')
        pumpShapeCB.addItem('Sech^2')
        pumpShapeLabel = QLabel()
        pumpShapeLabel.setText('Pump pulse shape')

        pumpShapeCorrectionFactorCheckBox = QCheckBox()
        pumpShapeCorrectionFactorLabel = QLabel()
        pumpShapeCorrectionFactorLabel.setText('Apply deconvolution factor?')


        layoutUpper.addWidget(fromLabel,             1, 2)
        layoutUpper.addWidget(singleLabel,           1, 3)
        layoutUpper.addWidget(toLabel,               1, 4)
        layoutUpper.addWidget(pumpwlLabel,           2, 1)
        layoutUpper.addWidget(pumpwlfromSB,          2, 2)
        layoutUpper.addWidget(pumpwlsingleSB,        2, 3)
        layoutUpper.addWidget(pumpwltoSB,            2, 4)
        layoutUpper.addWidget(pulsewidthLabel,       3, 1)
        layoutUpper.addWidget(pulsewidthfromSB,      3, 2)
        layoutUpper.addWidget(pulsewidthsingleSB,    3, 3)
        layoutUpper.addWidget(pulsewidthtoSB,        3, 4)

        layoutLower.addWidget(pumpShapeLabel,                       1, 1)
        layoutLower.addWidget(pumpShapeCB,                          1, 2)
        layoutLower.addWidget(pumpShapeCorrectionFactorLabel,       2, 1)
        layoutLower.addWidget(pumpShapeCorrectionFactorCheckBox,    2, 2)

        layoutPumpGroupBox.addLayout(layoutUpper)
        layoutPumpGroupBox.addLayout(layoutLower)

        PumpGroupBox.setLayout(layoutPumpGroupBox)
        layoutPump.addWidget(PumpGroupBox)

        return layoutPump

    def initLayoutCrystal(self):
        # crystal section
        layoutCrystal = QGridLayout()
        CrystalGroupBox = QGroupBox()
        layoutCrystalGroupBox = QGridLayout()

        CrystalPolingPeriodSpinBox = QDoubleSpinBox()
        CrystalPolingPeriodLabel = QLabel()
        CrystalMaterialComboBox = QComboBox()
        CrystalMaterialLabel = QLabel()
        CrystalNXComboBox = QComboBox()
        CrystalNYComboBox = QComboBox()
        CrystalNZComboBox = QComboBox()
        CrystalNXLabel = QLabel()
        CrystalNYLabel = QLabel()
        CrystalNZLabel = QLabel()

        CrystalPolingPeriodLabel.setText('Poling Period [µm]')
        CrystalMaterialLabel.setText('Material')
        CrystalNXLabel.setText('Refractive index (x)')
        CrystalNYLabel.setText('Refractive index (y)')
        CrystalNZLabel.setText('Refractive index (z)')
        CrystalGroupBox.setTitle('Crystal')

        layoutCrystalGroupBox.addWidget(CrystalPolingPeriodSpinBox, 1, 1)
        layoutCrystalGroupBox.addWidget(CrystalPolingPeriodLabel,   1, 2)
        layoutCrystalGroupBox.addWidget(CrystalMaterialComboBox,    2, 1)
        layoutCrystalGroupBox.addWidget(CrystalMaterialLabel,       2, 2)
        layoutCrystalGroupBox.addWidget(CrystalNXComboBox,          3, 1)
        layoutCrystalGroupBox.addWidget(CrystalNXLabel,             3, 2)
        layoutCrystalGroupBox.addWidget(CrystalNYComboBox,          4, 1)
        layoutCrystalGroupBox.addWidget(CrystalNYLabel,             4, 2)
        layoutCrystalGroupBox.addWidget(CrystalNZComboBox,          5, 1)
        layoutCrystalGroupBox.addWidget(CrystalNZLabel,             5, 2)

        for material in self.CrystalMaterials:
            CrystalMaterialComboBox.addItem(material)

        CrystalGroupBox.setLayout(layoutCrystalGroupBox)

        layoutCrystal.addWidget(CrystalGroupBox)

        return layoutCrystal

    def initLayoutTest(self):
        layoutTest = QGridLayout()
        Testbtn1 = QPushButton()
        Testbtn1.setText('Button 1')
        Testbtn2 = QPushButton()
        Testbtn2.setText('Button 2')
        TestLabel1 = QLabel()
        TestLabel1.setText('Label 1')
        TestLabel2 = QLabel()
        TestLabel2.setText('Label 2')

        layoutTest.addWidget(Testbtn1,      1, 2)
        layoutTest.addWidget(Testbtn2,      2, 2)
        layoutTest.addWidget(TestLabel1,    1, 1)
        layoutTest.addWidget(TestLabel2,    2, 1)

        return layoutTest

    def getProperties(self):
        self.CrystalMaterials = RefractiveIndex().materialList
        #self.findChild()
        #self.CrystalRefractiveIndices = RefractiveIndex().getAvailableRefractiveIndices()

    def showWindow(self):
        g = GUI()
        sys.exit(app.exec_())


if __name__ == '__main__':
    app = QApplication(sys.argv)
    g = GUI()
    sys.exit(app.exec_())