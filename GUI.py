#!/usr/bin/python3
import sys
from PyQt5.QtCore import QObject, QVariant
from PyQt5.QtWidgets import (QApplication, QWidget, QToolTip)
from PyQt5.QtWidgets import (QPushButton, QSpinBox, QLabel, QDoubleSpinBox, QGroupBox, QComboBox, QCheckBox, QMainWindow)
from PyQt5.QtWidgets import (QRadioButton)
from PyQt5.QtWidgets import (QGridLayout, QVBoxLayout, QHBoxLayout, QScrollArea, QSizePolicy)
from PyQt5.QtGui import QFont
from PyQt5 import QtCore, QtGui
import PyQt5.QtCore
from RefractiveIndex import RefractiveIndex
from PMC import PMC
import matplotlib.pyplot as plt
import numpy
import scipy
from pylab import *
from matplotlib.backends.backend_qt5agg import (FigureCanvasQTAgg as FigureCanvas, NavigationToolbar2QT as NavigationToolbar)

#from Multiphoton import MultiPhotonAnalysis
from Settings import Settings

class GUI(QMainWindow):
    def __init__(self, config, parent=None):
        super(GUI, self).__init__(parent)

        self.plotwindowcount=0
        self.pltwindowlist=[]

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
        self.ui_layoutPlotRefractiveIndex = self.initLayoutPlotRefractiveIndex()
        self.ui_layoutPlotPMC = self.initLayoutPlotPMC()
        self.ui_layoutPlotJSI = self.initLayoutPlotJSI()
        self.ui_layoutPurity = self.initLayoutPurity()

        #Group: set crystal properties
        self.ui_layout.addLayout(self.ui_layoutCrystal,             1, 1)
        #Group: set pump properties
        self.ui_layout.addLayout(self.ui_layoutPump,                2, 1)
        #Group: plot refractive indices
        self.ui_layout.addLayout(self.ui_layoutPlotRefractiveIndex, 1, 2)
        #Group: plot phase matching curve
        self.ui_layout.addLayout(self.ui_layoutPlotPMC,             2, 2)
        #Group: plot JSI
        self.ui_layout.addLayout(self.ui_layoutPlotJSI,             1, 3)
        #Group: Purity
        self.ui_layout.addLayout(self.ui_layoutPurity,              2, 3)

        self.centralWidget().setLayout(self.ui_layout)

        self.setProperties()

        self.initConnections()

        self.resize(self.sizeHint())

        self.show()

    def initLayoutPump(self):
        self.ui_layoutPump = QGridLayout()
        self.ui_PumpGroupBox = QGroupBox()
        self.ui_layoutPumpGroupBox = QVBoxLayout()
        self.ui_layoutUpper = QGridLayout()
        self.ui_layoutLower = QGridLayout()

        self.ui_PumpGroupBox.setTitle('Pump parameters')

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
        self.ui_pumpwlLabel.setText('Pump Wavelength [nm]')

        self.ui_pumpwlsingleSB.setRange(0, 10000)
        self.ui_pumpwlfromSB.setRange(0, 10000)
        self.ui_pumpwltoSB.setRange(0, 10000)

        self.ui_siwlsingleSB = QDoubleSpinBox()
        self.ui_siwlfromSB = QDoubleSpinBox()
        self.ui_siwltoSB = QDoubleSpinBox()
        self.ui_siwlLabel = QLabel()
        self.ui_siwlLabel.setText('S/I Wavelength [nm]')

        self.ui_siwlsingleSB.setRange(0, 10000)
        self.ui_siwlfromSB.setRange(0, 10000)
        self.ui_siwltoSB.setRange(0, 10000)

        self.ui_pulsewidthsingleSB = QDoubleSpinBox()
        self.ui_pulsewidthfromSB = QDoubleSpinBox()
        self.ui_pulsewidthtoSB = QDoubleSpinBox()
        self.ui_pulsewidthLabel = QLabel()
        self.ui_pulsewidthLabel.setText('Pulsewidth [ps]')

        self.ui_pulsewidthsingleSB.setRange(0, 1000)
        self.ui_pulsewidthfromSB.setRange(0, 1000)
        self.ui_pulsewidthtoSB.setRange(0, 1000)

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
        self.ui_layoutUpper.addWidget(self.ui_siwlLabel,             3, 1)
        self.ui_layoutUpper.addWidget(self.ui_siwlfromSB,            3, 2)
        self.ui_layoutUpper.addWidget(self.ui_siwlsingleSB,          3, 3)
        self.ui_layoutUpper.addWidget(self.ui_siwltoSB,              3, 4)
        self.ui_layoutUpper.addWidget(self.ui_pulsewidthLabel,       4, 1)
        self.ui_layoutUpper.addWidget(self.ui_pulsewidthfromSB,      4, 2)
        self.ui_layoutUpper.addWidget(self.ui_pulsewidthsingleSB,    4, 3)
        self.ui_layoutUpper.addWidget(self.ui_pulsewidthtoSB,        4, 4)

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
        self.ui_sublayoutCrystal = QGridLayout()

        self.ui_CrystalMaterialComboBox = QComboBox()
        self.ui_CrystalMaterialLabel = QLabel()
        self.ui_CrystalNXComboBox = QComboBox()
        self.ui_CrystalNYComboBox = QComboBox()
        self.ui_CrystalNZComboBox = QComboBox()
        self.ui_CrystalNXLabel = QLabel()
        self.ui_CrystalNYLabel = QLabel()
        self.ui_CrystalNZLabel = QLabel()


        self.ui_CrystalMaterialLabel.setText('Material')
        self.ui_CrystalNXLabel.setText('Refractive index (x)')
        self.ui_CrystalNYLabel.setText('Refractive index (y)')
        self.ui_CrystalNZLabel.setText('Refractive index (z)')
        self.ui_CrystalGroupBox.setTitle('Crystal parameters')

        #from/single/to sublayout
        self.ui_CrystalfromLabel = QLabel()
        self.ui_CrystaltoLabel = QLabel()
        self.ui_CrystalsingleLabel = QLabel()
        self.ui_CrystalfromLabel.setText('from')
        self.ui_CrystaltoLabel.setText('to')
        self.ui_CrystalsingleLabel.setText('single')

        self.ui_CrystalTsingleSB = QDoubleSpinBox()
        self.ui_CrystalTfromSB = QDoubleSpinBox()
        self.ui_CrystalTtoSB = QDoubleSpinBox()
        self.ui_CrystalTLabel = QLabel()
        self.ui_CrystalTLabel.setText('Temperature [°C]')

        self.ui_CrystalPolingPeriodsingleSB = QDoubleSpinBox()
        self.ui_CrystalPolingPeriodfromSB = QDoubleSpinBox()
        self.ui_CrystalPolingPeriodtoSB = QDoubleSpinBox()
        self.ui_CrystalPolingPeriodLabel = QLabel()
        self.ui_CrystalPolingPeriodLabel.setText('Poling Period [µm]')

        self.ui_CrystalTsingleSB.setRange(-273.15, 1000)
        self.ui_CrystalTfromSB.setRange(-273.15, 1000)
        self.ui_CrystalTtoSB.setRange(-273.15, 1000)

        self.ui_sublayoutCrystal.addWidget(self.ui_CrystalfromLabel,                1, 2)
        self.ui_sublayoutCrystal.addWidget(self.ui_CrystalsingleLabel,              1, 3)
        self.ui_sublayoutCrystal.addWidget(self.ui_CrystaltoLabel,                  1, 4)
        self.ui_sublayoutCrystal.addWidget(self.ui_CrystalTLabel,                   2, 1)
        self.ui_sublayoutCrystal.addWidget(self.ui_CrystalTfromSB,                  2, 2)
        self.ui_sublayoutCrystal.addWidget(self.ui_CrystalTsingleSB,                2, 3)
        self.ui_sublayoutCrystal.addWidget(self.ui_CrystalTtoSB,                    2, 4)
        self.ui_sublayoutCrystal.addWidget(self.ui_CrystalPolingPeriodLabel,        3, 1)
        self.ui_sublayoutCrystal.addWidget(self.ui_CrystalPolingPeriodfromSB,       3, 2)
        self.ui_sublayoutCrystal.addWidget(self.ui_CrystalPolingPeriodsingleSB,     3, 3)
        self.ui_sublayoutCrystal.addWidget(self.ui_CrystalPolingPeriodtoSB,         3, 4)


        #main  layout
        self.ui_layoutCrystalGroupBox.addWidget(self.ui_CrystalMaterialComboBox,    2, 1)
        self.ui_layoutCrystalGroupBox.addWidget(self.ui_CrystalMaterialLabel,       2, 2)
        self.ui_layoutCrystalGroupBox.addWidget(self.ui_CrystalNXComboBox,          3, 1)
        self.ui_layoutCrystalGroupBox.addWidget(self.ui_CrystalNXLabel,             3, 2)
        self.ui_layoutCrystalGroupBox.addWidget(self.ui_CrystalNYComboBox,          4, 1)
        self.ui_layoutCrystalGroupBox.addWidget(self.ui_CrystalNYLabel,             4, 2)
        self.ui_layoutCrystalGroupBox.addWidget(self.ui_CrystalNZComboBox,          5, 1)
        self.ui_layoutCrystalGroupBox.addWidget(self.ui_CrystalNZLabel,             5, 2)
        self.ui_layoutCrystalGroupBox.addLayout(self.ui_sublayoutCrystal,           6, 1, -1, 2)

        for material in self.CrystalMaterials:
            self.ui_CrystalMaterialComboBox.addItem(material)

        for refidx in self.CurrentAvailableRefractiveIndices[0]:
            self.ui_CrystalNXComboBox.addItem(refidx,refidx)
        for refidx in self.CurrentAvailableRefractiveIndices[1]:
            self.ui_CrystalNYComboBox.addItem(refidx,refidx)
        for refidx in self.CurrentAvailableRefractiveIndices[2]:
            self.ui_CrystalNZComboBox.addItem(refidx,refidx)

        #TODO: Set selected item
        #self.ui_CrystalNXComboBox.setCurrentText()
        idx = self.ui_CrystalNXComboBox.findData(self.lastRefractiveIndex[0])
        if (idx != -1):
            self.ui_CrystalNXComboBox.setCurrentIndex(idx)
        idx = self.ui_CrystalNYComboBox.findData(self.lastRefractiveIndex[1])
        if (idx != -1):
            self.ui_CrystalNYComboBox.setCurrentIndex(idx)
        idx = self.ui_CrystalNZComboBox.findData(self.lastRefractiveIndex[2])
        if (idx != -1):
            self.ui_CrystalNZComboBox.setCurrentIndex(idx)

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

    def initLayoutPlotRefractiveIndex(self):

        # |Grid Layout
        # |->GroupBox
        # ||->Grid Layout
        # |||->Scroll Area(layout?)
        # ||||->QWidget
        # |||||->Grid layout
        # ||||||->Checkboxes
        # |||->Button
        self.ui_layoutPlotRefractiveIndex = QGridLayout()
        self.ui_layoutPlotRefractiveIndexGroupBox = QGroupBox()
        self.ui_layoutPlotRefractiveIndexGroupBoxLayout = QGridLayout()
        self.ui_layoutPlotRefractiveIndexScrollArea = QScrollArea()
        self.ui_layoutPlotRefractiveIndexScrollAreaWidget = QWidget()
        self.ui_layoutPlotRefractiveIndexScrollAreaGrid = QGridLayout()

        self.ui_layoutPlotRefractiveIndexGroupBox.setTitle('Plot refractive indices')

        #TODO
        #vertical scroll area
        #containing checkboxes with all available refractive indices

        self.ui_PlotRefractiveIndex_Btn_Plot_T = QPushButton()
        self.ui_PlotRefractiveIndex_Btn_Plot_T.setText('Plot vs T')
        self.ui_PlotRefractiveIndex_Btn_Plot_wl = QPushButton()
        self.ui_PlotRefractiveIndex_Btn_Plot_wl.setText(r'Plot vs λp')
        #self.ui_PlotRefractiveIndex_Btn_Plot_wl.setText(r'Plot vs '+QChar(0xbb, 0x03))

        self.ui_layoutPlotRefractiveIndexScrollArea.setWidget(self.ui_layoutPlotRefractiveIndexScrollAreaWidget)
        self.ui_layoutPlotRefractiveIndexScrollAreaWidget.setLayout(self.ui_layoutPlotRefractiveIndexScrollAreaGrid)
        self.ui_layoutPlotRefractiveIndexScrollArea.setWidgetResizable(True)

        self.ui_layoutPlotRefractiveIndexGroupBoxLayout.addWidget(self.ui_layoutPlotRefractiveIndexScrollArea,  1, 1, 1, 2)
        self.ui_layoutPlotRefractiveIndexGroupBoxLayout.addWidget(self.ui_PlotRefractiveIndex_Btn_Plot_wl,      2, 1)
        self.ui_layoutPlotRefractiveIndexGroupBoxLayout.addWidget(self.ui_PlotRefractiveIndex_Btn_Plot_T,       2, 2)

        self.ui_layoutPlotRefractiveIndexGroupBox.setLayout(self.ui_layoutPlotRefractiveIndexGroupBoxLayout)

        self.ui_layoutPlotRefractiveIndex.addWidget(self.ui_layoutPlotRefractiveIndexGroupBox)

        #fill widget with checkboxes containing all available refractive indices
        i=0
        for material in self.CrystalMaterials:
            indices = RefractiveIndex().getAvailableRefractiveIndices(material)
            print(material, indices)
            for index in indices[0]:
                CBstring = '{0}:X:{1}'.format(material, index)
                self.ui_layoutPlotRefractiveIndexScrollAreaGrid.addWidget(QCheckBox(CBstring), i, 1)
                i=i+1
            for index in indices[1]:
                CBstring = '{0}:Y:{1}'.format(material, index)
                self.ui_layoutPlotRefractiveIndexScrollAreaGrid.addWidget(QCheckBox(CBstring), i, 1)
                i=i+1
            for index in indices[2]:
                CBstring = '{0}:Z:{1}'.format(material, index)
                self.ui_layoutPlotRefractiveIndexScrollAreaGrid.addWidget(QCheckBox(CBstring), i, 1)
                i=i+1
            self.ui_layoutPlotRefractiveIndexScrollArea.sizeHint()

        return self.ui_layoutPlotRefractiveIndex

    def initLayoutPlotPMC(self):
        self.ui_layoutPlotPMC = QGridLayout()
        self.ui_PlotPMCGroupBox= QGroupBox()
        self.ui_layoutPlotPMCGroupBox = QGridLayout()

        self.ui_PlotPMCGroupBox.setTitle('Plot PMC')

        self.ui_PlotPMClabelQPMorder = QLabel()
        self.ui_PlotPMClabelQPMorder.setText('QPM order')

        self.ui_PlotPMCSBQPMorder = QSpinBox()

        self.ui_PlotPMCvsT_Btn = QPushButton()
        self.ui_PlotPMCvsT_Btn.setText('Plot vs T')

        self.ui_PlotPMCvsPP_Btn = QPushButton()
        self.ui_PlotPMCvsPP_Btn.setText('Plot vs PP')

        # set limits
        self.ui_PlotPMCSBQPMorder.setRange(-10000, 10000)

        self.ui_layoutPlotPMCGroupBox.addWidget(self.ui_PlotPMClabelQPMorder,   1, 1)
        self.ui_layoutPlotPMCGroupBox.addWidget(self.ui_PlotPMCSBQPMorder,      1, 2)
        self.ui_layoutPlotPMCGroupBox.addWidget(self.ui_PlotPMCvsT_Btn,         2, 1, 1, 2)
        self.ui_layoutPlotPMCGroupBox.addWidget(self.ui_PlotPMCvsPP_Btn,        3, 1, 1, 2)

        self.ui_PlotPMCGroupBox.setLayout(self.ui_layoutPlotPMCGroupBox)
        self.ui_layoutPlotPMC.addWidget(self.ui_PlotPMCGroupBox)

        return self.ui_layoutPlotPMC

    def initLayoutPlotJSI(self):
        self.ui_layoutPlotJSI = QGridLayout()
        self.ui_PlotJSIGroupBox = QGroupBox()
        self.ui_layoutPlotJSIGroupBox = QGridLayout()

        self.ui_PlotJSIGroupBox.setTitle('TODO: Plot JSI/JSA')

        self.ui_PlotJSIorJSALayout = QGridLayout()
        self.ui_PlotJSI_plotJSIRadioButton = QRadioButton('Plot JSI')
        self.ui_PlotJSI_plotJSARadioButton = QRadioButton('Plot JSA')
        self.ui_PlotJSIorJSALayout.addWidget(self.ui_PlotJSI_plotJSIRadioButton, 1, 1)
        self.ui_PlotJSIorJSALayout.addWidget(self.ui_PlotJSI_plotJSARadioButton, 1, 2)
        self.ui_PlotJSI_plotJSIRadioButton.setChecked(True)

        self.ui_PlotJSI_plotBtn = QPushButton('TODO: Plot')

        ## Add controls
        self.ui_layoutPlotJSIGroupBox.addLayout(self.ui_PlotJSIorJSALayout,     1, 1)
        self.ui_layoutPlotJSIGroupBox.addWidget(self.ui_PlotJSI_plotBtn,        2, 1)

        self.ui_PlotJSIGroupBox.setLayout(self.ui_layoutPlotJSIGroupBox)
        self.ui_layoutPlotJSI.addWidget(self.ui_PlotJSIGroupBox)

        return self.ui_layoutPlotJSI

    def initLayoutPurity(self):
        self.ui_layoutPurity = QGridLayout()
        self.ui_PurityGroupBox = QGroupBox()
        self.ui_layoutPurityGroupBox = QGridLayout()

        self.ui_PurityGroupBox.setTitle('TODO: Purity')

        self.ui_Purity_plotvsTau_Btn = QPushButton('TODO: Plot vs τ')
        self.ui_Purity_plotvspwl_Btn = QPushButton('TODO: Plot vs λp')
        self.ui_Purity_plotvsL_Btn = QPushButton('TODO: Plot vs L')
        self.ui_Purity_plotvsTauandL_Btn = QPushButton('TODO: Plot vs L and τ')

        ## Add controls
        self.ui_layoutPurityGroupBox.addWidget(self.ui_Purity_plotvsTau_Btn,        1, 1, 1, 2)
        self.ui_layoutPurityGroupBox.addWidget(self.ui_Purity_plotvspwl_Btn,        2, 1, 1, 2)
        self.ui_layoutPurityGroupBox.addWidget(self.ui_Purity_plotvsL_Btn,          3, 1, 1, 2)
        self.ui_layoutPurityGroupBox.addWidget(self.ui_Purity_plotvsTauandL_Btn,    4, 1, 1, 2)

        self.ui_PurityGroupBox.setLayout(self.ui_layoutPurityGroupBox)
        self.ui_layoutPurity.addWidget(self.ui_PurityGroupBox)

        return self.ui_layoutPurity

    def getProperties(self):
        self.CrystalMaterials = RefractiveIndex().materialList
        self.CrystalMaterial = self.CrystalMaterials[0]#fallback

        lastMaterial = self.config.get('Crystal Material')
        if lastMaterial in self.CrystalMaterials:
            self.CrystalMaterial == lastMaterial
        self.CurrentAvailableRefractiveIndices = RefractiveIndex().getAvailableRefractiveIndices(self.CrystalMaterial)

        self.CrystalPolingPeriodFrom = self.config.get("Crystal Poling Period From")
        self.CrystalPolingPeriodSingle = self.config.get("Crystal Poling Period Single")
        self.CrystalPolingPeriodTo = self.config.get("Crystal Poling Period To")

        lastNX = self.config.get("Crystal Refractive Index X")
        if lastNX in self.CurrentAvailableRefractiveIndices[0]:
            self.CrystalNX = lastNX
        else:
            self.CrystalNX = self.CurrentAvailableRefractiveIndices[0]
        lastNY = self.config.get("Crystal Refractive Index Y")
        if lastNY in self.CurrentAvailableRefractiveIndices[1]:
            self.CrystalNY = lastNY
        else:
            self.CrystalNY = self.CurrentAvailableRefractiveIndices[1]
        lastNZ = self.config.get("Crystal Refractive Index Z")
        if lastNZ in self.CurrentAvailableRefractiveIndices[2]:
            self.CrystalNZ = lastNZ
        else:
            self.CrystalNZ = self.CurrentAvailableRefractiveIndices[2]

        self.lastRefractiveIndex = [self.CrystalNX,
                                    self.CrystalNY,
                                    self.CrystalNZ]
        self.PumpWlFrom=self.config.get("Pump wavelength from")
        self.PumpWlSingle=self.config.get("Pump wavelength single")
        self.PumpWlTo=self.config.get("Pump wavelength to")
        self.SIWlFrom = self.config.get("SI wavelength from")
        self.SIWlSingle = self.config.get("SI wavelength single")
        self.SIWlTo = self.config.get("SI wavelength to")
        self.PulsewidthFrom=self.config.get("Pump pulsewidth from")
        self.PulsewidthSingle=self.config.get("Pump pulsewidth single")
        self.PulsewidthTo=self.config.get("Pump pulsewidth to")
        self.PumpShape=self.config.get("Pump pulse shape")
        self.PumpShapeApplyDeconvolutionFactor=self.config.get("Pump pulsewidth apply deconvolution factor")
        self.CrystalTempFrom=self.config.get("Crystal Temperature from")
        self.CrystalTempSingle=self.config.get("Crystal Temperature single")
        self.CrystalTempTo=self.config.get("Crystal Temperature to")
        self.QPMOrder=self.config.get("QPM Order")

    def setProperties(self):
        self.ui_CrystalPolingPeriodsingleSB.setValue(self.CrystalPolingPeriodSingle*10**6)
        self.ui_CrystalPolingPeriodfromSB.setValue(self.CrystalPolingPeriodFrom * 10 ** 6)
        self.ui_CrystalPolingPeriodtoSB.setValue(self.CrystalPolingPeriodTo * 10 ** 6)
        self.ui_pumpwlfromSB.setValue(self.PumpWlFrom*10**9)
        self.ui_pumpwlsingleSB.setValue(self.PumpWlSingle*10**9)
        self.ui_pumpwltoSB.setValue(self.PumpWlTo*10**9)
        self.ui_siwlfromSB.setValue(self.SIWlFrom*10**9)
        self.ui_siwlsingleSB.setValue(self.SIWlSingle*10**9)
        self.ui_siwltoSB.setValue(self.SIWlTo*10**9)
        self.ui_pulsewidthfromSB.setValue(self.PulsewidthFrom*10**12)
        self.ui_pulsewidthsingleSB.setValue(self.PulsewidthSingle*10**12)
        self.ui_pulsewidthtoSB.setValue(self.PulsewidthTo*10**12)
        idx=0#Fallback
        idx=self.ui_pumpShapeCB.findText(self.PumpShape)
        self.ui_pumpShapeCB.setCurrentIndex(idx)
        self.ui_pumpShapeCorrectionFactorCheckBox.setChecked(self.PumpShapeApplyDeconvolutionFactor)
        self.ui_CrystalTfromSB.setValue(self.CrystalTempFrom)
        self.ui_CrystalTsingleSB.setValue(self.CrystalTempSingle)
        self.ui_CrystalTtoSB.setValue(self.CrystalTempTo)
        self.ui_PlotPMCSBQPMorder.setValue(self.QPMOrder)

    def initConnections(self):
        self.ui_CrystalPolingPeriodsingleSB.valueChanged.connect(self.getVarsFromGUI)
        self.ui_CrystalPolingPeriodfromSB.valueChanged.connect(self.getVarsFromGUI)
        self.ui_CrystalPolingPeriodtoSB.valueChanged.connect(self.getVarsFromGUI)
        self.ui_CrystalMaterialComboBox.currentIndexChanged.connect(self.getVarsFromGUI)
        self.ui_CrystalNXComboBox.currentIndexChanged.connect(self.getVarsFromGUI)
        self.ui_CrystalNYComboBox.currentIndexChanged.connect(self.getVarsFromGUI)
        self.ui_CrystalNZComboBox.currentIndexChanged.connect(self.getVarsFromGUI)
        self.ui_pumpwlsingleSB.valueChanged.connect(self.getVarsFromGUI)
        self.ui_pumpwlfromSB.valueChanged.connect(self.getVarsFromGUI)
        self.ui_pumpwltoSB.valueChanged.connect(self.getVarsFromGUI)
        self.ui_siwlfromSB.valueChanged.connect(self.getVarsFromGUI)
        self.ui_siwlsingleSB.valueChanged.connect(self.getVarsFromGUI)
        self.ui_siwltoSB.valueChanged.connect(self.getVarsFromGUI)
        self.ui_pulsewidthsingleSB.valueChanged.connect(self.getVarsFromGUI)
        self.ui_pulsewidthfromSB.valueChanged.connect(self.getVarsFromGUI)
        self.ui_pulsewidthtoSB.valueChanged.connect(self.getVarsFromGUI)
        self.ui_pumpShapeCB.currentIndexChanged.connect(self.getVarsFromGUI)
        self.ui_pumpShapeCorrectionFactorCheckBox.stateChanged.connect(self.getVarsFromGUI)
        self.ui_CrystalTsingleSB.valueChanged.connect(self.getVarsFromGUI)
        self.ui_CrystalTfromSB.valueChanged.connect(self.getVarsFromGUI)
        self.ui_CrystalTtoSB.valueChanged.connect(self.getVarsFromGUI)
        self.ui_PlotPMCSBQPMorder.valueChanged.connect(self.getVarsFromGUI)
        self.ui_PlotRefractiveIndex_Btn_Plot_T.pressed.connect(self.plot_RefIdx_vs_T)
        self.ui_PlotRefractiveIndex_Btn_Plot_wl.pressed.connect(self.plot_RefIdx_vs_wl)
        self.ui_PlotPMCvsT_Btn.pressed.connect(self.plot_pmc_wl_vs_T)
        self.ui_PlotPMCvsPP_Btn.pressed.connect(self.plot_pmc_wl_vs_PP)
        self.ui_Purity_plotvsTau_Btn.pressed.connect(self.plot_purity_vs_tau)
        self.ui_Purity_plotvspwl_Btn.pressed.connect(self.plot_purity_vs_pwl)
        self.ui_PlotJSI_plotBtn.pressed.connect(self.plot_jsi)
        self.ui_Purity_plotvsTauandL_Btn.pressed.connect(self.plot_purity_vs_Tau_and_L)
        self.ui_Purity_plotvsL_Btn.pressed.connect(self.plot_purity_vs_L)

    def plot_RefIdx_vs_T(self):
        pltwndidx=self.plotwindowcount
        self.open_new_plot_window()
        pltwnd=self.pltwindowlist[pltwndidx]
        
        Tmin = self.ui_CrystalTfromSB.value()
        Tmax = self.ui_CrystalTtoSB.value()
        wl = self.ui_pumpwlsingleSB.value()*10**(-9)
        RefIdxList=[]
        materialList=[]
        polList=[]
        paperList=[]
        for cb in self.ui_layoutPlotRefractiveIndexScrollAreaWidget.children():
            if isinstance(cb, QCheckBox):
                if cb.isChecked():
                    [material, pol, paper]=cb.text().replace('&','').split(':')
                    materialList.append(material)
                    polList.append(pol)
                    paperList.append(paper)
                    RefIdxList.append(RefractiveIndex().getSingleIDX(material, pol, paper))

        #NOTE: Let num (1000 here) be set in GUI
        plotrange=numpy.linspace(Tmin,Tmax,1000)

        for i in range(0,len(RefIdxList)):
            pltwnd.ax.plot(plotrange, RefIdxList[i](wl, plotrange),label='{0}, {1}, {2}'.format(materialList[i],polList[i],paperList[i]))

        pltwnd.ax.set_xlabel('Temperature [°C]')
        pltwnd.ax.set_ylabel('Refractive index')
        pltwnd.ax.set_title('Refractive indices vs Temperature')
        pltwnd.ax.annotate('Wavelength: {0:.2f}nm'.format(wl*10**9), xy=(0.55, 0.01), xycoords='axes fraction')
        pltwnd.ax.legend()
        pltwnd.canvas.draw()

    def plot_RefIdx_vs_wl(self):
        pltwndidx=self.plotwindowcount
        self.open_new_plot_window()
        pltwnd=self.pltwindowlist[pltwndidx]
        
        wlmin = self.ui_pumpwlfromSB.value()*10**(-9)
        wlmax = self.ui_pumpwltoSB.value()*10**(-9)
        T = self.ui_CrystalTsingleSB.value()
        RefIdxList = []
        materialList = []
        polList = []
        paperList = []
        for cb in self.ui_layoutPlotRefractiveIndexScrollAreaWidget.children():
            if isinstance(cb, QCheckBox):
                if cb.isChecked():
                    [material, pol, paper] = cb.text().replace('&', '').split(':')
                    materialList.append(material)
                    polList.append(pol)
                    paperList.append(paper)
                    RefIdxList.append(RefractiveIndex().getSingleIDX(material, pol, paper))

        # NOTE: Let num (1000 here) be set in GUI
        plotrange = numpy.linspace(wlmin, wlmax, 1000)

        for i in range(0, len(RefIdxList)):
            pltwnd.ax.plot(plotrange*10**9, RefIdxList[i](plotrange, T),
                     label='{0}, {1}, {2}'.format(materialList[i], polList[i], paperList[i]))

        pltwnd.ax.set_xlabel('Wavelength [nm]')
        pltwnd.ax.set_ylabel('Refractive index')
        pltwnd.ax.set_title('Refractive indices vs Wavelength')
        pltwnd.ax.annotate('Temperature: {0:.1f}°C'.format(T), xy=(0.55, 0.01), xycoords='axes fraction')
        pltwnd.ax.legend()
        pltwnd.canvas.draw()

    def plot_pmc_wl_vs_T(self):
        #init plot window
        pltwndidx=self.plotwindowcount
        self.open_new_plot_window()
        pltwnd=self.pltwindowlist[pltwndidx]

        #get variables from GUI
        lp = self.PumpWlSingle
        PP = self.CrystalPolingPeriodSingle
        Tmin = self.CrystalTempFrom
        Tmax = self.CrystalTempTo
        m = self.QPMOrder

        # get Ref indices
        nxfunc=RefractiveIndex().getSingleIDX(self.CrystalMaterial, "X", self.CrystalNX)
        nyfunc=RefractiveIndex().getSingleIDX(self.CrystalMaterial, "Y", self.CrystalNY)
        nzfunc=RefractiveIndex().getSingleIDX(self.CrystalMaterial, "Z", self.CrystalNZ)
        refidxfunc=[nxfunc,nyfunc,nzfunc]

        #prepare plotting
        plotrange = numpy.arange(Tmin,Tmax, (Tmax-Tmin)/250)

        #get PMC from PMC class
        pmc=PMC()

        #TODO: get refractive indices functions and give them to pmc.getSI_wl
        [siwl,idwl,Tcp]=pmc.getSI_wl_varT(lp,PP,plotrange,refidxfunc,m)

        #plot
        pltwnd.ax.plot(plotrange, siwl*10**9, lw=2, label='signal')
        pltwnd.ax.plot(plotrange, idwl*10**9, lw=2, label='idler')
        pltwnd.ax.set_xlabel('Temperature [°C]')
        pltwnd.ax.set_ylabel('S/I Wavelength [nm]')
        pltwnd.ax.set_title('Quasi Phase Matching Temperature Tuning Curve')
        pltwnd.ax.annotate('Crossing point temperature: {0:.1f}°C'.format(Tcp),xy=(0.55, 0.01), xycoords='axes fraction')
        pltwnd.ax.legend()
        pltwnd.canvas.draw()

    def plot_pmc_wl_vs_PP(self):
        #init plot window
        pltwndidx=self.plotwindowcount
        self.open_new_plot_window()
        pltwnd=self.pltwindowlist[pltwndidx]

        #get variables from GUI
        lp = self.PumpWlSingle
        PPmin = self.CrystalPolingPeriodFrom
        PPmax = self.CrystalPolingPeriodTo
        m = self.QPMOrder
        T = self.CrystalTempSingle

        # get Ref indices
        nxfunc=RefractiveIndex().getSingleIDX(self.CrystalMaterial, "X", self.CrystalNX)
        nyfunc=RefractiveIndex().getSingleIDX(self.CrystalMaterial, "Y", self.CrystalNY)
        nzfunc=RefractiveIndex().getSingleIDX(self.CrystalMaterial, "Z", self.CrystalNZ)
        refidxfunc=[nxfunc,nyfunc,nzfunc]

        #prepare plotting
        plotrange = numpy.arange(PPmin, PPmax, (PPmax-PPmin)/250)

        #get PMC from PMC class
        pmc=PMC()

        #TODO: get refractive indices functions and give them to pmc.getSI_wl
        [siwl,idwl,PPcp]=pmc.getSI_wl_varPP(lp,plotrange,T,refidxfunc,m)

        #plot
        pltwnd.ax.plot(plotrange*10**6, siwl*10**9, lw=2, label='signal')
        pltwnd.ax.plot(plotrange*10**6, idwl*10**9, lw=2, label='idler')
        pltwnd.ax.set_xlabel('Poling Period [µm]')
        pltwnd.ax.set_ylabel('S/I Wavelength [nm]')
        pltwnd.ax.set_title('Quasi Phase Matching Poling Period Tuning Curve')
        pltwnd.ax.annotate('Degenerate wavelengths (for {0:.1f}°C) at a poling period of {1:.2f} µm'.format(self.CrystalTempSingle,PPcp*10**6),xy=(0.01, 0.01), xycoords='axes fraction')
        pltwnd.ax.legend()
        pltwnd.canvas.draw()

    def plot_purity_vs_tau(self):
        pass

    def plot_purity_vs_pwl(self):
        pass

    def plot_purity_vs_L(self):
        pass

    def plot_purity_vs_Tau_and_L(self):
        pass

    def plot_jsi(self):
        pass

    def getVarsFromGUI(self):
        self.CrystalPolingPeriodSingle = self.ui_CrystalPolingPeriodsingleSB.value() * 10 ** (-6)
        self.CrystalPolingPeriodFrom = self.ui_CrystalPolingPeriodfromSB.value() * 10 ** (-6)
        self.CrystalPolingPeriodTo = self.ui_CrystalPolingPeriodtoSB.value() * 10 ** (-6)
        self.CrystalMaterial = self.ui_CrystalMaterialComboBox.currentText()
        self.CrystalNX = self.ui_CrystalNXComboBox.currentText()
        self.CrystalNY = self.ui_CrystalNYComboBox.currentText()
        self.CrystalNZ = self.ui_CrystalNZComboBox.currentText()
        self.CrystalTempSingle = self.ui_CrystalTsingleSB.value()
        self.CrystalTempFrom = self.ui_CrystalTfromSB.value()
        self.CrystalTempTo = self.ui_CrystalTtoSB.value()

        self.PumpWlSingle = self.ui_pumpwlsingleSB.value()*10**(-9)
        self.PumpWlFrom = self.ui_pumpwlfromSB.value()*10**(-9)
        self.PumpWlTo = self.ui_pumpwltoSB.value()*10**(-9)
        self.PulsewidthSingle = self.ui_pulsewidthsingleSB.value()*10**(-12)
        self.PulsewidthFrom = self.ui_pulsewidthfromSB.value()*10**(-12)
        self.PulsewidthTo = self.ui_pulsewidthtoSB.value()*10**(-12)
        self.PumpShape = self.ui_pumpShapeCB.currentText()
        self.PumpShapeApplyDeconvolutionFactor = self.ui_pumpShapeCorrectionFactorCheckBox.isChecked()

        self.SIWlFrom = self.ui_siwlfromSB.value()*10**(-9)
        self.SIWlSingle = self.ui_siwlsingleSB.value()*10**(-9)
        self.SIWlTo = self.ui_siwltoSB.value()*10**(-9)

        self.QPMOrder = self.ui_PlotPMCSBQPMorder.value()

        self.SaveSettings()

    def SaveSettings(self):
        self.config.set("Crystal Material", self.CrystalMaterial)
        self.config.set("Crystal Poling Period Single", self.CrystalPolingPeriodSingle)
        self.config.set("Crystal Poling Period From", self.CrystalPolingPeriodFrom)
        self.config.set("Crystal Poling Period To", self.CrystalPolingPeriodTo)
        self.config.set("Crystal Refractive Index X", self.CrystalNX)
        self.config.set("Crystal Refractive Index Y", self.CrystalNY)
        self.config.set("Crystal Refractive Index Z", self.CrystalNZ)

        self.config.set("QPM Order", self.QPMOrder)

        self.config.set("Pump wavelength from", self.PumpWlFrom)
        self.config.set("Pump wavelength single", self.PumpWlSingle)
        self.config.set("Pump wavelength to", self.PumpWlTo)

        self.config.set("SI wavelength from", self.SIWlFrom)
        self.config.set("SI wavelength single", self.SIWlSingle)
        self.config.set("SI wavelength to", self.SIWlTo)

        self.config.set("Pump pulsewidth from", self.PulsewidthFrom)
        self.config.set("Pump pulsewidth single", self.PulsewidthSingle)
        self.config.set("Pump pulsewidth to", self.PulsewidthTo)
        self.config.set("Pump pulse shape", self.PumpShape)
        self.config.set("Pump pulsewidth apply deconvolution factor", self.PumpShapeApplyDeconvolutionFactor)

        self.config.set("Crystal Temperature from", self.CrystalTempFrom)
        self.config.set("Crystal Temperature single", self.CrystalTempSingle)
        self.config.set("Crystal Temperature to", self.CrystalTempTo)

    def showWindow(self):
        g = GUI()
        sys.exit(app.exec_())

    def open_new_plot_window(self):
        print(self.plotwindowcount)
        self.pltwindowlist.append(PlotWindow())
        self.pltwindowlist[self.plotwindowcount].show()
        self.plotwindowcount = self.plotwindowcount+1


class PlotWindow(QWidget):
    def __init__(self):
        QWidget.__init__(self)
        #self.layout = QGridLayout()
        self.layout = QVBoxLayout()
        self.setLayout(self.layout)
        self.fig = figure(facecolor="white",figsize=(8.75*1.2, 5*1.2))
        self.ax = self.fig.add_subplot(111)
        self.ax.grid()
        self.canvas = FigureCanvas(self.fig)
        self.canvas.setParent(self)
        self.toolbar = NavigationToolbar(self.canvas, self)
        #self.addToolBar(self.toolbar)
        self.layout.addWidget(self.canvas)
        self.layout.addWidget(self.toolbar)

    def open(self):

        self.show(self)


if __name__ == '__main__':
    app = QApplication(sys.argv)
    g = GUI()
    sys.exit(app.exec_())
