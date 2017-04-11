#!/usr/bin/python3
import sys
from PyQt5.QtCore import QObject, QVariant
from PyQt5.QtWidgets import (QApplication, QWidget, QToolTip)
from PyQt5.QtWidgets import (QPushButton, QSpinBox, QLabel, QDoubleSpinBox, QGroupBox, QComboBox, QCheckBox, QMainWindow)
from PyQt5.QtWidgets import (QGridLayout, QVBoxLayout, QHBoxLayout, QScrollArea, QSizePolicy)
from PyQt5.QtGui import QFont
import PyQt5.QtCore
from RefractiveIndex import RefractiveIndex
import matplotlib.pyplot as plt
import numpy
import scipy

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
        self.ui_layoutPlotRefractiveIndex = self.initLayoutPlotRefractiveIndex()

        #Group: set crystal properties
        self.ui_layout.addLayout(self.ui_layoutCrystal,             1, 1)
        #Group: set pump properties
        self.ui_layout.addLayout(self.ui_layoutPump,                2, 1)
        #Group: set temperature
        self.ui_layout.addLayout(self.ui_layoutTemperature,         3, 1)
        #Group: plot refractive indices
        self.ui_layout.addLayout(self.ui_layoutPlotRefractiveIndex, 1, 2)

        self.centralWidget().setLayout(self.ui_layout)

        self.setProperties()

        self.initConnections()

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

        #set limits
        self.ui_TsingleSB.setRange(-273.15, 1000)
        self.ui_TfromSB.setRange(-273.15, 1000)
        self.ui_TtoSB.setRange(-273.15, 1000)

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

        self.ui_pumpwlsingleSB.setRange(0, 10000)
        self.ui_pumpwlfromSB.setRange(0, 10000)
        self.ui_pumpwltoSB.setRange(0, 10000)

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
        self.ui_PlotRefractiveIndex_Btn_Plot_wl.setText(r'Plot vs λ')
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

    def getProperties(self):
        self.CrystalMaterials = RefractiveIndex().materialList
        self.CrystalMaterial = self.CrystalMaterials[0]#fallback

        lastMaterial = self.config.get('Crystal Material')
        if lastMaterial in self.CrystalMaterials:
            self.CrystalMaterial == lastMaterial
        self.CurrentAvailableRefractiveIndices = RefractiveIndex().getAvailableRefractiveIndices(self.CrystalMaterial)


        self.CrystalPolingPeriod=self.config.get("Crystal Poling Period")


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
        self.PulsewidthFrom=self.config.get("Pump pulsewidth from")
        self.PulsewidthSingle=self.config.get("Pump pulsewidth single")
        self.PulsewidthTo=self.config.get("Pump pulsewidth to")
        self.PumpShape=self.config.get("Pump pulse shape")
        self.PumpShapeApplyDeconvolutionFactor=self.config.get("Pump pulsewidth apply deconvolution factor")
        self.TempFrom=self.config.get("Temperature from")
        self.TempSingle=self.config.get("Temperature single")
        self.TempTo=self.config.get("Temperature to")

    def setProperties(self):
        self.ui_CrystalPolingPeriodSpinBox.setValue(self.CrystalPolingPeriod)
        self.ui_pumpwlfromSB.setValue(self.PumpWlFrom)
        self.ui_pumpwlsingleSB.setValue(self.PumpWlSingle)
        self.ui_pumpwltoSB.setValue(self.PumpWlTo)
        self.ui_pulsewidthfromSB.setValue(self.PulsewidthFrom)
        self.ui_pulsewidthsingleSB.setValue(self.PulsewidthSingle)
        self.ui_pulsewidthtoSB.setValue(self.PulsewidthTo)
        idx=0#Fallback
        idx=self.ui_pumpShapeCB.findText(self.PumpShape)
        self.ui_pumpShapeCB.setCurrentIndex(idx)
        self.ui_pumpShapeCorrectionFactorCheckBox.setChecked(self.PumpShapeApplyDeconvolutionFactor)
        self.ui_TfromSB.setValue(self.TempFrom)
        self.ui_TsingleSB.setValue(self.TempSingle)
        self.ui_TtoSB.setValue(self.TempTo)

    def initConnections(self):
        self.ui_CrystalPolingPeriodSpinBox.valueChanged.connect(self.getVarsFromGUI)
        self.ui_CrystalMaterialComboBox.currentIndexChanged.connect(self.getVarsFromGUI)
        self.ui_CrystalNXComboBox.currentIndexChanged.connect(self.getVarsFromGUI)
        self.ui_CrystalNYComboBox.currentIndexChanged.connect(self.getVarsFromGUI)
        self.ui_CrystalNZComboBox.currentIndexChanged.connect(self.getVarsFromGUI)
        self.ui_pumpwlsingleSB.valueChanged.connect(self.getVarsFromGUI)
        self.ui_pumpwlfromSB.valueChanged.connect(self.getVarsFromGUI)
        self.ui_pumpwltoSB.valueChanged.connect(self.getVarsFromGUI)
        self.ui_pulsewidthsingleSB.valueChanged.connect(self.getVarsFromGUI)
        self.ui_pulsewidthfromSB.valueChanged.connect(self.getVarsFromGUI)
        self.ui_pulsewidthtoSB.valueChanged.connect(self.getVarsFromGUI)
        self.ui_pumpShapeCB.currentIndexChanged.connect(self.getVarsFromGUI)
        self.ui_pumpShapeCorrectionFactorCheckBox.stateChanged.connect(self.getVarsFromGUI)
        self.ui_TsingleSB.valueChanged.connect(self.getVarsFromGUI)
        self.ui_TfromSB.valueChanged.connect(self.getVarsFromGUI)
        self.ui_TtoSB.valueChanged.connect(self.getVarsFromGUI)
        self.ui_PlotRefractiveIndex_Btn_Plot_T.pressed.connect(self.plot_RefIdx_vs_T)
        self.ui_PlotRefractiveIndex_Btn_Plot_wl.pressed.connect(self.plot_RefIdx_vs_wl)

    def plot_RefIdx_vs_T(self):
        Tmin = self.ui_TfromSB.value()
        Tmax = self.ui_TtoSB.value()
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
                    print(RefIdxList)

        #NOTE: Let num (1000 here) be set in GUI
        plotrange=numpy.linspace(Tmin,Tmax,1000)

        for i in range(0,len(RefIdxList)):
            plt.plot(plotrange, RefIdxList[i](wl, plotrange),label='{0}, {1}, {2}'.format(materialList[i],polList[i],paperList[i]))

        plt.xlabel('Temperature [°C]')
        plt.ylabel('Refractive index')
        plt.title('Refractive indices vs Temperature')
        plt.annotate('Wavelength: {0:.2f}nm'.format(wl*10**9), xy=(0.55, 0.01), xycoords='axes fraction')
        plt.legend()
        plt.grid()
        plt.show()

    def plot_RefIdx_vs_wl(self):
        wlmin = self.ui_pumpwlfromSB.value()*10**(-9)
        wlmax = self.ui_pumpwltoSB.value()*10**(-9)
        T = self.ui_TsingleSB.value()
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
                    print(RefIdxList)

        # NOTE: Let num (1000 here) be set in GUI
        plotrange = numpy.linspace(wlmin, wlmax, 1000)

        for i in range(0, len(RefIdxList)):
            plt.plot(plotrange*10**9, RefIdxList[i](plotrange, T),
                     label='{0}, {1}, {2}'.format(materialList[i], polList[i], paperList[i]))

        plt.xlabel('Temperature [nm]')
        plt.ylabel('Refractive index')
        plt.title('Refractive indices vs Wavelength')
        plt.annotate('Temperature: {0:.1f}°C'.format(T), xy=(0.55, 0.01), xycoords='axes fraction')
        plt.legend()
        plt.grid()
        plt.show()

    def getVarsFromGUI(self):
        self.CrystalPolingPeriod = self.ui_CrystalPolingPeriodSpinBox.value()
        self.CrystalMaterial = self.ui_CrystalMaterialComboBox.currentText()
        self.CrystalNX = self.ui_CrystalNXComboBox.currentText()
        self.CrystalNY = self.ui_CrystalNYComboBox.currentText()
        self.CrystalNZ = self.ui_CrystalNZComboBox.currentText()
        self.PumpWlSingle = self.ui_pumpwlsingleSB.value()
        self.PumpWlFrom = self.ui_pumpwlfromSB.value()
        self.PumpWlTo = self.ui_pumpwltoSB.value()
        self.PulsewidthSingle = self.ui_pulsewidthsingleSB.value()
        self.PulsewidthFrom = self.ui_pulsewidthfromSB.value()
        self.PulsewidthTo = self.ui_pulsewidthtoSB.value()
        self.PumpShape = self.ui_pumpShapeCB.currentText()
        self.PumpShapeApplyDeconvolutionFactor = self.ui_pumpShapeCorrectionFactorCheckBox.isChecked()
        self.TempSingle = self.ui_TsingleSB.value()
        self.TempFrom = self.ui_TfromSB.value()
        self.TempTo = self.ui_TtoSB.value()

        self.SaveSettings()

    def SaveSettings(self):
        self.config.set("Crystal Material", self.CrystalMaterial)
        self.config.set("Crystal Poling Period", self.CrystalPolingPeriod)
        self.config.set("Crystal Refractive Index X", self.CrystalNX)
        self.config.set("Crystal Refractive Index Y", self.CrystalNY)
        self.config.set("Crystal Refractive Index Z", self.CrystalNZ)

        self.config.set("Pump wavelength from", self.PumpWlFrom)
        self.config.set("Pump wavelength single", self.PumpWlSingle)
        self.config.set("Pump wavelength to", self.PumpWlTo)
        self.config.set("Pump pulsewidth from", self.PulsewidthFrom)
        self.config.set("Pump pulsewidth single", self.PulsewidthSingle)
        self.config.set("Pump pulsewidth to", self.PulsewidthTo)
        self.config.set("Pump pulse shape", self.PumpShape)
        self.config.set("Pump pulsewidth apply deconvolution factor", self.PumpShapeApplyDeconvolutionFactor)

        self.config.set("Temperature from", self.TempFrom)
        self.config.set("Temperature single", self.TempSingle)
        self.config.set("Temperature to", self.TempTo)

    def showWindow(self):
        g = GUI()
        sys.exit(app.exec_())


if __name__ == '__main__':
    app = QApplication(sys.argv)
    g = GUI()
    sys.exit(app.exec_())