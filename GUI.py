#!/usr/bin/python3
from PyQt5.QtWidgets import (QApplication, QWidget, QToolTip, QPushButton, QSpinBox, QLabel,
                             QDoubleSpinBox, QGroupBox, QComboBox, QCheckBox, QMainWindow,
                             QRadioButton, QGridLayout, QVBoxLayout, QScrollArea, QMessageBox)
from PyQt5.QtGui import (QFont, QMouseEvent)
from RefractiveIndex import RefractiveIndex
from PMC import PMC
from JSI import JSI
from Filters import Filters
from Constants import Constants
from QTreimps import QHoverPushButton
import numpy
from pylab import *
from matplotlib.backends.backend_qt5agg import (FigureCanvasQTAgg as FigureCanvas,
                                                NavigationToolbar2QT as NavigationToolbar)

# noinspection PyAttributeOutsideInit
class GUI(QMainWindow):
    def __init__(self, config, parent=None):
        super(GUI, self).__init__(parent)

        self.plotwindowcount = 0
        self.pltwindowlist = []

        self.mainWidget = QWidget(self)
        self.setCentralWidget(self.mainWidget)

        self.config = config

        self.initUI()

        self.HighlightedSpinBox = "QSpinBox { background-color: alternate-background-color; }"
        self.HighlightedDoubleSpinBox = "QDoubleSpinBox { background-color: alternate-background-color; }"
        self.HighlightedComboBox = "QComboBox { background-color: alternate-background-color; }"
        self.HighlightedRadioButton = "QRadioButton { background-color: alternate-background-color; }"
        self.HighlightedWidget = "QWidget { background-color: alternate-background-color; }"

    def initUI(self):
        QToolTip.setFont(QFont('SansSerif', 10))
        self.setToolTip('Tooltip!')

        self.setWindowTitle('Multiphoton Analysis GUI')

        self.ui_layout = QGridLayout()

        self.getProperties()

        self.ui_layoutCrystal = self.initLayoutCrystal()
        self.ui_layoutPump = self.initLayoutPump()
        self.ui_layoutSIfilter = self.initLayoutSIfilter()
        self.ui_layoutPlotRefractiveIndex = self.initLayoutPlotRefractiveIndex()
        self.ui_layoutPlotPMC = self.initLayoutPlotPMC()
        self.ui_layoutPlotJSI = self.initLayoutPlotJSI()
        self.ui_layoutPurity = self.initLayoutPurity()
        self.ui_layoutGetEffPP = self.initLayoutGetEffPP()
        self.ui_layoutTcp = self.initLayoutTcp()
        self.ui_layoutHOM = self.initLayoutHOM()
        self.ui_layoutPlotFWHM = self.initLayoutPlotFWHM()

        # Group: set crystal properties
        self.ui_layout.addLayout(self.ui_layoutCrystal, 1, 1)
        # Group: set pump properties
        self.ui_layout.addLayout(self.ui_layoutPump, 2, 1)
        # Group: S/I spectral filtering
        self.ui_layout.addLayout(self.ui_layoutSIfilter, 3, 1)
        # Group: plot refractive indices
        self.ui_layout.addLayout(self.ui_layoutPlotRefractiveIndex, 1, 2)
        # Group: plot phase matching curve
        self.ui_layout.addLayout(self.ui_layoutPlotPMC, 2, 2)
        # Group: Get effective poling period
        self.ui_layout.addLayout(self.ui_layoutGetEffPP, 3, 2)
        # Group: Purity
        self.ui_layout.addLayout(self.ui_layoutPurity, 1, 3)
        # Group: plot JSI
        self.ui_layout.addLayout(self.ui_layoutPlotJSI, 2, 3)
        # Group: Crossing point temperature vs different variables
        self.ui_layout.addLayout(self.ui_layoutTcp, 3, 3)
        # Group: HOM
        self.ui_layout.addLayout(self.ui_layoutHOM, 2, 4)
        #Group: Plot FWHM of marginal spectral
        self.ui_layout.addLayout(self.ui_layoutPlotFWHM, 3, 4)

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
        self.ui_pumpShapeCorrectionFactorLabel.setText('TODO: Apply deconvolution factor?')
        # TODO: use deconvolution factor if set

        self.ui_layoutUpper.addWidget(self.ui_fromLabel, 1, 2)
        self.ui_layoutUpper.addWidget(self.ui_singleLabel, 1, 3)
        self.ui_layoutUpper.addWidget(self.ui_toLabel, 1, 4)
        self.ui_layoutUpper.addWidget(self.ui_pumpwlLabel, 2, 1)
        self.ui_layoutUpper.addWidget(self.ui_pumpwlfromSB, 2, 2)
        self.ui_layoutUpper.addWidget(self.ui_pumpwlsingleSB, 2, 3)
        self.ui_layoutUpper.addWidget(self.ui_pumpwltoSB, 2, 4)
        self.ui_layoutUpper.addWidget(self.ui_siwlLabel, 3, 1)
        self.ui_layoutUpper.addWidget(self.ui_siwlfromSB, 3, 2)
        self.ui_layoutUpper.addWidget(self.ui_siwlsingleSB, 3, 3)
        self.ui_layoutUpper.addWidget(self.ui_siwltoSB, 3, 4)
        self.ui_layoutUpper.addWidget(self.ui_pulsewidthLabel, 4, 1)
        self.ui_layoutUpper.addWidget(self.ui_pulsewidthfromSB, 4, 2)
        self.ui_layoutUpper.addWidget(self.ui_pulsewidthsingleSB, 4, 3)
        self.ui_layoutUpper.addWidget(self.ui_pulsewidthtoSB, 4, 4)

        self.ui_layoutLower.addWidget(self.ui_pumpShapeLabel, 1, 1)
        self.ui_layoutLower.addWidget(self.ui_pumpShapeCB, 1, 2)
        self.ui_layoutLower.addWidget(self.ui_pumpShapeCorrectionFactorLabel, 2, 1)
        self.ui_layoutLower.addWidget(self.ui_pumpShapeCorrectionFactorCheckBox, 2, 2)

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

        # from/single/to sublayout
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

        self.ui_CrystalTsingleSB.setRange(-273.15, 1000)
        self.ui_CrystalTfromSB.setRange(-273.15, 1000)
        self.ui_CrystalTtoSB.setRange(-273.15, 1000)

        self.ui_CrystalPolingPeriodsingleSB = QDoubleSpinBox()
        self.ui_CrystalPolingPeriodfromSB = QDoubleSpinBox()
        self.ui_CrystalPolingPeriodtoSB = QDoubleSpinBox()
        self.ui_CrystalPolingPeriodLabel = QLabel()
        self.ui_CrystalPolingPeriodLabel.setText('Poling Period [µm]')
        self.ui_CrystalPolingPeriodsingleSB.setDecimals(4)
        self.ui_CrystalPolingPeriodfromSB.setDecimals(4)
        self.ui_CrystalPolingPeriodtoSB.setDecimals(4)

        self.ui_CrystalLengthsingleSB = QDoubleSpinBox()
        self.ui_CrystalLengthfromSB = QDoubleSpinBox()
        self.ui_CrystalLengthtoSB = QDoubleSpinBox()
        self.ui_CrystalLengthLabel = QLabel()
        self.ui_CrystalLengthLabel.setText('Crystal length [mm]')


        self.ui_sublayoutCrystal.addWidget(self.ui_CrystalfromLabel,            1, 2)
        self.ui_sublayoutCrystal.addWidget(self.ui_CrystalsingleLabel,          1, 3)
        self.ui_sublayoutCrystal.addWidget(self.ui_CrystaltoLabel,              1, 4)
        self.ui_sublayoutCrystal.addWidget(self.ui_CrystalTLabel,               2, 1)
        self.ui_sublayoutCrystal.addWidget(self.ui_CrystalTfromSB,              2, 2)
        self.ui_sublayoutCrystal.addWidget(self.ui_CrystalTsingleSB,            2, 3)
        self.ui_sublayoutCrystal.addWidget(self.ui_CrystalTtoSB,                2, 4)
        self.ui_sublayoutCrystal.addWidget(self.ui_CrystalPolingPeriodLabel,    3, 1)
        self.ui_sublayoutCrystal.addWidget(self.ui_CrystalPolingPeriodfromSB,   3, 2)
        self.ui_sublayoutCrystal.addWidget(self.ui_CrystalPolingPeriodsingleSB, 3, 3)
        self.ui_sublayoutCrystal.addWidget(self.ui_CrystalPolingPeriodtoSB,     3, 4)
        self.ui_sublayoutCrystal.addWidget(self.ui_CrystalLengthLabel,          4, 1)
        self.ui_sublayoutCrystal.addWidget(self.ui_CrystalLengthfromSB,         4, 2)
        self.ui_sublayoutCrystal.addWidget(self.ui_CrystalLengthsingleSB,       4, 3)
        self.ui_sublayoutCrystal.addWidget(self.ui_CrystalLengthtoSB,           4, 4)

        # main  layout
        self.ui_layoutCrystalGroupBox.addWidget(self.ui_CrystalMaterialComboBox, 2, 1)
        self.ui_layoutCrystalGroupBox.addWidget(self.ui_CrystalMaterialLabel, 2, 2)
        self.ui_layoutCrystalGroupBox.addWidget(self.ui_CrystalNXComboBox, 3, 1)
        self.ui_layoutCrystalGroupBox.addWidget(self.ui_CrystalNXLabel, 3, 2)
        self.ui_layoutCrystalGroupBox.addWidget(self.ui_CrystalNYComboBox, 4, 1)
        self.ui_layoutCrystalGroupBox.addWidget(self.ui_CrystalNYLabel, 4, 2)
        self.ui_layoutCrystalGroupBox.addWidget(self.ui_CrystalNZComboBox, 5, 1)
        self.ui_layoutCrystalGroupBox.addWidget(self.ui_CrystalNZLabel, 5, 2)
        self.ui_layoutCrystalGroupBox.addLayout(self.ui_sublayoutCrystal, 6, 1, -1, 2)

        for material in self.CrystalMaterials:
            self.ui_CrystalMaterialComboBox.addItem(material)

        for refidx in self.CurrentAvailableRefractiveIndices[0]:
            self.ui_CrystalNXComboBox.addItem(refidx, refidx)
        for refidx in self.CurrentAvailableRefractiveIndices[1]:
            self.ui_CrystalNYComboBox.addItem(refidx, refidx)
        for refidx in self.CurrentAvailableRefractiveIndices[2]:
            self.ui_CrystalNZComboBox.addItem(refidx, refidx)

        # TODO: Set selected item
        # self.ui_CrystalNXComboBox.setCurrentText()
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

    def initLayoutHOM(self):
        self.ui_layoutHOM = QGridLayout()
        self.ui_HOMGroupBox = QGroupBox()
        self.ui_layoutHOMGroupBox = QGridLayout()
        self.ui_HOMGroupBox.setTitle('HOM')

        self.ui_HOM_Vis_DelayRange_SB = QDoubleSpinBox()
        self.ui_HOM_Vis_DelayRange_Label = QLabel('Delay range [ps]')

        self.ui_HOM_Vis_Resolution_SB = QSpinBox()
        self.ui_HOM_Vis_Resolution_Label = QLabel('Resolution')
        self.ui_HOM_Vis_Resolution_SB.setRange(7,10000)

        self.ui_HOM_PlotVis_Btn = QHoverPushButton()
        self.ui_HOM_PlotVis_Btn.setText('PlotHOMVis')
        self.ui_HOM_PlotVis_Btn.setObjectName('Plot HOM Visibility')

        self.ui_layoutHOMGroupBox.addWidget(self.ui_HOM_Vis_Resolution_Label, 1, 1)
        self.ui_layoutHOMGroupBox.addWidget(self.ui_HOM_Vis_Resolution_SB, 1, 2)
        self.ui_layoutHOMGroupBox.addWidget(self.ui_HOM_Vis_DelayRange_Label, 2, 1)
        self.ui_layoutHOMGroupBox.addWidget(self.ui_HOM_Vis_DelayRange_SB, 2, 2)
        self.ui_layoutHOMGroupBox.addWidget(self.ui_HOM_PlotVis_Btn, 3, 1, 1, -1)

        self.ui_HOMGroupBox.setLayout(self.ui_layoutHOMGroupBox)
        self.ui_layoutHOM.addWidget(self.ui_HOMGroupBox)

        return self.ui_layoutHOM

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

        # TODO
        # vertical scroll area
        # containing checkboxes with all available refractive indices

        self.ui_PlotRefractiveIndex_Btn_Plot_T = QHoverPushButton()
        self.ui_PlotRefractiveIndex_Btn_Plot_T.setText('Plot vs T')
        self.ui_PlotRefractiveIndex_Btn_Plot_T.setObjectName('Plot refractive indices vs T')
        self.ui_PlotRefractiveIndex_Btn_Plot_wl = QHoverPushButton()
        self.ui_PlotRefractiveIndex_Btn_Plot_wl.setText(r'Plot vs λp')
        self.ui_PlotRefractiveIndex_Btn_Plot_wl.setObjectName('Plot refractive indices vs wl')
        # self.ui_PlotRefractiveIndex_Btn_Plot_wl.setText(r'Plot vs '+QChar(0xbb, 0x03))

        self.ui_layoutPlotRefractiveIndexScrollArea.setWidget(self.ui_layoutPlotRefractiveIndexScrollAreaWidget)
        self.ui_layoutPlotRefractiveIndexScrollAreaWidget.setLayout(self.ui_layoutPlotRefractiveIndexScrollAreaGrid)
        self.ui_layoutPlotRefractiveIndexScrollArea.setWidgetResizable(True)

        self.ui_layoutPlotRefractiveIndexGroupBoxLayout.addWidget(self.ui_layoutPlotRefractiveIndexScrollArea, 1, 1, 1,
                                                                  2)
        self.ui_layoutPlotRefractiveIndexGroupBoxLayout.addWidget(self.ui_PlotRefractiveIndex_Btn_Plot_wl, 2, 1)
        self.ui_layoutPlotRefractiveIndexGroupBoxLayout.addWidget(self.ui_PlotRefractiveIndex_Btn_Plot_T, 2, 2)

        self.ui_layoutPlotRefractiveIndexGroupBox.setLayout(self.ui_layoutPlotRefractiveIndexGroupBoxLayout)

        self.ui_layoutPlotRefractiveIndex.addWidget(self.ui_layoutPlotRefractiveIndexGroupBox)

        # fill widget with checkboxes containing all available refractive indices
        i = 0
        for material in self.CrystalMaterials:
            indices = RefractiveIndex().getAvailableRefractiveIndices(material)
            print(material, indices)
            for index in indices[0]:
                CBstring = '{0}:X:{1}'.format(material, index)
                self.ui_layoutPlotRefractiveIndexScrollAreaGrid.addWidget(QCheckBox(CBstring), i, 1)
                i = i + 1
            for index in indices[1]:
                CBstring = '{0}:Y:{1}'.format(material, index)
                self.ui_layoutPlotRefractiveIndexScrollAreaGrid.addWidget(QCheckBox(CBstring), i, 1)
                i = i + 1
            for index in indices[2]:
                CBstring = '{0}:Z:{1}'.format(material, index)
                self.ui_layoutPlotRefractiveIndexScrollAreaGrid.addWidget(QCheckBox(CBstring), i, 1)
                i = i + 1
            self.ui_layoutPlotRefractiveIndexScrollArea.sizeHint()

        return self.ui_layoutPlotRefractiveIndex

    def initLayoutPlotPMC(self):
        self.ui_layoutPlotPMC = QGridLayout()
        self.ui_PlotPMCGroupBox = QGroupBox()
        self.ui_layoutPlotPMCGroupBox = QGridLayout()

        self.ui_PlotPMCGroupBox.setTitle('Plot PMC')

        self.ui_PlotPMClabelQPMorder = QLabel()
        self.ui_PlotPMClabelQPMorder.setText('QPM order')

        self.ui_PlotPMCSBQPMorder = QSpinBox()

        self.ui_PlotPMCvsT_Btn = QHoverPushButton()
        self.ui_PlotPMCvsT_Btn.setText('Plot vs T')
        self.ui_PlotPMCvsT_Btn.setObjectName('Plot PMC vs T')

        self.ui_PlotPMCvsPP_Btn = QHoverPushButton()
        self.ui_PlotPMCvsPP_Btn.setText('Plot vs PP')
        self.ui_PlotPMCvsPP_Btn.setObjectName('Plot PMC vs PP')

        # set limits
        self.ui_PlotPMCSBQPMorder.setRange(-10000, 10000)

        self.ui_layoutPlotPMCGroupBox.addWidget(self.ui_PlotPMClabelQPMorder, 1, 1)
        self.ui_layoutPlotPMCGroupBox.addWidget(self.ui_PlotPMCSBQPMorder, 1, 2)
        self.ui_layoutPlotPMCGroupBox.addWidget(self.ui_PlotPMCvsT_Btn, 2, 1, 1, 2)
        self.ui_layoutPlotPMCGroupBox.addWidget(self.ui_PlotPMCvsPP_Btn, 3, 1, 1, 2)

        self.ui_PlotPMCGroupBox.setLayout(self.ui_layoutPlotPMCGroupBox)
        self.ui_layoutPlotPMC.addWidget(self.ui_PlotPMCGroupBox)

        return self.ui_layoutPlotPMC

    def initLayoutSIfilter(self):
        self.ui_layoutSIfilter = QGridLayout()
        self.ui_SIfilterGroupBox = QGroupBox()
        self.ui_layoutSIfilterGroupBox = QVBoxLayout()
        self.ui_layoutSIfilterUpper = QGridLayout()

        self.ui_SIfilterGroupBox.setTitle('S/I spectral filtering')

        self.ui_SIfilterSignalLabel = QLabel('signal')
        self.ui_SIfilterIdlerLabel = QLabel('idler')
        self.ui_SIfilterTypeLabel = QLabel('Type')
        self.ui_SIfiltercenterWLLabel = QLabel('Center WL [nm]')
        self.ui_SIfilterBandwidthLabel = QLabel('FWHM [nm]')

        self.ui_SIfilterSignalCenterWL_SB = QDoubleSpinBox()
        self.ui_SIfilterIdlerCenterWL_SB = QDoubleSpinBox()
        self.ui_SIfilterSignalFWHM_SB = QDoubleSpinBox()
        self.ui_SIfilterIdlerFWHM_SB = QDoubleSpinBox()
        self.ui_SIfilterSignalCenterWL_SB.setRange(0,10000)
        self.ui_SIfilterIdlerCenterWL_SB.setRange(0,10000)
        self.ui_SIfilterSignalFWHM_SB.setRange(0,10000)
        self.ui_SIfilterIdlerFWHM_SB.setRange(0,10000)

        self.ui_SIfilterSignalType_CB = QComboBox()
        self.ui_SIfilterIdlerType_CB = QComboBox()

        for Type in self.SIfilterTypes:
            self.ui_SIfilterSignalType_CB.addItem(Type)
            self.ui_SIfilterIdlerType_CB.addItem(Type)

        self.ui_layoutSIfilterUpper.addWidget(self.ui_SIfilterTypeLabel,            1, 2)
        self.ui_layoutSIfilterUpper.addWidget(self.ui_SIfiltercenterWLLabel,        1, 3)
        self.ui_layoutSIfilterUpper.addWidget(self.ui_SIfilterBandwidthLabel,       1, 4)
        self.ui_layoutSIfilterUpper.addWidget(self.ui_SIfilterSignalLabel,          2, 1)
        self.ui_layoutSIfilterUpper.addWidget(self.ui_SIfilterSignalType_CB,        2, 2)
        self.ui_layoutSIfilterUpper.addWidget(self.ui_SIfilterSignalCenterWL_SB,    2, 3)
        self.ui_layoutSIfilterUpper.addWidget(self.ui_SIfilterSignalFWHM_SB,        2, 4)
        self.ui_layoutSIfilterUpper.addWidget(self.ui_SIfilterIdlerLabel,           3, 1)
        self.ui_layoutSIfilterUpper.addWidget(self.ui_SIfilterIdlerType_CB,         3, 2)
        self.ui_layoutSIfilterUpper.addWidget(self.ui_SIfilterIdlerCenterWL_SB,     3, 3)
        self.ui_layoutSIfilterUpper.addWidget(self.ui_SIfilterIdlerFWHM_SB,         3, 4)

        self.ui_layoutSIfilterGroupBox.addLayout(self.ui_layoutSIfilterUpper)

        self.ui_SIfilterGroupBox.setLayout(self.ui_layoutSIfilterGroupBox)
        self.ui_layoutSIfilter.addWidget(self.ui_SIfilterGroupBox)

        return self.ui_layoutSIfilter

    def initLayoutPlotJSI(self):
        self.ui_layoutPlotJSI = QGridLayout()
        self.ui_PlotJSIGroupBox = QGroupBox()
        self.ui_layoutPlotJSIGroupBox = QGridLayout()

        self.ui_PlotJSIGroupBox.setTitle('Plot JSI/JSA')

        self.ui_PlotJSI_Wlrange_SB = QDoubleSpinBox()
        self.ui_PlotJSI_Wlrange_label = QLabel('Wl range [nm]')

        self.ui_PlotJSIorJSALayout = QGridLayout()
        self.ui_PlotJSI_plotJSIRadioButton = QRadioButton('JSI')
        self.ui_PlotJSI_plotJSARadioButton = QRadioButton('JSA')
        self.ui_PlotJSIorJSALayout.addWidget(self.ui_PlotJSI_plotJSIRadioButton, 1, 1)
        self.ui_PlotJSIorJSALayout.addWidget(self.ui_PlotJSI_plotJSARadioButton, 2, 1)
        self.ui_PlotJSI_plotJSIRadioButton.setChecked(True)
        
        self.ui_PlotJSI_filterlossBtn = QHoverPushButton('TODO: Estimate filter losses')
        self.ui_PlotJSI_filterlossBtn.setObjectName('Esimate filter losses')

        self.ui_PlotJSI_plotBtn = QHoverPushButton('Plot')
        self.ui_PlotJSI_plotBtn.setObjectName('Plot JSA JSI')
        self.ui_PlotJSI_WLresolution_SB = QSpinBox()
        self.ui_PlotJSI_WLresolution_Label = QLabel('Resolution')
        self.ui_PlotJSI_WLresolution_SB.setMinimum(2)
        self.ui_PlotJSI_WLresolution_SB.setMaximum(10000)

        ## Add controls
        self.ui_layoutPlotJSIGroupBox.addLayout(self.ui_PlotJSIorJSALayout, 1, 1)
        self.ui_layoutPlotJSIGroupBox.addWidget(self.ui_PlotJSI_plotBtn, 1, 2)
        self.ui_layoutPlotJSIGroupBox.addWidget(self.ui_PlotJSI_filterlossBtn, 2, 1, 1, -1)
        self.ui_layoutPlotJSIGroupBox.addWidget(self.ui_PlotJSI_Wlrange_label, 3, 1)
        self.ui_layoutPlotJSIGroupBox.addWidget(self.ui_PlotJSI_Wlrange_SB, 3, 2)
        self.ui_layoutPlotJSIGroupBox.addWidget(self.ui_PlotJSI_WLresolution_Label,4,1)
        self.ui_layoutPlotJSIGroupBox.addWidget(self.ui_PlotJSI_WLresolution_SB,4,2)

        self.ui_PlotJSIGroupBox.setLayout(self.ui_layoutPlotJSIGroupBox)
        self.ui_layoutPlotJSI.addWidget(self.ui_PlotJSIGroupBox)

        return self.ui_layoutPlotJSI

    def initLayoutPurity(self):
        self.ui_layoutPurity = QGridLayout()
        self.ui_PurityGroupBox = QGroupBox()
        self.ui_layoutPurityGroupBox = QGridLayout()

        self.ui_PurityGroupBox.setTitle('Purity')

        self.ui_Purity_plotvsTau_Btn = QHoverPushButton('Plot vs τ')
        self.ui_Purity_plotvspwl_Btn = QHoverPushButton('TODO: Plot vs λp')
        self.ui_Purity_plotvsL_Btn = QHoverPushButton('Plot vs L')
        self.ui_Purity_plotvsTauandL_Btn = QHoverPushButton('Plot vs L and τ')
        self.ui_Purity_plotvsTau_Btn.setObjectName('Plot Purity vs tau')
        self.ui_Purity_plotvspwl_Btn.setObjectName('Plot Purity vs lp')
        self.ui_Purity_plotvsL_Btn.setObjectName('Plot Purity vs L')
        self.ui_Purity_plotvsTauandL_Btn.setObjectName('Plot Purity vs L and tau')

        self.ui_Purity_WLresolution_Label = QLabel('WL resolution')
        self.ui_Purity_WLrange_Label = QLabel('WL range [nm]')
        self.ui_Purity_Tauresolution_Label = QLabel('τ or L resolution')

        self.ui_Purity_WLresolution_SB = QSpinBox()
        self.ui_Purity_WLrange_SB = QSpinBox()
        self.ui_Purity_Tauresolution_SB = QSpinBox()

        self.ui_Purity_WLresolution_SB.setMinimum(2)
        self.ui_Purity_WLresolution_SB.setMaximum(10000)
        self.ui_Purity_WLrange_SB.setMinimum(2)
        self.ui_Purity_WLrange_SB.setMaximum(10000)
        self.ui_Purity_Tauresolution_SB.setMinimum(2)
        self.ui_Purity_Tauresolution_SB.setMaximum(10000)

        ## Add controls
        self.ui_layoutPurityGroupBox.addWidget(self.ui_Purity_Tauresolution_Label, 1, 1)
        self.ui_layoutPurityGroupBox.addWidget(self.ui_Purity_Tauresolution_SB, 1, 2)
        self.ui_layoutPurityGroupBox.addWidget(self.ui_Purity_WLresolution_Label, 2, 1)
        self.ui_layoutPurityGroupBox.addWidget(self.ui_Purity_WLresolution_SB, 2, 2)
        self.ui_layoutPurityGroupBox.addWidget(self.ui_Purity_WLrange_Label, 3, 1)
        self.ui_layoutPurityGroupBox.addWidget(self.ui_Purity_WLrange_SB, 3, 2)
        self.ui_layoutPurityGroupBox.addWidget(self.ui_Purity_plotvsTau_Btn, 4, 1, 1, 2)
        self.ui_layoutPurityGroupBox.addWidget(self.ui_Purity_plotvspwl_Btn, 5, 1, 1, 2)
        self.ui_layoutPurityGroupBox.addWidget(self.ui_Purity_plotvsL_Btn, 6, 1, 1, 2)
        self.ui_layoutPurityGroupBox.addWidget(self.ui_Purity_plotvsTauandL_Btn, 7, 1, 1, 2)

        self.ui_PurityGroupBox.setLayout(self.ui_layoutPurityGroupBox)
        self.ui_layoutPurity.addWidget(self.ui_PurityGroupBox)

        return self.ui_layoutPurity

    def initLayoutGetEffPP(self):
        self.ui_layoutGetEffPP = QGridLayout()
        self.ui_GetEffPPGroupBox = QGroupBox()
        self.ui_layoutGetEffPPGroupBox = QGridLayout()

        self.ui_GetEffPPGroupBox.setTitle('Get effective poling period')

        self.ui_GetEffPP_Btn = QHoverPushButton()
        self.ui_GetEffPP_Btn.setText('Get effective PP')
        self.ui_GetEffPP_Btn.setObjectName('Get effective PP')

        self.ui_layoutGetEffPPGroupBox.addWidget(self.ui_GetEffPP_Btn, 1, 1)

        self.ui_GetEffPPGroupBox.setLayout(self.ui_layoutGetEffPPGroupBox)
        self.ui_layoutGetEffPP.addWidget(self.ui_GetEffPPGroupBox)

        return self.ui_layoutGetEffPP

    def initLayoutPlotFWHM(self):
        self.ui_layoutPlotFWHM = QGridLayout()
        self.ui_PlotFWHMGroupBox = QGroupBox()
        self.ui_layoutPlotFWHMGroupBox = QGridLayout()

        self.ui_PlotFWHMGroupBox.setTitle('Plot marginal FWHM')

        self.ui_PlotFWHMvstau_Btn = QHoverPushButton()
        self.ui_PlotFWHMvstau_Btn.setText('Plot vs τ')
        self.ui_PlotFWHMvstau_Btn.setObjectName('Plot marginal FWHM vs tau')
        
        self.ui_PlotFWHMresolution_Label = QLabel('Resolution')
        self.ui_PlotFWHMresolution_SB = QSpinBox()
        self.ui_PlotFWHMresolution_SB.setRange(2,10000)
        
        self.ui_PlotFWHMprecision_Label = QLabel('Precision')
        self.ui_PlotFWHMprecision_SB = QSpinBox()
        self.ui_PlotFWHMprecision_SB.setRange(1,10)

        self.ui_layoutPlotFWHMGroupBox.addWidget(self.ui_PlotFWHMresolution_Label, 1, 1)
        self.ui_layoutPlotFWHMGroupBox.addWidget(self.ui_PlotFWHMresolution_SB, 1, 2)
        self.ui_layoutPlotFWHMGroupBox.addWidget(self.ui_PlotFWHMprecision_Label, 2, 1)
        self.ui_layoutPlotFWHMGroupBox.addWidget(self.ui_PlotFWHMprecision_SB, 2, 2)
        self.ui_layoutPlotFWHMGroupBox.addWidget(self.ui_PlotFWHMvstau_Btn, 3, 1, 1, -1)

        self.ui_PlotFWHMGroupBox.setLayout(self.ui_layoutPlotFWHMGroupBox)
        self.ui_layoutPlotFWHM.addWidget(self.ui_PlotFWHMGroupBox)

        return self.ui_layoutPlotFWHM

    def initLayoutTcp(self):
        self.ui_layoutTcp = QGridLayout()
        self.ui_TcpGroupBox = QGroupBox()
        self.ui_layoutTcpGroupBox = QGridLayout()

        self.ui_TcpGroupBox.setTitle('Crossing point temperature')

        self.ui_Tcp_vsPP_Btn = QHoverPushButton()
        self.ui_Tcp_vsPP_Btn.setText('Plot vs PP')
        self.ui_Tcp_vsPP_Btn.setObjectName('Plot Tcp vs PP')

        self.ui_Tcp_vslp_Btn = QHoverPushButton()
        self.ui_Tcp_vslp_Btn.setText('Plot vs λp')
        self.ui_Tcp_vslp_Btn.setObjectName('Plot Tcp vs pump wavelength')

        self.ui_layoutTcpGroupBox.addWidget(self.ui_Tcp_vsPP_Btn, 1, 1)
        self.ui_layoutTcpGroupBox.addWidget(self.ui_Tcp_vslp_Btn, 2, 1)

        self.ui_TcpGroupBox.setLayout(self.ui_layoutTcpGroupBox)
        self.ui_layoutTcp.addWidget(self.ui_TcpGroupBox)

        return self.ui_layoutTcp

    def getProperties(self):
        self.CrystalMaterials = RefractiveIndex().materialList
        self.CrystalMaterial = self.CrystalMaterials[0]  # fallback

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

        self.lastRefractiveIndex = [self.CrystalNX, self.CrystalNY, self.CrystalNZ]
        self.PumpWlFrom = self.config.get("Pump wavelength from")
        self.PumpWlSingle = self.config.get("Pump wavelength single")
        self.PumpWlTo = self.config.get("Pump wavelength to")
        self.SIWlFrom = self.config.get("SI wavelength from")
        self.SIWlSingle = self.config.get("SI wavelength single")
        self.SIWlTo = self.config.get("SI wavelength to")
        self.PulsewidthFrom = self.config.get("Pump pulsewidth from")
        self.PulsewidthSingle = self.config.get("Pump pulsewidth single")
        self.PulsewidthTo = self.config.get("Pump pulsewidth to")
        self.PumpShape = self.config.get("Pump pulse shape")
        self.PumpShapeApplyDeconvolutionFactor = self.config.get("Pump pulsewidth apply deconvolution factor")
        self.CrystalTempFrom = self.config.get("Crystal Temperature from")
        self.CrystalTempSingle = self.config.get("Crystal Temperature single")
        self.CrystalTempTo = self.config.get("Crystal Temperature to")
        self.CrystalLengthFrom = self.config.get("Crystal Length from")
        self.CrystalLengthSingle = self.config.get("Crystal Length single")
        self.CrystalLengthTo = self.config.get("Crystal Length to")
        self.QPMOrder = self.config.get("QPM Order")
        self.JSIwlRange = self.config.get("JSI wavelength range")
        self.JSIresolution = self.config.get("JSI resolution")

        self.PurityWLresolution = self.config.get("Purity wavelength resolution")
        self.PurityTauresolution = self.config.get("Purity tau resolution")
        self.PurityWLrange = self.config.get("Purity wavelength range")

        self.SIfilterIdlerType = self.config.get("SI filter Idler Type")
        self.SIfilterSignalType = self.config.get("SI filter Signal Type")
        self.SIfilterIdlerCenterWL = self.config.get("SI filter Idler center wavelength")
        self.SIfilterSignalCenterWL = self.config.get("SI filter Signal center wavelength")
        self.SIfilterIdlerFWHM = self.config.get("SI filter Idler FWHM")
        self.SIfilterSignalFWHM = self.config.get("SI filter Signal FWHM")
        self.SIfilterTypes = Filters().FilterList

        self.HOMresolution = self.config.get('HOM interference plot resolution')
        self.HOMdelayrange = self.config.get('HOM interference plot range')
        
        self.fwhmres = self.config.get('FWHM plot resolution')
        self.fwhmprecision = self.config.get('FWHM decimal precision')
        
    def setProperties(self):
        self.ui_CrystalPolingPeriodsingleSB.setValue(self.CrystalPolingPeriodSingle * 10 ** 6)
        self.ui_CrystalPolingPeriodfromSB.setValue(self.CrystalPolingPeriodFrom * 10 ** 6)
        self.ui_CrystalPolingPeriodtoSB.setValue(self.CrystalPolingPeriodTo * 10 ** 6)
        self.ui_pumpwlfromSB.setValue(self.PumpWlFrom * 10 ** 9)
        self.ui_pumpwlsingleSB.setValue(self.PumpWlSingle * 10 ** 9)
        self.ui_pumpwltoSB.setValue(self.PumpWlTo * 10 ** 9)
        self.ui_siwlfromSB.setValue(self.SIWlFrom * 10 ** 9)
        self.ui_siwlsingleSB.setValue(self.SIWlSingle * 10 ** 9)
        self.ui_siwltoSB.setValue(self.SIWlTo * 10 ** 9)
        self.ui_pulsewidthfromSB.setValue(self.PulsewidthFrom * 10 ** 12)
        self.ui_pulsewidthsingleSB.setValue(self.PulsewidthSingle * 10 ** 12)
        self.ui_pulsewidthtoSB.setValue(self.PulsewidthTo * 10 ** 12)
        idx = 0  # Fallback
        idx = self.ui_pumpShapeCB.findText(self.PumpShape)
        self.ui_pumpShapeCB.setCurrentIndex(idx)
        self.ui_pumpShapeCorrectionFactorCheckBox.setChecked(self.PumpShapeApplyDeconvolutionFactor)
        self.ui_CrystalTfromSB.setValue(self.CrystalTempFrom)
        self.ui_CrystalTsingleSB.setValue(self.CrystalTempSingle)
        self.ui_CrystalTtoSB.setValue(self.CrystalTempTo)
        self.ui_CrystalLengthfromSB.setValue(self.CrystalLengthFrom * 10 ** 3)
        self.ui_CrystalLengthsingleSB.setValue(self.CrystalLengthSingle * 10 ** 3)
        self.ui_CrystalLengthtoSB.setValue(self.CrystalLengthTo * 10 ** 3)
        self.ui_PlotPMCSBQPMorder.setValue(self.QPMOrder)
        self.ui_PlotJSI_Wlrange_SB.setValue(self.JSIwlRange * 10 ** 9)
        self.ui_PlotJSI_WLresolution_SB.setValue(self.JSIresolution)
        self.ui_Purity_WLresolution_SB.setValue(self.PurityWLresolution)
        self.ui_Purity_WLrange_SB.setValue(self.PurityWLrange * 10 ** 9)
        self.ui_Purity_Tauresolution_SB.setValue(self.PurityTauresolution)
        idx = 0
        idx = self.ui_SIfilterIdlerType_CB.findText(self.SIfilterIdlerType)
        self.ui_SIfilterIdlerType_CB.setCurrentIndex(idx)
        idx = 0
        idx = self.ui_SIfilterSignalType_CB.findText(self.SIfilterSignalType)
        self.ui_SIfilterSignalType_CB.setCurrentIndex(idx)
        self.ui_SIfilterIdlerCenterWL_SB.setValue(self.SIfilterIdlerCenterWL * 10 ** 9)
        self.ui_SIfilterSignalCenterWL_SB.setValue(self.SIfilterSignalCenterWL * 10 ** 9)
        self.ui_SIfilterIdlerFWHM_SB.setValue(self.SIfilterIdlerFWHM * 10 ** 9)
        self.ui_SIfilterSignalFWHM_SB.setValue(self.SIfilterSignalFWHM * 10 ** 9)

        self.ui_HOM_Vis_Resolution_SB.setValue(self.HOMresolution)
        self.ui_HOM_Vis_DelayRange_SB.setValue(self.HOMdelayrange*10**12)
        
        self.ui_PlotFWHMresolution_SB.setValue(self.fwhmres)
        self.ui_PlotFWHMprecision_SB.setValue(self.fwhmprecision)

    def initConnections(self):
        self.ui_PlotRefractiveIndex_Btn_Plot_T.pressed.connect(self.plot_RefIdx_vs_T)
        self.ui_PlotRefractiveIndex_Btn_Plot_wl.pressed.connect(self.plot_RefIdx_vs_wl)
        self.ui_PlotPMCvsT_Btn.pressed.connect(self.plot_pmc_wl_vs_T)
        self.ui_PlotPMCvsPP_Btn.pressed.connect(self.plot_pmc_wl_vs_PP)
        self.ui_Purity_plotvsTau_Btn.pressed.connect(self.plot_purity_vs_tau)
        self.ui_Purity_plotvspwl_Btn.pressed.connect(self.plot_purity_vs_pwl)
        self.ui_PlotJSI_plotBtn.pressed.connect(self.plot_jsi)
        self.ui_PlotJSI_filterlossBtn.pressed.connect(self.estimate_filter_losses)
        self.ui_Purity_plotvsTauandL_Btn.pressed.connect(self.plot_purity_vs_Tau_and_L)
        self.ui_Purity_plotvsL_Btn.pressed.connect(self.plot_purity_vs_L)
        self.ui_GetEffPP_Btn.pressed.connect(self.GetEffectivePolingPeriod)
        self.ui_Tcp_vslp_Btn.pressed.connect(self.plot_Tcp_vs_lp)
        self.ui_Tcp_vsPP_Btn.pressed.connect(self.plot_Tcp_vs_PP)
        self.ui_HOM_PlotVis_Btn.pressed.connect(self.plot_HOM_vis)
        self.ui_PlotFWHMvstau_Btn.pressed.connect(self.plot_FWHM_vs_tau)
        self.ui_pumpShapeCorrectionFactorCheckBox.stateChanged.connect(self.getVarsFromGUI)
        
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
        self.ui_CrystalTsingleSB.valueChanged.connect(self.getVarsFromGUI)
        self.ui_CrystalTfromSB.valueChanged.connect(self.getVarsFromGUI)
        self.ui_CrystalTtoSB.valueChanged.connect(self.getVarsFromGUI)
        self.ui_PlotPMCSBQPMorder.valueChanged.connect(self.getVarsFromGUI)
        self.ui_PlotJSI_Wlrange_SB.valueChanged.connect(self.getVarsFromGUI)
        self.ui_CrystalLengthfromSB.valueChanged.connect(self.getVarsFromGUI)
        self.ui_CrystalLengthsingleSB.valueChanged.connect(self.getVarsFromGUI)
        self.ui_CrystalLengthtoSB.valueChanged.connect(self.getVarsFromGUI)
        self.ui_PlotJSI_WLresolution_SB.valueChanged.connect(self.getVarsFromGUI)
        self.ui_Purity_WLresolution_SB.valueChanged.connect(self.getVarsFromGUI)
        self.ui_Purity_WLrange_SB.valueChanged.connect(self.getVarsFromGUI)
        self.ui_Purity_Tauresolution_SB.valueChanged.connect(self.getVarsFromGUI)
        self.ui_SIfilterIdlerType_CB.currentIndexChanged.connect(self.getVarsFromGUI)
        self.ui_SIfilterSignalType_CB.currentIndexChanged.connect(self.getVarsFromGUI)
        self.ui_SIfilterIdlerCenterWL_SB.valueChanged.connect(self.getVarsFromGUI)
        self.ui_SIfilterSignalCenterWL_SB.valueChanged.connect(self.getVarsFromGUI)
        self.ui_SIfilterIdlerFWHM_SB.valueChanged.connect(self.getVarsFromGUI)
        self.ui_SIfilterSignalFWHM_SB.valueChanged.connect(self.getVarsFromGUI)
        self.ui_HOM_Vis_DelayRange_SB.valueChanged.connect(self.getVarsFromGUI)
        self.ui_HOM_Vis_Resolution_SB.valueChanged.connect(self.getVarsFromGUI)
        self.ui_PlotFWHMresolution_SB.valueChanged.connect(self.getVarsFromGUI)
        self.ui_PlotFWHMprecision_SB.valueChanged.connect(self.getVarsFromGUI)

        self.ui_PlotPMCvsT_Btn.mouseentersignal.connect(self.MouseHoverEnter)
        self.ui_PlotRefractiveIndex_Btn_Plot_T.mouseentersignal.connect(self.MouseHoverEnter)
        self.ui_PlotRefractiveIndex_Btn_Plot_wl.mouseentersignal.connect(self.MouseHoverEnter)
        self.ui_PlotPMCvsPP_Btn.mouseentersignal.connect(self.MouseHoverEnter)
        self.ui_GetEffPP_Btn.mouseentersignal.connect(self.MouseHoverEnter)
        self.ui_PlotJSI_plotBtn.mouseentersignal.connect(self.MouseHoverEnter)
        self.ui_Purity_plotvsTau_Btn.mouseentersignal.connect(self.MouseHoverEnter)
        self.ui_Purity_plotvspwl_Btn.mouseentersignal.connect(self.MouseHoverEnter)
        self.ui_Purity_plotvsL_Btn.mouseentersignal.connect(self.MouseHoverEnter)
        self.ui_Purity_plotvsTauandL_Btn.mouseentersignal.connect(self.MouseHoverEnter)
        self.ui_Tcp_vslp_Btn.mouseentersignal.connect(self.MouseHoverEnter)
        self.ui_Tcp_vsPP_Btn.mouseentersignal.connect(self.MouseHoverEnter)
        self.ui_HOM_PlotVis_Btn.mouseentersignal.connect(self.MouseHoverEnter)
        self.ui_PlotFWHMvstau_Btn.mouseentersignal.connect(self.MouseHoverEnter)

        self.ui_PlotPMCvsT_Btn.mouseleavesignal.connect(self.MouseHoverLeave)
        self.ui_PlotRefractiveIndex_Btn_Plot_T.mouseleavesignal.connect(self.MouseHoverLeave)
        self.ui_PlotRefractiveIndex_Btn_Plot_wl.mouseleavesignal.connect(self.MouseHoverLeave)
        self.ui_PlotPMCvsPP_Btn.mouseleavesignal.connect(self.MouseHoverLeave)
        self.ui_GetEffPP_Btn.mouseleavesignal.connect(self.MouseHoverLeave)
        self.ui_PlotJSI_plotBtn.mouseleavesignal.connect(self.MouseHoverLeave)
        self.ui_Purity_plotvsTau_Btn.mouseleavesignal.connect(self.MouseHoverLeave)
        self.ui_Purity_plotvspwl_Btn.mouseleavesignal.connect(self.MouseHoverLeave)
        self.ui_Purity_plotvsL_Btn.mouseleavesignal.connect(self.MouseHoverLeave)
        self.ui_Purity_plotvsTauandL_Btn.mouseleavesignal.connect(self.MouseHoverLeave)
        self.ui_Tcp_vslp_Btn.mouseleavesignal.connect(self.MouseHoverLeave)
        self.ui_Tcp_vsPP_Btn.mouseleavesignal.connect(self.MouseHoverLeave)
        self.ui_HOM_PlotVis_Btn.mouseleavesignal.connect(self.MouseHoverLeave)
        self.ui_PlotFWHMvstau_Btn.mouseleavesignal.connect(self.MouseHoverLeave)

    def plot_RefIdx_vs_T(self):
        Tmin = self.ui_CrystalTfromSB.value()
        Tmax = self.ui_CrystalTtoSB.value()
        wl = self.ui_pumpwlsingleSB.value() * 10 ** (-9)
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
        plotrange = numpy.linspace(Tmin, Tmax, 1000)

        pltwndidx = self.plotwindowcount
        self.open_new_plot_window()
        pltwnd = self.pltwindowlist[pltwndidx]

        for i in range(0, len(RefIdxList)):
            pltwnd.ax.plot(plotrange, RefIdxList[i](wl, plotrange),
                           label='{0}, {1}, {2}'.format(materialList[i], polList[i], paperList[i]))

        pltwnd.ax.set_xlabel('Temperature [°C]')
        pltwnd.ax.set_ylabel('Refractive index')
        pltwnd.ax.set_title('Refractive indices vs Temperature')
        pltwnd.ax.annotate('Wavelength: {0:.2f}nm'.format(wl * 10 ** 9), xy=(0.55, 0.01), xycoords='axes fraction')
        pltwnd.ax.legend()
        pltwnd.canvas.draw()

    def plot_RefIdx_vs_wl(self):
        wlmin = self.ui_pumpwlfromSB.value() * 10 ** (-9)
        wlmax = self.ui_pumpwltoSB.value() * 10 ** (-9)
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

        pltwndidx = self.plotwindowcount
        self.open_new_plot_window()
        pltwnd = self.pltwindowlist[pltwndidx]

        for i in range(0, len(RefIdxList)):
            pltwnd.ax.plot(plotrange * 10 ** 9, RefIdxList[i](plotrange, T),
                           label='{0}, {1}, {2}'.format(materialList[i], polList[i], paperList[i]))

        pltwnd.ax.set_xlabel('Wavelength [nm]')
        pltwnd.ax.set_ylabel('Refractive index')
        pltwnd.ax.set_title('Refractive indices vs Wavelength')
        pltwnd.ax.annotate('Temperature: {0:.1f}°C'.format(T), xy=(0.55, 0.01), xycoords='axes fraction')
        pltwnd.ax.legend()
        pltwnd.canvas.draw()

    def plot_pmc_wl_vs_T(self):
        # get variables from GUI
        lp = self.PumpWlSingle
        PP = self.CrystalPolingPeriodSingle
        Tmin = self.CrystalTempFrom
        Tmax = self.CrystalTempTo
        m = self.QPMOrder

        # get Ref indices
        nxfunc = RefractiveIndex().getSingleIDX(self.CrystalMaterial, "X", self.CrystalNX)
        nyfunc = RefractiveIndex().getSingleIDX(self.CrystalMaterial, "Y", self.CrystalNY)
        nzfunc = RefractiveIndex().getSingleIDX(self.CrystalMaterial, "Z", self.CrystalNZ)
        refidxfunc = [nxfunc, nyfunc, nzfunc]

        # prepare plotting
        plotrange = numpy.arange(Tmin, Tmax, (Tmax - Tmin) / 250)

        [siwl, idwl, Tcp] = PMC().getSI_wl_varT(lp, PP, plotrange, refidxfunc, m)

        # plot
        # init plot window
        pltwndidx = self.plotwindowcount
        self.open_new_plot_window()
        pltwnd = self.pltwindowlist[pltwndidx]
        pltwnd.ax.plot(plotrange, siwl * 10 ** 9, lw=2, label='signal')
        pltwnd.ax.plot(plotrange, idwl * 10 ** 9, lw=2, label='idler')
        pltwnd.ax.set_xlabel('Temperature [°C]')
        pltwnd.ax.set_ylabel('S/I Wavelength [nm]')
        pltwnd.ax.set_title('Quasi Phase Matching Temperature Tuning Curve')
        pltwnd.ax.annotate('Crossing point temperature: {0:.1f}°C'.format(Tcp), xy=(0.55, 0.01),
                           xycoords='axes fraction')
        pltwnd.ax.legend()
        pltwnd.canvas.draw()

    def plot_pmc_wl_vs_PP(self):
        # get variables from GUI
        lp = self.PumpWlSingle
        PPmin = self.CrystalPolingPeriodFrom
        PPmax = self.CrystalPolingPeriodTo
        m = self.QPMOrder
        T = self.CrystalTempSingle

        # get Ref indices
        nxfunc = RefractiveIndex().getSingleIDX(self.CrystalMaterial, "X", self.CrystalNX)
        nyfunc = RefractiveIndex().getSingleIDX(self.CrystalMaterial, "Y", self.CrystalNY)
        nzfunc = RefractiveIndex().getSingleIDX(self.CrystalMaterial, "Z", self.CrystalNZ)
        refidxfunc = [nxfunc, nyfunc, nzfunc]

        # prepare plotting
        plotrange = numpy.arange(PPmin, PPmax, (PPmax - PPmin) / 250)

        [siwl, idwl, PPcp] = PMC().getSI_wl_varPP(lp, plotrange, T, refidxfunc, m)

        # plot
        # init plot window
        pltwndidx = self.plotwindowcount
        self.open_new_plot_window()
        pltwnd = self.pltwindowlist[pltwndidx]
        pltwnd.ax.plot(plotrange * 10 ** 6, siwl * 10 ** 9, lw=2, label='signal')
        pltwnd.ax.plot(plotrange * 10 ** 6, idwl * 10 ** 9, lw=2, label='idler')
        pltwnd.ax.set_xlabel('Poling Period [µm]')
        pltwnd.ax.set_ylabel('S/I Wavelength [nm]')
        pltwnd.ax.set_title('Quasi Phase Matching Poling Period Tuning Curve')
        pltwnd.ax.annotate(
            'Degenerate wavelengths (for {0:.1f}°C) at a poling period of {1:.2f} µm'.format(self.CrystalTempSingle,
                                                                                             PPcp * 10 ** 6),
            xy=(0.01, 0.01), xycoords='axes fraction')
        pltwnd.ax.legend()
        pltwnd.canvas.draw()

    def plot_purity_vs_tau(self):
        wlpts = self.PurityWLresolution
        taupts = self.PurityTauresolution
        pwl = self.PumpWlSingle
        PP = self.CrystalPolingPeriodSingle
        L = self.CrystalLengthSingle
        T = self.CrystalTempSingle
        m = self.QPMOrder
        taumin = self.PulsewidthFrom
        taumax = self.PulsewidthTo
        wlrange = self.PurityWLrange
        ffi = Filters().getFilterFunction(self.SIfilterIdlerType,self.SIfilterIdlerCenterWL,self.SIfilterIdlerFWHM)
        ffs = Filters().getFilterFunction(self.SIfilterSignalType,self.SIfilterSignalCenterWL,self.SIfilterSignalFWHM)
        spectralfilters = [ffs, ffi]
        FilterString = 'none'
        pumpshape = self.PumpShape
        calcGaussian = False
        calcSech = False
        if pumpshape == 'Gaussian':
            calcGaussian = True
        elif pumpshape == 'Sech^2':
            calcSech = True
        else:
            print('Error: pump shape unknown to JSA/JSI plot routine')

        # get Ref indices
        nxfunc = RefractiveIndex().getSingleIDX(self.CrystalMaterial, "X", self.CrystalNX)
        nyfunc = RefractiveIndex().getSingleIDX(self.CrystalMaterial, "Y", self.CrystalNY)
        nzfunc = RefractiveIndex().getSingleIDX(self.CrystalMaterial, "Z", self.CrystalNZ)
        refidxfunc = [nxfunc, nyfunc, nzfunc]

        taurange = numpy.arange(taumin, taumax, (taumax-taumin) / taupts)

        Tvec = numpy.arange(T, T + 1, 2)

        [ls, li, unused] = PMC().getSI_wl_varT(pwl, PP, Tvec, refidxfunc, m)

        signalrange = numpy.linspace(ls - wlrange / 2, ls + wlrange / 2, wlpts)
        idlerrange = numpy.linspace(li - wlrange / 2, li + wlrange / 2, wlpts)

        [purity, max, maxtau] = JSI().getpurity_vsTau(pwl, signalrange, idlerrange, taurange, T,
                                                      PP, L, refidxfunc, m, spectralfilters, pumpshape)

        AnnotateString = ''
        AnnotateString = AnnotateString + \
                         r'Maximum purity: {0:.3} at $\tau={1:.3}$ps'.format(max, maxtau * 10 ** 12)+'\n' \
                                                                                                     'Pump wavelength: {0:.2f}nm'.format(pwl*10**9)+'\n' \
                                                                                                                                                    'Poling period: {0:.4f}µm'.format(PP*10**6)+'\n' \
                                                                                                                                                                                                'Temperature: {0:.1f}°C'.format(T)+'\n' \
                                                                                                                                                                                                                                   'Crystal length: {0:.1f}mm'.format(L*10**3)+'\n'
        FilterString = ''
        print(self.SIfilterSignalType)
        if self.SIfilterSignalType.casefold() != 'none':
            FilterString = 'Spectral filters:\n' r'Signal: {0}, central wl {1:.2f}nm, FWHM {2:.2f}nm'.format(
                self.SIfilterSignalType, self.SIfilterSignalCenterWL * 10 ** 9,
                                         self.SIfilterSignalFWHM * 10 ** 9) + '\n'
        if self.SIfilterIdlerType.casefold() != 'none':
            FilterString = FilterString + r'Idler: {0}, central wl {1:.2f}nm, FWHM {2:.2f}nm'.format(
                self.SIfilterIdlerType, self.SIfilterIdlerCenterWL * 10 ** 9, self.SIfilterIdlerFWHM * 10 ** 9) + '\n'
        AnnotateString = AnnotateString + FilterString

        # plot
        # init plot window
        pltwndidx = self.plotwindowcount
        self.open_new_plot_window()
        pltwnd = self.pltwindowlist[pltwndidx]

        pltwnd.ax.plot(taurange * 10 ** 12, purity, lw=2, label='purity')
        pltwnd.ax.set_xlabel('Pulsewidth [ps]')
        pltwnd.ax.set_ylabel('Purity')
        pltwnd.ax.set_title('Purity')
        pltwnd.ax.annotate(AnnotateString, xy=(0.2, 0.1),
                           xycoords='axes fraction')
        pltwnd.ax.legend()
        pltwnd.canvas.draw()

    def plot_purity_vs_pwl(self):
        pass

    def plot_purity_vs_L(self):
        wlpts = self.PurityWLresolution
        pts = self.PurityTauresolution
        pwl = self.PumpWlSingle
        PP = self.CrystalPolingPeriodSingle
        Lmin = self.CrystalLengthFrom
        Lmax = self.CrystalLengthTo
        T = self.CrystalTempSingle
        m = self.QPMOrder
        tau = self.PulsewidthSingle
        wlrange = self.PurityWLrange
        ffi = Filters().getFilterFunction(self.SIfilterIdlerType, self.SIfilterIdlerCenterWL, self.SIfilterIdlerFWHM)
        ffs = Filters().getFilterFunction(self.SIfilterSignalType, self.SIfilterSignalCenterWL, self.SIfilterSignalFWHM)
        spectralfilters = [ffs, ffi]
        FilterString = 'none'
        pumpshape = self.PumpShape
        calcGaussian = False
        calcSech = False
        if pumpshape == 'Gaussian':
            calcGaussian = True
        elif pumpshape == 'Sech^2':
            calcSech = True
        else:
            print('Error: pump shape unknown to JSA/JSI plot routine')

        # get Ref indices
        nxfunc = RefractiveIndex().getSingleIDX(self.CrystalMaterial, "X", self.CrystalNX)
        nyfunc = RefractiveIndex().getSingleIDX(self.CrystalMaterial, "Y", self.CrystalNY)
        nzfunc = RefractiveIndex().getSingleIDX(self.CrystalMaterial, "Z", self.CrystalNZ)
        refidxfunc = [nxfunc, nyfunc, nzfunc]

        Lrange = numpy.arange(Lmin, Lmax, (Lmax- Lmin) / pts)

        Tvec = numpy.arange(T, T + 1, 2)

        [ls, li, unused] = PMC().getSI_wl_varT(pwl, PP, Tvec, refidxfunc, m)

        signalrange = numpy.linspace(ls - wlrange / 2, ls + wlrange / 2, wlpts)
        idlerrange = numpy.linspace(li - wlrange / 2, li + wlrange / 2, wlpts)

        [purity, max, maxL] = JSI().getpurity_vsL(pwl, signalrange, idlerrange, tau, T,
                                                  PP, Lrange, refidxfunc, m, spectralfilters, pumpshape)

        AnnotateString = ''
        AnnotateString = AnnotateString + \
                         r'Maximum purity: {0:.3f} at L={1:.2f}mm'.format(max, maxL*10**3) + '\n' \
                                                                                             'Pump wavelength: {0:.2f}nm'.format(
            pwl * 10 ** 9) + '\n' \
                             'Poling period: {0:.4f}µm'.format(PP * 10 ** 6) + '\n' \
                                                                               'Temperature: {0:.1f}°C'.format(T) + '\n' \
                                                                                                                    'Pulsewidth: {0:.2f}ps'.format(
            tau * 10 ** 12) + '\n'
        FilterString = ''
        if self.SIfilterSignalType.casefold() != 'none':
            FilterString = 'Spectral filters:\n' r'Signal: {0}, central wl {1:.2f}nm, FWHM {2:.2f}nm'.format(
                self.SIfilterSignalType, self.SIfilterSignalCenterWL * 10 ** 9,
                                         self.SIfilterSignalFWHM * 10 ** 9) + '\n'
        if self.SIfilterIdlerType.casefold() != 'none':
            FilterString = FilterString + r'Idler: {0}, central wl {1:.2f}nm, FWHM {2:.2f}nm'.format(
                self.SIfilterIdlerType, self.SIfilterIdlerCenterWL * 10 ** 9, self.SIfilterIdlerFWHM * 10 ** 9) + '\n'
        AnnotateString = AnnotateString + FilterString

        # plot
        # init plot window
        pltwndidx = self.plotwindowcount
        self.open_new_plot_window()
        pltwnd = self.pltwindowlist[pltwndidx]
        pltwnd.ax.plot(Lrange * 10 ** 3, purity, lw=2, label='purity')
        pltwnd.ax.set_xlabel('Crystal length [mm]')
        pltwnd.ax.set_ylabel('Purity')
        pltwnd.ax.set_title('Purity')
        pltwnd.ax.annotate(AnnotateString, xy=(0.2, 0.1),
                           xycoords='axes fraction')
        pltwnd.ax.legend()
        pltwnd.canvas.draw()

    def plot_purity_vs_Tau_and_L(self):
        wlpts = self.PurityWLresolution
        pts = self.PurityTauresolution
        pwl = self.PumpWlSingle
        PP = self.CrystalPolingPeriodSingle
        Lmin = self.CrystalLengthFrom
        Lmax = self.CrystalLengthTo
        T = self.CrystalTempSingle
        m = self.QPMOrder
        taumin = self.PulsewidthFrom
        taumax = self.PulsewidthTo
        wlrange = self.PurityWLrange
        ffi = Filters().getFilterFunction(self.SIfilterIdlerType, self.SIfilterIdlerCenterWL, self.SIfilterIdlerFWHM)
        ffs = Filters().getFilterFunction(self.SIfilterSignalType, self.SIfilterSignalCenterWL, self.SIfilterSignalFWHM)
        spectralfilters = [ffs, ffi]
        FilterString = 'none'
        pumpshape = self.PumpShape
        calcGaussian = False
        calcSech = False
        if pumpshape == 'Gaussian':
            calcGaussian = True
        elif pumpshape == 'Sech^2':
            calcSech = True
        else:
            print('Error: pump shape unknown to JSA/JSI plot routine')

        # get Ref indices
        nxfunc = RefractiveIndex().getSingleIDX(self.CrystalMaterial, "X", self.CrystalNX)
        nyfunc = RefractiveIndex().getSingleIDX(self.CrystalMaterial, "Y", self.CrystalNY)
        nzfunc = RefractiveIndex().getSingleIDX(self.CrystalMaterial, "Z", self.CrystalNZ)
        refidxfunc = [nxfunc, nyfunc, nzfunc]

        Lrange = numpy.arange(Lmin, Lmax, (Lmax - Lmin) / pts)
        Taurange = numpy.arange(taumin, taumax, (taumax - taumin) / pts)

        Tvec = numpy.arange(T, T + 1, 2)

        [ls, li, unused] = PMC().getSI_wl_varT(pwl, PP, Tvec, refidxfunc, m)

        signalrange = numpy.linspace(ls - wlrange / 2, ls + wlrange / 2, wlpts)
        idlerrange = numpy.linspace(li - wlrange / 2, li + wlrange / 2, wlpts)

        purity = JSI().getpurity_vsLandTau(pwl, signalrange, idlerrange, Taurange, T,
                                           PP, Lrange, refidxfunc, m, spectralfilters, pumpshape)

        AnnotateString = ''
        AnnotateString = AnnotateString + 'Pump wavelength: {0:.2f}nm'.format(pwl * 10 ** 9) + '\n'+ \
                         'Poling period: {0:.4f}µm'.format(PP * 10 ** 6) + '\n'
        FilterString = ''
        if self.SIfilterSignalType.casefold() != 'none':
            FilterString = 'Spectral filters:\n' r'Signal: {0}, central wl {1:.2f}nm, FWHM {2:.2f}nm'.format(
                self.SIfilterSignalType, self.SIfilterSignalCenterWL * 10 ** 9, self.SIfilterSignalFWHM * 10 ** 9) + '\n'
        if self.SIfilterIdlerType.casefold() != 'none':
            FilterString = FilterString + r'Idler: {0}, central wl {1:.2f}nm, FWHM {2:.2f}nm'.format(
                self.SIfilterIdlerType, self.SIfilterIdlerCenterWL * 10 ** 9, self.SIfilterIdlerFWHM * 10 ** 9) + '\n'
        AnnotateString = AnnotateString + FilterString

        # plot
        colormap = matplotlib.cm.jet
        xmin = numpy.min(Taurange) * 10 ** 12  # *taucfsech
        xmax = numpy.max(Taurange) * 10 ** 12  # *taucfsech
        ymin = numpy.min(Lrange) * 10 ** 3
        ymax = numpy.max(Lrange) * 10 ** 3

        # init plot window
        pltwndidx = self.plotwindowcount
        self.open_new_plot_window()
        pltwnd = self.pltwindowlist[pltwndidx]
        pltwnd.ax.grid('off')
        plot=pltwnd.ax.imshow(purity, cmap=colormap, vmin=abs(purity).min(), vmax=abs(purity).max(),aspect='auto',
                              origin='lower',interpolation='none', extent=[xmin,xmax,ymin,ymax])
        #(Lrange * 10 ** 3, purity, lw=2, label='purity')
        pltwnd.ax.set_xlabel('Pulse width [ps]')
        pltwnd.ax.set_ylabel('Crystal length [mm]')
        pltwnd.ax.set_title('Purity')
        pltwnd.ax.annotate(AnnotateString, xy=(0, 1),
                           xycoords='axes fraction')

        pltwnd.fig.subplots_adjust(bottom=0.2)
        pltwnd.cbar_ax = pltwnd.fig.add_axes([0.05, 0.1, 0.9, 0.025])
        pltwnd.cbar = pltwnd.fig.colorbar(plot, cax=pltwnd.cbar_ax, orientation='horizontal')
        pltwnd.cbar.set_label('a.u.', fontsize='medium', labelpad=-1)

        pltwnd.ax.legend()
        pltwnd.canvas.draw()

    def plot_jsi(self):
        numpts=self.JSIresolution
        pwl=self.PumpWlSingle
        PP=self.CrystalPolingPeriodSingle
        L=self.CrystalLengthSingle
        T=self.CrystalTempSingle
        m=self.QPMOrder
        tau=self.PulsewidthSingle
        wlrange=self.JSIwlRange
        ffi = Filters().getFilterFunction(self.SIfilterIdlerType, self.SIfilterIdlerCenterWL, self.SIfilterIdlerFWHM)
        ffs = Filters().getFilterFunction(self.SIfilterSignalType, self.SIfilterSignalCenterWL, self.SIfilterSignalFWHM)
        spectralfilters = [ffs, ffi]

        plotJSI=self.ui_PlotJSI_plotJSIRadioButton.isChecked()
        if plotJSI==True:
            calcJSA=False
            calcJSI=True
        else:
            calcJSA=True
            calcJSI=False

        pumpshape = self.PumpShape
        calcGaussian=False
        calcSech=False
        if pumpshape == 'Gaussian':
            calcGaussian=True
        elif pumpshape == 'Sech^2':
            calcSech=True
        else:
            print('Error: pump shape unknown to JSA/JSI plot routine')

        nxfunc = RefractiveIndex().getSingleIDX(self.CrystalMaterial, "X", self.CrystalNX)
        nyfunc = RefractiveIndex().getSingleIDX(self.CrystalMaterial, "Y", self.CrystalNY)
        nzfunc = RefractiveIndex().getSingleIDX(self.CrystalMaterial, "Z", self.CrystalNZ)
        refidxfunc = [nxfunc, nyfunc, nzfunc]

        Tvec=numpy.arange(T,T+1,2)

        [ls, li, unused] = PMC().getSI_wl_varT(pwl, PP, Tvec, refidxfunc, m)

        signalrange = numpy.linspace(ls - wlrange/2, ls + wlrange/2, numpts)
        idlerrange = numpy.linspace(li - wlrange/2, li + wlrange/2, numpts)

        [PE, PM, JS] = JSI().getplots(pwl, signalrange, idlerrange, tau, T, PP, L, refidxfunc,
                                      m, spectralfilters, plotJSI, pumpshape)

        #
        # plotting
        #

        print("plotting..")

        # init plot window
        pltwndidx = self.plotwindowcount
        self.open_new_plot_window()
        pltwnd = self.pltwindowlist[pltwndidx]

        # need customization for 3 plots in 1 window
        pltwnd.layout.removeWidget(pltwnd.canvas)
        pltwnd.layout.removeWidget(pltwnd.toolbar)
        pltwnd.fig = figure(facecolor="white")
        pltwnd.peplt = pltwnd.fig.add_subplot(131)
        pltwnd.pmplt = pltwnd.fig.add_subplot(132)
        pltwnd.jsplt = pltwnd.fig.add_subplot(133)
        pltwnd.peplt.set_aspect('equal')
        pltwnd.pmplt.set_aspect('equal')
        pltwnd.jsplt.set_aspect('equal')
        pltwnd.canvas = FigureCanvas(pltwnd.fig)
        pltwnd.canvas.setParent(pltwnd)
        pltwnd.toolbar = NavigationToolbar(pltwnd.canvas, pltwnd)
        # self.addToolBar(self.toolbar)
        pltwnd.layout.addWidget(pltwnd.canvas)
        pltwnd.layout.addWidget(pltwnd.toolbar)

        colormap = matplotlib.cm.jet
        # axes range
        xmin = numpy.min(signalrange) * 10 ** 9
        xmax = numpy.max(signalrange) * 10 ** 9
        ymin = numpy.min(idlerrange) * 10 ** 9
        ymax = numpy.max(idlerrange) * 10 ** 9
        # prepare subplots
        ppe = pltwnd.peplt.imshow(PE, cmap=colormap, vmin=PE.min(), vmax=PE.max(), aspect='auto',
                                  origin='lower', interpolation='none', extent=[xmin, xmax, ymin, ymax])
        ppm = pltwnd.pmplt.imshow(PM, cmap=colormap, vmin=PM.min(), vmax=PM.max(), aspect='auto',
                                  origin='lower', interpolation='none', extent=[xmin, xmax, ymin, ymax])
        pjs = pltwnd.jsplt.imshow(JS, cmap=colormap, vmin=JS.min(), vmax=JS.max(), aspect='auto',
                                  origin='lower', interpolation='none', extent=[xmin, xmax, ymin, ymax])
        pltwnd.peplt.set_xlabel(r'$\lambda_s$ [nm]')
        pltwnd.peplt.set_ylabel(r'$\lambda_i$ [nm]')
        # size of axes label
        for plot in [pltwnd.peplt, pltwnd.pmplt, pltwnd.jsplt]:
            plot.tick_params(axis='both', which='major', labelsize='medium')
            plot.tick_params(axis='both', which='minor', labelsize='medium')
        # label plot
        if calcJSA:
            pltwnd.peplt.set_title('PEA', fontsize=20)
            pltwnd.pmplt.set_title('PMA', fontsize=20)
            pltwnd.jsplt.set_title('JSA', fontsize=20)
            pltwnd.fig.suptitle('PEA, PMA and JSA', fontsize=24)
        elif calcJSI:
            pltwnd.peplt.set_title('PEI', fontsize=20)
            pltwnd.pmplt.set_title('PMI', fontsize=20)
            pltwnd.jsplt.set_title('JSI', fontsize=20)
            pltwnd.fig.suptitle('PEI, PMI and JSI', fontsize=24)
        # create legend
        pltwnd.fig.subplots_adjust(bottom=0.2)
        pltwnd.cbar_ax = pltwnd.fig.add_axes([0.05, 0.1, 0.9, 0.025])
        pltwnd.cbar = pltwnd.fig.colorbar(ppe, cax=pltwnd.cbar_ax, orientation='horizontal')
        pltwnd.cbar.set_label('a.u.', fontsize='medium', labelpad=-1)
        pad = 20
        parameterstring = 'Pump wavelength: {0:.2f}nm\n'.format(pwl * 10 ** 9) + 'Crystal Length: {0:.2f}mm\n'.format(
            L * 10 ** 3) + 'Poling period: {0:.2f}µm\n'.format(PP * 10 ** 6) + 'Temperature: {0:.2f}°C\n'.format(
            T) + 'Pulse duration: {0:.2f}ps'.format(tau * 10 ** 12)
        if calcGaussian:
            parameterstring = parameterstring + '\nGaussian beam shape'
        elif calcSech:
            parameterstring = parameterstring + '\nsech^2 beam shape'
        # state additional paramters on plot
        pltwnd.peplt.annotate(parameterstring, xy=(0.005, 0.83), xycoords='figure fraction', fontsize=9, color='r')

        pltwnd.peplt.set_aspect('equal')
        pltwnd.pmplt.set_aspect('equal')
        pltwnd.jsplt.set_aspect('equal')
        pltwnd.fig.set_size_inches(8,3)
        pltwnd.canvas.draw()
        pltwnd.resize(1200,600)
    
    def estimate_filter_losses(self):
        numpts=self.JSIresolution
        pwl=self.PumpWlSingle
        PP=self.CrystalPolingPeriodSingle
        L=self.CrystalLengthSingle
        T=self.CrystalTempSingle
        m=self.QPMOrder
        tau=self.PulsewidthSingle
        wlrange=self.JSIwlRange
        ffi = Filters().getFilterFunction(self.SIfilterIdlerType, self.SIfilterIdlerCenterWL, self.SIfilterIdlerFWHM)
        ffs = Filters().getFilterFunction(self.SIfilterSignalType, self.SIfilterSignalCenterWL, self.SIfilterSignalFWHM)
        spectralfilters = [ffs, ffi]
        
        #as reference: no filters:
        ffiref = Filters().getFilterFunction('None', 1, 1) #arguments: filtertype, central wl, bandwidth. only first argument does matter if type is 'None'
        ffsref = Filters().getFilterFunction('None', 1, 1)
        spectralfiltersref = [ffsref, ffiref]

        plotJSI=self.ui_PlotJSI_plotJSIRadioButton.isChecked()
        if plotJSI==True:
            calcJSA=False
            calcJSI=True
        else:
            calcJSA=True
            calcJSI=False

        pumpshape = self.PumpShape
        calcGaussian=False
        calcSech=False
        if pumpshape == 'Gaussian':
            calcGaussian=True
        elif pumpshape == 'Sech^2':
            calcSech=True
        else:
            print('Error: pump shape unknown to JSA/JSI plot routine')

        nxfunc = RefractiveIndex().getSingleIDX(self.CrystalMaterial, "X", self.CrystalNX)
        nyfunc = RefractiveIndex().getSingleIDX(self.CrystalMaterial, "Y", self.CrystalNY)
        nzfunc = RefractiveIndex().getSingleIDX(self.CrystalMaterial, "Z", self.CrystalNZ)
        refidxfunc = [nxfunc, nyfunc, nzfunc]

        Tvec=numpy.arange(T,T+1,2)

        [ls, li, unused] = PMC().getSI_wl_varT(pwl, PP, Tvec, refidxfunc, m)

        signalrange = numpy.linspace(ls - wlrange/2, ls + wlrange/2, numpts)
        idlerrange = numpy.linspace(li - wlrange/2, li + wlrange/2, numpts)

        [PE, PM, JS] = JSI().getplots(pwl, signalrange, idlerrange, tau, T, PP, L, refidxfunc,
                                      m, spectralfilters, plotJSI, pumpshape)
        [PEref,PMref,JSref] = JSI().getplots(pwl, signalrange, idlerrange, tau, T, PP, L, refidxfunc,
                                      m, spectralfiltersref, plotJSI, pumpshape)
        
        #
        # calculating losses
        #
        magnitude_withfilter = 0
        magnitude_wofilter = 0
        
        for i in range(0,len(JS)):
            for j in range(0,len(JS[i])):
                magnitude_withfilter = magnitude_withfilter + JS[i,j]
                magnitude_wofilter = magnitude_wofilter + JSref[i,j]
        
        filterlosses = (1- magnitude_withfilter/magnitude_wofilter)
        
        print('magnitude wo filter: ', magnitude_wofilter)
        print('magnitude w filter: ', magnitude_withfilter)
        print('filterlosses: ', filterlosses)
        
        #
        # plotting
        #

        print("plotting..")

        # init plot window
        pltwndidx = self.plotwindowcount
        self.open_new_plot_window()
        pltwnd = self.pltwindowlist[pltwndidx]

        # need customization for 3 plots in 1 window
        pltwnd.layout.removeWidget(pltwnd.canvas)
        pltwnd.layout.removeWidget(pltwnd.toolbar)
        pltwnd.fig = figure(facecolor="white")
        pltwnd.peplt = pltwnd.fig.add_subplot(131)
        pltwnd.pmplt = pltwnd.fig.add_subplot(132)
        #pltwnd.jsplt = pltwnd.fig.add_subplot(133)
        pltwnd.peplt.set_aspect('equal')
        pltwnd.pmplt.set_aspect('equal')
        #pltwnd.jsplt.set_aspect('equal')
        pltwnd.canvas = FigureCanvas(pltwnd.fig)
        pltwnd.canvas.setParent(pltwnd)
        pltwnd.toolbar = NavigationToolbar(pltwnd.canvas, pltwnd)
        # self.addToolBar(self.toolbar)
        pltwnd.layout.addWidget(pltwnd.canvas)
        pltwnd.layout.addWidget(pltwnd.toolbar)

        colormap = matplotlib.cm.jet
        # axes range
        xmin = numpy.min(signalrange) * 10 ** 9
        xmax = numpy.max(signalrange) * 10 ** 9
        ymin = numpy.min(idlerrange) * 10 ** 9
        ymax = numpy.max(idlerrange) * 10 ** 9
        # prepare subplots
        ppe = pltwnd.peplt.imshow(JS, cmap=colormap, vmin=PE.min(), vmax=PE.max(), aspect='auto',
                                  origin='lower', interpolation='none', extent=[xmin, xmax, ymin, ymax])
        ppm = pltwnd.pmplt.imshow(JSref, cmap=colormap, vmin=PM.min(), vmax=PM.max(), aspect='auto',
                                  origin='lower', interpolation='none', extent=[xmin, xmax, ymin, ymax])
        # pjs = pltwnd.jsplt.imshow(JSref, cmap=colormap, vmin=JS.min(), vmax=JS.max(), aspect='auto',
        #                          origin='lower', interpolation='none', extent=[xmin, xmax, ymin, ymax])
        pltwnd.peplt.set_xlabel(r'$\lambda_s$ [nm]')
        pltwnd.peplt.set_ylabel(r'$\lambda_i$ [nm]')
        # size of axes label
        for plot in [pltwnd.peplt, pltwnd.pmplt]: #, pltwnd.jsplt]:
            plot.tick_params(axis='both', which='major', labelsize='medium')
            plot.tick_params(axis='both', which='minor', labelsize='medium')
        # label plot
        if calcJSA:
            pltwnd.peplt.set_title('PEA', fontsize=20)
            pltwnd.pmplt.set_title('PMA', fontsize=20)
            #pltwnd.jsplt.set_title('JSA', fontsize=20)
            pltwnd.fig.suptitle('TODO: Plot JSA with and without filter.', fontsize=24)
        elif calcJSI:
            pltwnd.peplt.set_title('PEI', fontsize=20)
            pltwnd.pmplt.set_title('PMI', fontsize=20)
            #pltwnd.jsplt.set_title('JSI', fontsize=20)
            pltwnd.fig.suptitle('TODO: Plot JSI with and without filter.', fontsize=24)
        # create legend
        pltwnd.fig.subplots_adjust(bottom=0.2)
        pltwnd.cbar_ax = pltwnd.fig.add_axes([0.05, 0.1, 0.9, 0.025])
        pltwnd.cbar = pltwnd.fig.colorbar(ppe, cax=pltwnd.cbar_ax, orientation='horizontal')
        pltwnd.cbar.set_label('a.u.', fontsize='medium', labelpad=-1)
        pad = 20
        parameterstring = 'Pump wavelength: {0:.2f}nm\n'.format(pwl * 10 ** 9) + 'Crystal Length: {0:.2f}mm\n'.format(
            L * 10 ** 3) + 'Poling period: {0:.2f}µm\n'.format(PP * 10 ** 6) + 'Temperature: {0:.2f}°C\n'.format(
            T) + 'Pulse duration: {0:.2f}ps\n'.format(tau * 10 ** 12) + 'Filterlosses: {0:.2f}%'.format(filterlosses*100)
        if calcGaussian:
            parameterstring = parameterstring + '\nGaussian beam shape'
        elif calcSech:
            parameterstring = parameterstring + '\nsech^2 beam shape'
        # state additional paramters on plot
        pltwnd.peplt.annotate(parameterstring, xy=(0.005, 0.83), xycoords='figure fraction', fontsize=9, color='r')

        pltwnd.peplt.set_aspect('equal')
        pltwnd.pmplt.set_aspect('equal')
        #pltwnd.jsplt.set_aspect('equal')
        pltwnd.fig.set_size_inches(8,3)
        pltwnd.canvas.draw()
        pltwnd.resize(1200,600)

    def plot_Tcp_vs_PP(self):
        numpts = 100
        PPrange = numpy.linspace(self.CrystalPolingPeriodFrom, self.CrystalPolingPeriodTo, numpts)
        temp = self.CrystalTempSingle
        pwl = self.PumpWlSingle
        qpmorder = self.QPMOrder
        nxfunc = RefractiveIndex().getSingleIDX(self.CrystalMaterial, "X", self.CrystalNX)
        nyfunc = RefractiveIndex().getSingleIDX(self.CrystalMaterial, "Y", self.CrystalNY)
        nzfunc = RefractiveIndex().getSingleIDX(self.CrystalMaterial, "Z", self.CrystalNZ)
        refidxfunc = [nxfunc, nyfunc, nzfunc]
        Tcp = JSI().getTcpVsPP(PPrange, temp, pwl, refidxfunc, qpmorder)

        # plot
        # init plot window
        pltwndidx = self.plotwindowcount
        self.open_new_plot_window()
        pltwnd = self.pltwindowlist[pltwndidx]
        pltwnd.ax.plot(PPrange * 10 ** 6, Tcp, lw=2, label='Crossing point temperature')
        pltwnd.ax.set_xlabel('Poling period [µm]')
        pltwnd.ax.set_ylabel('Crossing point temperature [°C]')
        pltwnd.ax.set_title('Crossing point temperature vs crystal poling period')
        # pltwnd.ax.annotate(
        #    'Degenerate wavelengths (for {0:.1f}°C) at a poling period of {1:.2f} µm'.format(self.CrystalTempSingle,
        #                                                                                     PPcp * 10 ** 6),
        #    xy=(0.01, 0.01), xycoords='axes fraction')
        pltwnd.ax.legend()
        pltwnd.canvas.draw()

    def plot_Tcp_vs_lp(self):
        numpts=100
        pwlrange=numpy.linspace(self.PumpWlFrom,self.PumpWlTo,numpts)
        temp = self.CrystalTempSingle
        polingp = self.CrystalPolingPeriodSingle
        qpmorder=self.QPMOrder
        nxfunc = RefractiveIndex().getSingleIDX(self.CrystalMaterial, "X", self.CrystalNX)
        nyfunc = RefractiveIndex().getSingleIDX(self.CrystalMaterial, "Y", self.CrystalNY)
        nzfunc = RefractiveIndex().getSingleIDX(self.CrystalMaterial, "Z", self.CrystalNZ)
        refidxfunc = [nxfunc, nyfunc, nzfunc]
        Tcp = JSI().getTcpVslp(pwlrange, temp, polingp, refidxfunc, qpmorder)

        # plot
        # init plot window
        pltwndidx = self.plotwindowcount
        self.open_new_plot_window()
        pltwnd = self.pltwindowlist[pltwndidx]
        pltwnd.ax.plot(pwlrange * 10 ** 9, Tcp, lw=2, label='Crossing point temperature')
        pltwnd.ax.set_xlabel('Pump wavelength [nm]')
        pltwnd.ax.set_ylabel('Crossing point temperature [°C]')
        pltwnd.ax.set_title('Crossing point temperature vs pump wavelength')
        # pltwnd.ax.annotate(
        #    'Degenerate wavelengths (for {0:.1f}°C) at a poling period of {1:.2f} µm'.format(self.CrystalTempSingle,
        #                                                                                     PPcp * 10 ** 6),
        #    xy=(0.01, 0.01), xycoords='axes fraction')
        pltwnd.ax.legend()
        pltwnd.canvas.draw()

    def plot_HOM_vis(self):
        pwl = self.PumpWlSingle
        T = self.CrystalTempSingle
        PP = self.CrystalPolingPeriodSingle
        m = self.QPMOrder
        tau = self.PulsewidthSingle
        cl = self.CrystalLengthSingle
        pumpshape = self.PumpShape
        delayrange = numpy.linspace(-self.HOMdelayrange/2, self.HOMdelayrange/2,self.HOMresolution)
        JSIresolution = self.JSIresolution
        JSIwlrange = self.JSIwlRange

        ffi = Filters().getFilterFunction(self.SIfilterIdlerType, self.SIfilterIdlerCenterWL, self.SIfilterIdlerFWHM)
        ffs = Filters().getFilterFunction(self.SIfilterSignalType, self.SIfilterSignalCenterWL, self.SIfilterSignalFWHM)
        spectralfilters = [ffs, ffi]

        nxfunc = RefractiveIndex().getSingleIDX(self.CrystalMaterial, "X", self.CrystalNX)
        nyfunc = RefractiveIndex().getSingleIDX(self.CrystalMaterial, "Y", self.CrystalNY)
        nzfunc = RefractiveIndex().getSingleIDX(self.CrystalMaterial, "Z", self.CrystalNZ)
        refidxfunc = [nxfunc, nyfunc, nzfunc]

        Tvec = numpy.arange(T, T + 1, 2)
        [ls, li, unused] = PMC().getSI_wl_varT(pwl, PP, Tvec, refidxfunc, m)
        signalrange = numpy.linspace(ls - JSIwlrange / 2, ls + JSIwlrange / 2, JSIresolution)
        idlerrange = numpy.linspace(li - JSIwlrange / 2, li + JSIwlrange / 2, JSIresolution)

        [CoincProb,vis,fwhm] = JSI().getHOMinterference(pwl, T, PP, m, tau, cl, signalrange, idlerrange,
                                             JSIresolution, pumpshape, delayrange, refidxfunc, spectralfilters)

        # plot
        # init plot window
        pltwndidx = self.plotwindowcount
        self.open_new_plot_window()
        pltwnd = self.pltwindowlist[pltwndidx]
        pltwnd.ax.plot(delayrange * 10 ** 12, CoincProb, lw=2)
        pltwnd.ax.set_xlabel('Delay [ps]')
        pltwnd.ax.set_ylabel('Coincidence probability')
        pltwnd.ax.set_title('Hong-Ou-Mandel interference')
        pltwnd.ax.annotate('Visibility: {0:.3f} \n FWHM: {1:.3f}ps'.format(vis,fwhm*10**12), xy=(0.01, 0.01),
                           xycoords='axes fraction')
        pltwnd.ax.legend()
        pltwnd.canvas.draw()

    def plot_FWHM_vs_tau(self):
        pwl = self.PumpWlSingle
        T = self.CrystalTempSingle
        PP = self.CrystalPolingPeriodSingle
        m = self.QPMOrder
        cl = self.CrystalLengthSingle
        pumpshape = self.PumpShape
        JSIresolution = self.JSIresolution
        JSIwlrange = self.JSIwlRange
        taumin = self.PulsewidthFrom
        taumax = self.PulsewidthTo
        fwhmres = self.fwhmres
        decprec = self.fwhmprecision
        usetaucf = self.PumpShapeApplyDeconvolutionFactor

        ffi = Filters().getFilterFunction(self.SIfilterIdlerType, self.SIfilterIdlerCenterWL, self.SIfilterIdlerFWHM)
        ffs = Filters().getFilterFunction(self.SIfilterSignalType, self.SIfilterSignalCenterWL, self.SIfilterSignalFWHM)
        spectralfilters = [ffs, ffi]

        nxfunc = RefractiveIndex().getSingleIDX(self.CrystalMaterial, "X", self.CrystalNX)
        nyfunc = RefractiveIndex().getSingleIDX(self.CrystalMaterial, "Y", self.CrystalNY)
        nzfunc = RefractiveIndex().getSingleIDX(self.CrystalMaterial, "Z", self.CrystalNZ)
        refidxfunc = [nxfunc, nyfunc, nzfunc]

        Tvec = numpy.arange(T, T + 1, 2)
        ls=0
        li=0
        unused=0
        [ls, li, unused] = PMC().getSI_wl_varT(pwl, PP, Tvec, refidxfunc, m)
        signalrange = numpy.linspace(ls - JSIwlrange / 2, ls + JSIwlrange / 2, JSIresolution)
        idlerrange = numpy.linspace(li - JSIwlrange / 2, li + JSIwlrange / 2, JSIresolution)
        
        taurange = numpy.linspace(taumin, taumax, fwhmres)

        sigfwhm=[]
        idfwhm=[]
        [sigfwhm, idfwhm] = JSI().getFWHMvstau(pwl, signalrange, idlerrange, T, PP, m, cl, taurange, refidxfunc, spectralfilters, JSIresolution, pumpshape, decprec, usetaucf)
        
        
        # plot
        # init plot window
        pltwndidx = self.plotwindowcount
        self.open_new_plot_window()
        pltwnd = self.pltwindowlist[pltwndidx]
        ploterror=False
        if len(sigfwhm) == len(taurange):
            pltwnd.ax.plot(taurange * 10 ** 12, sigfwhm, lw=2, label='Signal FWHM')
        else:
            ploterror=True
        if len(idfwhm) == len(taurange):
            pltwnd.ax.plot(taurange * 10 ** 12, idfwhm, lw=2, label='Idler FWHM')
        else:
            ploterror=True
        
        if ploterror:
            pltwnd.close()
            msgbox=QMessageBox()
            #msgbox.setText('Error: could not find enough values for FWHM. Try increasing JSI resolution or JSI range.')
            msgbox.exec()
        else:
            pltwnd.ax.set_xlabel('Pump pulse width [ps]')
            pltwnd.ax.set_ylabel('FWHM [nm]')
            pltwnd.ax.set_title('FWHM of marginal sectra')
            #pltwnd.ax.annotate('Visibility: {0:.3f} \n FWHM: {1:.3f}ps'.format(vis,fwhm*10**12), xy=(0.01, 0.01),
                            #xycoords='axes fraction')
            pltwnd.ax.legend()
            pltwnd.canvas.draw()

    def GetEffectivePolingPeriod(self):
        nxfunc = RefractiveIndex().getSingleIDX(self.CrystalMaterial, "X", self.CrystalNX)
        nyfunc = RefractiveIndex().getSingleIDX(self.CrystalMaterial, "Y", self.CrystalNY)
        nzfunc = RefractiveIndex().getSingleIDX(self.CrystalMaterial, "Z", self.CrystalNZ)
        refidxfunc = [nxfunc, nyfunc, nzfunc]
        m=self.QPMOrder
        Tcp=self.CrystalTempSingle
        PPguess=self.CrystalPolingPeriodSingle
        lp=self.PumpWlSingle
        PP=JSI().GetEffectivePP(m, Tcp,PPguess,lp,refidxfunc)
        self.ui_CrystalPolingPeriodsingleSB.setValue(PP*10**6)

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
        self.CrystalLengthFrom = self.ui_CrystalLengthfromSB.value()*10**(-3)
        self.CrystalLengthSingle = self.ui_CrystalLengthsingleSB.value()*10**(-3)
        self.CrystalLengthTo = self.ui_CrystalLengthtoSB.value()*10**(-3)

        self.PumpWlSingle = self.ui_pumpwlsingleSB.value() * 10 ** (-9)
        self.PumpWlFrom = self.ui_pumpwlfromSB.value() * 10 ** (-9)
        self.PumpWlTo = self.ui_pumpwltoSB.value() * 10 ** (-9)
        self.PulsewidthSingle = self.ui_pulsewidthsingleSB.value() * 10 ** (-12)
        self.PulsewidthFrom = self.ui_pulsewidthfromSB.value() * 10 ** (-12)
        self.PulsewidthTo = self.ui_pulsewidthtoSB.value() * 10 ** (-12)
        self.PumpShape = self.ui_pumpShapeCB.currentText()
        self.PumpShapeApplyDeconvolutionFactor = self.ui_pumpShapeCorrectionFactorCheckBox.isChecked()

        self.SIWlFrom = self.ui_siwlfromSB.value() * 10 ** (-9)
        self.SIWlSingle = self.ui_siwlsingleSB.value() * 10 ** (-9)
        self.SIWlTo = self.ui_siwltoSB.value() * 10 ** (-9)

        self.QPMOrder = self.ui_PlotPMCSBQPMorder.value()

        self.JSIwlRange = self.ui_PlotJSI_Wlrange_SB.value() * 10 ** (-9)
        self.JSIresolution = self.ui_PlotJSI_WLresolution_SB.value()

        self.PurityWLresolution = self.ui_Purity_WLresolution_SB.value()
        self.PurityWLrange = self.ui_Purity_WLrange_SB.value() * 10 ** (-9)
        self.PurityTauresolution = self.ui_Purity_Tauresolution_SB.value()

        self.SIfilterIdlerType = self.ui_SIfilterIdlerType_CB.currentText()
        self.SIfilterSignalType = self.ui_SIfilterSignalType_CB.currentText()
        self.SIfilterIdlerCenterWL = self.ui_SIfilterIdlerCenterWL_SB.value() * 10 ** (-9)
        self.SIfilterSignalCenterWL = self.ui_SIfilterSignalCenterWL_SB.value() * 10 ** (-9)
        self.SIfilterIdlerFWHM = self.ui_SIfilterIdlerFWHM_SB.value() * 10 ** (-9)
        self.SIfilterSignalFWHM = self.ui_SIfilterSignalFWHM_SB.value() * 10 ** (-9)

        self.HOMresolution = self.ui_HOM_Vis_Resolution_SB.value()
        self.HOMdelayrange = self.ui_HOM_Vis_DelayRange_SB.value()*10**(-12)
        
        self.fwhmres = self.ui_PlotFWHMresolution_SB.value()
        self.fwhmprecision = self.ui_PlotFWHMprecision_SB.value()
        
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

        self.config.set("Crystal Length from", self.CrystalLengthFrom)
        self.config.set("Crystal Length single", self.CrystalLengthSingle)
        self.config.set("Crystal Length to", self.CrystalLengthTo)

        self.config.set("JSI wavelength range", self.JSIwlRange)
        self.config.set("JSI resolution", self.JSIresolution)

        self.config.set("Purity wavelength resolution", self.PurityWLresolution)
        self.config.set("Purity wavelength range", self.PurityWLrange)
        self.config.set("Purity tau resolution", self.PurityTauresolution)

        self.config.set("SI filter Idler Type", self.SIfilterIdlerType)
        self.config.set("SI filter Signal Type", self.SIfilterSignalType)
        self.config.set("SI filter Idler center wavelength", self.SIfilterIdlerCenterWL)
        self.config.set("SI filter Signal center wavelength", self.SIfilterSignalCenterWL)
        self.config.set("SI filter Idler FWHM", self.SIfilterIdlerFWHM)
        self.config.set("SI filter Signal FWHM", self.SIfilterSignalFWHM)

        self.config.set("HOM interference plot resolution", self.HOMresolution)
        self.config.set("HOM interference plot range", self.HOMdelayrange)
        
        self.config.set('FWHM plot resolution', self.fwhmres)
        self.config.set('FWHM decimal precision', self.fwhmprecision)

    def showWindow(self):
        g = GUI()
        sys.exit(app.exec_())

    def open_new_plot_window(self):
        print(self.plotwindowcount)
        self.pltwindowlist.append(PlotWindow())
        self.pltwindowlist[self.plotwindowcount].show()
        self.plotwindowcount = self.plotwindowcount + 1

    def toggle_usecf(self):
        self.PumpShapeApplyDeconvolutionFactor = self.ui_pumpShapeCorrectionFactorCheckBox.isChecked()
        Constants().usetaucf(self.PumpShapeApplyDeconvolutionFactor)

    def MouseHoverEnter(self, str):
        if str == 'Plot PMC vs T':
            self.ui_pumpwlsingleSB.setStyleSheet(self.HighlightedDoubleSpinBox)
            self.ui_CrystalPolingPeriodsingleSB.setStyleSheet(self.HighlightedDoubleSpinBox)
            self.ui_CrystalTtoSB.setStyleSheet(self.HighlightedDoubleSpinBox)
            self.ui_CrystalTfromSB.setStyleSheet(self.HighlightedDoubleSpinBox)
            self.ui_PlotPMCSBQPMorder.setStyleSheet(self.HighlightedSpinBox)
            self.ui_CrystalNXComboBox.setStyleSheet(self.HighlightedDoubleSpinBox)
            self.ui_CrystalNYComboBox.setStyleSheet(self.HighlightedDoubleSpinBox)
            self.ui_CrystalNZComboBox.setStyleSheet(self.HighlightedDoubleSpinBox)
            self.ui_CrystalMaterialComboBox.setStyleSheet(self.HighlightedComboBox)
            self.ui_CrystalNXComboBox.setStyleSheet(self.HighlightedComboBox)
            self.ui_CrystalNYComboBox.setStyleSheet(self.HighlightedComboBox)
            self.ui_CrystalNZComboBox.setStyleSheet(self.HighlightedComboBox)
        elif str == 'Plot PMC vs PP':
            self.ui_pumpwlsingleSB.setStyleSheet(self.HighlightedDoubleSpinBox)
            self.ui_CrystalPolingPeriodfromSB.setStyleSheet(self.HighlightedDoubleSpinBox)
            self.ui_CrystalPolingPeriodtoSB.setStyleSheet(self.HighlightedDoubleSpinBox)
            self.ui_CrystalTsingleSB.setStyleSheet(self.HighlightedDoubleSpinBox)
            self.ui_PlotPMCSBQPMorder.setStyleSheet(self.HighlightedSpinBox)
            self.ui_CrystalNXComboBox.setStyleSheet(self.HighlightedDoubleSpinBox)
            self.ui_CrystalNYComboBox.setStyleSheet(self.HighlightedDoubleSpinBox)
            self.ui_CrystalNZComboBox.setStyleSheet(self.HighlightedDoubleSpinBox)
            self.ui_CrystalMaterialComboBox.setStyleSheet(self.HighlightedComboBox)
            self.ui_CrystalNXComboBox.setStyleSheet(self.HighlightedComboBox)
            self.ui_CrystalNYComboBox.setStyleSheet(self.HighlightedComboBox)
            self.ui_CrystalNZComboBox.setStyleSheet(self.HighlightedComboBox)
        elif str == 'Plot refractive indices vs T':
            self.ui_CrystalTfromSB.setStyleSheet(self.HighlightedDoubleSpinBox)
            self.ui_CrystalTtoSB.setStyleSheet(self.HighlightedDoubleSpinBox)
            self.ui_pumpwlsingleSB.setStyleSheet(self.HighlightedDoubleSpinBox)
            self.ui_layoutPlotRefractiveIndexScrollAreaWidget.setStyleSheet(self.HighlightedWidget)
        elif str == 'Plot refractive indices vs wl':
            self.ui_CrystalTsingleSB.setStyleSheet(self.HighlightedDoubleSpinBox)
            self.ui_pumpwlfromSB.setStyleSheet(self.HighlightedDoubleSpinBox)
            self.ui_pumpwltoSB.setStyleSheet(self.HighlightedDoubleSpinBox)
            self.ui_layoutPlotRefractiveIndexScrollAreaWidget.setStyleSheet(self.HighlightedWidget)
        elif str == 'Get effective PP':
            self.ui_CrystalMaterialComboBox.setStyleSheet(self.HighlightedComboBox)
            self.ui_CrystalNXComboBox.setStyleSheet(self.HighlightedComboBox)
            self.ui_CrystalNYComboBox.setStyleSheet(self.HighlightedComboBox)
            self.ui_CrystalNZComboBox.setStyleSheet(self.HighlightedComboBox)
            self.ui_PlotPMCSBQPMorder.setStyleSheet(self.HighlightedSpinBox)
            self.ui_CrystalTsingleSB.setStyleSheet(self.HighlightedDoubleSpinBox)
            self.ui_CrystalPolingPeriodsingleSB.setStyleSheet(self.HighlightedDoubleSpinBox)
            self.ui_pumpwlsingleSB.setStyleSheet(self.HighlightedDoubleSpinBox)
        elif str == 'Plot JSA JSI':
            self.ui_pumpwlsingleSB.setStyleSheet(self.HighlightedDoubleSpinBox)
            self.ui_CrystalPolingPeriodsingleSB.setStyleSheet(self.HighlightedDoubleSpinBox)
            self.ui_CrystalLengthsingleSB.setStyleSheet(self.HighlightedDoubleSpinBox)
            self.ui_CrystalTsingleSB.setStyleSheet(self.HighlightedDoubleSpinBox)
            self.ui_PlotPMCSBQPMorder.setStyleSheet(self.HighlightedSpinBox)
            self.ui_pulsewidthsingleSB.setStyleSheet(self.HighlightedDoubleSpinBox)
            self.ui_PlotJSI_Wlrange_SB.setStyleSheet(self.HighlightedDoubleSpinBox)
            self.ui_SIfilterIdlerType_CB.setStyleSheet(self.HighlightedComboBox)
            self.ui_SIfilterIdlerCenterWL_SB.setStyleSheet(self.HighlightedDoubleSpinBox)
            self.ui_SIfilterIdlerFWHM_SB.setStyleSheet(self.HighlightedDoubleSpinBox)
            self.ui_SIfilterSignalType_CB.setStyleSheet(self.HighlightedComboBox)
            self.ui_SIfilterSignalCenterWL_SB.setStyleSheet(self.HighlightedDoubleSpinBox)
            self.ui_SIfilterSignalFWHM_SB.setStyleSheet(self.HighlightedDoubleSpinBox)
            self.ui_PlotJSI_plotJSARadioButton.setStyleSheet(self.HighlightedRadioButton)
            self.ui_PlotJSI_plotJSIRadioButton.setStyleSheet(self.HighlightedRadioButton)
            self.ui_pumpShapeCB.setStyleSheet(self.HighlightedComboBox)
            self.ui_CrystalMaterialComboBox.setStyleSheet(self.HighlightedComboBox)
            self.ui_CrystalNXComboBox.setStyleSheet(self.HighlightedComboBox)
            self.ui_CrystalNYComboBox.setStyleSheet(self.HighlightedComboBox)
            self.ui_CrystalNZComboBox.setStyleSheet(self.HighlightedComboBox)
            self.ui_PlotJSI_WLresolution_SB.setStyleSheet(self.HighlightedSpinBox)
        elif str == 'Plot Purity vs tau':
            self.ui_pumpwlsingleSB.setStyleSheet(self.HighlightedDoubleSpinBox)
            self.ui_CrystalPolingPeriodsingleSB.setStyleSheet(self.HighlightedDoubleSpinBox)
            self.ui_CrystalLengthsingleSB.setStyleSheet(self.HighlightedDoubleSpinBox)
            self.ui_CrystalTsingleSB.setStyleSheet(self.HighlightedDoubleSpinBox)
            self.ui_PlotPMCSBQPMorder.setStyleSheet(self.HighlightedSpinBox)
            self.ui_pulsewidthfromSB.setStyleSheet(self.HighlightedDoubleSpinBox)
            self.ui_pulsewidthtoSB.setStyleSheet(self.HighlightedDoubleSpinBox)
            self.ui_SIfilterIdlerType_CB.setStyleSheet(self.HighlightedComboBox)
            self.ui_SIfilterIdlerCenterWL_SB.setStyleSheet(self.HighlightedDoubleSpinBox)
            self.ui_SIfilterIdlerFWHM_SB.setStyleSheet(self.HighlightedDoubleSpinBox)
            self.ui_SIfilterSignalType_CB.setStyleSheet(self.HighlightedComboBox)
            self.ui_SIfilterSignalCenterWL_SB.setStyleSheet(self.HighlightedDoubleSpinBox)
            self.ui_SIfilterSignalFWHM_SB.setStyleSheet(self.HighlightedDoubleSpinBox)
            self.ui_pumpShapeCB.setStyleSheet(self.HighlightedComboBox)
            self.ui_CrystalMaterialComboBox.setStyleSheet(self.HighlightedComboBox)
            self.ui_CrystalNXComboBox.setStyleSheet(self.HighlightedComboBox)
            self.ui_CrystalNYComboBox.setStyleSheet(self.HighlightedComboBox)
            self.ui_CrystalNZComboBox.setStyleSheet(self.HighlightedComboBox)
            self.ui_Purity_WLresolution_SB.setStyleSheet(self.HighlightedSpinBox)
            self.ui_Purity_Tauresolution_SB.setStyleSheet(self.HighlightedSpinBox)
            self.ui_Purity_WLrange_SB.setStyleSheet(self.HighlightedSpinBox)
        elif str == 'Plot Purity vs lp':
            pass
        elif str == 'Plot Purity vs L':
            self.ui_pumpwlsingleSB.setStyleSheet(self.HighlightedDoubleSpinBox)
            self.ui_CrystalPolingPeriodsingleSB.setStyleSheet(self.HighlightedDoubleSpinBox)
            self.ui_CrystalLengthfromSB.setStyleSheet(self.HighlightedDoubleSpinBox)
            self.ui_CrystalLengthtoSB.setStyleSheet(self.HighlightedDoubleSpinBox)
            self.ui_CrystalTsingleSB.setStyleSheet(self.HighlightedDoubleSpinBox)
            self.ui_PlotPMCSBQPMorder.setStyleSheet(self.HighlightedSpinBox)
            self.ui_pulsewidthsingleSB.setStyleSheet(self.HighlightedDoubleSpinBox)
            self.ui_SIfilterIdlerType_CB.setStyleSheet(self.HighlightedComboBox)
            self.ui_SIfilterIdlerCenterWL_SB.setStyleSheet(self.HighlightedDoubleSpinBox)
            self.ui_SIfilterIdlerFWHM_SB.setStyleSheet(self.HighlightedDoubleSpinBox)
            self.ui_SIfilterSignalType_CB.setStyleSheet(self.HighlightedComboBox)
            self.ui_SIfilterSignalCenterWL_SB.setStyleSheet(self.HighlightedDoubleSpinBox)
            self.ui_SIfilterSignalFWHM_SB.setStyleSheet(self.HighlightedDoubleSpinBox)
            self.ui_pumpShapeCB.setStyleSheet(self.HighlightedComboBox)
            self.ui_CrystalMaterialComboBox.setStyleSheet(self.HighlightedComboBox)
            self.ui_CrystalNXComboBox.setStyleSheet(self.HighlightedComboBox)
            self.ui_CrystalNYComboBox.setStyleSheet(self.HighlightedComboBox)
            self.ui_CrystalNZComboBox.setStyleSheet(self.HighlightedComboBox)
            self.ui_Purity_WLresolution_SB.setStyleSheet(self.HighlightedSpinBox)
            self.ui_Purity_Tauresolution_SB.setStyleSheet(self.HighlightedSpinBox)
            self.ui_Purity_WLrange_SB.setStyleSheet(self.HighlightedSpinBox)
        elif str == 'Plot Purity vs L and tau':
            self.ui_pumpwlsingleSB.setStyleSheet(self.HighlightedDoubleSpinBox)
            self.ui_CrystalPolingPeriodsingleSB.setStyleSheet(self.HighlightedDoubleSpinBox)
            self.ui_CrystalLengthfromSB.setStyleSheet(self.HighlightedDoubleSpinBox)
            self.ui_CrystalLengthtoSB.setStyleSheet(self.HighlightedDoubleSpinBox)
            self.ui_CrystalTsingleSB.setStyleSheet(self.HighlightedDoubleSpinBox)
            self.ui_PlotPMCSBQPMorder.setStyleSheet(self.HighlightedSpinBox)
            self.ui_pulsewidthfromSB.setStyleSheet(self.HighlightedDoubleSpinBox)
            self.ui_pulsewidthtoSB.setStyleSheet(self.HighlightedDoubleSpinBox)
            self.ui_SIfilterIdlerType_CB.setStyleSheet(self.HighlightedComboBox)
            self.ui_SIfilterIdlerCenterWL_SB.setStyleSheet(self.HighlightedDoubleSpinBox)
            self.ui_SIfilterIdlerFWHM_SB.setStyleSheet(self.HighlightedDoubleSpinBox)
            self.ui_SIfilterSignalType_CB.setStyleSheet(self.HighlightedComboBox)
            self.ui_SIfilterSignalCenterWL_SB.setStyleSheet(self.HighlightedDoubleSpinBox)
            self.ui_SIfilterSignalFWHM_SB.setStyleSheet(self.HighlightedDoubleSpinBox)
            self.ui_pumpShapeCB.setStyleSheet(self.HighlightedComboBox)
            self.ui_CrystalMaterialComboBox.setStyleSheet(self.HighlightedComboBox)
            self.ui_CrystalNXComboBox.setStyleSheet(self.HighlightedComboBox)
            self.ui_CrystalNYComboBox.setStyleSheet(self.HighlightedComboBox)
            self.ui_CrystalNZComboBox.setStyleSheet(self.HighlightedComboBox)
            self.ui_Purity_WLresolution_SB.setStyleSheet(self.HighlightedSpinBox)
            self.ui_Purity_Tauresolution_SB.setStyleSheet(self.HighlightedSpinBox)
            self.ui_Purity_WLrange_SB.setStyleSheet(self.HighlightedSpinBox)
        elif str == 'Plot Tcp vs PP':
            self.ui_CrystalPolingPeriodfromSB.setStyleSheet(self.HighlightedDoubleSpinBox)
            self.ui_CrystalPolingPeriodtoSB.setStyleSheet(self.HighlightedDoubleSpinBox)
            self.ui_CrystalTsingleSB.setStyleSheet(self.HighlightedDoubleSpinBox)
            self.ui_pumpwlsingleSB.setStyleSheet(self.HighlightedDoubleSpinBox)
            self.ui_PlotPMCSBQPMorder.setStyleSheet(self.HighlightedSpinBox)
            self.ui_CrystalMaterialComboBox.setStyleSheet(self.HighlightedComboBox)
            self.ui_CrystalNXComboBox.setStyleSheet(self.HighlightedComboBox)
            self.ui_CrystalNYComboBox.setStyleSheet(self.HighlightedComboBox)
            self.ui_CrystalNZComboBox.setStyleSheet(self.HighlightedComboBox)
        elif str == 'Plot Tcp vs pump wavelength':
            self.ui_CrystalPolingPeriodsingleSB.setStyleSheet(self.HighlightedDoubleSpinBox)
            self.ui_CrystalTsingleSB.setStyleSheet(self.HighlightedDoubleSpinBox)
            self.ui_pumpwlfromSB.setStyleSheet(self.HighlightedDoubleSpinBox)
            self.ui_pumpwltoSB.setStyleSheet(self.HighlightedDoubleSpinBox)
            self.ui_PlotPMCSBQPMorder.setStyleSheet(self.HighlightedSpinBox)
            self.ui_CrystalMaterialComboBox.setStyleSheet(self.HighlightedComboBox)
            self.ui_CrystalNXComboBox.setStyleSheet(self.HighlightedComboBox)
            self.ui_CrystalNYComboBox.setStyleSheet(self.HighlightedComboBox)
            self.ui_CrystalNZComboBox.setStyleSheet(self.HighlightedComboBox)
        elif str == 'Plot HOM Visibility':
            self.ui_pumpwlsingleSB.setStyleSheet(self.HighlightedDoubleSpinBox)
            self.ui_CrystalTsingleSB.setStyleSheet(self.HighlightedDoubleSpinBox)
            self.ui_CrystalPolingPeriodsingleSB.setStyleSheet(self.HighlightedDoubleSpinBox)
            self.ui_PlotPMCSBQPMorder.setStyleSheet(self.HighlightedSpinBox)
            self.ui_pulsewidthsingleSB.setStyleSheet(self.HighlightedDoubleSpinBox)
            self.ui_CrystalLengthsingleSB.setStyleSheet(self.HighlightedDoubleSpinBox)
            self.ui_pumpShapeCB.setStyleSheet(self.HighlightedComboBox)
            self.ui_HOM_Vis_DelayRange_SB.setStyleSheet(self.HighlightedDoubleSpinBox)
            self.ui_HOM_Vis_Resolution_SB.setStyleSheet(self.HighlightedSpinBox)
            self.ui_PlotJSI_WLresolution_SB.setStyleSheet(self.HighlightedSpinBox)
            self.ui_PlotJSI_Wlrange_SB.setStyleSheet(self.HighlightedDoubleSpinBox)
            self.ui_SIfilterIdlerType_CB.setStyleSheet(self.HighlightedComboBox)
            self.ui_SIfilterIdlerCenterWL_SB.setStyleSheet(self.HighlightedDoubleSpinBox)
            self.ui_SIfilterIdlerFWHM_SB.setStyleSheet(self.HighlightedDoubleSpinBox)
            self.ui_SIfilterSignalType_CB.setStyleSheet(self.HighlightedComboBox)
            self.ui_SIfilterSignalCenterWL_SB.setStyleSheet(self.HighlightedDoubleSpinBox)
            self.ui_SIfilterSignalFWHM_SB.setStyleSheet(self.HighlightedDoubleSpinBox)
            self.ui_CrystalMaterialComboBox.setStyleSheet(self.HighlightedComboBox)
            self.ui_CrystalNXComboBox.setStyleSheet(self.HighlightedComboBox)
            self.ui_CrystalNYComboBox.setStyleSheet(self.HighlightedComboBox)
            self.ui_CrystalNZComboBox.setStyleSheet(self.HighlightedComboBox)
        elif str == 'Plot marginal FWHM vs tau':
            self.ui_pulsewidthfromSB.setStyleSheet(self.HighlightedDoubleSpinBox)
            self.ui_pulsewidthtoSB.setStyleSheet(self.HighlightedDoubleSpinBox)
            self.ui_pumpwlsingleSB.setStyleSheet(self.HighlightedDoubleSpinBox)
            self.ui_CrystalTsingleSB.setStyleSheet(self.HighlightedDoubleSpinBox)
            self.ui_PlotJSI_Wlrange_SB.setStyleSheet(self.HighlightedDoubleSpinBox)
            self.ui_PlotJSI_WLresolution_SB.setStyleSheet(self.HighlightedSpinBox)
            self.ui_PlotFWHMresolution_SB.setStyleSheet(self.HighlightedSpinBox)
            self.ui_PlotFWHMprecision_SB.setStyleSheet(self.HighlightedSpinBox)
            self.ui_SIfilterIdlerType_CB.setStyleSheet(self.HighlightedComboBox)
            self.ui_SIfilterIdlerCenterWL_SB.setStyleSheet(self.HighlightedDoubleSpinBox)
            self.ui_SIfilterIdlerFWHM_SB.setStyleSheet(self.HighlightedDoubleSpinBox)
            self.ui_SIfilterSignalType_CB.setStyleSheet(self.HighlightedComboBox)
            self.ui_SIfilterSignalCenterWL_SB.setStyleSheet(self.HighlightedDoubleSpinBox)
            self.ui_SIfilterSignalFWHM_SB.setStyleSheet(self.HighlightedDoubleSpinBox)
            self.ui_CrystalMaterialComboBox.setStyleSheet(self.HighlightedComboBox)
            self.ui_CrystalNXComboBox.setStyleSheet(self.HighlightedComboBox)
            self.ui_CrystalNYComboBox.setStyleSheet(self.HighlightedComboBox)
            self.ui_CrystalNZComboBox.setStyleSheet(self.HighlightedComboBox)
            self.ui_CrystalPolingPeriodsingleSB.setStyleSheet(self.HighlightedSpinBox)
            self.ui_PlotPMCSBQPMorder.setStyleSheet(self.HighlightedSpinBox)
            self.ui_CrystalLengthsingleSB.setStyleSheet(self.HighlightedSpinBox)
            self.ui_pumpShapeCB.setStyleSheet(self.HighlightedComboBox)
        elif str == 'Esimate filter losses':
            self.ui_pumpwlsingleSB.setStyleSheet(self.HighlightedDoubleSpinBox)
            self.ui_CrystalPolingPeriodsingleSB.setStyleSheet(self.HighlightedDoubleSpinBox)
            self.ui_CrystalLengthsingleSB.setStyleSheet(self.HighlightedDoubleSpinBox)
            self.ui_CrystalTsingleSB.setStyleSheet(self.HighlightedDoubleSpinBox)
            self.ui_PlotPMCSBQPMorder.setStyleSheet(self.HighlightedSpinBox)
            self.ui_pulsewidthsingleSB.setStyleSheet(self.HighlightedDoubleSpinBox)
            self.ui_PlotJSI_Wlrange_SB.setStyleSheet(self.HighlightedDoubleSpinBox)
            self.ui_SIfilterIdlerType_CB.setStyleSheet(self.HighlightedComboBox)
            self.ui_SIfilterIdlerCenterWL_SB.setStyleSheet(self.HighlightedDoubleSpinBox)
            self.ui_SIfilterIdlerFWHM_SB.setStyleSheet(self.HighlightedDoubleSpinBox)
            self.ui_SIfilterSignalType_CB.setStyleSheet(self.HighlightedComboBox)
            self.ui_SIfilterSignalCenterWL_SB.setStyleSheet(self.HighlightedDoubleSpinBox)
            self.ui_SIfilterSignalFWHM_SB.setStyleSheet(self.HighlightedDoubleSpinBox)
            self.ui_PlotJSI_plotJSARadioButton.setStyleSheet(self.HighlightedRadioButton)
            self.ui_PlotJSI_plotJSIRadioButton.setStyleSheet(self.HighlightedRadioButton)
            self.ui_pumpShapeCB.setStyleSheet(self.HighlightedComboBox)
            self.ui_CrystalMaterialComboBox.setStyleSheet(self.HighlightedComboBox)
            self.ui_CrystalNXComboBox.setStyleSheet(self.HighlightedComboBox)
            self.ui_CrystalNYComboBox.setStyleSheet(self.HighlightedComboBox)
            self.ui_CrystalNZComboBox.setStyleSheet(self.HighlightedComboBox)
            self.ui_PlotJSI_WLresolution_SB.setStyleSheet(self.HighlightedSpinBox)
        else:
            print('No Object name specified in QHoverPushButton!')

    def MouseHoverLeave(self, str):
        self.ui_pumpwlfromSB.setStyleSheet("")
        self.ui_pumpwlsingleSB.setStyleSheet("")
        self.ui_pumpwltoSB.setStyleSheet("")
        self.ui_pulsewidthfromSB.setStyleSheet("")
        self.ui_pulsewidthsingleSB.setStyleSheet("")
        self.ui_pulsewidthtoSB.setStyleSheet("")
        self.ui_CrystalPolingPeriodfromSB.setStyleSheet("")
        self.ui_CrystalPolingPeriodsingleSB.setStyleSheet("")
        self.ui_CrystalPolingPeriodtoSB.setStyleSheet("")
        self.ui_CrystalLengthfromSB.setStyleSheet("")
        self.ui_CrystalLengthsingleSB.setStyleSheet("")
        self.ui_CrystalLengthtoSB.setStyleSheet("")
        self.ui_CrystalTtoSB.setStyleSheet("")
        self.ui_CrystalTsingleSB.setStyleSheet("")
        self.ui_CrystalTfromSB.setStyleSheet("")
        self.ui_PlotPMCSBQPMorder.setStyleSheet("")
        self.ui_CrystalNXComboBox.setStyleSheet("")
        self.ui_CrystalNYComboBox.setStyleSheet("")
        self.ui_CrystalNZComboBox.setStyleSheet("")
        self.ui_CrystalMaterialComboBox.setStyleSheet("")
        self.ui_pumpShapeCB.setStyleSheet("")
        self.ui_SIfilterIdlerType_CB.setStyleSheet("")
        self.ui_SIfilterIdlerCenterWL_SB.setStyleSheet("")
        self.ui_SIfilterIdlerFWHM_SB.setStyleSheet("")
        self.ui_SIfilterSignalType_CB.setStyleSheet("")
        self.ui_SIfilterSignalCenterWL_SB.setStyleSheet("")
        self.ui_SIfilterSignalFWHM_SB.setStyleSheet("")
        self.ui_PlotJSI_plotJSARadioButton.setStyleSheet("")
        self.ui_PlotJSI_plotJSIRadioButton.setStyleSheet("")
        self.ui_PlotJSI_Wlrange_SB.setStyleSheet("")
        self.ui_PlotJSI_WLresolution_SB.setStyleSheet("")
        self.ui_Purity_WLresolution_SB.setStyleSheet("")
        self.ui_Purity_Tauresolution_SB.setStyleSheet("")
        self.ui_Purity_WLrange_SB.setStyleSheet("")
        self.ui_layoutPlotRefractiveIndexScrollAreaWidget.setStyleSheet("")
        self.ui_HOM_Vis_DelayRange_SB.setStyleSheet("")
        self.ui_HOM_Vis_Resolution_SB.setStyleSheet("")
        self.ui_PlotFWHMresolution_SB.setStyleSheet("")
        self.ui_PlotFWHMprecision_SB.setStyleSheet("")

class PlotWindow(QWidget):
    def __init__(self):
        QWidget.__init__(self)
        # self.layout = QGridLayout()
        self.layout = QVBoxLayout()
        self.setLayout(self.layout)
        self.fig = figure(facecolor="white", figsize=(8.75 * 1.2, 5 * 1.2))
        self.ax = self.fig.add_subplot(111)
        self.ax.grid()
        self.canvas = FigureCanvas(self.fig)
        self.canvas.setParent(self)
        self.toolbar = NavigationToolbar(self.canvas, self)
        # self.addToolBar(self.toolbar)
        self.layout.addWidget(self.canvas)
        self.layout.addWidget(self.toolbar)

    def open(self):
        self.show(self)


if __name__ == '__main__':
    app = QApplication(sys.argv)
    g = GUI()
    sys.exit(app.exec_())
