#!/usr/bin/python3

from PyQt5.QtWidgets import QPushButton
from PyQt5 import QtCore
from PyQt5.QtCore import QObject, pyqtSignal

#reimplementation of the QPushButton
#change: emit signal when moush hovers over button
class QHoverPushButton(QPushButton):
    mouseentersignal = pyqtSignal(str)
    mouseleavesignal = pyqtSignal(str)

    def __init__(self, parent=None):
        QPushButton.__init__(self, parent)
        self.setMouseTracking(True)


    def enterEvent(self, a0: QtCore.QEvent):
        self.mouseentersignal.emit(self.objectName())


    def leaveEvent(self, a0: QtCore.QEvent):
        self.mouseleavesignal.emit(self.objectName())
