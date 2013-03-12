# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'explorer.ui'
#
# Created: Tue Mar 12 12:23:16 2013
#      by: PyQt4 UI code generator 4.9.1
#
# WARNING! All changes made in this file will be lost!

from PyQt4 import QtCore, QtGui

try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    _fromUtf8 = lambda s: s

class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        MainWindow.setObjectName(_fromUtf8("MainWindow"))
        MainWindow.resize(800, 600)
        self.centralwidget = QtGui.QWidget(MainWindow)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Maximum, QtGui.QSizePolicy.Maximum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.centralwidget.sizePolicy().hasHeightForWidth())
        self.centralwidget.setSizePolicy(sizePolicy)
        self.centralwidget.setObjectName(_fromUtf8("centralwidget"))
        self.horizontalLayoutWidget = QtGui.QWidget(self.centralwidget)
        self.horizontalLayoutWidget.setGeometry(QtCore.QRect(10, 0, 781, 551))
        self.horizontalLayoutWidget.setObjectName(_fromUtf8("horizontalLayoutWidget"))
        self.horizontalLayout = QtGui.QHBoxLayout(self.horizontalLayoutWidget)
        self.horizontalLayout.setMargin(0)
        self.horizontalLayout.setObjectName(_fromUtf8("horizontalLayout"))
        self.show3d = Show3DWithSlider(self.horizontalLayoutWidget)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Preferred, QtGui.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.show3d.sizePolicy().hasHeightForWidth())
        self.show3d.setSizePolicy(sizePolicy)
        self.show3d.setObjectName(_fromUtf8("show3d"))
        self.horizontalLayout.addWidget(self.show3d)
        self.gridLayout = QtGui.QGridLayout()
        self.gridLayout.setObjectName(_fromUtf8("gridLayout"))
        self.btn_MC_chain = QtGui.QPushButton(self.horizontalLayoutWidget)
        self.btn_MC_chain.setObjectName(_fromUtf8("btn_MC_chain"))
        self.gridLayout.addWidget(self.btn_MC_chain, 3, 0, 1, 1)
        self.btn_db_sample_min = QtGui.QPushButton(self.horizontalLayoutWidget)
        self.btn_db_sample_min.setObjectName(_fromUtf8("btn_db_sample_min"))
        self.gridLayout.addWidget(self.btn_db_sample_min, 2, 0, 1, 1)
        spacerItem = QtGui.QSpacerItem(20, 40, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Expanding)
        self.gridLayout.addItem(spacerItem, 1, 0, 1, 1)
        self.list_replicas = QtGui.QListWidget(self.horizontalLayoutWidget)
        self.list_replicas.setObjectName(_fromUtf8("list_replicas"))
        self.gridLayout.addWidget(self.list_replicas, 0, 0, 1, 1)
        self.horizontalLayout.addLayout(self.gridLayout)
        MainWindow.setCentralWidget(self.centralwidget)
        self.menubar = QtGui.QMenuBar(MainWindow)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 800, 25))
        self.menubar.setObjectName(_fromUtf8("menubar"))
        MainWindow.setMenuBar(self.menubar)
        self.statusbar = QtGui.QStatusBar(MainWindow)
        self.statusbar.setObjectName(_fromUtf8("statusbar"))
        MainWindow.setStatusBar(self.statusbar)

        self.retranslateUi(MainWindow)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

    def retranslateUi(self, MainWindow):
        MainWindow.setWindowTitle(QtGui.QApplication.translate("MainWindow", "MainWindow", None, QtGui.QApplication.UnicodeUTF8))
        self.btn_MC_chain.setText(QtGui.QApplication.translate("MainWindow", "MC chain", None, QtGui.QApplication.UnicodeUTF8))
        self.btn_db_sample_min.setToolTip(QtGui.QApplication.translate("MainWindow", "<html><head/><body><p>select a starting point for the monte carlo chain randomly</p></body></html>", None, QtGui.QApplication.UnicodeUTF8))
        self.btn_db_sample_min.setText(QtGui.QApplication.translate("MainWindow", "get coords from minima", None, QtGui.QApplication.UnicodeUTF8))

from pygmin.gui.show3d import Show3DWithSlider
