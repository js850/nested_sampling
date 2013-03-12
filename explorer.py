# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'explorer.ui'
#
# Created: Tue Mar 12 15:21:05 2013
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
        self.verticalLayout_4 = QtGui.QVBoxLayout()
        self.verticalLayout_4.setContentsMargins(10, -1, -1, -1)
        self.verticalLayout_4.setObjectName(_fromUtf8("verticalLayout_4"))
        self.list_replicas = QtGui.QListWidget(self.horizontalLayoutWidget)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.list_replicas.sizePolicy().hasHeightForWidth())
        self.list_replicas.setSizePolicy(sizePolicy)
        self.list_replicas.setObjectName(_fromUtf8("list_replicas"))
        self.verticalLayout_4.addWidget(self.list_replicas)
        spacerItem = QtGui.QSpacerItem(20, 40, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Preferred)
        self.verticalLayout_4.addItem(spacerItem)
        self.formLayout_6 = QtGui.QFormLayout()
        self.formLayout_6.setFieldGrowthPolicy(QtGui.QFormLayout.AllNonFixedFieldsGrow)
        self.formLayout_6.setContentsMargins(-1, 10, -1, -1)
        self.formLayout_6.setObjectName(_fromUtf8("formLayout_6"))
        self.label = QtGui.QLabel(self.horizontalLayoutWidget)
        self.label.setObjectName(_fromUtf8("label"))
        self.formLayout_6.setWidget(2, QtGui.QFormLayout.LabelRole, self.label)
        self.lineEdit_Emax = QtGui.QLineEdit(self.horizontalLayoutWidget)
        self.lineEdit_Emax.setObjectName(_fromUtf8("lineEdit_Emax"))
        self.formLayout_6.setWidget(2, QtGui.QFormLayout.FieldRole, self.lineEdit_Emax)
        self.btn_db_sample_min = QtGui.QPushButton(self.horizontalLayoutWidget)
        self.btn_db_sample_min.setObjectName(_fromUtf8("btn_db_sample_min"))
        self.formLayout_6.setWidget(4, QtGui.QFormLayout.FieldRole, self.btn_db_sample_min)
        self.btn_MC_chain = QtGui.QPushButton(self.horizontalLayoutWidget)
        self.btn_MC_chain.setObjectName(_fromUtf8("btn_MC_chain"))
        self.formLayout_6.setWidget(5, QtGui.QFormLayout.FieldRole, self.btn_MC_chain)
        self.lineEdit_stepsize = QtGui.QLineEdit(self.horizontalLayoutWidget)
        self.lineEdit_stepsize.setObjectName(_fromUtf8("lineEdit_stepsize"))
        self.formLayout_6.setWidget(1, QtGui.QFormLayout.FieldRole, self.lineEdit_stepsize)
        self.label_2 = QtGui.QLabel(self.horizontalLayoutWidget)
        self.label_2.setObjectName(_fromUtf8("label_2"))
        self.formLayout_6.setWidget(1, QtGui.QFormLayout.LabelRole, self.label_2)
        self.label_3 = QtGui.QLabel(self.horizontalLayoutWidget)
        self.label_3.setObjectName(_fromUtf8("label_3"))
        self.formLayout_6.setWidget(0, QtGui.QFormLayout.LabelRole, self.label_3)
        self.lineEdit_mciter = QtGui.QLineEdit(self.horizontalLayoutWidget)
        self.lineEdit_mciter.setObjectName(_fromUtf8("lineEdit_mciter"))
        self.formLayout_6.setWidget(0, QtGui.QFormLayout.FieldRole, self.lineEdit_mciter)
        self.verticalLayout_4.addLayout(self.formLayout_6)
        self.horizontalLayout.addLayout(self.verticalLayout_4)
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
        self.label.setToolTip(QtGui.QApplication.translate("MainWindow", "<html><head/><body><p>set Emax manually</p></body></html>", None, QtGui.QApplication.UnicodeUTF8))
        self.label.setText(QtGui.QApplication.translate("MainWindow", "Emax", None, QtGui.QApplication.UnicodeUTF8))
        self.btn_db_sample_min.setToolTip(QtGui.QApplication.translate("MainWindow", "<html><head/><body><p>select a starting point for the monte carlo chain randomly</p></body></html>", None, QtGui.QApplication.UnicodeUTF8))
        self.btn_db_sample_min.setText(QtGui.QApplication.translate("MainWindow", "get coords from minima", None, QtGui.QApplication.UnicodeUTF8))
        self.btn_MC_chain.setText(QtGui.QApplication.translate("MainWindow", "MC chain", None, QtGui.QApplication.UnicodeUTF8))
        self.label_2.setToolTip(QtGui.QApplication.translate("MainWindow", "<html><head/><body><p>step size for the Monte Carlo chain</p></body></html>", None, QtGui.QApplication.UnicodeUTF8))
        self.label_2.setText(QtGui.QApplication.translate("MainWindow", "stepsize", None, QtGui.QApplication.UnicodeUTF8))
        self.label_3.setToolTip(QtGui.QApplication.translate("MainWindow", "<html><head/><body><p>number of iterations in the Monte Carlo chain</p></body></html>", None, QtGui.QApplication.UnicodeUTF8))
        self.label_3.setText(QtGui.QApplication.translate("MainWindow", "mciter", None, QtGui.QApplication.UnicodeUTF8))

from pygmin.gui.show3d import Show3DWithSlider
