import sys
import numpy as np

from bh_sampling import NestedSamplingBS

from explorer import Ui_MainWindow
from nested_sampling import Replica

from PyQt4 import QtCore, QtGui

class QReplicaInList(QtGui.QListWidgetItem):
    
    def __init__(self, replica):
        text="%.4f"%(replica.energy)
        QtGui.QListWidgetItem.__init__(self, text)
        self.r = replica
        
    def __lt__(self, item2):
        #sort the energies in the list highest to lowest
        return self.r.energy > item2.r.energy


class NSGUI(QtGui.QMainWindow):
    def __init__(self, app, system, nested_sampling):
        QtGui.QWidget.__init__(self)
        self.ui = Ui_MainWindow()
        self.ui.setupUi(self)
        
        
        self.app = app
        self.system = system    
        self.nested_sampling = nested_sampling
        
        self.show3d = self.ui.show3d
        self.ui.show3d.setSystem(system)
        
        self.update_list()
        
        from OpenGL.GLUT import glutInit
        glutInit()
        
        self.selected = None
        self.Emax = max([r.energy for r in self.nested_sampling.replicas])
    
    def update_list(self):
        self.ui.list_replicas.clear()
        for r in self.nested_sampling.replicas:
            item = QReplicaInList(r)
            self.ui.list_replicas.addItem(item)
    
    def set_selected(self, x, energy):
        self.selected = Replica(x, energy)
        self.show3d.setCoords(x, index=1)
        self.show3d.setCoords(None, index=2)

    
    def get_Emax(self):
        val = self.ui.lineEdit_Emax.text()
        try:
            self.Emax = float(val)
        except:
            print "Emax: can't turn", val, "into a float"
        if self.Emax is None:
            print "Emax is not set"
            self.Emax = self.nested_sampling.replicas[-1].energy
        return self.Emax
    
    def on_btn_db_sample_min_clicked(self, clicked=None):
        if clicked is None: return
        print clicked
        Emax = self.get_Emax()
        x, energy = self.nested_sampling.get_starting_configuration_minima(Emax)
        self.set_selected(x, energy)

    def mc_event(self, x=None, energy=None, **kwargs):
        self.mc_path.append(x)
        self.mc_path_energy.append(energy)

    def finish_cycle(self, x, energy):
        self.nested_sampling.pop_replica()
        self.nested_sampling.add_new_replica(x, energy)
        self.update_list()
        
        self.selected = None
        self.Emax = None

    def load_stepsize(self):
        stepsize = None
        val = self.ui.lineEdit_stepsize.text()
        try:
            stepsize = float(val)
            self.nested_sampling.takestep.stepsize = stepsize
        except:
            print "stepsize: can't convert", val, "into a float"
            return None
        
        return stepsize
    
    def load_mciter(self):
        mciter = None
        val = self.ui.lineEdit_mciter.text()
        try:
            mciter = int(val)
            self.nested_sampling.mciter = mciter
        except:
            print "stepsize: can't convert", val, "into an int"
            return None
        
        return mciter
    
    def on_btn_db_sample_coords_clicked(self, clicked=None):
        if clicked is None: return
        Emax = self.get_Emax()
        m = self.nested_sampling.bh_sampler.sample_minimum(Emax)        
        x = self.nested_sampling.bh_sampler.sample_coords_from_basin(m, Emax)
        pot = self.system.get_potential()
        e = pot.getEnergy(x)
        self.set_selected(x, e)
        self.show3d.setCoords(m.coords, index=2)
        print "energy of minimum        ", m.energy
        print "energy of starting config", e
        print "Emax                     ", Emax


    
    def on_btn_MC_chain_clicked(self, clicked=None):
        if clicked is None: return
        if self.selected is None:
            print "select a starting point first" 
            return
        Emax = self.get_Emax()
        if self.selected.energy >= self.Emax:
            print "energy of selected must be less than Emax"
            return
        x, energy = self.selected.x, self.selected.energy

        self.mc_path = [x]
        self.mc_path_energy = [energy]
        events = [self.mc_event]
        
        self.load_stepsize()
        self.load_mciter()
        
        
        mc = self.nested_sampling.do_monte_carlo_chain(self.selected.x, Emax, self.selected.energy, events=events)
        
        path = np.array(self.mc_path)
        labels = ["step %d: E=%s" % (i, str(e)) for i, e in enumerate(self.mc_path_energy)]
        self.show3d.setCoordsPath(path, labels=labels)
        self.finish_cycle(mc.x, mc.energy)
        

    def on_list_replicas_itemClicked(self, item):
        if item is None: return
        self.set_selected(item.r.x, item.r.energy)
        




def run_gui(system, nested_sampling):
    """
    The top level function that will launch the gui for a given system
    
    Parameters
    ----------
    system : System class
        A pele system, derrived from BaseSystemClass.  All information 
        about the system is in this class.
    db : str, optional
        connect to the database at this file location
        
    """
    app = QtGui.QApplication(sys.argv)
    
#    sys.excepthook = excepthook

    myapp = NSGUI(app, system, nested_sampling)
#    if db is not None:
#        myapp.connect_db(db)
        
#    refresh_timer = QtCore.QTimer()
#    refresh_timer.timeout.connect(refresh_pl)
#    refresh_timer.start(0.)
    myapp.show()
    sys.exit(app.exec_()) 
       
#def run_gui(systemtype):
#    app = QtGui.QApplication(sys.argv)
#    import pylab as pl
#    myapp = MyForm(systemtype)
#    refresh_timer = QtCore.QTimer()
#    refresh_timer.timeout.connect(refresh_pl)
#    refresh_timer.start(0.)
#    
#    myapp.show()
#    sys.exit(app.exec_())

def build_lj_nested_sampling(system, db):
    from pele.takestep import RandomDisplacement
    nreplicas = 20
    mciter = 10000
    nminima = 10000
    minima = db.minima()
    if len(minima) > nminima:
        minima = minima[:nminima]
    takestep = RandomDisplacement(stepsize=0.5)
    accept_tests = system.get_config_tests()
    ns = NestedSamplingBS(system, nreplicas, takestep, minima, mciter=mciter, accept_tests=accept_tests)
    return ns

if __name__ == "__main__":
    from lj_run import LJClusterNew
    natoms = 31
    system = LJClusterNew(31)

    dbname = "lj%d.db" % (natoms)
    db = system.create_database(dbname)

    ns = build_lj_nested_sampling(system, db)

    run_gui(system, ns)
