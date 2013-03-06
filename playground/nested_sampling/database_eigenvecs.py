"""
attach eigenvalues and eigenvectors to the database
"""
import numpy as np

from sqlalchemy.orm import relationship, backref, deferred
from sqlalchemy import Column, Integer, Float, PickleType, ForeignKey

from pygmin.storage.database import Base
from pygmin.storage import Database

class EigenPair(Base):
    """an eigenvalue and associated eigenvectors
    
    actually frequencies are stored rather than eigenvalues.  Frequencies are
    eigenvalues with the appropriate weighting for (from the mass or from rigid bodies)
    """
    __tablename__ = "tbl_eigenpair"
    _id = Column(Integer, primary_key=True)

    vector = deferred(Column(PickleType))
    freq = Column(Float)

    _minimum_id = Column(Integer, ForeignKey('tbl_minima._id'))
    minimum = relationship("Minimum",
                            primaryjoin="Minimum._id==EigenPair._minimum_id",
                            backref='eigenpair')
    def __init__(self, m, freq, vector):
        self.minimum = m
        self.freq = freq
        self.vector = np.copy(vector)

if __name__ == "__main__":
    from pygmin.systems import LJCluster
    from pygmin.utils.hessian import get_sorted_eig
    system = LJCluster(20)
    db = system.create_database("test.sqlite")
    coords, energy = system.get_random_minimized_configuration()[:2]
    print energy
    m = db.addMinimum(energy, coords)
    
    pot = system.get_potential()
    e, g, hess = pot.getEnergyGradientHessian(coords)
    eval, evec = get_sorted_eig(hess)
    if False:
        epair = EigenPair(m, eval[0], evec[:,0])
        print m.eigenpair 
        print m.eigenpair[0].freq, eval[0]
        epair = EigenPair(m, eval[1], evec[:,1])
        print m.eigenpair 
        print m.eigenpair[0].freq, eval[0]
        print m.eigenpair[1].freq, eval[1]
    else:
        for i in range(len(eval)):
            epair = EigenPair(m, eval[i], evec[:,i])
#            print m.eigenpair
        db.session.commit()
        print m.eigenpair[-1].freq
        print len(m.eigenpair)
        db.session.delete(m.eigenpair[0])
        db.session.commit()
        print len(m.eigenpair)
        
