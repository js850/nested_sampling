import sys
import Pyro4
import Pyro4.util
import socket

#from nested_sampling._mc_walker import MCWalkerParallelWrapper
try:
    import queue
except ImportError:
    import Queue as queue
from nested_sampling import Forwarditem
import uuid

class pyro_worker(object):
    """ the worker starts a demon that picks work times from the queue"""
    
    def __init__(self, dispatcher_URI, mc_runner, worker_name=None, host=None, port=0, serializer='pickle', server_type='multiplex'):
        
        if host==None:
            hostname = socket.gethostname()
            self.host = Pyro4.socketutil.getIpAddress(hostname, workaround127=True)
            print "host IP address was found to be {0}".format(self.host)
        else:
            self.host = host
        
        if worker_name == None:
            self.worker_name = "n.s.worker.{0}@{1}".format(uuid.uuid4(),self.host)
        else:
            self.worker_name = worker_name
         
        self.port = port
        self.dispatcher_URI = dispatcher_URI
        self.serializer = serializer
        self.server_type = server_type
        self.mc_runner = mc_runner
        
    def _run(self,item):
        mc = self.mc_runner(item.x, item.stepsize, item.Emax, item.energy, item.seed)
        mc.processedBy = self.worker_name
        mc.Id = item.Id
        return mc
    
    def _start_worker(self):
        
        sys.excepthook = Pyro4.util.excepthook
        Pyro4.config.SERIALIZER = self.serializer
        Pyro4.config.SERIALIZERS_ACCEPTED.add(self.serializer)
        Pyro4.config.SERVERTYPE=self.server_type
        
        self.dispatcher = Pyro4.Proxy(self.dispatcher_URI)
        print "This is worker {0}".format(self.worker_name)
        print("getting work from dispatcher.")
        
        while True:
            try:
                item = self.dispatcher.getWork()
            except queue.Empty:
                pass
            else:
                result = self._run(item)
                self.dispatcher.putResult(result)
        