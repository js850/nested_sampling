import sys
import Pyro4
import Pyro4.util
#from nested_sampling._mc_walker import MCWalkerParallelWrapper
import argparse
import uuid
try:
    import queue
except ImportError:
    import Queue as queue
from workitem import Workitem

class pyro_worker(object):
    """ the worker starts a demon that picks work times from the queue"""
    
    def __init__(self, dispatcher_URI, worker_name=None, host=None, port=0, serializer='pickle', server_type='multiplex'):
        
        if host==None:
            self.host = Pyro4.socketutil.getIpAddress(None, workaround127=True)
            print "host IP address was found to be {0}".format(host)
        else:
            self.host = host
        
        if worker_name == None:
            self.worker_name = "nested.sampling.worker@{0}".format(host)
        else:
            self.worker_name = worker_name
         
        self.port = port
        self.dispatcher_URI = dispatcher_URI
        self.serializer = serializer
        self.server_type = server_type
                 
    def _run(self,item):
        result = item.mc_runner(item.x, item.stepsize, item.Emax, item.energy, item.seed)
        item.x = result.x
        item.energy = result.energy
        item.processedBy = self.worker_name
    
    def _start_worker(self):
        
        sys.excepthook = Pyro4.util.excepthook
        Pyro4.config.SERIALIZER = self.serializer
        Pyro4.config.SERIALIZERS_ACCEPTED.add(self.serializer)
        Pyro4.config.SERVERTYPE=self.server_type
        
        self.dispatcher = Pyro4.core.Proxy(self.dispatcher_URI)
        print "This is worker: {0}".format(self.worker_name)
        print("getting work from dispatcher.")
        
        while True:
            try:
                item = self.dispatcher.getWork()
            except queue.Empty:
                pass
            else:
                self._run(item)
                self.dispatcher.putResult(item)

    

def main():   
    parser = argparse.ArgumentParser(description="must pass the URI of the dispatcher")
    parser.add_argument("dispatcher_URI", type=str, help="name for the worker")
    parser.add_argument("--worker-name", type=str, help="name for the worker",default=None)
    parser.add_argument("--host", type=str, help="address of the host (node on which the worker is started)",default=None)
    parser.add_argument("--port", type=int, help="port number on which the worker is started)",default=0)
    parser.add_argument("--server-type", type=str, help="multiplex or threaded",default="multiplex")
    args = parser.parse_args()
    
    if args.worker_name != None:
        worker_name = args.worker_name
    else:
        worker_name = "{0}".format(uuid.uuid4())
    
    dispatcher_URI = args.dispatcher_URI
    host = args.host
    port = args.port
    server_type = args.server_type
    
    worker = pyro_worker(dispatcher_URI, worker_name=worker_name, host=host, port=port, server_type=server_type)
    worker._start_worker()
       
if __name__ == "__main__":
    main()
        