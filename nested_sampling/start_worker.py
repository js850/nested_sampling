import sys
import Pyro4
import Pyro4.util
from nested_sampling._mc_walker import MCWalkerParallelWrapper
import argparse
import uuid

class pyro_worker(object):
    """ the worker starts a demon and registers its uri with the name server (passed to the )"""
    
    def __init__(self, worker_name, job_name, nsIP, host=None, port=0, serializer='pickle'):
        
        self.worker_name = "nested.sampling.{0}.{1}".format(job_name, worker_name)
        self.host = host
        self.port = port
        self.nsIP = nsIP
        self.serializer = serializer
        
        sys.excepthook = Pyro4.util.excepthook
        Pyro4.config.BROADCAST_ADDRS = self.nsIP
        Pyro4.config.SERIALIZER = self.serializer
        Pyro4.config.SERIALIZERS_ACCEPTED.add(self.serializer)
        self.ns = Pyro4.locateNS()
    
    def _start_core(self):
        worker = MCWalkerParallelWrapper()
        self.daemon = Pyro4.Daemon(host=self.host,port=self.port)
        self.worker_uri = self.daemon.register(worker)
        self.ns.register(self.worker_name, self.worker_uri)
        print "{0} is listening".format(self.worker_name)
        self.daemon.requestLoop()
    
    def name_and_uri(self):
        return self.worker_name, self.worker_uri

def main():   
    parser = argparse.ArgumentParser(description="must pass a name for the worker to be registered with the name server", 
                                                epilog="and the IP address where the name server is kept")
    parser.add_argument("job_name", type=str, help="name of the job")
    parser.add_argument("name_server_IP", type=str, help="IP address of the machine hosting the Name Server")
    parser.add_argument("--worker-name", type=str, help="name for the worker",default=None)
    parser.add_argument("--host", type=str, help="address of the host (node on which the worker is started)",default=None)
    parser.add_argument("--port", type=int, help="port number on which the worker is started)",default=0)
    args = parser.parse_args()
    
    if args.worker_name != None:
        worker_name = args.worker_name
    else:
        worker_name = "{0}".format(uuid.uuid4())
    job_name = args.job_name
    nsIP = args.name_server_IP
    host = args.host
    port = args.port
    
    worker = pyro_worker(worker_name, job_name, nsIP, host=host, port=port)
    worker._start_core()
       
if __name__ == "__main__":
    main()
        