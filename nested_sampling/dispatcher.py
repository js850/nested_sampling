#===============================================================================
# code adapted from distributed computing and event loop example of Pyro4
#===============================================================================
from __future__ import print_function
import sys
import argparse
from Pyro4 import core
import logging
try:
    import queue
except ImportError:
    import Queue as queue
import Pyro4
from Pyro4.errors import PyroError, NamingError

log=logging.getLogger("Pyro4.naming")

class DispatcherQueue(object):
    def __init__(self):
        self.workqueue = queue.Queue()
        self.resultqueue = queue.Queue()
    def putWork(self, item):
        self.workqueue.put(item)
    def getWork(self, timeout=5):
        return self.workqueue.get(block=True, timeout=timeout)
    def putResult(self, item):
        self.resultqueue.put(item)
    def getResult(self, timeout=5):
        return self.resultqueue.get(block=True, timeout=timeout)
    def workQueueSize(self):
        return self.workqueue.qsize()
    def resultQueueSize(self):
        return self.resultqueue.qsize()

def main():
    parser = argparse.ArgumentParser(description="dispatcher queue")
    parser.add_argument("--host", type=str, help="address of the host (node on which the worker is started)",default=None)
    parser.add_argument("--server-type", type=str, help="multiplex or threaded",default="multiplex")
    args = parser.parse_args()
    
    sys.excepthook = Pyro4.util.excepthook
    Pyro4.config.SERIALIZER = 'pickle'
    Pyro4.config.SERIALIZERS_ACCEPTED.add('pickle')
    Pyro4.config.SERVERTYPE=args.server_type
    
    if args.host == None:
        host = Pyro4.socketutil.getIpAddress(None, workaround127=True)
    else:
        host = args.host
    
    daemon = Pyro4.Daemon(host=host)
    dispatcher_uri = daemon.register(DispatcherQueue())
    print(dispatcher_uri)
    out_dispatcher_uri = open('dispatcher_uri.dat','w+')
    out_dispatcher_uri.write(str(dispatcher_uri)) 
    daemon.requestLoop()

#    hostname=socket.gethostname()
#    ns_command = "python -m Pyro4.naming " + "-n " + str(my_ip) + " &" 
#    print(ns_command)
#    os.system(ns_command)

if __name__ == "__main__":
    main()