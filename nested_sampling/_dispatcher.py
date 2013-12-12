#===============================================================================
# code adapted from distributed computing example of Pyro4
#===============================================================================

from __future__ import print_function
try:
    import queue
except ImportError:
    import Queue as queue
import Pyro4

class DispatcherQueue(object):
    def __init__(self):
        self.workqueue = queue.Queue()
        self.resultqueue = queue.Queue()
    def putWork(self, item):
        self.workqueue.put(item)
    def getWork(self, timeout=0):
        return self.workqueue.get(block=True, timeout=timeout)
    def putResult(self, item):
        self.resultqueue.put(item)
    def getResult(self, timeout=0):
        return self.resultqueue.get(block=True, timeout=timeout)
    def workQueueSize(self):
        return self.workqueue.qsize()
    def resultQueueSize(self):
        return self.resultqueue.qsize()