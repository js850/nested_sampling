from nested_sampling import DispatcherQueue
import argparse
import sys
import Pyro4
import socket

def main():
    parser = argparse.ArgumentParser(description="dispatcher queue")
    parser.add_argument("--host", type=str, help="address of the host (node on which the worker is started)",default=None)
    parser.add_argument("--server-type", type=str, help="multiplex or threaded",default="multiplex")
    args = parser.parse_args()
    
    if args.host == None:
        hostname=socket.gethostname()
        host = Pyro4.socketutil.getIpAddress(hostname, workaround127=True)
    else:
        host = args.host
    
    sys.excepthook = Pyro4.util.excepthook
    Pyro4.config.SERIALIZER = 'pickle'
    Pyro4.config.SERIALIZERS_ACCEPTED.add('pickle')
    Pyro4.config.SERVERTYPE=args.server_type
    
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