from nested_sampling import pyro_worker
import argparse
from sens.models._lj_tools import LJClusterSENS

def main():   
    parser = argparse.ArgumentParser(description="must pass the URI of the dispatcher")
    parser.add_argument("natoms", type=int, help="number of atoms")
    parser.add_argument("dispatcher_URI", type=str, help="name for the worker")
    parser.add_argument("-n","--mciter", type=int, default=1000, help="number of steps in the monte carlo walk")
    parser.add_argument("--radius", type=float, default=2.5, help="maintain atoms in a sphere of this radius")
    parser.add_argument("--worker-name", type=str, help="name for the worker",default=None)
    parser.add_argument("--host", type=str, help="address of the host (node on which the worker is started)",default=None)
    parser.add_argument("--port", type=int, help="port number on which the worker is started)",default=0)
    parser.add_argument("--server-type", type=str, help="multiplex or threaded",default="multiplex")
    args = parser.parse_args()
    
    #===========================================================================
    # MUST SET THE SYSTEM
    #===========================================================================
    system = LJClusterSENS(args.natoms, args.radius)
    
    dispatcher_URI = args.dispatcher_URI
    worker_name = args.worker_name
    host = args.host
    port = args.port
    server_type = args.server_type
    mciter = args.mciter
    
    worker = pyro_worker(dispatcher_URI, system, mciter, worker_name=worker_name, host=host, port=port, server_type=server_type)
    worker._start_worker()
       
if __name__ == "__main__":
    main()