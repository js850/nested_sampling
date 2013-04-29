import multiprocessing as mp


class MCRunner(mp.Process):
    def __init__(self, conn, mc_runner):
        mp.Process.__init__(self)
        self.conn = conn
        self.mc_runner = mc_runner

    def do_MC(self, x0, mciter, stepsize, Emax, seed):
        return self.mc_runner(x0, mciter, stepsize, Emax, seed) 
     
    def run(self):
        while 1:
            message = self.conn.recv()
            #print "message", message
            if message == "kill":
                #print "terminating", self.name
                return
            elif message[0] == "do mc":
                #print "received message: calculating gradient"
#                print "recieved message", message[0]
                x0, mciter, stepsize, Emax, seed = message[1:]
                res = self.do_MC(x0, mciter, stepsize, Emax, seed)
#                print "sending results back"
                self.conn.send(res)
            else:
                print "error: unknown message: %s\n%s" % (self.name, message)
