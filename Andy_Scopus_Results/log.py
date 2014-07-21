import os
import datetime

class Logger(object):

    def __init__(self, filename, loc='logs', init=True):
        if not os.path.isdir(loc) and init:
            os.makedirs(loc)

        self.logfile = os.path.join(loc, filename)
        open(self.logfile, 'w')

    def log(self, message, status="INFO"):
        with open(self.logfile, 'a') as logfile:
            logfile.write("[" + status + "]" + str(datetime.datetime.now()) + "\n")
            logfile.write(message + "\n")

