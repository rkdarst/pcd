# Richard Darst, October 2013

import io
import os
import select
#import subprocess
import time
import threading


class FDSplitter(object):  #io.IOBase):
    """Class for splitting a file descriptor (for subprocess).

    The subprocess module (and subprocesses in general, for any unix
    process) use file descriptors (integers), not file objects, as the
    standard output.  For subprocess, if you want to capture that
    output and direct it multiple places, you can not use a class with
    a .write() method.  You must use one with a .fileno() method.
    This class creates a pipe (os.pipe(), a unix concept different
    from shell pipes) and uses that pipe as the target for subprocess
    stdout.  There is a daemon thread that reads data from that pipe,
    and writes it to all of the children file descriptors.

    Base classes:

    http://stackoverflow.com/questions/8684091/how-can-i-implement-a-posix-file-descriptor-in-python-3
    http://stackoverflow.com/questions/4713932/decorate-delegate-a-file-object-to-add-functionality/

    """
    bufsize = 1024
    def __init__(self, fds, flush_every=4096):
        """Class init.

        fds: list of file objects or (integer) file descriptors.

        flush_every: buffers are flushed every time this much data has accumulated."""
        self._run = None
        self._done = False
        self.fds = fds
        self.pipe = os.pipe()
        self.thread = threading.Thread(target=self._flusher)
        self.thread.setDaemon(True)
        self.thread.start()
        self.flush_every = flush_every
    def _flusher(self):
        self._run = True
        buf = b''
        unflushed_data = 0
        #while self._run:
        while True:
            #print "="*20+"flusher loop"
            #for fh in select.select([self.pipe[0]], [], [], 0)[0]:
            #    new = os.read(fh, self.bufsize)
            #    buf += new
            #    while b'\n' in buf:
            #        data, buf = buf.split(b'\n', 1)
            #        self.write(data+'\n')
            #time.sleep(.1)
            #select.select([self.pipe[0]], [], [])
            new = os.read(self.pipe[0], self.bufsize)
            #print "new data: %r"%new
            if len(new) == 0:
                break
            unflushed_data += len(new)
            buf += new
            while b'\n' in buf:
                data, buf = buf.split(b'\n', 1)
                self.write(data+'\n')
            if unflushed_data > self.flush_every:
                self.flush()
                unflushed_data = 0
            #time.sleep(1)
        #print "="*20+"Breaking from flusher loop"
        self.write(buf)
        self.flush()
        self._done = True

    def write(self, data):
        #print "data:", data.strip()
        for fd in self.fds:
            #print "writing to fd:", fd
            if isinstance(fd, int):
                os.write(fd, data)
            else:
                fd.write(data)
    def flush(self):
        for fd in self.fds:
            if isinstance(fd, int):
                os.fsync(fd)
            else:
                fd.flush()
    def fileno(self):
        return self.pipe[1]
    def wait(self):
        """Waits for completion (all data to be read and the pipe to
        close) and return."""
        #print select.select([], [self.pipe[1]], [])
        while not self._done:
            time.sleep(.1)

    def close(self, force=False):
        """Close internal pipes.

        You MUST call close() on this object eventually, or else the process can hang."""
        if self._run or force:
            #print '='*20, "Closing in", os.getpid()
            self._run = False
            #while self._run is not None:
            #    time.sleep(1)
            #os.close(self.pipe[0])
            os.close(self.pipe[1])
        #self._run = False
    def __del__(self):
        self.close()
        os.close(self.pipe[0])
