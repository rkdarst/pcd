#!/usr/bin/env python

# Richard Darst, November 2013

import sys

import pcd.util
import pcd.ioutil

def run():
    stack = [ ]
    for cmd in sys.argv[1:]:
        cmd = pcd.util.leval(cmd)

        if cmd == 'read':
            pass
        else:
            stack.append(cmd)


if __name__ == "__main__":
    run()
