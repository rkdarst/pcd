#!/usr/bin/env python

# Richard Darst, November 2013

import sys

import pcd.cmty
import pcd.util
import pcd.ioutil

def run():
    stack = [ ]
    for cmd in sys.argv[1:]:
        cmd = pcd.util.leval(cmd)

        if cmd == 'open_cmty':
            cmtys = pcd.cmty.CommunityFile(stack[-1])
        elif cmd == 'open_graph':
            g = pcd.ioutil.read_any(stack[-1])
        elif cmd == 'interact':
            from fitz import interact ; interact.interact()
        else:
            stack.append(cmd)


if __name__ == "__main__":
    run()
