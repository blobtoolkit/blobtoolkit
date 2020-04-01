#!/usr/bin/env python3

# pylint: disable=no-member, too-many-branches, too-many-statements, too-many-locals, W0603, W0703

"""
Host a collection of BlobDirs.

Usage:
    blobtools host [--port INT]  [--api-port INT]
                   [--hostname STRING] [--viewer DIRECTORY] DIRECTORY

Arguments:
    DIRECTORY             Directory containing one or more BlobDirs.

Options:
    --port INT            HTTP port number. [Default: 8080]
    --api-port INT        API port number. [Default: 8000]
    --hostname STRING     Hostname used to connect to API. [Default: localhost]
    --viewer DIRECTORY    Path to BlobToolKit viewer. [Default: ../viewer]
"""

import os
import shlex
import signal
import socket
import time
from pathlib import Path
from subprocess import PIPE, Popen

from docopt import docopt

PIDS = []


def test_port(port, service):
    """Exit if port is already in use."""
    port = int(port)
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as skt:
        try:
            skt.bind(('', port))
        except OSError:
            if service == 'test':
                return False
            print("ERROR: Port %d already in use, unable to host %s." % (port, service))
            print("       Use: `lsof -nP -iTCP:%d | grep LISTEN` to find the associated process." % port)
            print("       It may take ~30s for this port to become available when restarting %s." % service)
            exit(1)
    return True


def start_api(cwd, port, api_port, hostname, directory):
    """Start BlobToolKit API."""
    cmd = 'npm run api'
    origins = "http://localhost:%d http://localhost null" % int(port)
    if hostname != 'localhost':
        origins += " http://%s:%d http://%s" % (hostname, int(port), hostname)
    process = Popen(shlex.split(cmd),
                    stdout=PIPE,
                    stderr=PIPE,
                    encoding='ascii',
                    cwd=cwd,
                    env=dict(os.environ,
                             BTK_API_PORT=api_port,
                             BTK_FILE_PATH=directory,
                             BTK_ORIGINS=origins)
                    )
    return process


def start_viewer(cwd, port, api_port, hostname):
    """Start BlobToolKit viewer."""
    cmd = 'npm run client'
    api_url = "http://%s:%d/api/v1" % (hostname, int(api_port))
    process = Popen(shlex.split(cmd),
                    stdout=PIPE,
                    stderr=PIPE,
                    encoding='ascii',
                    cwd=cwd,
                    env=dict(os.environ,
                             BTK_HOST=hostname,
                             BTK_CLIENT_PORT=port,
                             BTK_API_PORT=api_port,
                             BTK_API_URL=api_url)
                    )
    return process


def main():
    """Entrypoint for blobtools host."""
    global PIDS
    args = docopt(__doc__)
    viewer_dir = Path(__file__).resolve().parent.parent / args['--viewer']
    viewer_dir.resolve().absolute()
    if not viewer_dir.exists():
        print("ERROR: Viewer could not be found at '%s'" % args['--viewer'])
        exit(1)
    elif not (viewer_dir / 'package.json').exists():
        print("ERROR: Directory '%s' does not appear to contain the viewer code" % args['--viewer'])
        exit(1)
    path = Path(args['DIRECTORY'])
    if not path.exists():
        print("ERROR: Directory '%s' does not exist" % args['DIRECTORY'])
        exit(1)
    if (path / 'meta.json').exists():
        print("WARNING: Directory '%s' appears to be a BlobDir." % args['DIRECTORY'])
        print("         Hosting the parent directory instead.")
        path = path.resolve().parent
    test_port(args['--api-port'], 'BlobtoolKit API')
    test_port(args['--port'], 'BlobtoolKit viewer')
    api = start_api(viewer_dir,
                    args['--port'],
                    args['--api-port'],
                    args['--hostname'],
                    path.absolute())
    PIDS.append(api.pid)
    print("Starting BlobToolKit API on port %d (pid: %d)" % (int(args['--api-port']), api.pid))
    time.sleep(2)
    viewer = start_viewer(viewer_dir,
                          args['--port'],
                          args['--api-port'],
                          args['--hostname'])
    PIDS.append(viewer.pid)
    print("Starting BlobToolKit viewer on port %d (pid: %d)" % (int(args['--port']), viewer.pid))
    time.sleep(2)
    ready = False
    url = "http://%s:%d" % (args['--hostname'], int(args['--port']))
    while True:
        time.sleep(1)
        if api.poll() is not None:
            for line in api.stdout.readlines():
                print(line.strip())
            for line in api.stderr.readlines():
                print(line.strip())
            if viewer.poll() is not None:
                for line in viewer.stdout.readlines():
                    print(line.strip())
                for line in viewer.stderr.readlines():
                    print(line.strip())
                try:
                    os.kill(viewer.pid, signal.SIGTERM)
                except ProcessLookupError:
                    pass
            break
        elif viewer.poll() is not None:
            for line in viewer.stdout.readlines():
                print(line.strip())
            for line in viewer.stderr.readlines():
                print(line.strip())
            try:
                os.kill(api.pid, signal.SIGTERM)
            except ProcessLookupError:
                pass
            break
        elif not ready:
            print("Visit %s to use the interactive BlobToolKit Viewer." % url)
            ready = True
        time.sleep(1)


if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        pass
    finally:
        for pid in PIDS:
            print(pid)
            try:
                os.kill(pid, signal.SIGTERM)
            except ProcessLookupError:
                pass
