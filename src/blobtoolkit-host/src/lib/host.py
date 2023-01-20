#!/usr/bin/env python3

# pylint: disable=no-member, too-many-branches, too-many-statements, too-many-locals, W0603, W0703

"""
Host a collection of BlobDirs.

Usage:
    blobtoolkit-host [--port INT]  [--api-port INT]
                     [-h|--help] [-v|--version]
                     [--hostname STRING] DIRECTORY

Arguments:
    DIRECTORY             Directory containing one or more BlobDirs.

Options:
    --port INT            HTTP port number. [Default: 8080]
    --api-port INT        API port number. [Default: 8000]
    --hostname STRING     Hostname used to connect to API. [Default: localhost]
    -h, --help            Show this
    -v, --version         Show version number

"""


import contextlib
import os
import platform
import shlex
import signal
import socket
import sys
import time
from pathlib import Path
from shutil import which
from subprocess import PIPE
from subprocess import Popen

import psutil
from docopt import docopt

from .version import __version__

PIDS = []


def iter_user_procs(process):
    """Iterate through processes owned by the current user."""
    username = psutil.Process.username(
        process
    )  # get username in format used by processes
    for proc in psutil.process_iter(["username"]):
        if proc.info["username"] == username:
            yield proc


def close_ports(process, args):
    """Close processes running on chosen ports."""
    for proc in iter_user_procs(process):
        try:
            for conns in proc.connections(kind="inet"):
                if conns.laddr.port == args["--api-port"]:
                    proc.send_signal(signal.SIGTERM)
                elif conns.laddr.port == args["--port"]:
                    proc.send_signal(signal.SIGTERM)
        except psutil.NoSuchProcess:
            continue
        except psutil.ZombieProcess:
            continue
        except psutil.AccessDenied:
            continue


def kill_child_processes(parent_pid, args=None, sig=signal.SIGTERM):
    """Kill all child processes."""
    try:
        parent = psutil.Process(parent_pid)
    except psutil.NoSuchProcess:
        return
    # if args:
    #     close_ports(parent, args)
    children = parent.children(recursive=True)
    for process in children:
        try:
            process.send_signal(sig)
        except psutil.NoSuchProcess:
            continue
        except psutil.ZombieProcess:
            continue
    parent.send_signal(sig)


def test_port(port, service):
    """Exit if port is already in use."""
    port = int(port)
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as skt:
        try:
            skt.bind(("", port))
        except OSError:
            if service == "test":
                return False
            print("ERROR: Port %d already in use, unable to host %s." % (port, service))
            print(
                "       Use: `lsof -nP -iTCP:%d | grep LISTEN` to find the associated process."
                % port
            )
            print(
                f"       It may take ~30s for this port to become available when restarting {service}."
            )

            sys.exit(1)
    return True


def find_binary(tool):
    """Find a binary executable for the blobtoolkit viewer or API."""
    script_dir = os.path.dirname(os.path.realpath(__file__))
    system = platform.system()
    # arch, _ = platform.architecture()
    default_binaries = {"api": "blobtoolkit-api", "viewer": "blobtoolkit-viewer"}
    if system == "Darwin":
        binaries = {
            "api": "blobtoolkit-api-macos",
            "viewer": "blobtoolkit-viewer-macos",
        }
    elif system == "Linux":
        binaries = {
            "api": "blobtoolkit-api-linux",
            "viewer": "blobtoolkit-viewer-linux",
        }
    elif system == "Windows":
        binaries = {
            "api": "blobtoolkit-api-win.exe",
            "viewer": "blobtoolkit-viewer-win.exe",
        }
        default_binaries = {
            "api": "blobtoolkit-api.exe",
            "viewer": "blobtoolkit-viewer.exe",
        }
    if which(default_binaries[tool]) is not None:
        return default_binaries[tool]
    if which(binaries[tool]) is not None:
        return default_binaries[tool]
    executable = binaries[tool]
    executable_path = os.path.join(
        os.path.dirname(script_dir), "data", "bin", executable
    )
    if os.path.isfile(executable_path):
        return executable_path
    print(
        f"ERROR: {default_binaries[tool]} executable was not found. Please add {default_binaries[tool]} to your PATH."
    )

    sys.exit(1)


def start_api(port, api_port, hostname, directory):
    """Start BlobToolKit API."""
    cmd = find_binary("api")
    # cmd = "blobtoolkit-api"
    origins = "http://localhost:%d http://localhost null" % int(port)
    if hostname != "localhost":
        origins += " http://%s:%d http://%s" % (hostname, int(port), hostname)
    if directory == "_":
        env = dict(
            os.environ,
            BTK_API_PORT=api_port,
            BTK_ORIGINS=origins,
        )
    else:
        env = dict(
            os.environ,
            BTK_API_PORT=api_port,
            BTK_FILE_PATH=directory,
            BTK_ORIGINS=origins,
        )
    return Popen(
        shlex.split(cmd),
        stdout=PIPE,
        stderr=PIPE,
        encoding="ascii",
        env=env,
    )


def start_viewer(port, api_port, hostname):
    """Start BlobToolKit viewer."""
    cmd = find_binary("viewer")
    # cmd = "blobtoolkit-viewer"
    api_url = "http://%s:%d/api/v1" % (hostname, int(api_port))
    return Popen(
        shlex.split(cmd),
        stdout=PIPE,
        stderr=PIPE,
        encoding="ascii",
        env=dict(
            os.environ,
            BTK_HOST=hostname,
            BTK_CLIENT_PORT=port,
            BTK_API_PORT=api_port,
            BTK_API_URL=api_url,
        ),
    )


def main(args):
    """Entrypoint for blobtools host."""
    global PIDS
    directory = args["DIRECTORY"]
    if directory != "_":
        path = Path(directory)
        if not path.exists():
            print("ERROR: Directory '%s' does not exist" % directory)
            sys.exit(1)
        if (path / "meta.json").exists():
            print("WARNING: Directory '%s' appears to be a BlobDir." % directory)
            print("         Hosting the parent directory instead.")
            path = path.resolve().parent
        directory = path.absolute()
    test_port(args["--api-port"], "BlobtoolKit API")
    test_port(args["--port"], "BlobtoolKit viewer")
    api = start_api(
        args["--port"],
        args["--api-port"],
        args["--hostname"],
        directory,
    )
    PIDS.append(api.pid)
    print(
        "Starting BlobToolKit API on port %d (pid: %d)"
        % (int(args["--api-port"]), api.pid)
    )
    time.sleep(2)
    viewer = start_viewer(args["--port"], args["--api-port"], args["--hostname"])
    PIDS.append(viewer.pid)
    print(
        "Starting BlobToolKit viewer on port %d (pid: %d)"
        % (int(args["--port"]), viewer.pid)
    )
    time.sleep(2)
    ready = False
    url = "http://%s:%d" % (args["--hostname"], int(args["--port"]))
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
                with contextlib.suppress(ProcessLookupError):
                    os.kill(viewer.pid, signal.SIGTERM)
            break
        if viewer.poll() is not None:
            for line in viewer.stdout.readlines():
                print(line.strip())
            for line in viewer.stderr.readlines():
                print(line.strip())
            with contextlib.suppress(ProcessLookupError):
                os.kill(api.pid, signal.SIGTERM)
            break
        if not ready:
            print(f"Visit {url} to use the interactive BlobToolKit Viewer.")
            ready = True
        time.sleep(1)


def cli(rename=None):
    """Entry point."""
    if sys.argv[0].endswith("-host"):
        args = docopt(__doc__, version=__version__)
    else:
        docs = __doc__
        if rename is not None:
            docs = docs.replace("blobtoolkit-host", rename)
        if len(sys.argv) == sys.argv.index(__name__.split(".")[-1]) + 1:
            args = docopt(docs, argv=[])
        else:
            args = docopt(docs, version=__version__)
    try:
        main(args)
    except KeyboardInterrupt:
        pass
    finally:
        for pid in PIDS:
            kill_child_processes(pid, args=args, sig=signal.SIGKILL)
            time.sleep(2)


if __name__ == "__main__":
    cli()
