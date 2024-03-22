#!/usr/bin/env python3

# pylint: disable=no-member, too-many-branches, too-many-locals, too-many-statements, broad-except

"""
Generate plots using BlobToolKit Viewer.

Usage:
  blobtools view [--format STRING...] [--host STRING] [--interactive]
  [--out PATH] [--param STRING...] [--ports RANGE] [--prefix STRING]
  [--preview STRING...] [--driver STRING] [--driver-log PATH]
  [--local] [--remote] [--plot] [--timeout INT] [--view STRING...] DIRECTORY

Options:
      --format STRING         Image format (svg|png). [Default: png]
      --host STRING           Hostname. [Default: http://localhost]
      --interactive           Start interactive session (opens dataset in Firefox/Chromium). [Default: False]
      --out PATH              Directory for outfiles. [Default: .]
      --param key=value       Query string parameter.
      --ports RANGE           Port range for viewer and API. [Default: 8000-8099]
      --prefix STRING         URL prefix. [Default: view]
      --preview STRING        Field name.
      --driver STRING         Webdriver to use (chromium or firefox). [Default: firefox]
      --driver-log PATH       Path to driver logfile for debugging. [Default: /dev/null]
      --local                 Start viewer for local session. [Default: False]
      --remote                Start viewer for remote session. [Default: False]
      --plot                  Use blobtk plot to generate plots. [Default: False]
      --timeout INT           Time to wait for page load in seconds. Default (0) is no timeout. [Default: 0]
      --view STRING           Plot type (blob|cumulative|snail). [Default: blob]

Examples:
    # 1. Create a snail plot:
    blobtools view --plot --view snail /path/to/BlobDir

    # 2. Start an interactive session:
    blobtools view --interactive /path/to/BlobDir

    # 3. Start an remote session:
    blobtools view --remote /path/to/BlobDir

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
from site import getsitepackages
from subprocess import PIPE
from subprocess import Popen
from traceback import format_exc

from blobtk import plot
from docopt import docopt
from pyvirtualdisplay import Display
from pyvirtualdisplay.abstractdisplay import XStartError
from selenium import webdriver
from selenium.common.exceptions import WebDriverException
from selenium.webdriver.chrome.service import Service
from selenium.webdriver.common.by import By
from selenium.webdriver.firefox.options import Options
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver.support.ui import WebDriverWait
from tolkein import tolog

from .version import __version__

LOGGER = tolog.logger(__name__)


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


def file_ready(file_path, timeout, callback):
    """Check if file is ready."""
    start_time = time.time()
    while not os.path.exists(file_path):
        elapsed_time = time.time() - start_time
        if timeout and elapsed_time > timeout:
            callback("Timeout waiting for file")
            return False
        # flush nfs cache by chowning parent to current owner
        parent = os.path.dirname(os.path.abspath(file_path))
        os.chown(parent, os.stat(parent).st_uid, os.stat(parent).st_gid)
        time.sleep(1)
    if os.path.isfile(file_path):
        print(f" - found {file_path}")
        return True
    callback(f"{file_path} is not a file")
    return False


def test_loc(args):
    """See if dataset needs to be hosted and, if so, find an empty port."""
    info = args["--host"].split(":")
    dataset = Path(args["DIRECTORY"]).name
    level = "dataset"
    if len(info) >= 2 and info[1] != "//localhost":
        loc = f'{args["--host"]}/{args["--prefix"]}/{dataset}/dataset/{dataset}'
        return loc, None, None, None, level
    if len(info) == 1 and info[0] != "localhost":
        # need to add test for http vs https
        loc = f'http://{args["--host"]}/{args["--prefix"]}/{dataset}/dataset/{dataset}'
        return loc, None, None, None, level
    if len(info) == 3:
        port = info[2]
        for i in range(10):
            if available := test_port(port, "test"):
                if i == 9:
                    print(f"ERROR: No service running on port {port}")
                    print(f'       Unable to connect to {args["--host"]}')
                    sys.exit(1)
                time.sleep(1)
            else:
                loc = f'{args["--host"]}/{args["--prefix"]}/{dataset}/dataset/{dataset}'
                return loc, None, None, None, level
    if args["DIRECTORY"] == "_":
        parent = "_"
        level = "blobdir"
    else:
        if not Path(args["DIRECTORY"]).exists():
            print(
                "ERROR: DIRECTORY '%s' must be a valid path to begin hosting."
                % args["DIRECTORY"]
            )
            sys.exit(1)
        dataset = Path(args["DIRECTORY"]).name
        if (
            Path(f'{args["DIRECTORY"]}/meta.json').is_file()
            or Path(f'{args["DIRECTORY"]}/meta.json.gz').is_file()
        ):
            parent = Path(args["DIRECTORY"]).resolve().absolute().parent
        else:
            level = "blobdir"
            parent = Path(args["DIRECTORY"]).resolve().absolute()
    port_range = args["--ports"].split("-")
    api_port = False
    port = False
    for i in range(int(port_range[0]), int(port_range[1])):
        if test_port(i, "test"):
            if not api_port:
                api_port = i
                continue
            if not port:
                port = i
            break
    # directory = Path(__file__).resolve().parent.parent
    cmd = "blobtools host --port %d --api-port %d %s" % (port, api_port, parent)
    process = Popen(shlex.split(cmd), stdout=PIPE, stderr=PIPE, encoding="ascii")
    loc = "%s:%d/%s" % (args["--host"], port, args["--prefix"])
    loc += f"/{dataset}/dataset/{dataset}" if level == "dataset" else "/all"
    print("Initializing viewer")
    time.sleep(2)
    for i in range(4):
        poll = process.poll()
        if poll is None:
            time.sleep(1)
            if i > 1 and test_port(port, "test"):
                break
        else:
            print(process.stdout.read(), file=sys.stderr)
            print(process.stderr.read(), file=sys.stderr)
            print("ERROR: Viewer quit unexpectedly", file=sys.stderr)
            print(f"Unable to run: {cmd}", file=sys.stderr)
            sys.exit(1)
    return loc, process, port, api_port, level


def check_permission(file_path):
    """Check if current user has delete permission on a file."""
    file_dirname = os.path.dirname(file_path)
    if os.path.isfile(file_path):
        if os.access(file_path, os.W_OK):
            if os.access(file_dirname, os.W_OK | os.X_OK):
                return True
    return False


def start_display():
    """Start Xvfb display."""
    try:
        display = Display(visible=0, size=(800, 600))
        display.start()
    except XStartError:
        LOGGER.error("Unable to start X server")
        if platform.system() == "Darwin":
            LOGGER.error(
                "XQuartz (https://www.xquartz.org) must be installed and running to use this command"
            )
        sys.exit(1)
    except FileNotFoundError as err:
        if "Xvfb" in str(err):
            if "/usr/X11/bin" not in os.environ["PATH"]:
                try:
                    os.environ["PATH"] += os.pathsep + os.path.join(
                        "/usr", "x11", "bin"
                    )
                    LOGGER.info("Xvfb not found in PATH, trying /usr/X11/bin")
                    display = Display(visible=0, size=(800, 600))
                    display.start()
                except Exception:
                    LOGGER.error("Xvfb not found in PATH or /usr/X11/bin")
                    sys.exit(1)
            else:
                LOGGER.error("Xvfb not found in PATH")
                sys.exit(1)
        else:
            raise err
    return display


def firefox_driver(args):
    """Start firefox."""
    if not (which("firefox") or which("geckodriver")):
        try:
            LOGGER.info("Firefox and geckodriver must be available for blobtools view")
            LOGGER.info("Attempting to install the appropriate geckodriver version")
            import geckodriver_autoinstaller

            geckodriver_autoinstaller.install()
            LOGGER.info("Successfully installed geckodriver")

        except Exception:
            LOGGER.error(
                "Unable to install automatically. Try `conda install -c conda-forge firefox geckodriver` or `blobtools view --driver chromium` to use chrome browser."
            )
            sys.exit(1)

    outdir = os.path.abspath(args["--out"])
    os.makedirs(Path(outdir), exist_ok=True)

    options = Options()
    options.set_preference("browser.download.folderList", 2)
    options.set_preference("browser.download.manager.showWhenStarting", False)
    options.set_preference("browser.download.dir", outdir)
    options.set_preference("browser.download.lastDir", args["--out"])
    options.set_preference(
        "browser.helperApps.neverAsk.saveToDisk",
        "image/png, image/svg+xml, text/csv, text/plain, application/json",
    )
    options.headless = not args["--interactive"]

    display = start_display()
    service = Service(
        service_log_path=args["--driver-log"],
    )
    try:
        driver = webdriver.Firefox(
            options=options,
            service=service,
        )
    except WebDriverException:
        LOGGER.error(
            "Unable to start Firefox using geckodriver, possible version mismatch."
        )
        LOGGER.error(
            "Try `conda install -c conda-forge firefox geckodriver` or use `blobtools view --driver chromium` to use chrome browser."
        )
        if check_permission(which("geckodriver")):
            LOGGER.error(
                "Or try deleting the existing geckodriver executable from %s then rerun this command to install a compatible version.",
                which("geckodriver"),
            )
        else:
            LOGGER.error(
                "Or ask your administrator to ensure that the versions of Firefox and geckodriver are compatible"
            )
        display.stop()
        sys.exit(1)
    time.sleep(2)
    return driver, display


def chromium_driver(args):
    """Start chromium browser."""
    os.environ["WDM_LOG_LEVEL"] = "0"
    try:
        LOGGER.info("Chrome and chromedriver must be available for blobtools view")
        LOGGER.info("Attempting to install the appropriate chromedriver version")
        import chromedriver_autoinstaller

        chromedriver_autoinstaller.install()
        # ChromeDriverManager()

        # service_object = Service(ChromeDriverManager().install())
        LOGGER.info("Successfully installed chromedriver")
    except Exception as err:
        # if "LATEST_RELEASE_" in str(err):
        #     try:
        #         _, version = str(err).split("LATEST_RELEASE_")
        #         from selenium.webdriver.chrome.service import Service as ChromeService

        #         service_object = ChromeService(
        #             ChromeDriverManager(
        #                 latest_release_url="https://googlechromelabs.github.io/chrome-for-testing/last-known-good-versions-with-downloads.json",
        #                 version=version,
        #             ).install()
        #         )
        #     except Exception as err:
        #         print(str(err))
        #         LOGGER.error(
        #             "Unable to install automatically. Check Chrome is available or try `blobtools view --driver firefox` to use firefox browser."
        #         )
        #         sys.exit(1)
        # else:
        LOGGER.error(
            "Unable to install automatically. Check Chrome is available or try `blobtools view --driver firefox` to use firefox browser."
        )
        sys.exit(1)
    outdir = os.path.abspath(args["--out"])
    os.makedirs(Path(outdir), exist_ok=True)

    options = webdriver.ChromeOptions()
    if not args["--interactive"]:
        options.add_argument("headless")
    prefs = {
        "profile.default_content_settings.popups": 0,
        "download.default_directory": outdir,
    }

    options.add_experimental_option("prefs", prefs)
    display = start_display()
    chromedriver_path = which("chromedriver")
    if chromedriver_path is None:
        chromedriver_path = os.path.join(
            getsitepackages()[0], "chromedriver_binary", "chromedriver"
        )
    service = Service(
        executable_path=chromedriver_path,
        service_log_path=args["--driver-log"],
    )
    try:
        driver = webdriver.Chrome(
            options=options,
            # service=service_object,
            # executable_path=add option to set binary location,
            service=service,
        )
    except WebDriverException:
        LOGGER.error(
            "Unable to start Chrome using chromedriver, possible version mismatch."
        )
        LOGGER.error("Try `blobtools view --driver firefox` to use firefox browser.")
        if check_permission(which("chromedriver")):
            LOGGER.error(
                "Or try deleting the existing chromedriver executable from %s then rerun this command to install a compatible version.",
                which("chromedriver"),
            )
        else:
            LOGGER.error(
                "Or ask your administrator to ensure that the versions of Chrome and chromedriver are compatible"
            )
        display.stop()
        sys.exit(1)
    time.sleep(2)
    return driver, display


def stop_driver(driver, display, viewer):
    driver.quit()
    display.stop()
    if viewer is not None:
        viewer.send_signal(signal.SIGINT)


def static_view(args, loc, viewer):
    """Generate static images."""
    qstr = "staticThreshold=Infinity" + "&nohitThreshold=Infinity"
    qstr += "&plotGraphics=svg"
    file_stem = Path(args["DIRECTORY"]).name
    if file_stem == "_":
        file_stem = "FXWY01"
    if args["--format"] == "svg":
        qstr += "&svgThreshold=Infinity"
    shape = "circle"
    for param in args["--param"]:
        qstr += f"&{str(param)}"
        key, value = param.split("=")
        if key == "plotShape":
            shape = value

    timeout = int(args["--timeout"])
    outdir = os.path.abspath(args["--out"])
    driver = None

    def handle_error(err):
        """Release resources before quitting."""
        if viewer is not None:
            viewer.send_signal(signal.SIGINT)
        if driver is not None:
            driver.quit()
            display.stop()
        print(err)
        sys.exit(1)

    if args["--driver"] == "firefox":
        """View dataset in Firefox."""
        try:
            driver, display = firefox_driver(args)
        except Exception as err:
            format_exc(err)
            LOGGER.error(
                "Unable to start Firefox. Try `conda install -c conda-forge firefox geckodriver` or `blobtools view --driver chromium` to use Chrome browser."
            )
            sys.exit(1)

    elif args["--driver"] == "chromium":
        """View dataset in Chromium browser."""
        try:
            driver, display = chromium_driver(args)
        except Exception:
            LOGGER.error(
                "Unable to start Chrome. Try `blobtools view --driver firefox` to use Firefox browser."
            )
            sys.exit(1)
    else:
        handle_error(f'{args["--driver"]} is not a valid driver')

    try:
        view = args["--view"][0]
        if args["--preview"]:
            qstr += "#Filters"
        url = f"{loc}/{view}?{qstr}"
        print(f"Loading {url}")
        for i in range(10):
            try:
                driver.get(url)
            except Exception as err:
                if i < 5:
                    time.sleep(1)
                    continue
                handle_error(err)
            break

        for next_view in args["--view"]:
            if next_view != view:
                view = next_view
                url = f"{loc}/{view}?{qstr}"
                print(f"Navigating to {url}")
                try:
                    driver.get(url)
                except Exception as err:
                    handle_error(err)
            for fmt in args["--format"]:
                file = f"{file_stem}.{view}"
                if view == "blob":
                    file += f".{shape}"
                elif view == "busco":
                    view = f"all_{view}"
                    if fmt not in ("csv", "json"):
                        fmt = "json"
                file += f".{fmt}"
                print(f"Fetching {file}")
                el_id = f"{view}_save_{fmt}"
                print(f" - waiting for element {el_id}")
                unstable = True
                start_time = time.time()
                while unstable:
                    elapsed_time = time.time() - start_time
                    if timeout and elapsed_time > timeout:
                        handle_error("Timeout waiting for file")
                    try:
                        element = WebDriverWait(driver, timeout).until(
                            EC.visibility_of_element_located((By.ID, el_id))
                        )
                        element.click()
                        unstable = False
                        file_name = f"{outdir}/{file}"
                        print(" - waiting for file '%s'" % file_name)
                        file_ready(file_name, timeout, handle_error)
                    except Exception as err:
                        unstable = True
                        time.sleep(1)

        for preview in args["--preview"]:
            print(f"Creating {preview} preview")
            for fmt in args["--format"]:
                el_id = f"{preview}_preview_save_{fmt}"
                file = f'{Path(args["DIRECTORY"]).name}.{preview}.preview.{fmt}'
                try:
                    element = WebDriverWait(driver, timeout).until(
                        EC.visibility_of_element_located((By.ID, el_id))
                    )
                    element.click()
                    file_name = f"{outdir}/{file}"
                    print("waiting for file '%s'" % file_name)
                    file_ready(file_name, timeout, handle_error)
                except Exception as err:
                    handle_error(err)
        if viewer is not None:
            viewer.send_signal(signal.SIGINT)
        driver.quit()
        display.stop()
    except Exception as err:
        handle_error(err)
    return True


def plot_view(args):
    """Generate static images using blobtk plot."""
    param_names = {
        "segments": "segments",
        "circumferenceScale": "max_span",
        "radiusScale": "max_scaffold",
        "maxSpan": "max_span",
        "maxScaffold": "max_scaffold",
        "xField": "x_field",
        "yField": "y_field",
        "zField": "z_field",
        "catField": "cat_field",
        "resolution": "resolution",
        "histHeight": "hist_height",
        "zReducer": "reducer_function",
        "zScale": "scale_function",
        "plotScale": "scale_factor",
        "reducerFunction": "reducer_function",
        "scaleFunction": "scale_function",
        "scaleFactor": "scale_factor",
        "xLimit": "x_limit",
        "yLimit": "y_limit",
        "catCount": "cat_count",
        "showLegend": "show_legend",
        "legend": "show_legend",
        "catOrder": "cat_order",
        "origin": "origin",
        "palette": "palette",
        "shape": "shape",
        # color,
    }
    file_stem = Path(args["DIRECTORY"]).name
    outdir = os.path.abspath(args["--out"])
    if args["--format"] == "svg":
        qstr = "&svgThreshold=Infinity"
    shape = "circle"
    suffix = ["Min", "Max", "Inv", "Key", "Keys"]
    filters = []
    flags = {}
    for param in args["--param"]:
        key, value = param.split("=")
        if key == "plotShape":
            shape = value
        parts = key.split("--")
        if len(parts) == 2:
            if any(parts[1] in s for s in suffix):
                filters.append(param)
            elif parts[1] == "Order":
                flags["cat_order"] = value
        elif key in param_names:
            flags[param_names[key]] = value
        elif key in param_names.values():
            flags[key] = value
        elif key.startswith("color"):
            if "color" not in flags:
                flags["color"] = []
            flags["color"].append(f'{key.replace("color", "")}={value}')
    for view in args["--view"]:
        for image_format in args["--format"]:
            if view == "blob":
                # if shape != "circle":
                outfile = f"{outdir}/{file_stem}.{view}.{shape}.{image_format}"
            else:
                outfile = f"{outdir}/{file_stem}.{view}.{image_format}"
            print(outfile)

            plot.plot(
                blobdir=args["DIRECTORY"],
                output=outfile,
                view=view,
                filter=filters,
                **flags,
            )
    quit()


def interactive_view(args, loc, viewer, level):
    if args["--driver"] == "firefox":
        """View dataset in Firefox."""
        driver, display = firefox_driver(args)
    elif args["--driver"] == "chromium":
        """View dataset in Chromium browser."""
        driver, display = chromium_driver(args)
    qstr = "".join(f"&{str(param)}" for param in args["--param"])
    try:
        view = args["--view"][0]
        if args["--preview"]:
            qstr += "#Filters"
        if level == "dataset":
            url = f"{loc}/{view}"
            if qstr:
                url += f"?{qstr}"
        else:
            url = loc if loc.endswith("all") else f"{loc}/all"
        print(f"Loading {url}")
        try:
            driver.get(url)
        except Exception as err:
            print(err)
        poll = viewer.poll()
        while poll is None:
            time.sleep(5)
            poll = viewer.poll()
        stop_driver(driver, display, viewer)
    except Exception as err:
        print(err)
        stop_driver(driver, display, viewer)
    return True


def remote_view(args, loc, viewer, port, api_port, level, remote):
    """View dataset remotely."""
    qstr = "".join(f"&{str(param)}" for param in args["--param"])
    try:
        view = args["--view"][0]
        if args["--preview"]:
            qstr += "#Filters"
        if level == "dataset":
            url = f"{loc}/{view}"
            if qstr:
                url += f"?{qstr}"
            print(f"View dataset at {url}")
        else:
            print(f"View datasets at {loc}")
        if remote:
            print("For remote access use:")
            print(
                "    ssh -L %d:127.0.0.1:%d -L %d:127.0.0.1:%d username@remote_host"
                % (port, port, api_port, api_port)
            )
        while True:
            time.sleep(5)
    except Exception as err:
        print("remote exception")
        print(err)
        if viewer is not None:
            viewer.send_signal(signal.SIGINT)
    return True


def main(args):
    """Entrypoint for blobtools view."""
    if not args["--plot"]:
        loc, viewer, port, api_port, level = test_loc(args)
        with contextlib.suppress(KeyboardInterrupt):
            if args["--interactive"]:
                interactive_view(args, loc, viewer, level)
            elif args["--remote"]:
                remote_view(args, loc, viewer, port, api_port, level, True)
            elif args["--local"]:
                remote_view(args, loc, viewer, port, api_port, level, False)
            static_view(args, loc, viewer)
        time.sleep(1)
        if viewer is not None:
            viewer.send_signal(signal.SIGINT)
            time.sleep(1)
    else:
        plot_view(args)


def cli():
    """Entry point."""
    if len(sys.argv) == sys.argv.index(__name__.split(".")[-1]) + 1:
        args = docopt(__doc__, argv=[])
    else:
        args = docopt(__doc__, version=__version__)
    if not os.path.exists(os.environ["HOME"]):
        os.mkdir(os.environ["HOME"])
    main(args)


if __name__ == "__main__":
    cli()
