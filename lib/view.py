#!/usr/bin/env python3

# pylint: disable=no-member, too-many-branches, too-many-locals, too-many-statements, broad-except

"""
Generate plots using BlobToolKit Viewer.

Usage:
  blobtools view [--format STRING...] [--host STRING] [--out PATH]
  [--param STRING...] [--ports RANGE] [--prefix STRING] [--preview STRING...]
  [--timeout INT] [--view STRING...] DATASET

Options:
      --format STRING     Image format (svg|png). [Default: png]
      --host STRING       Hostname. [Default: http://localhost]
      --ports RANGE       Port range for viewer and API. [Default: 8000-8099]
      --out PATH          Directory for outfiles. [Default: .]
      --param key=value   Query string parameter.
      --prefix STRING     URL prefix. [Default: view]
      --preview STRING    Field name.
      --timeout INT       Time to wait for page load in seconds. Default (0) is no timeout. [Default: 0]
      --view STRING       Plot type (blob|cumulative|snail). [Default: blob]
"""
import os
import shlex
import time
from pathlib import Path
from subprocess import PIPE, Popen
from docopt import docopt
from selenium import webdriver
from selenium.webdriver.firefox.options import Options
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from pyvirtualdisplay import Display
from host import test_port


def file_ready(file_path):
    """Check if file is ready."""
    while not os.path.exists(file_path):
        time.sleep(1)
        if os.path.isfile(file_path):
            return True
        raise ValueError("%s isn't a file!" % file_path)


def test_loc(args):
    """See if dataset needs to be hosted and, if so, find an empty port."""
    info = args['--host'].split(':')
    if len(info) >= 2 and info[1] != '//localhost':
        return "%s/%s/dataset/%s" % (args['--host'], args['--prefix'], args['DATASET'])
    if len(info) == 1 and info[0] != 'localhost':
        # need to add test for http vs https
        return "http://%s/%s/dataset/%s" % (args['--host'], args['--prefix'], args['DATASET'])
    if len(info) == 3:
        port = info[2]
        available = test_port(port, 'test')
        if available:
            print("ERROR: No service running on port %s" % port)
            print("       Unable to connect to %s" % args['--host'])
            exit(1)
        else:
            return "%s/%s/dataset/%s" % (args['--host'], args['--prefix'], args['DATASET'])
    if not Path(args['DATASET']).exists():
        print("ERROR: DATASET '%s' must be a valid path to begin hosting.")
        exit(1)
    dataset = Path(args['DATASET']).name
    parent = Path(args['DATASET']).resolve().absolute().parent
    port_range = args['--ports'].split('-')
    api_port = False
    port = False
    for i in range(int(port_range[0]), int(port_range[1])):
        if test_port(i, 'test'):
            if not api_port:
                api_port = i
                continue
            if not port:
                port = i
            break
    print("api: %d, viewer: %d" % (api_port, port))
    directory = Path(__file__).resolve().parent.parent
    cmd = "%s/blobtools host --port %d --api-port %d %s" % (directory, port, api_port, parent)
    process = Popen(shlex.split(cmd),
                    stdout=PIPE,
                    stderr=PIPE,
                    encoding='ascii')
    loc = "%s:%d/%s/dataset/%s" % (args['--host'], port, args['--prefix'], dataset)
    print("Waiting 15 seconds for Viewer to start...")
    for i in range(1, 15):
        poll = process.poll()
        if poll is None:
            time.sleep(1)
        else:
            print("ERROR: Viewer quit unexpectedly")
            print("Unable to run: %s" % cmd)
            exit(1)
    return loc, process


def main():
    """Entrypoint for blobtools view."""
    args = docopt(__doc__)
    loc, viewer = test_loc(args)
    try:  # pylint: disable=too-many-nested-blocks
        qstr = "staticThreshold=Infinity"
        qstr += "&nohitThreshold=Infinity"
        qstr += "&plotGraphics=svg"
        if args['--format'] == 'svg':
            qstr += "&svgThreshold=Infinity"
        shape = 'square'
        for param in args['--param']:
            qstr += "&%s" % str(param)
            key, value = param.split('=')
            if key == 'plotShape':
                shape = value

        timeout = int(args['--timeout'])

        outdir = os.path.abspath(args['--out'])
        os.makedirs(Path(outdir), exist_ok=True)

        profile = webdriver.FirefoxProfile()
        profile.set_preference('browser.download.folderList', 2)
        profile.set_preference('browser.download.manager.showWhenStarting', False)
        profile.set_preference('browser.download.dir', outdir)
        profile.set_preference('browser.download.lastDir', args['--out'])
        profile.set_preference('browser.helperApps.neverAsk.saveToDisk',
                               'image/png, image/svg+xml, text/csv, text/plain, application/json')

        options = Options()
        options.set_headless(headless=True)
        display = Display(visible=0, size=(800, 600))
        display.start()
        driver = webdriver.Firefox(options=options, firefox_profile=profile)
        try:
            view = args['--view'][0]
            if args['--preview']:
                qstr += '#Filters'
            url = "%s/%s?%s" % (loc, view, qstr)
            print("Loading %s" % url)
            try:
                driver.get(url)
            except Exception as err:
                print(err)

            for next_view in args['--view']:
                if next_view != view:
                    view = next_view
                    url = "%s/%s?%s" % (loc, view, qstr)
                    print("Navigating to %s" % url)
                    try:
                        driver.get(url)
                    except Exception as err:
                        print(err)
                for fmt in args['--format']:
                    file = "%s.%s" % (args['DATASET'], view)
                    if view == 'blob':
                        file += ".%s" % shape
                    elif view == 'busco':
                        view = "all_%s" % view
                        if fmt not in ('csv', 'json'):
                            fmt = 'json'
                    file += ".%s" % fmt
                    print("Fetching %s" % file)
                    el_id = "%s_save_%s" % (view, fmt)
                    print("waiting for element %s" % el_id)
                    unstable = True
                    while unstable:
                        try:
                            element = WebDriverWait(driver, timeout).until(
                                EC.visibility_of_element_located((By.ID, el_id))
                            )
                            element.click()
                            unstable = False
                            file_name = "%s/%s" % (outdir, file)
                            print("waiting for file '%s'" % file_name)
                            file_ready(file_name)
                        except Exception as err:
                            time.sleep(1)

            for preview in args['--preview']:
                print("Creating %s preview" % preview)
                for fmt in args['--format']:
                    el_id = "%s_preview_save_%s" % (preview, fmt)
                    file = "%s.%s.preview.%s" % (args['DATASET'], preview, fmt)
                    try:
                        element = WebDriverWait(driver, timeout).until(
                            EC.visibility_of_element_located((By.ID, el_id))
                        )
                        element.click()
                        file_name = "%s/%s" % (outdir, file)
                        print("waiting for file '%s'" % file_name)
                        file_ready(file_name)
                    except Exception as err:
                        print(err)
            driver.quit()
            display.popen.terminate()
            os.killpg(os.getpgid(viewer.pid), 15)
        except Exception as err:
            print(err)
            driver.quit()
            display.popen.terminate()
            os.killpg(os.getpgid(viewer.pid), 15)
    except Exception as err:
        os.killpg(os.getpgid(viewer.pid), 15)


if __name__ == '__main__':
    main()
