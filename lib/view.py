#!/usr/bin/env python3

# pylint: disable=no-member, too-many-branches, too-many-locals, too-many-statements, broad-except

"""
Generate plots using BlobToolKit Viewer.

Usage:
  blobtools view [--format STRING...] [--host STRING] [--out PATH]
  [--param STRING...] [--prefix STRING] [--preview STRING...]
  [--timeout INT] [--view STRING...] DATASET

Options:
      --format STRING     Image format (svg|png). [Default: png]
      --host STRING       Hostname. [Default: http://localhost:8080]
      --out PATH          Directory for outfiles. [Default: .]
      --param key=value   Query string parameter.
      --prefix STRING     URL prefix. [Default: view]
      --preview STRING    Field name.
      --timeout INT       Time to wait for page load in seconds. Default (0) is no timeout. [Default: 0]
      --view STRING       Plot type (blob|cumulative|snail). [Default: blob]
"""
import os
import time
from pathlib import Path
from docopt import docopt
from selenium import webdriver
from selenium.webdriver.firefox.options import Options
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from pyvirtualdisplay import Display


def file_ready(file_path):
    """Check if file is ready."""
    while not os.path.exists(file_path):
        time.sleep(1)
        if os.path.isfile(file_path):
            return True
        raise ValueError("%s isn't a file!" % file_path)


def main():
    """Entrypoint for blobtools view."""
    args = docopt(__doc__)
    loc = "%s/%s/dataset/%s" % (args['--host'], args['--prefix'], args['DATASET'])
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
    options.set_headless(headless=False)
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
    except Exception as err:
        print(err)
        driver.quit()
        display.popen.terminate()


if __name__ == '__main__':
    main()
