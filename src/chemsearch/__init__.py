# -*- coding: utf-8 -*-
import os
import sys
import logging
import argparse
from pkg_resources import get_distribution, DistributionNotFound


__author__ = "Stephen Gaffney"
__copyright__ = "Stephen Gaffney"
__license__ = "gpl3"

_logger = logging.getLogger(__name__)


def main(args):
    """Main entry point allowing external calls

    Args:
      args ([str]): command line parameter list
    """
    from . import admin, drive
    args = parse_args(args)
    setup_logging(args.loglevel)
    local_archive_path = args.path
    admin.reload_env()

    meta = drive.Meta().build()
    drive.create_local_archive(meta.molfiles,
                               local_root=local_archive_path,
                               files_resource=meta.files_resource)
    summary_df = admin.assemble_archive_metadata(local_archive_path)
    _logger.info(f"Info generated for {len(summary_df)} molecules.")
    _logger.debug("Script complete.")


def run():
    """Entry point for console_scripts
    """
    main(sys.argv[1:])


def setup_logging(loglevel):
    """Setup basic logging

    Args:
      loglevel (int): minimum loglevel for emitting messages
    """
    log_format = "[%(asctime)s] %(levelname)s:%(name)s:%(message)s"
    logging.basicConfig(level=loglevel, stream=sys.stdout,
                        format=log_format, datefmt="%Y-%mol-%d %H:%M:%S")


def parse_args(args):
    """Parse command line parameters

    Args:
      args ([str]): command line parameters as list of strings

    Returns:
      :obj:`argparse.Namespace`: command line parameters namespace
    """
    parser = argparse.ArgumentParser(
        description="Create local molecules repository.")
    parser.add_argument(
        "--version",
        action="version",
        version="chemsearch {ver}".format(ver=__version__))
    parser.add_argument(
        dest="path",
        help="path to local db directory (e.g. local_db)",
        type=str)
    parser.add_argument(
        "-v",
        "--verbose",
        dest="loglevel",
        help="set loglevel to INFO",
        action="store_const",
        const=logging.INFO)
    parser.add_argument(
        "-vv",
        "--very-verbose",
        dest="loglevel",
        help="set loglevel to DEBUG",
        action="store_const",
        const=logging.DEBUG)
    return parser.parse_args(args)


try:
    # Change here if project is renamed and does not equal the package name
    dist_name = __name__
    __version__ = get_distribution(dist_name).version
except DistributionNotFound:
    __version__ = 'unknown'
finally:
    del get_distribution, DistributionNotFound
