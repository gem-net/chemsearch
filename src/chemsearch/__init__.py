# -*- coding: utf-8 -*-
import os
import sys
import logging
import argparse
from typing import Union
from pkg_resources import get_distribution, DistributionNotFound


__author__ = "Stephen Gaffney"
__copyright__ = "Stephen Gaffney"
__license__ = "gpl3"


def _create_logger(log_level: Union[str, None] = None):
    if log_level is None:
        log_level = os.environ.get('LOG_LEVEL', 'INFO').upper()
    log_level_int = getattr(logging, log_level)

    ch = logging.StreamHandler()
    ch.setLevel(log_level_int)
    formatter = logging.Formatter(
        fmt='%(asctime)s - %(message)s',
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    ch.setFormatter(formatter)

    _logger = logging.getLogger(__name__)
    _logger.addHandler(ch)
    return _logger


logger = _create_logger()


def main(args):
    """Main entry point allowing external calls

    Args:
      args ([str]): command line parameter list
    """
    from . import admin, drive
    args = parse_args(args)
    local_archive_path = args.path
    admin.reload_env()

    meta = drive.Meta().build()
    drive.create_local_archive(meta.molfiles,
                               local_root=local_archive_path,
                               files_resource=meta.files_resource)
    summary_df = admin.assemble_archive_metadata(local_archive_path)
    logger.info(f"Info generated for {len(summary_df)} molecules.")
    logger.debug("Script complete.")


def run():
    """Entry point for console_scripts
    """
    main(sys.argv[1:])


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


if sys.version_info[:2] >= (3, 8):
    # TODO: Import directly (no need for conditional) when `python_requires = >= 3.8`
    from importlib.metadata import PackageNotFoundError, version  # pragma: no cover
else:
    from importlib_metadata import PackageNotFoundError, version  # pragma: no cover

try:
    # Change here if project is renamed and does not equal the package name
    dist_name = __name__
    __version__ = version(dist_name)
except PackageNotFoundError:  # pragma: no cover
    __version__ = "unknown"
finally:
    del version, PackageNotFoundError
