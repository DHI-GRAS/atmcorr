import sys
import logging


def set_logger(
        set_exception_handler=True, name=None, level='INFO',
        external_loggers=[]):
    """Create a logger for an entire processing session"""
    myname = name or __name__
    rootlogger = logging.getLogger(myname)
    all_loggers = external_loggers + [myname]

    for lname in all_loggers:
        logger = logging.getLogger(lname)
        logger.setLevel(level)

    logformatter = logging.Formatter(
            "%(asctime)s %(levelname)s %(message)s",
            datefmt='%Y-%m-%d %H:%M:%S')

    consolehandler = logging.StreamHandler()
    consolehandler.setFormatter(logformatter)
    for lname in all_loggers:
        logger = logging.getLogger(lname)
        logger.addHandler(consolehandler)

    if set_exception_handler:
        # define handler function
        def log_uncaught_exception(exc_type, exc_value, exc_traceback):
            if issubclass(exc_type, KeyboardInterrupt):
                sys.__excepthook__(exc_type, exc_value, exc_traceback)
                return
            rootlogger.critical(
                    "Uncaught exception!", exc_info=(exc_type, exc_value, exc_traceback))

        # install exception handler
        sys.excepthook = log_uncaught_exception

    return rootlogger


def set_cli_logger(**kwargs):
    logging.captureWarnings(True)
    set_logger(name='atmospheric_correction', **kwargs)
