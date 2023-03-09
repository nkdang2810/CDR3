import logging
from pathlib import Path
from datetime import datetime
import time

class CustomFormatter(logging.Formatter):

    grey = "\x1b[38;20m"
    yellow = "\x1b[33;20m"
    red = "\x1b[31;20m"
    bold_red = "\x1b[31;1m"
    reset = "\x1b[0m"
    format = "%(asctime)s - %(levelname)s - %(message)s"

    FORMATS = {
        logging.DEBUG: grey + format + reset,
        logging.INFO: grey + format + reset,
        logging.WARNING: yellow + format + reset,
        logging.ERROR: red + format + reset,
        logging.CRITICAL: bold_red + format + reset
    }

    def format(self, record):
        log_fmt = self.FORMATS.get(record.levelno)
        formatter = logging.Formatter(log_fmt)
        return formatter.format(record)

class Log():
    def __init__(self, path):
        # logging set up
        Path.mkdir(path, parents=True, exist_ok=True)

        self.rootLogger = logging.getLogger()
        self.rootLogger.setLevel(logging.DEBUG)
        # logging.getLogger('matplotlib').setLevel(logging.ERROR)

        self.fileHandler = logging.FileHandler(path.joinpath(
            f'{datetime.today().strftime("%m-%d--%H-%M-%S")}.log'))
        self.fileHandler.setFormatter(logging.Formatter("%(asctime)s - %(levelname)s - %(message)s"))
        self.rootLogger.addHandler(self.fileHandler)

        self.consoleHandler = logging.StreamHandler()
        self.consoleHandler.setFormatter(CustomFormatter())
        self.rootLogger.addHandler(self.consoleHandler)
        self.start_time_entire = time.time()

    def detach(self):
        running_time = time.time() - self.start_time_entire
        self.info(
            f"Total running time: {int(running_time//60)}:{int(running_time%60):0>2d}\n")
        self.rootLogger.removeHandler(self.fileHandler)
        self.rootLogger.removeHandler(self.consoleHandler)

    def info(self, mess):
        logging.info(mess)

    def error(self, mess):
        logging.error(mess)

    def warning(self,mess):
        logging.warning(mess)
        