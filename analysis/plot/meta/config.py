import logging
import os

# Log Levels
LOG_LEVEL = {
    "debug": logging.DEBUG,
    "info": logging.INFO,
    "warning": logging.WARNING,
    "error": logging.ERROR,
}

# Cohort groupings for GLM and Plots
EXPOSED_COHORTS = ["RADAR", "CRU", "PILOT"]
LEVELS_BASE = ["INOVA", "RADAR", "CRU", "PILOT"]
CATEGORY_ORDER_BASE = LEVELS_BASE

LEVELS_SUBGROUPS = [
    "INOVA::NO_EXPOSURE",
    "RADAR::NO_EXPOSURE",
    "RADAR::EXPOSED",
    "CRU::EXPOSED",
    "CRU::NO_EXPOSURE",
    "PILOT::EXPOSED",
    "PILOT::NO_EXPOSURE",
]
CATEGORY_ORDER_SUBGROUPS = LEVELS_SUBGROUPS

# Cache Dir
# Cache dir is taken as environment table at load time, before the actual
# arguments are parsed.
CACHE_DIR = os.path.abspath(os.getenv("CACHE_DIR", "./output"))
if not os.path.isdir(CACHE_DIR):
    os.makedirs(CACHE_DIR, exist_ok=True)
