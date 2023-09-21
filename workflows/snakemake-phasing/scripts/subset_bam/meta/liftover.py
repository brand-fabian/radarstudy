from liftover import get_lifter
import typing
import os
import logging
from collections import namedtuple
import pandas

logger = logging.getLogger("subset-bam")

LIFTER_REFERENCES = {
    "grch37": "hg19",
    "hg19": "hg19",
    "grch38": "hg38",
    "hg38": "hg38",
}

class Forklift:
    ColumnInfo = namedtuple("ColumnInfo", ["site"])
    DefaultColumnInfo = ColumnInfo(site="locus")
    """Liftover locations from source to target reference"""
    def __init__(self, source: str, target: str):
        self._lifter = get_lifter(
            LIFTER_REFERENCES[source],
            LIFTER_REFERENCES[target]
        )

    def liftover(self, sites: pandas.DataFrame, col_info: ColumnInfo = DefaultColumnInfo):
        def _lift(site: str):
            try:
                return ":".join(map(str, self._lifter.query(*[
                    x[1] if x[0] != 1 else int(x[1]) for x in enumerate(site.split(":")) 
                ])[0][:-1]))
            except KeyError as e:
                logger.warning("failed liftover for {}".format(site))
                return None
            except IndexError as e:
                logger.warning("no location found for {}".format(site))
                return None
        remapped_col = list(sites[col_info.site].apply(_lift))
        return remapped_col