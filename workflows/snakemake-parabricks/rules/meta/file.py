from posixpath import abspath
import typing
import pandas
import os
from enum import Enum


class FILE_TYPE(Enum):
    RAW_READS = 'raw-reads'
    ALIGNMENT = 'alignment'
    VARIANTS  = 'variants'
    UNKNOWN   = 'unknown'


INDEX_SUFFIX = {
    FILE_TYPE.ALIGNMENT: [ 'bai', 'bam.bai' ],
    FILE_TYPE.VARIANTS: [ 'idx', 'tbi' ],
}


def ensure_file(type: FILE_TYPE, files: pandas.DataFrame, samples: pandas.DataFrame,
                is_paired_end: bool = True, ensure_idx: bool = True) -> typing.Dict[str, typing.List[typing.Union[str, typing.Tuple[str, str]]]]:
    """Ensure the files exist and have the correct setup.

    We assume a paired_end library by default for `FILE_TYPE.RAW_READS`, and
    look for index file for bam files by default (e.g. len(fastqs for one sample) % 2 == 0).
    
    Parameters
    ----------
    files: pandas.DataFrame
           DataFrame containing file paths for files concerning the workflow.
    samples: pandas.DataFrame
             DataFrame samplesheet with sample descriptions for the workflow.
    is_paired_end: bool
                   Check for paired end fastq (raw-reads) files.
    ensure_idx: bool
                Also verify that the index file for the corresponding file exists.
    
    Returns
    -------
    dict
        Dictionary mapping sample_id's to a set of fastq files. Returns values
        as tuples if ensure_idx is set. Tuple values are the paths to the file
        and the retrieved index respectively.

    Raises
    ------
    Error
        If a file could not be found or the paired_end assumption is violated.
    """
    ret_val = {}
    known_samples = set(samples['sample_id']) if samples is not None else None
    for _, row in files[files['type'] == type.value].iterrows():
        assert row['sample_id'] in known_samples if known_samples is not None else True
        if row['sample_id'] not in ret_val:
            ret_val[row['sample_id']] = []
        assert os.path.isfile(row['path'])
        if ensure_idx:
            idx = None
            for x in INDEX_SUFFIX[type]:
                if os.path.isfile("{}.{}".format(row['path'], x)):
                    idx = "{}.{}".format(row['path'], x)
            assert idx is not None
            ret_val[row['sample_id']].append((os.path.abspath(row['path']), os.path.abspath(idx)))
        else:
            ret_val[row['sample_id']].append(os.path.abspath(row['path']))
    
    assert all(len(arr) % 2 == 0 for arr in ret_val.values()) if is_paired_end else True
    return ret_val


def ensure_fastqs(files: pandas.DataFrame, samples: pandas.DataFrame,
                  is_paired_end: bool = True) -> typing.Dict[str, typing.List[typing.Union[str, typing.Tuple[str, str]]]]:
    return ensure_file(FILE_TYPE.RAW_READS, files, samples, is_paired_end=is_paired_end, ensure_idx=False)


def ensure_bams(files: pandas.DataFrame, samples: pandas.DataFrame,
                ensure_idx: bool = True) -> typing.Dict[str, typing.List[typing.Union[str, typing.Tuple[str, str]]]]:
    return ensure_file(FILE_TYPE.ALIGNMENT, files, samples, ensure_idx=ensure_idx, is_paired_end=False)


def ensure_vcfs(files: pandas.DataFrame, samples: pandas.DataFrame,
                ensure_idx: bool = True) -> typing.Dict[str, typing.List[typing.Union[str, typing.Tuple[str, str]]]]:
    return ensure_file(FILE_TYPE.VARIANTS, files, samples, ensure_idx=ensure_idx, is_paired_end=False)