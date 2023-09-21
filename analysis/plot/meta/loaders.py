import logging
import os
import pickle
import sys
import typing

import hail as hl
import pandas

from meta.checkpoint import NamedCheckpoint
from meta.config import CACHE_DIR
from meta.helpers import downsample_age_match
from meta.serializers import FILE_TYPE

logger = logging.getLogger(__name__)

###############################################################################
# Input data loader                                                           #
###############################################################################
class Loader:
    pass


class MsdnLoader(Loader):
    @staticmethod
    @NamedCheckpoint.checkpoint(
        prefix=CACHE_DIR, name="msdn_final.ht", output_type=FILE_TYPE.HAIL_TABLE
    )
    def load_input(
        input_ht: str,
        excluded_samples: str = r"R_(0048|0158|0290|0385|0180|0214|0066|0146|0273).*",
    ) -> hl.Table:
        """Load a MSDN hail table.

        Load the msdn hail table from the path `input_ht`. Samples matched by the
        `excluded_samples` regex are excluded from the results table.
        """
        msdn_ht = hl.read_table(input_ht)

        for ref, alt in [
            (
                "C",
                "A",
            )
        ]:
            msdn_ht = msdn_ht.filter(
                (
                    (msdn_ht.alleles[0] == ref)
                    & (msdn_ht.alleles[1] == alt)
                    & (msdn_ht.AF > 0.45)
                    & (msdn_ht.AF < 0.55)
                )
                | ((msdn_ht.alleles[0] != ref) | (msdn_ht.alleles[1] != alt))
            )
            msdn_ht = msdn_ht.filter(
                (
                    (msdn_ht.previous.alleles[0] == ref)
                    & (msdn_ht.previous.alleles[1] == alt)
                    & (msdn_ht.previous.AF > 0.45)
                    & (msdn_ht.previous.AF < 0.55)
                )
                | (
                    (msdn_ht.previous.alleles[0] != ref)
                    | (msdn_ht.previous.alleles[1] != alt)
                )
            )

        for ref, alt, cohort, cutoff in [
            ("A", "T", "CRU", 0.45),
            ("G", "T", "CRU", 0.45),
            ("C", "T", "CRU", 0.35),
            ("G", "A", "CRU", 0.35),
        ]:
            msdn_ht = msdn_ht.filter(
                (
                    (msdn_ht.alleles[0] == ref)
                    & (msdn_ht.alleles[1] == alt)
                    & (msdn_ht.AF > cutoff)
                    & (msdn_ht.AF < 1 - cutoff)
                    & (msdn_ht.cohort == cohort)
                )
                | (
                    (msdn_ht.alleles[0] != ref)
                    | (msdn_ht.alleles[1] != alt)
                    | (msdn_ht.cohort != cohort)
                )
            )
            msdn_ht = msdn_ht.filter(
                (
                    (msdn_ht.previous.alleles[0] == ref)
                    & (msdn_ht.previous.alleles[1] == alt)
                    & (msdn_ht.previous.AF > cutoff)
                    & (msdn_ht.previous.AF < 1 - cutoff)
                    & (msdn_ht.cohort == cohort)
                )
                | (
                    (msdn_ht.previous.alleles[0] != ref)
                    | (msdn_ht.previous.alleles[1] != alt)
                    | (msdn_ht.cohort != cohort)
                )
            )

        msdn_ht = msdn_ht.filter(~msdn_ht.s.matches(excluded_samples))
        msdn_ht = msdn_ht.annotate(
            cohort=(
                hl.case()
                .when(msdn_ht.s.startswith("R_"), "RADAR")
                .when(msdn_ht.s.startswith("UTRI"), "CRU")
                .when(msdn_ht.s.matches(r"\d{2}-\d{4}"), "PILOT")
                .default("INOVA")
            )
        )
        return msdn_ht

    @staticmethod
    @NamedCheckpoint.checkpoint(
        prefix=CACHE_DIR, name="msdn_df.pickle", output_type=FILE_TYPE.PICKLE
    )
    def load_dataframe(msdn_ht: hl.Table) -> pandas.DataFrame:
        """Convert a msdn ht to a pandas DataFrame"""
        return msdn_ht.to_pandas()

    @staticmethod
    @NamedCheckpoint.checkpoint(
        prefix=CACHE_DIR, name="metadata.pickle", output_type=FILE_TYPE.PICKLE
    )
    def load_metadata(
        radar: str,
        inova: str,
        cru: str,
        pilot: str = None,
        excluded_samples: str = r"R_(0048|0158|0290|0385|0180|0214|0066|0146|0273).*",
    ) -> pandas.DataFrame:
        """Load study metadata from three given input paths.

        Parameters
        ----------
        radar : str
                Path to the radar metadata table (.csv)
        inova : str
                Path to the inova metadata table (.xlsx)
        cru : str
              Path to the CRU metadata table (.csv)
        excluded_samples : str
                           Regex with sample names to exclude by match statements.

        Returns
        -------
        pandas.DataFrame
            A DataFrame containing metadata of all three studies sorted by
            the cohort key.
        """
        radarstudy_meta = pandas.read_csv(radar)
        radarstudy_meta["cohort"] = "RADAR"

        cru_meta = pandas.read_csv(cru)
        cru_meta["cohort"] = "CRU"

        pilot_meta = None
        if pilot is not None:
            pilot_meta = pandas.read_csv(pilot)
            pilot_meta["cohort"] = "PILOT"

        inova_meta = pandas.read_excel(inova)
        inova_data = []
        for _, row in inova_meta.iterrows():
            fam_id = "-".join(row.BabyID.split("-")[0:2])
            inova_data.append(
                [
                    0,
                    row.Baby_sample,
                    row.Father_Age,
                    row.Mother_Age,
                    "UNKNOWN",
                    0,
                    0,
                    "NO_EXPOSURE",  # exposure_status
                    None,  # device
                    0,  # monthly_dosage
                    0,  # monthly_dosage_gonads
                    0,  # sum_dosage
                    0,  # sum_dosage_gonads
                    "NA",  # no_exposure_remark
                    "INOVA",  # group|cohort
                ]
            )

        md_dfs = [
            radarstudy_meta,
            pandas.DataFrame(inova_data, columns=radarstudy_meta.columns),
            cru_meta,
        ]
        if pilot_meta is not None:
            md_dfs.append(pilot_meta)

        age_data = pandas.concat(md_dfs).sort_values(["cohort", "rid"])
        for col in [
            "monthly_dosage",
            "monthly_dosage_gonads",
            "sum_dosage",
            "sum_dosage_gonads",
            "father_age",
            "mother_age",
            "service_begin",
            "service_end",
        ]:
            age_data[col] = age_data[col].fillna(0)
        for col in ["no_exposure_remark", "device"]:
            age_data[col] = age_data[col].fillna("UNKNOWN")

        age_data = age_data[~age_data.rid.str.match(excluded_samples)]
        return age_data

    @staticmethod
    @NamedCheckpoint.checkpoint(
        prefix=CACHE_DIR,
        name="msdn_downsample_age_match.pickle",
        output_type=FILE_TYPE.PICKLE,
    )
    def load_downsample_age_match(
        msdn_df: pandas.DataFrame, age_data: pandas.DataFrame, factor: int
    ):
        downsampled_metadata_df, matching = downsample_age_match(
            age_data, factor=factor
        )
        downsample_msdn_df = msdn_df[
            msdn_df["s"].isin(set(downsampled_metadata_df["rid"]))
        ]
        return downsample_msdn_df, downsampled_metadata_df, matching

    @staticmethod
    @NamedCheckpoint.checkpoint(
        prefix=CACHE_DIR,
        name="graphtyper_msdn_filtered.pickle",
        output_type=FILE_TYPE.PICKLE,
    )
    def load_graphtyper_filter(
        msdn_df: pandas.DataFrame,
        graphtyper_col: str = "graphtyper_evidence",
        expected_values: typing.Set[str] = set(
            ["CONFIRMATION", "PARTIAL_CONFIRMATION"]
        ),
    ):
        msdn_df = msdn_df[msdn_df[graphtyper_col].isin(expected_values)]
        return msdn_df


class DnmLoader(Loader):
    @staticmethod
    @NamedCheckpoint.checkpoint(
        prefix=CACHE_DIR, name="dnm_final.ht", output_type=FILE_TYPE.HAIL_TABLE
    )
    def load_input(
        input_mt: str,
        excluded_samples: str = r"R_(0048|0158|0290|0385|0180|0214|0066|0146|0273).*",
    ):
        dnm_mt = hl.read_matrix_table(input_mt)
        dnm_mt = dnm_mt.annotate_cols(
            cohort=(
                hl.case()
                .when(dnm_mt.s.startswith("R_"), "RADAR")
                .when(dnm_mt.s.startswith("UTRI"), "CRU")
                .when(dnm_mt.s.matches(r"\d{2}-\d{4}"), "PILOT")
                .default("INOVA")
            )
        )
        dnm_mt = dnm_mt.filter_rows(dnm_mt.locus.contig != "X")

        for ref, alt in [
            (
                "C",
                "A",
            )
        ]:
            dnm_mt = dnm_mt.filter_entries(
                (
                    (dnm_mt.alleles[0] == ref)
                    & (dnm_mt.alleles[1] == alt)
                    & (dnm_mt.AF > 0.45)
                    & (dnm_mt.AF < 0.55)
                )
                | ((dnm_mt.alleles[0] != ref) | (dnm_mt.alleles[1] != alt))
            )

        for ref, alt, cohort in [
            ("A", "T", "CRU"),
            ("G", "T", "CRU"),
            ("C", "T", "CRU"),
            ("G", "A", "CRU"),
        ]:
            dnm_mt = dnm_mt.filter_entries(
                (
                    (dnm_mt.alleles[0] == ref)
                    & (dnm_mt.alleles[1] == alt)
                    & (dnm_mt.AF > 0.45)
                    & (dnm_mt.AF < 0.55)
                    & (dnm_mt.cohort == cohort)
                )
                | (
                    (dnm_mt.alleles[0] != ref)
                    | (dnm_mt.alleles[1] != alt)
                    | (dnm_mt.cohort != cohort)
                )
            )

        dnm_mt = dnm_mt.annotate_rows(
            is_transversion=hl.is_transversion(dnm_mt.alleles[0], dnm_mt.alleles[1]),
            context=dnm_mt.locus.sequence_context(after=1),
            isCpG=(dnm_mt.locus.sequence_context(after=1) == "CG"),
        )
        dnm_mt = dnm_mt.annotate_rows(
            ref=hl.if_else(dnm_mt.isCpG, "CpG", dnm_mt.alleles[0]),
            alt=hl.if_else(
                (dnm_mt.isCpG & (dnm_mt.alleles[1] == "T")), "TpG", dnm_mt.alleles[1]
            ),
        )

        dnm_mt = dnm_mt.filter_cols(~dnm_mt.s.matches(excluded_samples))

        dnm_ht = dnm_mt.entries()
        dnm_ht = dnm_ht.filter(hl.is_defined(dnm_ht.GT))
        return dnm_ht

    @staticmethod
    @NamedCheckpoint.checkpoint(
        prefix=CACHE_DIR, name="dnm_df.pickle", output_type=FILE_TYPE.PICKLE
    )
    def load_dataframe(dnm_ht: hl.Table, has_graphtyper: bool = False):
        cols = [
            "iAF",
            "is_transversion",
            "isCpG",
            "ref",
            "alt",
            "AD",
            "DP",
            "GT",
            "p_de_novo",
            "mother",
            "father",
            "AF",
            "cohort",
        ]
        if has_graphtyper:
            cols.extend([
                "graphtyper",
                "graphtyper_evidence",
            ])
        dnm_df = dnm_ht.select(*cols).to_pandas()
        dnm_df["str"] = dnm_df["ref"] + ">" + dnm_df["alt"]
        return dnm_df

    @staticmethod
    @NamedCheckpoint.checkpoint(
        prefix=CACHE_DIR, name="metadata.pickle", output_type=FILE_TYPE.PICKLE
    )
    def load_metadata(
        radar: str,
        inova: str,
        cru: str,
        pilot: str = None,
        excluded_samples: str = r"R_(0048|0158|0290|0385|0180|0214|0066|0146|0273).*",
    ):
        return MsdnLoader.load_metadata(
            radar, inova, cru, pilot=pilot, excluded_samples=excluded_samples
        )

    @staticmethod
    @NamedCheckpoint.checkpoint(
        prefix=CACHE_DIR,
        name="downsample_age_match.pickle",
        output_type=FILE_TYPE.PICKLE,
    )
    def load_downsample_age_match(
        dnm_df: pandas.DataFrame, age_data: pandas.DataFrame, factor: int
    ):
        downsampled_metadata_df, matching = downsample_age_match(
            age_data, factor=factor, reference_data=dnm_df
        )
        downsampled_dnm_df = dnm_df[
            dnm_df["s"].isin(set(downsampled_metadata_df["rid"]))
        ]
        return downsampled_dnm_df, downsampled_metadata_df, matching

    @staticmethod
    @NamedCheckpoint.checkpoint(
        prefix=CACHE_DIR, name="isolated_dnm_df.pickle", output_type=FILE_TYPE.PICKLE
    )
    def load_isolated(
        dnm_df: pandas.DataFrame, msdn_df: pandas.DataFrame
    ) -> pandas.DataFrame:
        non_isolated_loci = set(msdn_df["locus"]) | set(msdn_df["previous.locus"])
        isolated_dnm_df = dnm_df[~dnm_df["locus"].isin(non_isolated_loci)]
        return isolated_dnm_df

    @staticmethod
    @NamedCheckpoint.checkpoint(
        prefix=CACHE_DIR,
        name="graphtyper_dnm_filtered.pickle",
        output_type=FILE_TYPE.PICKLE,
    )
    def load_graphtyper_filter(
        dnm_df: pandas.DataFrame,
        graphtyper_col: str = "graphtyper_evidence",
        expected_values: typing.Set[str] = set(["CONFIRMATION"]),
    ):
        dnm_df = dnm_df[dnm_df[graphtyper_col].isin(expected_values)]
        return dnm_df


class PhasingLoader(Loader):
    @staticmethod
    @NamedCheckpoint.checkpoint(
        prefix=CACHE_DIR, name="phasing.ht", output_type=FILE_TYPE.HAIL_TABLE
    )
    def load_phasing(path: str) -> hl.Table:
        phase_df = None
        with open(os.path.abspath(path), "rb") as pickle_f:
            phase_df = pickle.load(pickle_f)

        if phase_df is None:
            logger.error("failed to load phasing data from {}".format(path))
            sys.exit(1)

        phase_ht = hl.Table.from_pandas(phase_df)
        phase_ht = phase_ht.annotate(
            locus=hl.parse_locus(phase_ht.locus),
            alleles=phase_ht.alleles,
        )
        phase_ht = phase_ht.key_by("locus", "alleles", "sample")
        return phase_ht

    @staticmethod
    @NamedCheckpoint.checkpoint(
        prefix=CACHE_DIR,
        name="dnm_phase.pickle",
        output_type=FILE_TYPE.PICKLE,
    )
    def annotate_dnm_phase(dnm_ht: hl.Table, phase_ht: hl.Table) -> pandas.DataFrame:
        dnm_phase_ht = dnm_ht.annotate(
            is_phased=hl.if_else(
                hl.is_defined(phase_ht[(dnm_ht.locus, dnm_ht.alleles, dnm_ht.s)]),
                phase_ht[(dnm_ht.locus, dnm_ht.alleles, dnm_ht.s)].is_phased,
                hl.missing("bool"),
            ),
            is_paternal=hl.if_else(
                hl.is_defined(phase_ht[(dnm_ht.locus, dnm_ht.alleles, dnm_ht.s)]),
                phase_ht[(dnm_ht.locus, dnm_ht.alleles, dnm_ht.s)].is_paternal,
                hl.missing("bool"),
            ),
        )
        dnm_phase_ht = dnm_phase_ht.annotate(
            phasing=(
                hl.case()
                .when(dnm_phase_ht.is_phased & dnm_phase_ht.is_paternal, "PATERNAL")
                .when(dnm_phase_ht.is_phased, "MATERNAL")
                .default("UNPHASED")
            ),
        )
        dnm_phase_ht = dnm_phase_ht.persist(storage_level="MEMORY_ONLY")
        return dnm_phase_ht.to_pandas()

    @staticmethod
    @NamedCheckpoint.checkpoint(
        prefix=CACHE_DIR,
        name="msdn_phase.pickle",
        output_type=FILE_TYPE.PICKLE,
    )
    def annotate_msdn_phase(msdn_ht: hl.Table, phase_ht: hl.Table) -> pandas.DataFrame:
        msdn_ht = msdn_ht.key_by("locus", "alleles", "s")
        msdn_phase_ht = msdn_ht.annotate(
            is_phased=hl.if_else(
                hl.is_defined(phase_ht[(msdn_ht.locus, msdn_ht.alleles, msdn_ht.s)]),
                phase_ht[(msdn_ht.locus, msdn_ht.alleles, msdn_ht.s)].is_phased,
                hl.missing("bool"),
            ),
            is_paternal=hl.if_else(
                hl.is_defined(phase_ht[(msdn_ht.locus, msdn_ht.alleles, msdn_ht.s)]),
                phase_ht[(msdn_ht.locus, msdn_ht.alleles, msdn_ht.s)].is_paternal,
                hl.missing("bool"),
            ),
            previous=hl.Struct(
                **{
                    **msdn_ht.previous,
                    "is_phased": phase_ht[
                        (
                            msdn_ht.previous.locus,
                            msdn_ht.previous.alleles,
                            msdn_ht.previous.s,
                        )
                    ].is_phased,
                    "is_paternal": phase_ht[
                        (
                            msdn_ht.previous.locus,
                            msdn_ht.previous.alleles,
                            msdn_ht.previous.s,
                        )
                    ].is_paternal,
                }
            ),
        )
        msdn_phase_ht = msdn_phase_ht.annotate(
            phasing=(
                hl.case(missing_false=True)
                .when(
                    (msdn_phase_ht.is_phased & msdn_phase_ht.is_paternal)
                    | (
                        msdn_phase_ht.previous.is_phased
                        & msdn_phase_ht.previous.is_paternal
                    ),
                    "PATERNAL",
                )
                .when(
                    msdn_phase_ht.is_phased | msdn_phase_ht.previous.is_phased,
                    "MATERNAL",
                )
                .default("UNPHASED")
            )
        )
        msdn_phase_ht = msdn_phase_ht.persist(storage_level="MEMORY_ONLY")
        return msdn_phase_ht.to_pandas()

    @staticmethod
    @NamedCheckpoint.checkpoint(
        prefix=CACHE_DIR,
        name="phasing_downsample_age_match.pickle",
        output_type=FILE_TYPE.PICKLE,
    )
    def load_downsample_age_match(
        dnm_df: pandas.DataFrame,
        msdn_df: pandas.DataFrame,
        age_data: pandas.DataFrame,
        factor: int,
    ):
        downsampled_metadata_df, matching = downsample_age_match(
            age_data, factor=factor, reference_data=dnm_df
        )

        downsampled_dnm_df = dnm_df[
            dnm_df["s"].isin(set(downsampled_metadata_df["rid"]))
        ]
        downsampled_msdn_df = msdn_df[
            msdn_df["s"].isin(set(downsampled_metadata_df["rid"]))
        ]

        return (
            downsampled_dnm_df,
            downsampled_msdn_df,
            downsampled_metadata_df,
            matching,
        )


class GraphtyperLoader(Loader):
    @staticmethod
    @NamedCheckpoint.checkpoint(
        prefix=CACHE_DIR,
        name="graphtyper.pickle",
        output_type=FILE_TYPE.PICKLE,
    )
    def load_graphtyper(
        paths: typing.List[str], graphtyper_is_called_col: str = "graphtyper_called"
    ) -> pandas.DataFrame:
        def _load_pickle(path: str, is_called_col=graphtyper_is_called_col):
            with open(path, "rb") as f:
                df = pickle.load(f)
                return df[df[is_called_col]]

        df = pandas.concat(map(_load_pickle, paths))
        return df

    @staticmethod
    @NamedCheckpoint.checkpoint(
        prefix=CACHE_DIR,
        name="graphtyper.ht",
        output_type=FILE_TYPE.HAIL_TABLE,
    )
    def load_graphtyper_hail(
        graphtyper_df: pandas.DataFrame,
    ):
        graphtyper_ht = hl.Table.from_pandas(
            graphtyper_df[
                [
                    "locus",
                    "alleles",
                    "s",
                    "graphtyper_called",
                    "graphtyper_het",
                    "graphtyper_phased",
                ]
            ]
        )
        graphtyper_ht = graphtyper_ht.annotate(
            locus=hl.parse_locus(graphtyper_ht.locus),
        ).key_by("locus", "alleles", "s")
        return graphtyper_ht

    @staticmethod
    @NamedCheckpoint.checkpoint(
        prefix=CACHE_DIR,
        name="graphtyper_dnm.ht",
        output_type=FILE_TYPE.HAIL_TABLE,
    )
    def annotate_graphtyper_dnm(
        graphtyper_ht: hl.Table,
        dnm_ht: hl.Table,
    ) -> hl.Table:
        dnm_gt_ht = dnm_ht.annotate(
            graphtyper=hl.Struct(
                **{
                    "has_call": hl.if_else(
                        hl.is_defined(
                            graphtyper_ht[(dnm_ht.locus, dnm_ht.alleles, dnm_ht.s)]
                        ),
                        graphtyper_ht[
                            (dnm_ht.locus, dnm_ht.alleles, dnm_ht.s)
                        ].graphtyper_called,
                        hl.missing("bool"),
                    ),
                    "is_het": hl.if_else(
                        hl.is_defined(
                            graphtyper_ht[(dnm_ht.locus, dnm_ht.alleles, dnm_ht.s)]
                        ),
                        graphtyper_ht[
                            (dnm_ht.locus, dnm_ht.alleles, dnm_ht.s)
                        ].graphtyper_het,
                        hl.missing("bool"),
                    ),
                }
            ),
        )

        dnm_gt_ht = dnm_gt_ht.annotate(
            graphtyper_evidence=(
                hl.case(missing_false=True)
                .when(
                    dnm_gt_ht.graphtyper.has_call & dnm_gt_ht.graphtyper.is_het,
                    "CONFIRMATION",
                )
                .when(
                    dnm_gt_ht.graphtyper.has_call & ~dnm_gt_ht.graphtyper.is_het,
                    "CONTRADICTION",
                )
                .default("NO EVIDENCE")
            ),
        )
        return dnm_gt_ht

    @staticmethod
    @NamedCheckpoint.checkpoint(
        prefix=CACHE_DIR,
        name="graphtyper_msdn.ht",
        output_type=FILE_TYPE.HAIL_TABLE,
    )
    def annotate_graphtyper_msdn(
        graphtyper_ht: hl.Table,
        msdn_ht: hl.Table,
    ):
        msdn_ht = msdn_ht.key_by("locus", "alleles", "s")
        msdn_gt_ht = msdn_ht.annotate(
            graphtyper=hl.Struct(
                has_call=hl.if_else(
                    hl.is_defined(
                        graphtyper_ht[(msdn_ht.locus, msdn_ht.alleles, msdn_ht.s)]
                    ),
                    graphtyper_ht[
                        (msdn_ht.locus, msdn_ht.alleles, msdn_ht.s)
                    ].graphtyper_called,
                    hl.missing("bool"),
                ),
                is_het=hl.if_else(
                    hl.is_defined(
                        graphtyper_ht[(msdn_ht.locus, msdn_ht.alleles, msdn_ht.s)]
                    ),
                    graphtyper_ht[
                        (msdn_ht.locus, msdn_ht.alleles, msdn_ht.s)
                    ].graphtyper_het,
                    hl.missing("bool"),
                ),
            ),
            previous=hl.Struct(
                **{
                    **msdn_ht.previous,
                    "graphtyper": hl.Struct(
                        has_call=graphtyper_ht[
                            (
                                msdn_ht.previous.locus,
                                msdn_ht.previous.alleles,
                                msdn_ht.previous.s,
                            )
                        ].graphtyper_called,
                        is_het=graphtyper_ht[
                            (
                                msdn_ht.previous.locus,
                                msdn_ht.previous.alleles,
                                msdn_ht.previous.s,
                            )
                        ].graphtyper_het,
                    ),
                }
            ),
        )

        msdn_gt_ht = msdn_gt_ht.annotate(
            graphtyper_evidence=(
                hl.case(missing_false=True)
                .when(
                    (
                        msdn_gt_ht.graphtyper.has_call
                        & msdn_gt_ht.graphtyper.is_het
                        & msdn_gt_ht.previous.graphtyper.has_call
                        & msdn_gt_ht.previous.graphtyper.is_het
                    ),
                    "CONFIRMATION",
                )
                .when(
                    (msdn_gt_ht.graphtyper.has_call & msdn_gt_ht.graphtyper.is_het)
                    | (
                        msdn_gt_ht.previous.graphtyper.has_call
                        & msdn_gt_ht.previous.graphtyper.is_het
                    ),
                    "PARTIAL_CONFIRMATION",
                )
                .when(
                    (msdn_gt_ht.graphtyper.has_call & ~msdn_gt_ht.graphtyper.is_het)
                    | (
                        msdn_gt_ht.previous.graphtyper.has_call
                        & ~msdn_gt_ht.previous.graphtyper.is_het
                    ),
                    "CONTRADICTION",
                )
                .default("NO EVIDENCE")
            )
        )

        return msdn_gt_ht
