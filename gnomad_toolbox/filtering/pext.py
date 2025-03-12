"""Functions to filter the gnomAD pext HTs."""

from typing import Optional

import hail as hl

from gnomad_toolbox.filtering.variant import get_single_variant


def get_pext_for_variant(
    variant: Optional[str] = None,
    contig: Optional[str] = None,
    position: Optional[int] = None,
    ref: Optional[str] = None,
    alt: Optional[str] = None,
    **kwargs,
) -> hl.Table:
    """
    Get annotation-level pext score for a single variant from the gnomAD HT.

    .. note::

        The annotation-level pext score only sums the expression of transcripts on which a variant has a given annotation.
        For more information, please see https://gnomad.broadinstitute.org/help/pext.

    :param variant: Variant string in the format "chr12-235245-A-C" or
        "chr12:235245:A:C". If provided, `contig`, `position`, `ref`, and `alt` are
        ignored.
    :param contig: Chromosome of the variant. Required if `variant` is not provided.
    :param position: Variant position. Required if `variant` is not provided.
    :param ref: Reference allele. Required if `variant` is not provided.
    :param alt: Alternate allele. Required if `variant` is not provided.
    :param kwargs: Additional arguments to pass to `_get_dataset`.
    :return: Table with the single variant.
    """
    return get_single_variant(
        variant=variant,
        contig=contig,
        position=position,
        ref=ref,
        alt=alt,
        dataset="pext",
        data_type="annotation_level",
        **kwargs,
    )
