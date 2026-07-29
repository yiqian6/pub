#!/usr/bin/env python3
"""
Convert a phased VCF to a dated ARG (ancestral recombination graph) using tsinfer and tsdate.
Additionally, generate a sample mapping file (tree_node_id -> sample) for TRACE.

This script performs the following steps:
  1. Extract sample IDs from the VCF header.
  2. Write a sample map TSV file (optional, default alongside the output .trees file).
  3. Convert the VCF to Zarr format using vcf2zarr (if not already present).
  4. Add the REF_allele array to the Zarr store if missing.
  5. Run tsinfer to infer an ARG.
  6. Run tsdate to date the ARG and output a final .trees file.

The output .trees file is the primary input for TRACE. The sample map is a separate
auxiliary file that can be passed to TRACE via the TRACE_SAMPLE_MAP environment variable.
"""

import argparse
import gzip
from pathlib import Path
import shutil
import subprocess


def extract_sample_ids(vcf_path: str) -> list:
    """
    Extract sample IDs from the VCF header (#CHROM line).

    Parameters
    ----------
    vcf_path : str
        Path to the input VCF file (can be .vcf or .vcf.gz).

    Returns
    -------
    list of str
        List of sample IDs in the order they appear in the VCF.

    Raises
    ------
    SystemExit
        If no sample IDs are found in the VCF header.
    """
    samples = []
    # Choose the appropriate file opener based on file extension.
    open_func = gzip.open if str(vcf_path).endswith(".gz") else open
    with open_func(vcf_path, "rt") as f:
        for line in f:
            if line.startswith("#CHROM"):
                parts = line.strip().split("\t")
                # The first 9 columns are fixed (chrom, pos, id, ref, alt, qual, filter, info, format)
                samples = parts[9:]
                break
    if not samples:
        raise SystemExit(f"No sample IDs found in VCF header: {vcf_path}")
    return samples


def write_sample_map(sample_ids: list, output_path: str) -> None:
    """
    Write a sample mapping file (TSV) that maps tree_node_id (0-based) to sample ID.

    The mapping is created by assigning each sample an integer index starting from 0,
    which corresponds to the order of samples in the VCF. This order is preserved in
    the ARG's sample nodes, making this mapping valid for TRACE.

    Parameters
    ----------
    sample_ids : list of str
        List of sample IDs in VCF order.
    output_path : str
        Path where the sample map TSV file will be written.

    Returns
    -------
    None
    """
    with open(output_path, "w") as f:
        f.write("tree_node_id\tsample\n")
        for idx, samp in enumerate(sample_ids):
            f.write(f"{idx}\t{samp}\n")
    print(f"Sample map written to {output_path}")


def main():
    """
    Main entry point: parse arguments, generate sample map, and run ARG inference.

    The script processes a phased VCF file to produce a dated .trees file for TRACE.
    It also optionally generates a sample map file to ensure accurate sample mapping
    in downstream TRACE analyses.
    """
    parser = argparse.ArgumentParser(
        description=(
            "Prepare a dated .trees input for TRACE from a phased VCF. "
            "Also generates a sample_map.tsv for TRACE sample mapping."
        )
    )
    parser.add_argument(
        "--vcf",
        required=True,
        help="Input phased VCF file (can be .vcf or .vcf.gz)",
    )
    parser.add_argument(
        "--out",
        required=True,
        help="Output .trees file (dated ARG)",
    )
    parser.add_argument(
        "--sample-map",
        default=None,
        help=(
            "Output sample map file (default: <out_prefix>.sample_map.tsv). "
            "If not specified, the sample map is written alongside the .trees file."
        ),
    )
    parser.add_argument(
        "--ancestral-state",
        default="REF_allele",
        help="Name of the ancestral state array in Zarr (default: REF_allele)",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=8,
        help="Number of threads for tsinfer inference (default: 8)",
    )
    parser.add_argument(
        "--mutation-rate",
        type=float,
        default=1.25e-8,
        help="Mutation rate per generation per bp (default: 1.25e-8)",
    )
    parser.add_argument(
        "--population-size",
        type=float,
        default=30000,
        help="Effective population size for tsdate (default: 30000)",
    )
    args = parser.parse_args()

    # Import required packages (may be heavy, so import after argument parsing)
    import numpy as np
    import tsdate
    import tsinfer
    import zarr

    # -------------------- Step 1: Extract sample IDs from VCF --------------------
    sample_ids = extract_sample_ids(args.vcf)

    # -------------------- Step 2: Generate sample map --------------------
    # Determine the output path for the sample map.
    if args.sample_map is None:
        out_path = Path(args.out)
        sample_map_path = out_path.parent / (out_path.stem + ".sample_map.tsv")
    else:
        sample_map_path = Path(args.sample_map)
    # Ensure the directory exists.
    sample_map_path.parent.mkdir(parents=True, exist_ok=True)
    write_sample_map(sample_ids, sample_map_path)

    # -------------------- Step 3: Prepare output paths for ARG --------------------
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    vcz_path = out.with_suffix(".vcz")          # Zarr store for tsinfer
    inferred_path = out.with_suffix(".inferred.trees")  # Un-dated ARG

    # -------------------- Step 4: Convert VCF to Zarr if needed --------------------
    if not vcz_path.exists():
        converter = shutil.which("vcf2zarr")
        if converter is None:
            raise SystemExit(
                "vcf2zarr is missing; install bio2zarr in conda environment gu"
            )
        subprocess.run(
            [converter, "convert", args.vcf, str(vcz_path)],
            check=True,
        )

    # -------------------- Step 5: Add REF_allele if missing --------------------
    # The REF_allele dataset is required by tsinfer as the ancestral state.
    root = zarr.open(str(vcz_path), mode="r+")
    if args.ancestral_state not in root:
        alleles = root["variant_allele"][:]
        refs = np.asarray([str(row[0]) for row in alleles], dtype="U50")
        root.create_dataset(
            args.ancestral_state,
            data=refs,
            shape=refs.shape,
            dtype=refs.dtype,
        )

    # -------------------- Step 6: Run tsinfer to infer the ARG --------------------
    # VariantData uses the VCF data and the ancestral state we just added.
    data = tsinfer.VariantData(str(vcz_path), ancestral_state=args.ancestral_state)
    inferred_ts = tsinfer.infer(data, num_threads=args.threads)
    inferred_ts.dump(str(inferred_path))

    # -------------------- Step 7: Simplify and date the ARG --------------------
    # Remove unary nodes (nodes with only one child) to simplify the tree.
    clean_ts = inferred_ts.simplify(keep_unary=False)
    # Date the tree using the inside-outside method.
    dated_ts = tsdate.date(
        clean_ts,
        mutation_rate=args.mutation_rate,
        population_size=args.population_size,
        method="inside_outside",
    )
    dated_ts.dump(str(out))

    # -------------------- Final output messages --------------------
    print(f"Dated ARG written to {out}")
    print(f"Sample map written to {sample_map_path}")


if __name__ == "__main__":
    main()
