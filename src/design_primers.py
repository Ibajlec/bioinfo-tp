import argparse
import json
import os
import sys

from Bio import SeqIO
import primer3

def load_config(path: str) -> dict:
    """Load configuration from a JSON file."""
    with open(path) as fh:
        return json.load(fh)

def load_seq(fastafile: str) -> tuple[str, str]:
    """Return the sequence and its identifier from a FASTA file."""
    record = next(SeqIO.parse(fastafile, "fasta"))
    return str(record.seq), record.id

def design(cfg: dict, seq: str, sid: str) -> list[dict]:
    """Return a list of primers for the given sequence."""

    seq_args = {"SEQUENCE_ID": sid, "SEQUENCE_TEMPLATE": seq}
    global_args = {
        "PRIMER_TASK": "generic",
        "PRIMER_NUM_RETURN": cfg["PRIMER_NUM_RETURN"],
        "PRIMER_MIN_SIZE": cfg["PRIMER_MIN_SIZE"],
        "PRIMER_OPT_SIZE": (cfg["PRIMER_MIN_SIZE"] + cfg["PRIMER_MAX_SIZE"]) // 2,
        "PRIMER_MAX_SIZE": cfg["PRIMER_MAX_SIZE"],
        "PRIMER_MIN_GC": cfg["PRIMER_MIN_GC_PERCENT"],
        "PRIMER_MAX_GC": cfg["PRIMER_MAX_GC_PERCENT"],
        "PRIMER_MAX_TM": cfg["PRIMER_MAX_TM"],
    }
    if cfg.get("AVOID_GC_AT_3_PRIME"):
        global_args["PRIMER_3_END_GC_CLAMP"] = 0

    res = primer3.bindings.designPrimers(seq_args, global_args)

    primers = []
    for i in range(cfg["PRIMER_NUM_RETURN"]):
        primers.append({
            "id": sid,
            "forward": res[f"PRIMER_LEFT_{i}_SEQUENCE"],
            "tm_fwd": res[f"PRIMER_LEFT_{i}_TM"],
            "reverse": res[f"PRIMER_RIGHT_{i}_SEQUENCE"],
            "tm_rev": res[f"PRIMER_RIGHT_{i}_TM"],
            "product": res[f"PRIMER_PAIR_{i}_PRODUCT_SIZE"],
        })
    return primers

def main(argv: list[str] | None = None) -> None:
    """CLI entry point."""
    parser = argparse.ArgumentParser(description="Design PCR primers")
    parser.add_argument("fasta", help="Input FASTA file")
    parser.add_argument("params", help="JSON file with primer3 parameters")
    parser.add_argument("output", help="Directory to store the results")
    args = parser.parse_args(argv)

    seq, sid = load_seq(args.fasta)
    cfg = load_config(args.params)
    primers = design(cfg, seq, sid)

    os.makedirs(args.output, exist_ok=True)
    base = os.path.splitext(os.path.basename(args.fasta))[0]
    out_fn = os.path.join(args.output, f"resultados_{base}.json")

    with open(out_fn, "w") as fh:
        json.dump({sid: primers}, fh, indent=2)

    print(f"→ {len(primers)} cebadores guardados en {out_fn}")


if __name__ == "__main__":
    try:
        main()
    except FileNotFoundError as exc:
        sys.exit(f"Error: {exc}")
    except Exception as exc:  # pragma: no cover - unexpected errors
        sys.exit(f"Unexpected error: {exc}")

  
