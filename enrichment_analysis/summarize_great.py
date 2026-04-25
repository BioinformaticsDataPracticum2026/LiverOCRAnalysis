import pandas as pd
from pathlib import Path

GREAT_DIR = Path("results/enrichment/great")
OUTPUT = Path("results/enrichment/summary/great_summary.tsv")


def main():
    all_dfs = []

    for d in GREAT_DIR.iterdir():
        if not d.is_dir():
            continue

        file = d / "gobp.csv"
        if not file.exists():
            continue

        print(f"[PROCESS] {d.name}")

        df = pd.read_csv(file)

        # use adjusted p-value
        df = df.sort_values("Binom_Adjp_BH").head(10)

        df["dataset"] = d.name

        all_dfs.append(
            df[["dataset", "name", "Binom_Adjp_BH", "Binom_Fold_Enrichment"]]
        )

    final = pd.concat(all_dfs)

    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    final.to_csv(OUTPUT, sep="\t", index=False)

    print(f"[DONE] saved → {OUTPUT}")


if __name__ == "__main__":
    main()
