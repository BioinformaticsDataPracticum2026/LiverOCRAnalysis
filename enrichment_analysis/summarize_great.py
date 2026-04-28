import pandas as pd
from pathlib import Path

# folders for the GREAT results and the final summary table
GREAT_DIR = Path("results/enrichment/great")
OUTPUT = Path("results/enrichment/summary/great_summary.tsv")


def main():
    # store the top GO terms from each GREAT result folder
    all_dfs = []

    # go through each dataset folder inside the GREAT results directory
    for d in GREAT_DIR.iterdir():
        if not d.is_dir():
            continue

        # GREAT GO Biological Process result file
        file = d / "gobp.csv"
        if not file.exists():
            continue

        print(f"[PROCESS] {d.name}")

        # read the GREAT output table
        df = pd.read_csv(file)

        # use adjusted p-value
        df = df.sort_values("Binom_Adjp_BH").head(10)

        # keep track of which dataset this result came from
        df["dataset"] = d.name

        # keep only the columns needed for the summary
        all_dfs.append(
            df[["dataset", "name", "Binom_Adjp_BH", "Binom_Fold_Enrichment"]]
        )

    # combine all dataset summaries into one table
    final = pd.concat(all_dfs)

    # make the output folder if it does not exist
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)

    # save the combined summary table
    final.to_csv(OUTPUT, sep="\t", index=False)

    print(f"[DONE] saved → {OUTPUT}")


if __name__ == "__main__":
    main()