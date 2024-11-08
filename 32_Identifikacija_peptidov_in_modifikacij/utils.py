import numpy as np
import pandas as pd
import re


def read_files(filenames, sep="\t"):
    return [
        pd.read_csv(filename, sep=sep).assign(sample=filename.split('\\')[-1].split('.')[0]) for filename in filenames
    ]


mods_dict = {
    "cC": "C_carbamidomethylation",
    "C(57.0214)": "C_carbamidomethylation",
    "[Common Fixed:Carbamidomethyl on C]": "C_carbamidomethylation",
    "(Oxidation (M))": "M_oxidation",
    "(Carbamidomethylation (C))": "C_carbamidomethylation",
    "M(15.9949)": "M_oxidation",
    "oxM": "M_oxidation",
    "[Common Variable:Oxidation on M]": "M_oxidation",
}


# define function for removing duplicates and sorting
def _sort_and_remove_duplicates(protein):
    proteins = protein.split(", ")
    unique_proteins = sorted(set(proteins))
    return ", ".join(unique_proteins)


# define MSF parser
def MSF_parser(df, mapper=mods_dict):

    # select and rename cols
    df = df[["sample", "Peptide", "PeptideProphet Probability", "Assigned Modifications", "Protein", "Mapped Proteins"]]
    df = df.rename(
        columns={
            "Peptide": "peptide",
            "PeptideProphet Probability": "MSF_probability",
            "Assigned Modifications": "mods_type",
            "Protein": "proteinID",
            "Mapped Proteins": 'mapped_proteins'
        }
    )

    # add mods col
    df["mods"] = df["mods_type"].apply(lambda x: "modified" if isinstance(x, str) else "unmodified")

    # combine proteinID and mapped proteins
    df['mapped_proteins'] = df['mapped_proteins'].fillna('')
    df['proteinID'] = df['proteinID'].str.cat(df['mapped_proteins'], sep=', ', na_rep='').str.strip(', ')

    # remove '-p1', remove duplicates and sort proteinID col
    df["proteinID"] = df["proteinID"].str.replace("-p1", "").apply(_sort_and_remove_duplicates)

    # use mapper to get wished output
    for key, val in mapper.items():
        df["mods_type"] = df["mods_type"].str.replace(key, val)
    df["mods_type"] = df["mods_type"].astype(str)

    # order mods_type
    def _mods_type_sorter(val: str):
        return ", ".join(
            sorted(val.split(", "), key=lambda el: None if el == "nan" else int(re.search("^\d+", el).group(0)))
        )

    df["mods_type"] = df["mods_type"].apply(_mods_type_sorter)

    # select max values for those with matching keys
    df = df.groupby(["sample", "peptide", "proteinID", "mods", "mods_type"], as_index=False)["MSF_probability"].max()

    return df.reset_index(drop=True)


# define MM parser
def MM_parser(df, df_all, mapper=mods_dict):

    # add new col from df_all, select and rename cols
    df["MM_score"] = df["Sequence"].map(df_all.set_index("Full Sequence")["Score"])
    df = df[["sample", "MM_score", "Base Sequence", "Sequence", "Protein Groups"]]
    df = df.rename(columns={"Base Sequence": "peptide", "Sequence": "mods_type", "Protein Groups": "proteinID"})

    # add mods col
    df["mods"] = df.apply(lambda row: "modified" if row["mods_type"] != row["peptide"] else "unmodified", axis=1)

    # remove '-p1', remove duplicates and sort proteinID col
    df["proteinID"] = (
        df["proteinID"]
        .str.replace("-p1", "")
        .str.replace(" | ", ", ")
        .str.replace("|", ", ")
        .str.replace(";", ", ")
        .apply(_sort_and_remove_duplicates)
    )

    # order mods type
    def _mods_type_extract(val):
        res = []
        while re.search("(\[[^\[\]]+\])", val):
            m = re.search("(\[[^\[\]]+\])", val)
            res.append(f"{m.start()}{mapper[m.group(1)]}")
            val = val[: m.start()] + val[m.end() :]

        if not res:
            return "nan"

        return ", ".join(res)

    df["mods_type"] = df["mods_type"].apply(_mods_type_extract)

    # select max values for those with matching keys
    df = df.groupby(["sample", "peptide", "proteinID", "mods_type", "mods"], as_index=False)["MM_score"].max()

    return df


# define MQ parser
def MQ_parser(df, mapper=mods_dict):

    # select and rename cols
    df = df[["sample", "Sequence", "Modified sequence", "Score", "Proteins"]]
    df = df.rename(
        columns={"Sequence": "peptide", "Score": "MQ_score", "Modified sequence": "mods_type", "Proteins": "proteinID"}
    )

    # remove contaminants
    df = df.drop(df[df["proteinID"].str.contains("CON__")].index)

    # add mods col
    df["mods_type"] = df["mods_type"].str.replace("_", "")
    df["mods"] = df.apply(lambda row: "unmodified" if row["mods_type"] == row["peptide"] else "modified", axis=1)

    # remove '-p1', remove duplicates and sort proteinID col
    df["proteinID"] = df["proteinID"].str.replace("-p1", "").str.replace(";", ", ").apply(_sort_and_remove_duplicates)

    # order mods type
    def _mods_type_extract(val):
        res = []
        while re.search("(\(\w+\s+\([A-Z]\)\))", val):
            m = re.search("(\(\w+\s+\([A-Z]\)\))", val)
            res.append(f"{m.start()}{mapper[m.group(0)]}")
            val = val[: m.start()] + val[m.end() :]

        if not res:
            return "nan"

        return ", ".join(res)

    df["mods_type"] = df["mods_type"].apply(_mods_type_extract)

    # select max values for those with matching keys
    df = df.groupby(["sample", "peptide", "proteinID", "mods", "mods_type"], as_index=False)["MQ_score"].max()
    return df


# define AP parser
def AP_parser(df, mapper=mods_dict):

    # select and rename cols
    df = df[["sample", "sequence", "sequence_naked", "score", "protein"]]
    df = df.rename(
        columns={"sequence_naked": "peptide", "score": "AP_score", "sequence": "mods_type", "protein": "proteinID"}
    )

    # add mods col
    df["mods"] = df.apply(lambda row: "modified" if row["mods_type"] != row["peptide"] else "unmodified", axis=1)

    # remove '-p1', remove duplicates and sort proteinID col
    df["proteinID"] = df["proteinID"].str.replace("-p1", "").str.replace(",", ", ").apply(_sort_and_remove_duplicates)

    # order mods type
    def _mods_type_extract(val):
        res = []
        while re.search("([a-z]+[A-Z])", val):
            m = re.search("([a-z]+[A-Z])", val)
            res.append(f"{m.start()+1}{mapper[m.group(1)]}")
            val = val[: m.start()] + val[m.end() - 1 :]

        if not res:
            return "nan"

        return ", ".join(res)

    df["mods_type"] = df.mods_type.apply(_mods_type_extract)

    # select max values for those with matching keys
    df = df.groupby(["sample", "peptide", "proteinID", "mods", "mods_type"], as_index=False)["AP_score"].max()
    return df
