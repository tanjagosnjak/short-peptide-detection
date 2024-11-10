import numpy as np
import pandas as pd
import re

print('okk')

def read_files(filenames, sep="\t"):
    return [
        pd.read_csv(filename, sep=sep).assign(sample=filename.split('\\')[-1].split('.')[0]) for filename in filenames
    ]


mods_dict = {
    "T(-14.0156)": "T->S",
    '[1 nucleotide substitution:T->S on T]': 'T->S',
    '(Thr->Ser)': 'T->S',
    "cC": "C_carbamidomethylation",
    "C(57.0214)": "C_carbamidomethylation",
    "[Common Fixed:Carbamidomethyl on C]": "C_carbamidomethylation",
    "(Oxidation (M))": "M_oxidation",
    "(Carbamidomethylation (C))": "C_carbamidomethylation",
    "M(15.9949)": "M_oxidation",
    "oxM": "M_oxidation",
    "[Common Variable:Oxidation on M]": "M_oxidation",
    "[Common Biological:Acetylation on X]": 'Acetylation'
}


# define function for removing duplicates and sorting
def _sort_and_remove_duplicates(protein):
    proteins = protein.split(", ")
    unique_proteins = sorted(set(proteins))
    return ", ".join(unique_proteins)


# Define MSF parser
def MSF_parser(df, mapper=mods_dict):
    # Select and rename cols
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

    # Add mods col
    df["mods"] = df["mods_type"].apply(lambda x: "modified" if isinstance(x, str) and x != "nan" else "unmodified")

    # Combine proteinID and mapped proteins
    df['mapped_proteins'] = df['mapped_proteins'].fillna('')
    df['proteinID'] = df['proteinID'].str.cat(df['mapped_proteins'], sep=', ', na_rep='').str.strip(', ')

    # Remove '-p1', remove duplicates and sort proteinID col
    df["proteinID"] = df["proteinID"].str.replace("-p1", "").apply(_sort_and_remove_duplicates)

    # Use mapper to get wished output
    for key, val in mapper.items():
        df["mods_type"] = df["mods_type"].str.replace(key, val)
    df["mods_type"] = df["mods_type"].astype(str)

    # Order mods_type
    def _mods_type_sorter(val: str):
        def extract_number(el):
            match = re.search("^\d+", el)
            return int(match.group(0)) if match else float('inf')

        return ", ".join(sorted(val.split(", "), key=extract_number))

    df["mods_type"] = df["mods_type"].apply(_mods_type_sorter)

    return df.reset_index(drop=True)



# define MM parser
def MM_parser(df, mapper=mods_dict):

    # add new col from df_all, select and rename cols
    df["MM_score"] = 'mm_score'
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
    
    # Apply the extraction function
    def _mods_type_extract(val):
        res = []
        while re.search(r"(\(\w+\s+\([A-Z]\)\))|(\(\w+->\w+\))", val):
            m = re.search(r"(\(\w+\s+\([A-Z]\)\))|(\(\w+->\w+\))", val)
            res.append(f"{m.start()}{mapper[m.group(0)]}")
            val = val[:m.start()] + val[m.end():]

        if not res:
            return "nan"

        return ", ".join(res)

    df["mods_type"] = df["mods_type"].apply(_mods_type_extract)

    # select max values for those with matching keys
    df = df.groupby(["sample", "peptide", "proteinID", "mods", "mods_type"], as_index=False)["MQ_score"].max()
    return df