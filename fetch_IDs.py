import requests
import pandas as pd
import time

# ---- Your accession list ----
accessions = [
    "A1A4S6", "A1X283", "A7KAX9", "O00160", "O00305", "O00459", "O00499", "O14559", "O14936", "O15068",
    "O15117", "O15259", "O43150", "O43281", "O43295", "O43307", "O43586", "O43639", "O60229", "O60437",
    "O60504", "O60861", "O75044", "O75563", "O75791", "O75886", "O75962", "O75995", "O76041", "O94868",
    "O94875", "O94885", "O95049", "O95153", "P00519", "P02549", "P06239", "P06241", "P07947", "P07948",
    "P08631", "P09769", "P12931", "P14317", "P14598", "P15498", "P15924", "P16333", "P16885", "P19174",
    "P19878", "P20929", "P20936", "P27986", "P41240", "P42679", "P42680", "P42681", "P42684", "P42685",
    "P46108", "P46109", "P49418", "P51451", "P52735", "P54284", "P55345", "P56945", "P57075", "P62993",
    "P78352", "P80192", "P98171", "Q00013", "Q02641", "Q02779", "Q03001", "Q06187", "Q07157", "Q07912",
    "Q08289", "Q08881", "Q12774", "Q12929", "Q12959", "Q12965", "Q13239", "Q13368", "Q13387", "Q13402",
    "Q13470", "Q13588", "Q13625", "Q13813", "Q13882", "Q14155", "Q14168", "Q14185", "Q14247", "Q14511",
    "Q14847", "Q15052", "Q15080", "Q15149", "Q15642", "Q15700", "Q15811", "Q16584", "Q16674", "Q5HYK7",
    "Q5JRA6", "Q5T0N5", "Q5T2T1", "Q5TCX8", "Q5TCZ1", "Q5VST9", "Q5VV41", "Q5VWT5", "Q6PIF6", "Q6UXY1",
    "Q6XZF7", "Q6ZMT1", "Q6ZUM4", "Q7Z6B7", "Q7Z6J0", "Q86UR1", "Q86WN1", "Q86WV1", "Q8IVI9", "Q8IWW6",
    "Q8IZD9", "Q8IZP0", "Q8N157", "Q8N1I0", "Q8N2Y8", "Q8N3R9", "Q8N5V2", "Q8NFA2", "Q8TDM6", "Q8TE67",
    "Q8TE68", "Q8TEC5", "Q8TEJ3", "Q8TF17", "Q8TF42", "Q8WUF5", "Q8WV41", "Q92608", "Q92783", "Q92796",
    "Q92817", "Q92882", "Q92968", "Q96B97", "Q96DR7", "Q96HU1", "Q96JB8", "Q96JP2", "Q96KQ4", "Q96MF2",
    "Q96N96", "Q96PC5", "Q96QH2", "Q96RF0", "Q96RU3", "Q99469", "Q99961", "Q99962", "Q99963", "Q9BRR9",
    "Q9BVN2", "Q9BX66", "Q9BY11", "Q9BYB0", "Q9BYC5", "Q9H3Y6", "Q9H6Q3", "Q9H6R6", "Q9H6S3", "Q9H7D0",
    "Q9NQ75", "Q9NR46", "Q9NR80", "Q9NSI8", "Q9NYB9", "Q9NZM3", "Q9NZQ3", "Q9P0V3", "Q9P2A4", "Q9UDY2",
    "Q9UFD9", "Q9UHR4", "Q9UJU6", "Q9UKN7", "Q9UKS6", "Q9UKW4", "Q9ULH1", "Q9UNA1", "Q9UNF0", "Q9UPN3",
    "Q9UPX8", "Q9UQB8", "Q9UQF2", "Q9Y371", "Q9Y566", "Q9Y5K6", "Q9Y5X1", "A6NI28", "O15034", "Q6ZN28",
    "Q8IW93", "Q8WXD9", "Q96HL8", "Q9NXL2", "Q9NZW5", "A1IGU5", "A4FU49", "A6NI72", "A6NJZ7", "A6NNM3",
    "A8MVU1", "Q8TC17", "Q8WXE0", "Q9NRC9", "Q8TE82"
]

# ---- Function to fetch data in batches ----


def fetch_uniprot_data(accessions, batch_size=100):
    base_url = "https://rest.uniprot.org/uniprotkb/stream"
    all_rows = []

    for i in range(0, len(accessions), batch_size):
        batch = accessions[i:i+batch_size]
        query = " OR ".join([f"accession:{acc}" for acc in batch])

        params = {
            "format": "tsv",
            "fields": "accession,xref_interpro",
            "query": query
        }

        print(f"Fetching batch {i//batch_size + 1}...")
        response = requests.get(base_url, params=params)

        if response.status_code != 200:
            raise Exception(f"Request failed: {response.status_code}")

        lines = response.text.strip().split("\n")
        header = lines[0].split("\t")

        for line in lines[1:]:
            values = line.split("\t")
            row = dict(zip(header, values))
            all_rows.append(row)

        time.sleep(1)

    return pd.DataFrame(all_rows)


# ---- Fetch data ----
df = fetch_uniprot_data(accessions)

# ---- Clean + expand InterPro ----
df = df.rename(columns={
    "Entry": "ACCESSION",
    "InterPro": "INTERPRO"
})

# Drop rows with no InterPro annotation
df = df.dropna(subset=["INTERPRO"])

# Split multiple InterPro IDs into separate rows
df["INTERPRO"] = df["INTERPRO"].str.split(";")
df = df.explode("INTERPRO")

# Remove whitespace
df["INTERPRO"] = df["INTERPRO"].str.strip()

# ---- Count unique proteins per InterPro ----
result = (
    df.groupby("INTERPRO")["ACCESSION"]
    .nunique()
    .reset_index()
)

result = result.rename(columns={"ACCESSION": "PROTEIN_COUNT"})

# ---- Save output ----
result.to_csv("db_results.csv", index=False)

print("Done! File saved as db_results.csv")
