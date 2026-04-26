import pandas as pd

def reverse_complement(seq):
    comp = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
    return ''.join(comp.get(b, 'N') for b in reversed(seq.upper()))

df = pd.read_csv("data/processed/splice_sites_10k.csv")

print(f"Before fix - GT at donors: {(df[df['site_type']=='donor']['sequence'].str[200:202]=='GT').mean()*100:.1f}%")

# Reverse complement minus strand sequences
mask = df["strand"] == "-"
df.loc[mask, "sequence"] = df.loc[mask, "sequence"].apply(reverse_complement)

print(f"After fix  - GT at donors: {(df[df['site_type']=='donor']['sequence'].str[200:202]=='GT').mean()*100:.1f}%")
print(f"After fix  - AG at acceptors: {(df[df['site_type']=='acceptor']['sequence'].str[198:200]=='AG').mean()*100:.1f}%")

df.to_csv("data/processed/splice_sites_10k.csv", index=False)
print("Saved fixed sequences")
