import pandas as pd

df = pd.read_csv("SNP.csv", header=None)

print(df.iloc[:5, :10])