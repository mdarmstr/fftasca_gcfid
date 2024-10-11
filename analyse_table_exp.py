import pandas as pd

def create_vector(filename):
    vector = [0, 0, 0, 0]  # Initialize vector with zeros
    
    # Condition for "Blank"
    if "Blank" in filename:
        return vector
    
    # Condition for "24h" and "72h"
    if "24" in filename:
        vector[0] = 1
    elif "72" in filename:
        vector[0] = 2
    
    # Condition for 5th entry from the end
    if filename[-6] == 'N':
        vector[1] = 1
    elif filename[-6] == 'P':
        vector[1] = 2
    elif filename[-6] == 'W':
        vector[1] = 3
    
    # Condition for 4th entry from the end
    if filename[-5] == 'm':
        vector[2] = 1

    vector[3] = int(filename[-4])
    
    return vector

csv_file = 'table_exp.csv'

# Read the CSV file without a header
df = pd.read_csv(csv_file, header=None)

# Convert the first column to a list
filenames = df[0].tolist()

print(filenames)

vectors = {}

for file in filenames:
    vectors[file] = create_vector(file)

vectors_df = pd.DataFrame.from_dict(vectors, orient='index', columns=['Position 1', 'Position 2', 'Position 3', 'Position 4'])

print(vectors_df)

vectors_df.to_csv('Y_peak.csv',index=False, header=False)
