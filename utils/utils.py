def categorize(row):
    if row['Is_DEG_Paired'] and row['Is_DEG_Indep']:
        return 'Both (Concordant)'
    elif row['Is_DEG_Paired']:
        return 'Paired-Only'
    elif row['Is_DEG_Indep']:
        return 'Independent-Only'
    else:
        return 'Neither'