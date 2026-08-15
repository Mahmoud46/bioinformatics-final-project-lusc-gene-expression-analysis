def categorize_overlap(is_paired, is_independent):
    if is_paired and is_independent:
        return 'Both (Concordant)'
    elif is_paired:
        return 'Paired-Only'
    elif is_independent:
        return 'Independent-Only'
    else:
        return 'Neither'